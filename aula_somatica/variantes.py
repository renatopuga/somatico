"""Entrada VCF, checagem de referência e métricas por alelo e amostra."""

from collections import Counter
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory


CONTIGS_PRIMARIOS = {str(chrom) for chrom in range(1, 23)} | {"X", "Y", "M", "MT"}
TAMANHOS = {
    "GRCh37": {"1": 249250621, "9": 141213431},
    "GRCh38": {"1": 248956422, "9": 138394717},
}
SEQUENCIA = re.compile(r"[ACGT]+")
FORMATOS_USADOS = {"DP", "AD", "AF"}
FORMATO_HEADER = re.compile(r"##FORMAT=<ID=([^,]+),Number=([^,]+)")
AS_FILTER = "AULA_AS_FILTER"


def contig_canonico(contig):
    sem_prefixo = contig.removeprefix("chr")
    if sem_prefixo in {"M", "MT"}:
        return "MT"
    return sem_prefixo if sem_prefixo in CONTIGS_PRIMARIOS else contig


def montagem_canonica(montagem):
    valor = montagem.strip().lower()
    if valor in {"grch37", "hg19", "b37"} or re.fullmatch(r"grch37\.p\d+", valor):
        return "GRCh37"
    if valor in {"grch38", "hg38"} or re.fullmatch(r"grch38\.p\d+", valor):
        return "GRCh38"
    raise ValueError(f"Montagem não suportada: {montagem!r}. Use GRCh37 ou GRCh38.")


@dataclass(frozen=True)
class Alelo:
    montagem: str
    cromossomo: str
    posicao: int
    ref: str
    alt: str

    @property
    def chave(self):
        return f"{self.montagem}:{self.cromossomo}:{self.posicao}:{self.ref}:{self.alt}"

    @property
    def sem_ancora(self):
        """Coordenada da primeira base alterada ou do flanco direito de uma inserção."""
        pos, ref, alt = self.posicao, self.ref, self.alt
        while ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        while ref and alt and ref[0] == alt[0]:
            pos, ref, alt = pos + 1, ref[1:], alt[1:]
        return self.montagem, self.cromossomo, pos, ref, alt


def abrir_vcf(caminho):
    caminho = Path(caminho)
    with caminho.open("rb") as f:
        compactado = f.read(2) == b"\x1f\x8b"
    return gzip.open(caminho, "rt") if compactado else caminho.open()


def ler_header(stream):
    metadados = []
    for linha in stream:
        if linha.startswith("##"):
            metadados.append(linha.strip())
        elif linha.startswith("#CHROM\t"):
            return metadados, linha.rstrip("\r\n").split("\t")[9:]
        else:
            break
    raise ValueError("O arquivo precisa ser um VCF com cabeçalho #CHROM.")


def validar_montagem(montagem, contigs, origem):
    for nome, tamanho in contigs.items():
        for build, comprimentos in TAMANHOS.items():
            if build != montagem and comprimentos.get(contig_canonico(nome)) == tamanho:
                raise ValueError(f"{origem}: comprimento de {nome} indica {build}, mas foi informado {montagem}.")


def verificar_vcf(vcf, montagem, tumor, normal):
    with abrir_vcf(vcf) as f:
        metadados, amostras = ler_header(f)
    if not tumor or tumor not in amostras or (normal and normal not in amostras) or tumor == normal:
        raise ValueError(f"Confira os nomes de tumor e normal no VCF. Amostras encontradas: {', '.join(amostras)}.")
    contigs = {}
    for linha in metadados:
        if linha.startswith("##contig=<"):
            nome = re.match(r"##contig=<ID=([^,>]+)", linha)
            tamanho = re.search(r",length=(\d+)([,>])", linha, re.IGNORECASE)
            if nome and tamanho:
                contigs[nome[1]] = int(tamanho[1])
        if linha.startswith(("##reference=", "##assembly=")):
            build = re.search(r"(?<![A-Za-z0-9])(GRCh37|hg19|b37|GRCh38|hg38)(?![A-Za-z0-9])", linha, re.IGNORECASE)
            if build and montagem_canonica(build[1]) != montagem:
                raise ValueError(f"A referência declarada no cabeçalho do VCF indica {build[1]}, mas foi informado {montagem}.")
    validar_montagem(montagem, contigs, "VCF")
    return metadados, contigs


def valor_info(info, tag):
    campo = re.search(rf"(?:^|;){tag}=([^;]+)", info)
    return campo[1] if campo else ""


def preparar_vcf_para_norm(vcf, destino):
    # O VCF de trabalho usa somente os campos da aula. O original conserva as
    # demais anotações, inclusive tags específicas de cada versão do caller.
    with abrir_vcf(vcf) as origem, Path(destino).open("w", encoding="utf-8") as trabalho:
        for linha in origem:
            if linha.startswith("##INFO="):
                continue
            if linha.startswith("##FORMAT="):
                campo = FORMATO_HEADER.match(linha)
                if not campo or campo[1] not in FORMATOS_USADOS:
                    continue
            if linha.startswith("#CHROM\t"):
                trabalho.write(f'##INFO=<ID={AS_FILTER},Number=A,Type=String,Description="Filtros por alelo do caller">\n')
            if linha.startswith("#"):
                trabalho.write(linha)
                continue
            if not linha.strip():
                continue
            campos = linha.rstrip("\r\n").split("\t")
            as_filter = valor_info(campos[7], "AS_FilterStatus")
            campos[7] = "."
            if as_filter and as_filter != ".":
                n_allelos = len(campos[4].split(","))
                filtros = as_filter.split("|")
                if len(filtros) != n_allelos:
                    alternativo = as_filter.split(",")
                    filtros = alternativo if len(alternativo) == n_allelos else ["indeterminado"] * n_allelos
                campos[7] = f"{AS_FILTER}=" + ",".join(f.strip().replace(",", "&") or "." for f in filtros)
            formato = campos[8].split(":")
            indices = [i for i, tag in enumerate(formato) if tag in FORMATOS_USADOS]
            campos[8] = ":".join(formato[i] for i in indices) or "."
            for j in range(9, len(campos)):
                valores = campos[j].split(":")
                campos[j] = ":".join(valores[i] if i < len(valores) else "." for i in indices) or "."
            trabalho.write("\t".join(campos) + "\n")


def normalizar(vcf, referencia, montagem, tumor, normal, saida):
    montagem = montagem_canonica(montagem)
    _, contigs_vcf = verificar_vcf(vcf, montagem, tumor, normal)
    referencia = Path(referencia)
    fai = Path(f"{referencia}.fai")
    if not fai.exists() or fai.stat().st_mtime < referencia.stat().st_mtime:
        subprocess.run(["samtools", "faidx", str(referencia)], check=True)
    with fai.open() as f:
        contigs_fasta = {campos[0]: int(campos[1]) for linha in f if (campos := linha.split())}
    validar_montagem(montagem, contigs_fasta, "FASTA")
    for nome in contigs_vcf.keys() & contigs_fasta.keys():
        if contigs_vcf[nome] != contigs_fasta[nome]:
            raise ValueError(f"VCF e FASTA declaram comprimentos diferentes para {nome}.")
    with TemporaryDirectory(dir=Path(saida).parent) as pasta:
        entrada = Path(pasta) / "metricas.vcf"
        preparar_vcf_para_norm(vcf, entrada)
        comando = ["bcftools", "norm", "--no-version", "-m", "-any", "--check-ref", "e", "-f", str(referencia),
                   "--old-rec-tag", "AULA_ORIG", "-Ov", "-o", str(saida), str(entrada)]
        try:
            subprocess.run(comando, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as erro:
            Path(saida).unlink(missing_ok=True)
            raise ValueError(f"Falha ao normalizar. Confira REF e nomes dos contigs no VCF e no FASTA.\n{erro.stderr.strip()}") from erro
    return Path(saida)


def inteiro(texto):
    return int(texto) if texto.isascii() and texto.isdecimal() else None


def metricas(formato, amostra, numeros, multialelico):
    valores = dict(zip(formato.split(":"), amostra.split(":")))
    dp = inteiro(valores.get("DP", ""))
    ad = valores.get("AD", "").split(",")
    ref, alt = (inteiro(ad[0]), inteiro(ad[1])) if numeros.get("AD") == "R" and len(ad) == 2 else (None, None)
    af, fonte = None, ""
    if numeros.get("AF") == "A" or (not multialelico and numeros.get("AF") in {"1", "."}):
        try:
            valor = float(valores.get("AF", ""))
        except ValueError:
            valor = None
        if valor is not None and 0 <= valor <= 1:
            af, fonte = valor, "FORMAT/AF"
    if af is None and dp and alt is not None and alt <= dp:
        af, fonte = alt / dp, "FORMAT/AD_ALT ÷ FORMAT/DP"
    return {"dp": dp, "ad_ref": ref, "ad_alt": alt, "vaf": af, "vaf_fonte": fonte}


def ler_chamadas(vcf, montagem, tumor, normal=None):
    """Ler um VCF bialélico normalizado, mantendo os filtros originais."""
    montagem = montagem_canonica(montagem)
    metadados, _ = verificar_vcf(vcf, montagem, tumor, normal)
    numeros = dict(FORMATO_HEADER.findall("\n".join(metadados)))
    fontes = [linha for linha in metadados if linha.startswith(("##source=", "##DRAGENVersion=", "##gatkVersion="))]
    origens = {"tumor": tumor, "normal": normal}
    contagem = Counter()
    chamadas = {}
    with abrir_vcf(vcf) as f:
        _, amostras = ler_header(f)
        indices = {origem: amostras.index(nome) + 9 for origem, nome in origens.items() if nome}
        for linha in f:
            if not linha.strip():
                continue
            campos = linha.rstrip("\r\n").split("\t")
            if len(campos) != 9 + len(amostras):
                raise ValueError("Linha do VCF com número inesperado de colunas.")
            chrom, pos, _, ref, alt, qual, filtros, info, formato = campos[:9]
            contagem["registros"] += 1
            if "," in alt:
                raise ValueError("O VCF ainda contém registros multialélicos. Execute a normalização.")
            if not SEQUENCIA.fullmatch(ref) or not SEQUENCIA.fullmatch(alt) or ref == alt or len(ref) == len(alt) > 1:
                contagem["fora_escopo_snv_indel"] += 1
                continue
            original = valor_info(info, "AULA_ORIG")
            original_campos = original.split("|")
            multialelico = len(original_campos) >= 4 and "," in original_campos[3]
            filtro_alelo = valor_info(info, AS_FILTER)
            passa_qc = filtros == "PASS" and filtro_alelo in {"", ".", "PASS", "SITE"}
            alelo = Alelo(montagem, contig_canonico(chrom), int(pos), ref, alt)
            if alelo.chave in chamadas:
                raise ValueError(f"Alelo duplicado após normalização: {alelo.chave}. Confira o VCF de entrada.")
            linha_saida = {"variant_key": alelo.chave, "assembly": montagem, "chrom": alelo.cromossomo,
                           "pos": alelo.posicao, "ref": ref, "alt": alt, "filter": filtros,
                           "allele_filter": filtro_alelo, "passa_qc": passa_qc,
                           "qual": None if qual == "." else qual, "original_record": original}
            for origem in origens:
                amostra = campos[indices[origem]] if origem in indices else ""
                medidas = metricas(formato, amostra, numeros, multialelico)
                linha_saida.update({f"{origem}_{campo}": valor for campo, valor in medidas.items()})
            chamadas[alelo.chave] = (alelo, linha_saida)
            contagem["pass" if passa_qc else "outros_filtros"] += 1
    return chamadas, dict(contagem), fontes
