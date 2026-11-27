"""Prepara a referência e executa Mutect2 no par JAK2 da aula."""

import csv
import gzip
from hashlib import file_digest
import math
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import NamedTemporaryFile
from urllib.request import Request, urlopen
from zipfile import ZipFile


GATK_VERSAO = "4.6.2.0"
GATK_URL = f"https://github.com/broadinstitute/gatk/releases/download/{GATK_VERSAO}/gatk-{GATK_VERSAO}.zip"
FASTA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz"
REGIAO_JAK2 = "9:5021937-5126899"


def checksum(caminho):
    with Path(caminho).open("rb") as f:
        return file_digest(f, "sha256").hexdigest()


def baixar(url, destino, timeout=30):
    """Guardar no cache só após o download completo."""
    destino = Path(destino)
    if destino.is_file() and destino.stat().st_size:
        return destino
    destino.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=destino.parent, delete=False) as parcial:
        temporario = Path(parcial.name)
        try:
            with urlopen(Request(url, headers={"User-Agent": "aula-somatica/1"}), timeout=timeout) as resposta:
                shutil.copyfileobj(resposta, parcial)
            parcial.close()
            if not temporario.stat().st_size:
                raise ValueError(f"Download vazio: {url}")
            temporario.replace(destino)
        finally:
            temporario.unlink(missing_ok=True)
    return destino


def java17():
    candidatos = sorted(Path("/usr/lib/jvm").glob("java-17-openjdk-*/bin/java"))
    java_no_path = shutil.which("java")
    if java_no_path:
        candidatos.append(Path(java_no_path))
    for candidato in candidatos:
        versao = subprocess.run([str(candidato), "-version"], check=True, capture_output=True, text=True)
        if 'version "17.' in versao.stderr + versao.stdout:
            return str(candidato)
    return None


def preparar_ferramentas(chamar_mutect2):
    pacotes = [programa for programa in ("samtools", "bcftools") if not shutil.which(programa)]
    if chamar_mutect2 and not java17():
        pacotes.append("openjdk-17-jre-headless")
    if pacotes:
        prefixo = [] if os.geteuid() == 0 else ["sudo"]
        subprocess.run([*prefixo, "apt-get", "update", "-qq"], check=True)
        subprocess.run([*prefixo, "apt-get", "install", "-y", "-qq", *pacotes], check=True)
    return {programa: subprocess.check_output([programa, "--version"], text=True).splitlines()[0]
            for programa in ("samtools", "bcftools")}


def fasta_aula(cache):
    cache = Path(cache)
    fasta = cache / "hg19-9.fa"
    if fasta.is_file() and fasta.stat().st_size:
        return fasta
    compactado = baixar(FASTA_URL, cache / "hg19-chr9.fa.gz")
    with gzip.open(compactado, "rb") as origem, NamedTemporaryFile(dir=cache, delete=False) as parcial:
        temporario = Path(parcial.name)
        try:
            if origem.readline().strip() != b">chr9":
                raise ValueError("O FASTA da UCSC não contém o cabeçalho esperado >chr9.")
            parcial.write(b">9\n")
            shutil.copyfileobj(origem, parcial)
            parcial.close()
            temporario.replace(fasta)
        finally:
            temporario.unlink(missing_ok=True)
    return fasta


def gatk_jar(cache):
    cache = Path(cache)
    nome = f"gatk-package-{GATK_VERSAO}-local.jar"
    jar = cache / nome
    if jar.is_file() and jar.stat().st_size:
        return jar
    pacote = baixar(GATK_URL, cache / f"gatk-{GATK_VERSAO}.zip", timeout=60)
    with ZipFile(pacote) as zip_gatk, zip_gatk.open(f"gatk-{GATK_VERSAO}/{nome}") as origem:
        with NamedTemporaryFile(dir=cache, delete=False) as parcial:
            temporario = Path(parcial.name)
            try:
                shutil.copyfileobj(origem, parcial)
                parcial.close()
                temporario.replace(jar)
            finally:
                temporario.unlink(missing_ok=True)
    return jar


def amostra_bam(bam):
    subprocess.run(["samtools", "quickcheck", str(bam)], check=True)
    cabecalho = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
    amostras = {campo[3:] for linha in cabecalho.splitlines() if linha.startswith("@RG\t")
                for campo in linha.split("\t")[1:] if campo.startswith("SM:") and campo[3:]}
    if len(amostras) != 1:
        raise ValueError(f"Esperava um ID SM no BAM {bam.name}. Encontrei {sorted(amostras)}.")
    return amostras.pop()


def sitios_cobertos(caminho):
    with caminho.open() as f:
        linhas = (linha for linha in f if not linha.startswith("#"))
        dados = csv.DictReader(linhas, delimiter="\t")
        return sum(sum(int(linha[campo]) for campo in ("ref_count", "alt_count", "other_alt_count")) > 10
                   for linha in dados)


def estimativa_valida(caminho, amostra):
    with caminho.open() as f:
        for linha in f:
            campos = linha.strip().split("\t")
            if len(campos) >= 3 and campos[0] == amostra:
                fracao, erro = map(float, campos[1:3])
                # A estimativa 0 ± 1 é o retorno do GATK quando o modelo não estima a fração.
                if math.isfinite(fracao) and 0 <= fracao <= 1 and math.isfinite(erro) and 0 <= erro < 1:
                    return {"fracao": fracao, "erro_padrao": erro}
    return None


def mutect2_aula(repo, cache, saida, referencia=None, jar=None):
    repo, saida = Path(repo), Path(saida)
    tumor_bam, normal_bam = repo / "tumor_JAK2.bam", repo / "normal_JAK2.bam"
    gnomad = repo / "af-only-gnomad-chr9.vcf.gz"
    tumor, normal = amostra_bam(tumor_bam), amostra_bam(normal_bam)
    if (tumor, normal) != ("WP043", "WP044"):
        raise ValueError(f"O par esperado é WP043/WP044. Encontrei {tumor}/{normal}.")
    referencia = Path(referencia) if referencia else fasta_aula(cache)
    jar = Path(jar) if jar else gatk_jar(cache)
    java = java17()
    if not java:
        raise RuntimeError("Java 17 é necessário para executar o GATK desta aula.")
    prefixo = [java, "-Xmx4g", "-jar", str(jar)]
    subprocess.run(["samtools", "faidx", str(referencia)], check=True)
    with Path(f"{referencia}.fai").open() as f:
        if not any(linha.split("\t")[:2] == ["9", "141213431"] for linha in f):
            raise ValueError("A referência da chamada precisa ter o contig 9 de GRCh37 com 141213431 bases.")
    dicionario = referencia.with_suffix(".dict")
    if not dicionario.exists() or dicionario.stat().st_mtime < referencia.stat().st_mtime:
        subprocess.run([*prefixo, "CreateSequenceDictionary", "-R", str(referencia), "-O", str(dicionario)], check=True)
    com_log = saida / "gatk.log"
    with com_log.open("w") as log:
        def gatk(ferramenta, *argumentos):
            try:
                subprocess.run([*prefixo, ferramenta, *map(str, argumentos)], check=True,
                               stdout=log, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as erro:
                trecho = com_log.read_text(errors="replace").splitlines()[-20:]
                raise RuntimeError(f"{ferramenta} falhou. Log em {com_log}.\n" + "\n".join(trecho)) from erro

        bruto = saida / "somatic.vcf.gz"
        gatk("Mutect2", "-R", referencia, "-I", tumor_bam, "-I", normal_bam, "-normal", normal,
             "--germline-resource", gnomad, "--native-pair-hmm-threads", 2, "-L", REGIAO_JAK2, "-O", bruto)
        tabela_tumor, tabela_normal = saida / "tumor.pileups.tsv", saida / "normal.pileups.tsv"
        for bam, tabela in ((tumor_bam, tabela_tumor), (normal_bam, tabela_normal)):
            gatk("GetPileupSummaries", "-I", bam, "-V", gnomad, "-L", REGIAO_JAK2, "-O", tabela)
        contagem = {"tumor": sitios_cobertos(tabela_tumor), "normal": sitios_cobertos(tabela_normal)}
        contaminacao = {"sitios_com_mais_de_10_leituras": contagem, "uso": "sem_estimativa_utilizavel_no_recorte"}
        parametros_filtro = []
        if min(contagem.values()) > 0:
            tabela = saida / "contamination.tsv"
            gatk("CalculateContamination", "-I", tabela_tumor, "-matched", tabela_normal, "-O", tabela)
            estimativa = estimativa_valida(tabela, tumor)
            if estimativa is not None:
                parametros_filtro = ["--contamination-table", tabela]
                contaminacao.update(estimativa, uso="estimativa_exploratoria_aplicada")
        vcf = saida / "filtered.vcf.gz"
        gatk("FilterMutectCalls", "-R", referencia, "-V", bruto, *parametros_filtro, "-O", vcf)
    versao = subprocess.check_output([*prefixo, "--version"], text=True, stderr=subprocess.STDOUT)
    metadados = {"gatk": versao.strip(), "gatk_jar_sha256": checksum(jar), "intervalo": REGIAO_JAK2,
                 "contaminacao": contaminacao, "germline_resource_sha256": checksum(gnomad),
                 "bam_sha256": {"tumor": checksum(tumor_bam), "normal": checksum(normal_bam)}}
    return vcf, referencia, tumor, normal, metadados
