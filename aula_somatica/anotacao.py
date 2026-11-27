"""Anotação VEP e cruzamento de alelos com as tabelas do CIViC."""

from collections import defaultdict
import csv
from datetime import date
import json
from pathlib import Path
import re
from time import sleep
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urlencode
from urllib.request import Request, urlopen

from .chamada import baixar, checksum
from .variantes import Alelo, SEQUENCIA, contig_canonico, montagem_canonica


VEP_SERVIDORES = {"GRCh37": "https://grch37.rest.ensembl.org", "GRCh38": "https://rest.ensembl.org"}
VEP_CAMINHO = "/vep/homo_sapiens/region"
MESES = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()
CIVIC_TABELAS = {
    "variantes": ("VariantSummaries.tsv", {"variant_id", "chromosome", "start", "stop",
                    "reference_build", "reference_bases", "variant_bases", "single_variant_molecular_profile_id"}),
    "perfis": ("MolecularProfileSummaries.tsv", {"molecular_profile_id", "name", "variant_ids"}),
    "evidencias": ("ClinicalEvidenceSummaries.tsv", {"molecular_profile_id", "evidence_id", "disease",
                   "doid", "evidence_type", "evidence_status", "evidence_level", "significance"}),
    "assercoes": ("AssertionSummaries.tsv", {"molecular_profile_id", "assertion_id", "disease",
                   "doid", "assertion_type", "amp_category", "significance"}),
}
VEP_COLUNAS = ("symbol", "gene_id", "transcript", "consequence", "hgvsc", "hgvsp",
               "amino_acids", "canonical", "mane_select")
CIVIC_COLUNAS = (
    "variant_key", "civic_variant_id", "civic_variant_url", "match", "molecular_profile_id",
    "molecular_profile", "aplicabilidade_perfil", "registro", "registro_id", "registro_url",
    "doenca", "doid", "finalidade", "direcao", "significancia", "terapias", "status_civic",
    "flag_civic", "nivel_evidencia_civic", "amp_category_civic", "evidence_ids",
    "fonte_tipo", "fonte_id", "fonte_citacao", "fonte_url", "comparacao_pv_simulada",
)
JAK2_V617F = {"GRCh37:9:5073770:G:T", "GRCh38:9:5073770:G:T"}


def json_http(url, payload=None):
    dados = None if payload is None else json.dumps(payload).encode()
    request = Request(url, data=dados, headers={"Accept": "application/json", "Content-Type": "application/json"})
    for tentativa in range(3):
        try:
            with urlopen(request, timeout=20) as resposta:
                return json.load(resposta)
        except HTTPError as erro:
            if erro.code not in {429, 502, 503, 504} or tentativa == 2:
                raise
            espera = erro.headers.get("Retry-After", "")
            sleep(min(10, float(espera)) if espera.isdecimal() else tentativa + 1)


def anotar_vep(alelos, montagem, fonte, recorte, arquivo=None, fallback=False):
    chaves = {alelo.chave: alelo for alelo in alelos}
    pacote = {"assembly": montagem, "fonte_solicitada": fonte, "registros": {}, "avisos": []}
    if fonte == "arquivo":
        if not arquivo:
            raise ValueError("Informe VEP_ARQUIVO, um vep_response.json exportado pelo notebook.")
        anterior = json.loads(Path(arquivo).read_text(encoding="utf-8"))
        if anterior.get("assembly") != montagem or not isinstance(anterior.get("registros"), dict):
            raise ValueError("O arquivo VEP precisa ter a mesma montagem e o formato vep_response.json.")
        pacote["registros"] = {chave: anterior["registros"][chave] for chave in chaves if chave in anterior["registros"]}
        pacote["metadados_arquivo"] = {"sha256": checksum(arquivo), "ensembl": anterior.get("ensembl")}
    elif fonte == "rest":
        servidor = VEP_SERVIDORES[montagem]
        parametros = {"merged": 1, "hgvs": 1, "canonical": 1, "transcript_version": 1}
        if montagem == "GRCh38":
            parametros["mane"] = 1
        endpoint = f"{servidor}{VEP_CAMINHO}?{urlencode(parametros)}"
        pacote["ensembl"] = {"endpoint": endpoint}
        pedidos = list(chaves.values())
        try:
            for inicio in range(0, len(pedidos), 200):
                lote = pedidos[inicio:inicio + 200]
                enviados = {f"aula{inicio + i}": alelo for i, alelo in enumerate(lote)}
                payload = {"variants": [f"{a.cromossomo} {a.posicao} {id_} {a.ref} {a.alt} . . ."
                                        for id_, a in enviados.items()]}
                resposta = json_http(endpoint, payload)
                if not isinstance(resposta, list):
                    raise ValueError("Resposta inesperada do VEP para uma consulta em lote.")
                for resultado in resposta:
                    campos = resultado.get("input", "").split()
                    if len(campos) < 5 or campos[2] not in enviados:
                        raise ValueError("A resposta do VEP não pode ser ligada ao alelo solicitado.")
                    alelo = enviados[campos[2]]
                    if (resultado.get("assembly_name") != montagem or contig_canonico(campos[0]) != alelo.cromossomo
                            or campos[1] != str(alelo.posicao) or campos[3:5] != [alelo.ref, alelo.alt]):
                        raise ValueError(f"Montagem ou alelo divergente na resposta do VEP: {alelo.chave}.")
                    entrada = pacote["registros"].setdefault(alelo.chave, {"fonte": "vep_rest", "respostas": []})
                    entrada["respostas"].append(resultado)
            if pedidos:
                try:
                    pacote["ensembl"]["dados"] = json_http(f"{servidor}/info/data")
                except (URLError, TimeoutError, ValueError) as erro:
                    pacote["avisos"].append(f"O VEP respondeu, mas a versão dos dados ficou indisponível ({type(erro).__name__}).")
        except (URLError, TimeoutError) as erro:
            if not fallback:
                raise RuntimeError("VEP REST indisponível. Use VEP_FONTE='arquivo' com um resultado salvo.") from erro
            pacote["avisos"].append(f"VEP REST indisponível ({type(erro).__name__}). Usando o recorte onde houver anotação.")
    elif fonte != "recorte":
        raise ValueError("VEP_FONTE deve ser rest, arquivo ou recorte.")
    if fonte == "recorte" or (fonte == "rest" and fallback):
        for chave in chaves:
            if chave in recorte:
                pacote["registros"].setdefault(chave, recorte[chave])
    return pacote


def linhas_vep(registro):
    if not registro:
        return []
    if "tabela" in registro:
        return registro["tabela"]
    linhas = []
    for resposta in registro.get("respostas", []):
        transcritos = resposta.get("transcript_consequences", [])
        if not transcritos:
            linhas.append({"consequence": resposta.get("most_severe_consequence", "")})
        for transcrito in transcritos:
            linhas.append({
                "symbol": transcrito.get("gene_symbol", ""), "gene_id": transcrito.get("gene_id", ""),
                "transcript": transcrito.get("transcript_id", ""),
                "consequence": "&".join(transcrito.get("consequence_terms", [])),
                "hgvsc": unquote(transcrito.get("hgvsc", "")), "hgvsp": unquote(transcrito.get("hgvsp", "")),
                "amino_acids": transcrito.get("amino_acids", ""), "canonical": transcrito.get("canonical", ""),
                "mane_select": transcrito.get("mane_select", ""),
            })
    return linhas


def ler_tsv(caminho, obrigatorios):
    with Path(caminho).open(encoding="utf-8-sig", newline="") as f:
        leitor = csv.DictReader(f, delimiter="\t")
        ausentes = obrigatorios - set(leitor.fieldnames or ())
        if ausentes:
            raise ValueError(f"Colunas ausentes em {Path(caminho).name}: {', '.join(sorted(ausentes))}.")
        return list(leitor)


def civic_mensal(release, cache, pasta=None):
    dia = date.fromisoformat(release)
    if dia.day != 1:
        raise ValueError("CIVIC_RELEASE deve ser o primeiro dia do mês em AAAA-MM-DD.")
    etiqueta = f"01-{MESES[dia.month - 1]}-{dia.year}"
    fontes, tabelas = {}, {}
    for categoria, (nome, obrigatorios) in CIVIC_TABELAS.items():
        online = f"https://civicdb.org/downloads/{etiqueta}/{etiqueta}-{nome}"
        espelho = f"https://civic-aws-opendata.s3.us-west-2.amazonaws.com/{nome[:-4]}/date={etiqueta}/{nome}"
        if pasta:
            alternativos = [f"{etiqueta}-{nome}", nome]
            if categoria in {"evidencias", "assercoes"}:
                alternativos.extend([f"{etiqueta}-Accepted{nome}", f"Accepted{nome}"])
            caminho = next((Path(pasta) / f for f in alternativos if (Path(pasta) / f).is_file()), None)
            if caminho is None:
                raise FileNotFoundError(f"Arquivo ausente em CIVIC_PASTA: {nome}.")
            fonte = "arquivo_local"
        else:
            pasta_cache = Path(cache) / "civic" / etiqueta
            destino, destino_aws = pasta_cache / "civicdb" / nome, pasta_cache / "aws" / nome
            fonte = online
            if not destino.exists() and destino_aws.exists():
                fonte, caminho = espelho, destino_aws
            else:
                try:
                    caminho = baixar(online, destino, timeout=12)
                except (URLError, TimeoutError):
                    fonte = espelho
                    caminho = baixar(espelho, destino_aws, timeout=12)
        tabelas[categoria] = ler_tsv(caminho, obrigatorios)
        fontes[categoria] = {"origem": fonte, "arquivo": caminho.name, "sha256": checksum(caminho)}
    return tabelas, {"escopo": "release_mensal", "release": release, "arquivos": fontes}


def obter_civic(fonte, release, cache, recorte, pasta=None, fallback=False):
    if fonte not in {"mensal", "recorte"}:
        raise ValueError("CIVIC_FONTE deve ser mensal ou recorte.")
    metadados = {"escopo": "recorte_didatico", "origem": recorte["origem"]}
    if fonte == "mensal":
        try:
            return civic_mensal(release, cache, pasta)
        except (URLError, TimeoutError) as erro:
            if not fallback or pasta:
                raise RuntimeError("Release do CIViC indisponível. Informe CIVIC_PASTA ou selecione o recorte didático.") from erro
            metadados.update(release_solicitada=release, aviso=f"Release indisponível ({type(erro).__name__}).")
    return recorte["tabelas"], metadados


def chave_coordenadas_civic(variante):
    try:
        montagem = montagem_canonica(variante["reference_build"])
        inicio, fim = int(variante["start"]), int(variante["stop"])
    except (TypeError, ValueError):
        return None
    cromossomo = contig_canonico(variante["chromosome"])
    ref, alt = (variante[campo].strip().upper() for campo in ("reference_bases", "variant_bases"))
    ref, alt = ("" if base == "-" else base for base in (ref, alt))
    if (not cromossomo or inicio < 1 or ref == alt or
            any(base and not SEQUENCIA.fullmatch(base) for base in (ref, alt))):
        return None
    if ref and fim == inicio + len(ref) - 1:
        return Alelo(montagem, cromossomo, inicio, ref, alt).sem_ancora
    if not ref and fim == inicio + 1:
        return montagem, cromossomo, fim, ref, alt
    return None


def hgvs_codificante(texto):
    return texto if re.fullmatch(r"(?:NM_\d+|ENST\d+)\.\d+:c\.\S+", texto) else ""


def sinalizado(*registros):
    return any(r.get("is_flagged", "").lower() in {"true", "1", "t"} for r in registros)


def indice_civic(tabelas):
    coordenadas, hgvs, perfis, registros = (defaultdict(list) for _ in range(4))
    for variante in tabelas["variantes"]:
        chave = chave_coordenadas_civic(variante)
        if chave:
            coordenadas[chave].append(variante)
        for descricao in variante.get("hgvs_descriptions", "").split(","):
            if hgvsc := hgvs_codificante(descricao.strip()):
                hgvs[hgvsc].append(variante)
    for perfil in tabelas["perfis"]:
        for variante in perfil["variant_ids"].replace(" ", "").split(","):
            perfis[variante].append(perfil)
    for tipo, categoria in (("EID", "evidencias"), ("AID", "assercoes")):
        for registro in tabelas[categoria]:
            if registro.get("evidence_status", registro.get("assertion_status", "accepted")).lower() == "accepted":
                registros[registro["molecular_profile_id"]].append((tipo, registro))
    return coordenadas, hgvs, perfis, registros


def evidencias_civic(alelo, transcritos, indice, escopo):
    coordenadas, hgvs, perfis, registros = indice
    correspondencias = coordenadas.get(alelo.sem_ancora)
    criterio = "coordenada_exata"
    if not correspondencias:
        candidatas = (v for t in transcritos for v in hgvs.get(hgvs_codificante(t.get("hgvsc", "")), []))
        correspondencias = list({v["variant_id"]: v for v in candidatas}.values())
        criterio = "hgvs_requer_revisao"
    if not correspondencias:
        return ("sem_correspondencia_na_release" if escopo == "release_mensal" else "sem_correspondencia_no_recorte"), []
    if len(correspondencias) > 1:
        criterio = "correspondencia_ambigua"
    linhas = []
    for variante in correspondencias:
        variante_id = variante["variant_id"]
        for perfil in perfis[variante_id]:
            perfil_id = perfil["molecular_profile_id"]
            perfil_unitario = (perfil_id == variante["single_variant_molecular_profile_id"] and
                               perfil["variant_ids"].strip() == variante_id)
            for tipo, registro in registros[perfil_id]:
                id_prefixo = "evidence" if tipo == "EID" else "assertion"
                finalidade = registro[f"{id_prefixo}_type"]
                doid = registro["doid"].strip().removeprefix("DOID:")
                flag = sinalizado(variante, perfil, registro)
                escopo_perfil = "variante_isolada" if perfil_unitario else "perfil_requer_revisao"
                if flag:
                    escopo_perfil = "curadoria_sinalizada"
                comparacao = ""
                if alelo.chave in JAK2_V617F:
                    if criterio != "coordenada_exata" or not perfil_unitario or flag or not doid:
                        comparacao = "requer_revisao"
                    elif doid != "8997":
                        comparacao = "outra_doenca"
                    else:
                        comparacao = "mesma_doenca_e_finalidade" if finalidade.lower() == "diagnostic" else "outra_finalidade"
                fonte_tipo, fonte_id = registro.get("source_type", ""), registro.get("citation_id", "")
                linhas.append({
                    "variant_key": alelo.chave, "civic_variant_id": variante_id,
                    "civic_variant_url": variante.get("variant_civic_url", ""), "match": criterio,
                    "molecular_profile_id": perfil_id, "molecular_profile": perfil["name"],
                    "aplicabilidade_perfil": escopo_perfil, "registro": tipo, "registro_id": registro[f"{id_prefixo}_id"],
                    "registro_url": registro.get(f"{id_prefixo}_civic_url", ""), "doenca": registro["disease"],
                    "doid": f"DOID:{doid}" if doid else "", "finalidade": finalidade,
                    "direcao": registro.get(f"{id_prefixo}_direction", ""), "significancia": registro["significance"],
                    "terapias": registro.get("therapies", ""), "status_civic": registro.get("evidence_status", "accepted"),
                    "flag_civic": flag, "nivel_evidencia_civic": registro.get("evidence_level", ""),
                    "amp_category_civic": registro.get("amp_category", ""), "evidence_ids": registro.get("evidence_item_ids", ""),
                    "fonte_tipo": fonte_tipo, "fonte_id": fonte_id, "fonte_citacao": registro.get("citation", ""),
                    "fonte_url": f"https://pubmed.ncbi.nlm.nih.gov/{fonte_id}/" if fonte_tipo.lower() == "pubmed" else "",
                    "comparacao_pv_simulada": comparacao,
                })
    return criterio if linhas else f"{criterio}_sem_registros_aceitos", linhas


def reunir_anotacoes(chamadas, pacote_vep, tabelas_civic, escopo_civic):
    indice = indice_civic(tabelas_civic)
    variantes, evidencias = [], []
    for chave, (alelo, chamada) in chamadas.items():
        registro = pacote_vep["registros"].get(chave)
        transcritos = linhas_vep(registro)
        status_vep = registro["fonte"] if registro else "sem_anotacao_disponivel"
        if chamada["passa_qc"]:
            status_civic, itens = evidencias_civic(alelo, transcritos, indice, escopo_civic)
            evidencias.extend(itens)
        else:
            status_vep, status_civic = "filtrada", "filtrada"
        for transcrito in transcritos or [{}]:
            variantes.append({**chamada, "vep_status": status_vep, "civic_status": status_civic,
                              **{coluna: transcrito.get(coluna, "") for coluna in VEP_COLUNAS}})
    return variantes, evidencias
