"""Exporta tabelas e registra a proposta do aluno e a decisão do revisor."""

import csv
from datetime import date
from pathlib import Path
from tempfile import NamedTemporaryFile

from .anotacao import JAK2_V617F, VEP_COLUNAS


METRICAS = ("dp", "ad_ref", "ad_alt", "vaf", "vaf_fonte")
VARIANTES_COLUNAS = (
    "variant_key", "assembly", "chrom", "pos", "ref", "alt", "filter", "allele_filter", "passa_qc", "qual", "original_record",
    *(f"{amostra}_{campo}" for amostra in ("tumor", "normal") for campo in METRICAS),
    "vep_status", "civic_status", *VEP_COLUNAS,
)
REVISAO_COLUNAS = (
    "variant_key", "cenario", "contexto", "finalidade", "tier_proposto", "nivel_proposto",
    "fonte_proposta", "justificativa_proposta", "tier_revisado", "nivel_revisado", "fonte_revisada",
    "justificativa_revisao", "revisor", "data_revisao", "pendencias", "decisao",
)
TIER_NIVEIS = {"I": {"A", "B"}, "II": {"C", "D"}, "III": {""}, "IV": {""}, "": {""}}
ICC_PV = "https://pmc.ncbi.nlm.nih.gov/articles/PMC9479031/"
AMP_ASCO_CAP = "https://pubmed.ncbi.nlm.nih.gov/27993330/"


def escrever_tsv(caminho, colunas, linhas):
    caminho = Path(caminho)
    with NamedTemporaryFile(dir=caminho.parent, mode="w", encoding="utf-8", newline="", delete=False) as f:
        temporario = Path(f.name)
        try:
            escritor = csv.DictWriter(f, fieldnames=colunas, delimiter="\t", lineterminator="\n")
            escritor.writeheader()
            escritor.writerows(linhas)
            f.close()
            temporario.replace(caminho)
        finally:
            temporario.unlink(missing_ok=True)


def ficha_revisao(chamadas):
    fichas = []
    for chave, (_, chamada) in chamadas.items():
        if not chamada["passa_qc"]:
            continue
        comum = {"variant_key": chave, "decisao": "pendente"}
        fichas.append({**comum, "cenario": "observado", "contexto": "",
                       "pendencias": "Contexto clínico a documentar."})
        if chave in JAK2_V617F:
            fichas.append({**comum, "cenario": "pv_simulada", "contexto": "Policitemia vera (cenário fictício)",
                           "finalidade": "diagnostica", "pendencias": "Conferir critérios do ICC e aplicabilidade das fontes."})
    return fichas


def com_revisoes(fichas, alteracoes, hoje=None):
    resultado = {(ficha["variant_key"], ficha["cenario"]): dict(ficha) for ficha in fichas}
    hoje = hoje or date.today().isoformat()
    campos_editaveis = set(REVISAO_COLUNAS) - {"variant_key", "cenario", "decisao", "data_revisao"}
    chaves_editadas = set()
    for alteracao in alteracoes:
        chave = (alteracao.get("variant_key"), alteracao.get("cenario"))
        if chave not in resultado or chave in chaves_editadas:
            raise ValueError(f"Ficha inexistente ou repetida: {chave}.")
        chaves_editadas.add(chave)
        desconhecidos = alteracao.keys() - campos_editaveis - {"variant_key", "cenario"}
        if desconhecidos:
            raise ValueError(f"Campos não editáveis: {', '.join(sorted(desconhecidos))}.")
        ficha = resultado[chave]
        ficha.update({k: str(v).strip() for k, v in alteracao.items() if k in campos_editaveis and v is not None})
        tier_proposto, nivel_proposto = ficha.get("tier_proposto", ""), ficha.get("nivel_proposto", "")
        tier_revisado, nivel_revisado = ficha.get("tier_revisado", ""), ficha.get("nivel_revisado", "")
        for tier, nivel in ((tier_proposto, nivel_proposto), (tier_revisado, nivel_revisado)):
            if nivel not in TIER_NIVEIS.get(tier, set()):
                raise ValueError(f"Combinação inválida de tier e nível na ficha {chave}.")
        if tier_proposto:
            exigidos = ("contexto", "finalidade", "fonte_proposta", "justificativa_proposta")
            if not all(ficha.get(campo) for campo in exigidos):
                raise ValueError(f"Preencha contexto, finalidade, fonte e justificativa em {chave}.")
            ficha["decisao"] = "proposta_registrada"
        if tier_revisado:
            if not tier_proposto or not all(ficha.get(campo) for campo in ("revisor", "fonte_revisada", "justificativa_revisao")):
                raise ValueError(f"A decisão revisada requer proposta, revisor, fonte e justificativa em {chave}.")
            ficha.update(decisao="revista", data_revisao=hoje)
    return list(resultado.values())
