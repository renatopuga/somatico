"""Regressões para os pontos em que a ligação por alelo ou contexto pode falhar."""

import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from urllib.error import URLError

from aula_somatica import anotacao, chamada, revisao, variantes


PASTA = Path(__file__).parent
RECORTE = json.loads((PASTA / "recorte_aula.json").read_text(encoding="utf-8"))
CHAVE_JAK2 = "GRCh37:9:5073770:G:T"


class TestAula(unittest.TestCase):
    def setUp(self):
        self.chamadas, self.contagens, _ = variantes.ler_chamadas(
            PASTA / "exemplo_sintetico.vcf", "GRCh37", "TUMOR_DEMO", "NORMAL_DEMO",
        )
        self.jak2 = self.chamadas[CHAVE_JAK2][0]

    def test_exemplo_recorte_nao_transformado_em_tier(self):
        pacote = anotacao.anotar_vep([self.jak2], "GRCh37", "recorte", RECORTE["vep"])
        tabela, evidencias = anotacao.reunir_anotacoes(
            self.chamadas, pacote, RECORTE["civic"]["tabelas"], "recorte_didatico",
        )
        por_chave = {linha["variant_key"]: linha for linha in tabela}
        self.assertEqual((self.contagens["pass"], self.contagens["outros_filtros"]), (2, 1))
        self.assertEqual(por_chave[CHAVE_JAK2]["tumor_vaf"], 0.28)
        self.assertEqual(por_chave[CHAVE_JAK2]["normal_vaf"], 0)
        self.assertEqual(por_chave[CHAVE_JAK2]["transcript"], "NM_001322194.2")
        self.assertEqual(por_chave["GRCh37:9:5073681:C:CT"]["civic_status"], "sem_correspondencia_no_recorte")
        self.assertEqual(por_chave["GRCh37:9:5073770:G:A"]["civic_status"], "filtrada")
        self.assertEqual([(e["registro_id"], e["comparacao_pv_simulada"]) for e in evidencias], [("1", "outra_doenca")])
        fichas = revisao.ficha_revisao(self.chamadas)
        self.assertEqual(len(fichas), 3)
        self.assertTrue(all(not f.get("tier_proposto") and f["decisao"] == "pendente" for f in fichas))

    def test_ficha_controle_nivel_e_revisor(self):
        fichas = revisao.ficha_revisao(self.chamadas)
        formulario = {"variant_key": CHAVE_JAK2, "cenario": "pv_simulada", "tier_proposto": "I",
                      "nivel_proposto": "C", "fonte_proposta": revisao.ICC_PV,
                      "justificativa_proposta": "Critério molecular no ICC."}
        with self.assertRaisesRegex(ValueError, "Combinação inválida"):
            revisao.com_revisoes(fichas, [formulario])
        formulario.update(nivel_proposto="A", tier_revisado="I", nivel_revisado="A")
        with self.assertRaisesRegex(ValueError, "revisor"):
            revisao.com_revisoes(fichas, [formulario])
        formulario.update(revisor="Revisor A", fonte_revisada=revisao.ICC_PV, justificativa_revisao="Fonte conferida.")
        revista = revisao.com_revisoes(fichas, [formulario], hoje="2026-11-27")
        linha = next(r for r in revista if r["cenario"] == "pv_simulada")
        self.assertEqual((linha["decisao"], linha["data_revisao"], linha["tier_revisado"]), ("revista", "2026-11-27", "I"))

    def test_coordenadas_civic_indel_e_build(self):
        base = {"reference_build": "GRCh37", "chromosome": "chr9", "start": "5", "stop": "6"}
        insercao = {**base, "reference_bases": "", "variant_bases": "CT"}
        delecao = {**base, "reference_bases": "CT", "variant_bases": "-"}
        self.assertEqual(anotacao.chave_coordenadas_civic(insercao), variantes.Alelo("GRCh37", "9", 5, "A", "ACT").sem_ancora)
        self.assertEqual(anotacao.chave_coordenadas_civic(delecao), variantes.Alelo("GRCh37", "9", 4, "ACT", "A").sem_ancora)
        self.assertNotEqual(anotacao.chave_coordenadas_civic({**insercao, "reference_build": "GRCh38"}),
                            variantes.Alelo("GRCh37", "9", 5, "A", "ACT").sem_ancora)

    def test_civic_release_tsv_filtra_submetidos_e_marca_perfil_composto(self):
        tabelas = {k: [dict(v) for v in valores] for k, valores in RECORTE["civic"]["tabelas"].items()}
        item = tabelas["evidencias"][0]
        tabelas["evidencias"].append({**item, "evidence_id": "2", "evidence_status": "submitted"})
        tabelas["perfis"].append({"molecular_profile_id": "teste-composto", "name": "Teste", "variant_ids": "64, 9"})
        tabelas["assercoes"].append({"molecular_profile_id": "teste-composto", "assertion_id": "3", "disease": "PV",
                                    "doid": "8997", "assertion_type": "Diagnostic", "amp_category": "I", "significance": "Positive"})
        with TemporaryDirectory() as pasta:
            for chave, (nome, colunas) in anotacao.CIVIC_TABELAS.items():
                campos = colunas | {c for linha in tabelas[chave] for c in linha}
                revisao.escrever_tsv(Path(pasta) / nome, sorted(campos), tabelas[chave])
            dados, meta = anotacao.civic_mensal("2026-11-01", pasta, pasta=pasta)
        criterio, resultados = anotacao.evidencias_civic(
            self.jak2, [], anotacao.indice_civic(dados), meta["escopo"],
        )
        self.assertEqual(criterio, "coordenada_exata")
        self.assertEqual({r["registro_id"] for r in resultados}, {"1", "3"})
        registro = next(r for r in resultados if r["registro"] == "AID")
        self.assertEqual(registro["comparacao_pv_simulada"], "requer_revisao")
        self.assertEqual(registro["aplicabilidade_perfil"], "perfil_requer_revisao")
        self.assertEqual(len(meta["arquivos"]), 4)

    def test_hgvs_com_build_distinto_e_apenas_candidato(self):
        alelo = variantes.Alelo("GRCh38", "9", 5073770, "G", "T")
        transcritos = [{"hgvsc": "NM_004972.4:c.1849G>T"}]
        criterio, resultados = anotacao.evidencias_civic(
            alelo, transcritos, anotacao.indice_civic(RECORTE["civic"]["tabelas"]), "recorte_didatico",
        )
        self.assertEqual(criterio, "hgvs_requer_revisao")
        self.assertEqual(resultados[0]["comparacao_pv_simulada"], "requer_revisao")

    def test_vep_restaura_ligacao_por_alelo_e_montagem(self):
        requests = []

        def resposta(url, payload=None):
            requests.append(url)
            if payload is None:
                return {"releases": [104]}
            self.assertEqual(payload, {"variants": ["9 5073770 aula0 G T . . ."]})
            return [{"input": payload["variants"][0], "assembly_name": "GRCh37", "transcript_consequences": [
                {"gene_symbol": "JAK2", "hgvsc": "NM_004972.4%3Ac.1849G%3ET"},
            ]}]

        with patch.object(anotacao, "json_http", side_effect=resposta):
            pacote = anotacao.anotar_vep([self.jak2], "GRCh37", "rest", RECORTE["vep"])
        self.assertIn("grch37.rest.ensembl.org", requests[0])
        self.assertEqual(anotacao.linhas_vep(pacote["registros"][CHAVE_JAK2])[0]["hgvsc"], "NM_004972.4:c.1849G>T")
        with patch.object(anotacao, "json_http", side_effect=URLError("offline")):
            fallback = anotacao.anotar_vep([self.jak2], "GRCh37", "rest", RECORTE["vep"], fallback=True)
        self.assertEqual(fallback["registros"][CHAVE_JAK2]["fonte"], "vep_104_3_historico")
        self.assertTrue(fallback["avisos"])

    def test_contaminacao_indisponivel_e_distinta_de_zero_estimado(self):
        with TemporaryDirectory() as pasta:
            tabela = Path(pasta) / "contamination.tsv"
            tabela.write_text("sample\tcontamination\terror\nWP043\t0.0\t1.0\n")
            self.assertIsNone(chamada.estimativa_valida(tabela, "WP043"))
            tabela.write_text("sample\tcontamination\terror\nWP043\t0.0\t0.02\n")
            self.assertEqual(chamada.estimativa_valida(tabela, "WP043"), {"fracao": 0, "erro_padrao": 0.02})

    @unittest.skipUnless(shutil.which("samtools") and shutil.which("bcftools"), "requer samtools e bcftools")
    def test_norm_vcf_multialelico_respeita_filtro_alelo_e_metricas(self):
        with TemporaryDirectory() as pasta:
            raiz = Path(pasta)
            fasta, vcf, saida = (raiz / n for n in ("mini.fa", "entrada.vcf", "norm.vcf"))
            fasta.write_text(">chr9\nCCGTT\n")
            vcf.write_text("\n".join([
                "##fileformat=VCFv4.2", "##reference=GRCh38", "##contig=<ID=chr9,length=5>",
                '##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description="Filtros por alelo">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
                '##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">',
                '##FORMAT=<ID=MMQ,Number=R,Type=Integer,Description="Campo fora do exercício">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT",
                "chr9\t3\t.\tG\tT,A\t.\tPASS\tAS_FilterStatus=SITE|weak_evidence,strand_bias\tDP:AD:AF:MMQ\t100:60,25,10:0.25,.:60,60,60,60\n",
            ]))
            variantes.normalizar(vcf, fasta, "GRCh38", "T", None, saida)
            chamadas, _, _ = variantes.ler_chamadas(saida, "GRCh38", "T")
            aprovada = chamadas["GRCh38:9:3:G:T"][1]
            filtrada = chamadas["GRCh38:9:3:G:A"][1]
            self.assertTrue(aprovada["passa_qc"])
            self.assertFalse(filtrada["passa_qc"])
            self.assertEqual((aprovada["tumor_ad_alt"], aprovada["tumor_vaf"]), (25, 0.25))
            self.assertEqual((filtrada["tumor_ad_alt"], filtrada["tumor_vaf"]), (10, 0.10))
            self.assertEqual(filtrada["allele_filter"], "weak_evidence&strand_bias")
            self.assertIsNone(aprovada["normal_vaf"])
            self.assertEqual(len(revisao.ficha_revisao(chamadas)), 1)
            with self.assertRaisesRegex(ValueError, "GRCh38"):
                variantes.normalizar(vcf, fasta, "GRCh37", "T", None, saida)
            fasta.write_text(">chr9\nCCATT\n")
            with self.assertRaisesRegex(ValueError, "Confira REF"):
                variantes.normalizar(vcf, fasta, "GRCh38", "T", None, saida)


if __name__ == "__main__":
    unittest.main()
