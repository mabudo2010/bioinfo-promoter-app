import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import io
import gzip
import re
import os

from collections import defaultdict
from Bio import Seq, SeqIO
from scipy.stats import ttest_ind, pearsonr, spearmanr, levene, mannwhitneyu, chi2_contingency, fisher_exact

try:
    from statistics import calcular_importancia_e_associacao
except ImportError:
    pass

try:
    from regulatory_analysis import analisar_impacto_regulatorio, mapear_motivos_nativos_referencia
except ImportError:
    def analisar_impacto_regulatorio(*args, **kwargs): return [], []
    def mapear_motivos_nativos_referencia(*args, **kwargs): return []

st.set_page_config(page_title="Bioinfo Promoter Analyzer Pro", layout="wide")

# ==============================================================================
# SISTEMA DE INTERNACIONALIZAÇÃO (i18n)
# ==============================================================================
if 'lang' not in st.session_state:
    st.session_state['lang'] = 'PT'

def tr(pt, en):
    """Devolve a string em PT ou EN consoante a linguagem selecionada no session_state."""
    return en if st.session_state['lang'] == 'EN' else pt

def atualizar_idioma():
    """Callback para alterar a língua sem abortar a execução (evita perda de ficheiros)."""
    escolha = st.session_state.lang_selector
    st.session_state['lang'] = 'PT' if 'Português' in escolha else 'EN'

# ==============================================================================
# SELETOR DE LÍNGUA (SIDEBAR)
# ==============================================================================
st.sidebar.markdown("### 🌐 Idioma / Language")
st.sidebar.radio(
    "Selecione / Select:", 
    ["🇵🇹 Português", "🇬🇧 English"], 
    index=0 if st.session_state['lang'] == 'PT' else 1,
    horizontal=True,
    label_visibility="collapsed",
    key="lang_selector",
    on_change=atualizar_idioma
)

st.sidebar.markdown("---")


st.title("🧬 Bioinfo Promoter Analyser Pro")
st.markdown(tr(
    "Compare o perfil mutacional, impacto regulatório e correlação fenotípica real das amostras baseado nos dados clínicos consolidados.",
    "Compare the mutational profile, regulatory impact, and real phenotypic correlation of the samples based on consolidated clinical data."
))

MUTACOES_FENOTIPO = {
    -124: tr("Associação com Hiperexpressão de CDR1 e Resistência a Azóis (Fluconazol)", "Association with CDR1 Overexpression and Azole Resistance (Fluconazole)"),
    -342: tr("Ganho de função do promotor induzido por stress xenobiótico", "Promoter gain of function induced by xenobiotic stress"),
    -450: tr("Mutação reportada em isolados clínicos resistentes em ambiente hospitalar", "Mutation reported in resistant clinical isolates in hospital environments")
}

config_plotly_export = {
    'toImageButtonOptions': {
        'format': 'svg',
        'filename': 'promoter_physical_mapping_pro',
        'height': 800,
        'width': 1300,
        'scale': 3
    },
    'displaylogo': False,
    'modeBarButtonsToRemove': ['sendDataToCloud', 'hoverClosestCartesian', 'hoverCompareCartesian']
}

# ==============================================================================
# CONFIGURAÇÃO: FATORES DE TRANSCRIÇÃO PARA MATRIZ GLOBAL
# ==============================================================================
TF_LIST = [
    ("Amt1p", "WTWHNGCTGW"), ("Ap1p", "TTACAAA"), ("Ap1p", "TTAGTAA"), ("Ap1p", "VTTACWAAB"),
    ("Ap1p", "VDTASTAA"), ("Ap1p", "GWWSTASWWCT"), ("Ap5p", "TTACKAA"), ("Ap5p", "ATTACTA"),
    ("Ap5p", "CASACAGWAS"), ("Ap5p", "AYAGTGGTTCAT"), ("Atf2p", "ACACCC"), ("CAGL0B03421g", "CCCCTCC"),
    ("Cst6p", "TGACGTCA"), ("Hap5p", "CCAAT"), ("Msn2p", "CCYYCCYYYW"), ("Msn4p", "CCYYCCYYYW"),
    ("Pdr1p", "BCCRYYRGD"), ("Pdr1p", "TCCRYGGA"), ("Pdr1p", "TCCACGGA"), ("Pdr1p", "HYCCGTGGR"),
    ("Pdr1p", "TCCACGGAA"), ("Pdr1p", "RTTCCGTGGAA"), ("Pho4p", "CCCGTG"), ("Pho4p", "CACGTG"),
    ("Rpn4p", "TTGCAAA"), ("Yap3p", "ATTACGTAAT"), ("Yap6p", "AATKACV"), ("Yap6p", "WATKACGW"),
    ("Yap6p", "ATGACGT"), ("Yap6p", "ATKACGTMAT"), ("Yap6p", "RTKACGTMAY"), ("Yap6p", "GCACAGA"),
    ("Yap6p", "ATTACGTW"), ("Yap6p", "YCGWCGS"), ("Yap7p", "TTACTAA"), ("Yap7p", "WGMTTACTRR"),
    ("Yap7p", "TTAGTAA"), ("Yap7p", "ATTACTAA"), ("Yap7p", "MTKASTMA"), ("Yap7p", "DATTACTAAD"),
    ("Yap7p", "GMTKASTMAK"), ("Yap7p", "TTASTAAK"), ("Yap7p", "AGTCATM"), ("Yap7p", "KTACTAAD"),
]

IUPAC_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

def seq_to_regex(seq):
    return "".join(IUPAC_REGEX.get(c, c) for c in seq.upper())


# ==============================================================================
# FUNÇÕES UTILITÁRIAS
# ==============================================================================
def normalizar_chrom(chrom_nome):
    if pd.isna(chrom_nome): return ""
    c_str = str(chrom_nome).strip()
    c_str = c_str.split('_')[0]
    c_str = re.sub(r'^[Cc]hr', '', c_str)
    return c_str

def normalizar_id_amostra(id_str):
    if pd.isna(id_str): return ""
    s = str(id_str).strip()
    s = re.sub(r'^[Ss]ample[_-]?', '', s)
    s = re.sub(r'^[Ii]solate[_-]?', '', s)
    s = re.sub(r'\.(vcf|fasta|fa|bam|sam)(\.gz)?$', '', s, flags=re.IGNORECASE)
    s = s.replace('-', '_').replace(' ', '_').upper()
    s = re.sub(r'_+', '_', s).strip('_')
    return s

def corresponder_amostra(paciente_vcf, mapa_clinico):
    norm_vcf = normalizar_id_amostra(paciente_vcf)
    if norm_vcf in mapa_clinico:
        return mapa_clinico[norm_vcf]
        
    def limpar_id(s):
        s_limpo = re.sub(r'\W+', '', str(s)).upper()
        for prefix in ['SAMPLE', 'ISOLATE', 'AMOSTRA', 'ISO', 'CEPA', 'STRAIN', 'CG', 'CAN', 'GLABRATA', 'PACIENTE', 'PATIENT', 'ID']:
            if s_limpo.startswith(prefix) and len(s_limpo) > len(prefix):
                s_limpo = s_limpo[len(prefix):]
        return re.sub(r'^0+', '', s_limpo)
        
    id_vcf_limpo = limpar_id(paciente_vcf)
    
    for k_clin, row_c in mapa_clinico.items():
        id_clin_limpo = limpar_id(k_clin)
        if id_vcf_limpo == id_clin_limpo and len(id_vcf_limpo) > 0:
            return row_c
            
    for k_clin, row_c in mapa_clinico.items():
        id_clin_limpo = re.sub(r'\W+', '', str(k_clin)).upper()
        vcf_simples = re.sub(r'\W+', '', str(paciente_vcf)).upper()
        if (len(vcf_simples) > 2 and vcf_simples in id_clin_limpo) or (len(id_clin_limpo) > 2 and id_clin_limpo in vcf_simples):
            return row_c
            
    return None

def extrair_valor_mic(val):
    if pd.isna(val): return 0.0
    val_str = str(val).strip().replace(',', '.')
    match = re.search(r'(\d+(?:\.\d+)?)', val_str)
    if match:
        try: return float(match.group(1))
        except ValueError: return 0.0
    return 0.0

def extrair_e_traduzir_gt(info_campo, ref, alt):
    if pd.isna(info_campo): return "./."
    info_str = str(info_campo).strip()
    if info_str.startswith('.') or not info_str or info_str == "./.": return "./."
    
    gt = info_str.split(':')[0].replace('|', '/')
    sub_gts = gt.split('/')
    
    alelos_alt = [a.strip() for a in str(alt).replace('/', ',').split(',') if a.strip()]
    alelos_todos = [str(ref).strip()] + alelos_alt
    
    alelos_mutados = []
    for a_idx in sub_gts:
        if a_idx != '.' and a_idx.isdigit():
            idx = int(a_idx)
            if idx > 0:
                if idx < len(alelos_todos): alelos_mutados.append(alelos_todos[idx])
                else: alelos_mutados.append(alelos_todos[-1])
                    
    if not alelos_mutados: return "WT"
        
    alelos_mutados_unicos = sorted(list(set(alelos_mutados)))
    return " / ".join([f"{str(ref).strip()}->{a_mut}" for a_mut in alelos_mutados_unicos])

def formatar_p_value(p_val):
    if pd.isna(p_val): return "N/A"
    return "< 0.001" if p_val < 0.001 else f"{p_val:.4f}"

def formatar_p_texto(p_val):
    if pd.isna(p_val): return "N/A"
    return "p < 0.001" if p_val < 0.001 else f"p = {p_val:.4f}"

def agrupar_sitios_fisicos(lista_motivos, margem_bp=10, footprint_alvo=12):
    if not lista_motivos:
        return []
    
    motivos_por_ft = {}
    for m in lista_motivos:
        ft = str(m.get("FT", "")).strip()
        if ft not in motivos_por_ft:
            motivos_por_ft[ft] = []
        motivos_por_ft[ft].append(m)
        
    sitios_agrupados = []
    
    for ft, motivos in motivos_por_ft.items():
        motivos_ordenados = sorted(motivos, key=lambda x: min(int(x.get('Pos_5p', 0)), int(x.get('Pos_3p', 0))))
        
        sitios_unicos = []
        for m in motivos_ordenados:
            p_min_raw = min(int(m.get('Pos_5p', 0)), int(m.get('Pos_3p', 0)))
            p_max_raw = max(int(m.get('Pos_5p', 0)), int(m.get('Pos_3p', 0)))
            
            # Ancoragem biológica exata de footprint de 12 bp
            centro = (p_min_raw + p_max_raw) / 2.0
            p_min = int(centro - (footprint_alvo / 2.0))
            p_max = int(centro + (footprint_alvo / 2.0))
            
            if not sitios_unicos:
                sitios_unicos.append({"inicio": p_min, "fim": p_max, "motivos": [m], "FT": ft})
            else:
                ultimo_sitio = sitios_unicos[-1]
                if p_min <= ultimo_sitio["fim"] + margem_bp:
                    ultimo_sitio["fim"] = max(ultimo_sitio["fim"], p_max)
                    ultimo_sitio["motivos"].append(m)
                else:
                    sitios_unicos.append({"inicio": p_min, "fim": p_max, "motivos": [m], "FT": ft})
        
        sitios_agrupados.extend(sitios_unicos)
        
    return sitios_agrupados

def resolver_caminho_fasta(file_uploader_obj=None):
    if file_uploader_obj is not None:
        nome_temp = "temp_user_fasta.gz" if file_uploader_obj.name.endswith(".gz") else "temp_user_fasta.fasta"
        with open(nome_temp, "wb") as f: f.write(file_uploader_obj.getvalue())
        return nome_temp
        
    candidatos = [
        "C_glabrata_CBS138_current_orf_genomic_1000.fasta.gz", "C_glabrata_CBS138_current_orf_genomic_1000.fasta",
        "C_glabrata_CBS138_current_1000bp_upstream.fasta.gz", "C_glabrata_CBS138_current_1000bp_upstream.fasta",
        "genoma_referencia.fasta", "genoma_referencia.fasta.gz"
    ]
    for c in candidatos:
        if os.path.exists(c): return c
            
    for ficheiro in os.listdir("."):
        if ficheiro.endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")): return ficheiro
            
    return "genoma_referencia.fasta"

def exibir_botoes_exportacao(df, nome_ficheiro):
    st.markdown(tr("**Opções de Exportação:**", "**Export Options:**"))
    col1, col2, _ = st.columns([1, 1, 3])
    with col1:
        formato = st.selectbox(tr("Selecione o formato:", "Select format:"), options=["CSV", "XLSX"], label_visibility="collapsed", key=f"sel_{nome_ficheiro}")
    with col2:
        if formato == "CSV":
            dados = df.to_csv(index=False).encode('utf-8')
            mime_type = "text/csv"; extensao = "csv"
        else:
            try:
                buffer = io.BytesIO()
                with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
                    df.to_excel(writer, index=False, sheet_name='Data' if st.session_state['lang'] == 'EN' else 'Dados')
                dados = buffer.getvalue()
                mime_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"; extensao = "xlsx"
            except Exception:
                st.error(tr("Erro: Executa 'pip install xlsxwriter' no terminal.", "Error: Run 'pip install xlsxwriter' in the terminal."))
                return 
        st.download_button(label=tr(f"📥 Exportar como {formato}", f"📥 Export as {formato}"), data=dados, file_name=f"{nome_ficheiro}.{extensao}", mime=mime_type, key=f"btn_{formato}_{nome_ficheiro}")

# ==============================================================================
# FUNÇÕES DE CACHE (PROCESSAMENTO GLOBAL)
# ==============================================================================
@st.cache_data(show_spinner=False)
def extrair_matrizes_globais(gff_content_bytes, caminho_fasta, prom_len):
    tf_patterns = {}
    for tf, consensus in TF_LIST:
        fwd_regex = seq_to_regex(consensus)
        rev_cons = str(Seq.Seq(consensus).reverse_complement())
        rev_regex = seq_to_regex(rev_cons)
        if fwd_regex == rev_regex: regex_obj = re.compile(fwd_regex, re.IGNORECASE)
        else: regex_obj = re.compile(f"{fwd_regex}|{rev_regex}", re.IGNORECASE)
        tf_patterns[(tf, consensus)] = regex_obj

    gene_coords = {}
    if gff_content_bytes.startswith(b"\x1f\x8b"):
        gff_str = gzip.decompress(gff_content_bytes).decode('utf-8', errors="ignore")
    else:
        gff_str = gff_content_bytes.decode('utf-8', errors="ignore")
        
    for line in gff_str.splitlines():
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        if len(parts) < 9: continue
        chrom, source, feature_type, start, end, score, strand, phase, attr = parts

        if feature_type in ["gene", "CDS"]:
            id_match = re.search(r"(?:ID|Name)=([^;]+)", attr)
            if not id_match: continue
            gene_id = id_match.group(1)
            if gene_id in gene_coords: continue

            chrom_limpo = chrom.split("_")[0]
            start, end = int(start), int(end)
            if strand == "+": prom_start, prom_end = max(1, start - prom_len), start - 1
            else: prom_start, prom_end = end + 1, end + prom_len

            gene_coords[gene_id] = {
                "Cromossoma": chrom_limpo, "Strand": strand, "Start": start, "End": end,
                "Promotor_Start": prom_start, "Promotor_End": prom_end, "Tamanho_Promotor": (prom_end - prom_start) + 1,
            }

    gene_tf_counts = defaultdict(lambda: defaultdict(int))
    gene_consensus_counts = defaultdict(lambda: defaultdict(int))
    consensus_totals = defaultdict(int)
    genes_processados = []

    def open_file_smart_local(filepath):
        with open(filepath, "rb") as test_f: magic_bytes = test_f.read(2)
        if magic_bytes == b"\x1f\x8b": return gzip.open(filepath, "rt", encoding="utf-8", errors="ignore")
        else: return open(filepath, "r", encoding="utf-8", errors="ignore")

    with open_file_smart_local(caminho_fasta) as f_handle:
        for record in SeqIO.parse(f_handle, "fasta"):
            gene_id = record.id; genes_processados.append(gene_id)
            promoter_seq = str(record.seq[:prom_len])
            motivos_encontrados = []

            for (tf, consensus), regex_obj in tf_patterns.items():
                for match in regex_obj.finditer(promoter_seq):
                    motivos_encontrados.append({"FT": tf, "Consenso": consensus, "Pos_5p": match.start(), "Pos_3p": match.end() - 1})
                    gene_consensus_counts[gene_id][f"{tf} ({consensus})"] += 1
                    consensus_totals[(tf, consensus)] += 1

            sitios_reais = agrupar_sitios_fisicos(motivos_encontrados, margem_bp=10, footprint_alvo=12)
            for sitio in sitios_reais: gene_tf_counts[gene_id][sitio["FT"]] += 1

    all_unique_tfs = sorted(list(set(tf for tf, _ in TF_LIST)))
    all_detailed_cols = [f"{tf} ({cons})" for tf, cons in TF_LIST]
    matriz_fts_rows, matriz_detalhada_rows = [], []

    for gene_id in sorted(genes_processados):
        coords = gene_coords.get(gene_id, {"Cromossoma": "N/A", "Strand": "N/A", "Start": "N/A", "End": "N/A", "Promotor_Start": "N/A", "Promotor_End": "N/A", "Tamanho_Promotor": prom_len})
        base_info = {"Gene_ID": gene_id, "Cromossoma": coords["Cromossoma"], "Start": coords["Start"], "End": coords["End"], "Promotor_Start": coords["Promotor_Start"], "Promotor_End": coords["Promotor_End"], "Strand": coords["Strand"], "Tamanho_Promotor": coords["Tamanho_Promotor"]}
        
        row_fts = base_info.copy()
        for tf in all_unique_tfs: row_fts[tf] = gene_tf_counts[gene_id][tf]
        matriz_fts_rows.append(row_fts)

        row_detalhada = base_info.copy()
        for col_name in all_detailed_cols: row_detalhada[col_name] = gene_consensus_counts[gene_id][col_name]
        matriz_detalhada_rows.append(row_detalhada)

    resumo_rows = [{"Fator de Transcrição": tf, "Consenso (IUPAC)": consensus, "Total de Binding Sites no Genoma": consensus_totals[(tf, consensus)]} for tf, consensus in TF_LIST]
    
    return pd.DataFrame(resumo_rows), pd.DataFrame(matriz_fts_rows), pd.DataFrame(matriz_detalhada_rows)

@st.cache_data(show_spinner=False)
def executar_mutagenese_in_silico(vcf_bytes, gff_bytes, caminho_fasta, prom_len):
    # 1. Preparar Regex dos TFs
    tf_patterns = {}
    for tf, consensus in TF_LIST:
        fwd_regex = seq_to_regex(consensus)
        rev_cons = str(Seq.Seq(consensus).reverse_complement())
        rev_regex = seq_to_regex(rev_cons)
        if fwd_regex == rev_regex:
            regex_obj = re.compile(fwd_regex, re.IGNORECASE)
        else:
            regex_obj = re.compile(f"{fwd_regex}|{rev_regex}", re.IGNORECASE)
        tf_patterns[(tf, consensus)] = regex_obj
        
    # 2. Processar GFF
    gene_coords = {}
    promoters_by_chrom = defaultdict(list)
    
    if gff_bytes.startswith(b"\x1f\x8b"):
        gff_str = gzip.decompress(gff_bytes).decode('utf-8', errors="ignore")
    else:
        gff_str = gff_bytes.decode('utf-8', errors="ignore")
        
    for line in gff_str.splitlines():
        if line.startswith("#"): continue
        parts = line.strip().split("\t")
        if len(parts) < 9: continue
        
        chrom, source, feature_type, start, end, score, strand, phase, attr = parts
        
        if feature_type in ["gene", "CDS"]:
            id_match = re.search(r"(?:ID|Name)=([^;]+)", attr)
            if not id_match: continue
            gene_id = id_match.group(1)
            
            if gene_id in gene_coords: continue
            
            chrom_norm = normalizar_chrom(chrom)
            start, end = int(start), int(end)
            
            if strand == "+":
                prom_start = max(1, start - prom_len)
                prom_end = start - 1
            else:
                prom_start = end + 1
                prom_end = end + prom_len
                
            info = {
                "Gene_ID": gene_id,
                "Cromossoma_Raw": chrom,
                "Cromossoma_Norm": chrom_norm,
                "Strand": strand,
                "Start": start,
                "End": end,
                "Promotor_Start": min(prom_start, prom_end),
                "Promotor_End": max(prom_start, prom_end)
            }
            
            gene_coords[gene_id] = info
            promoters_by_chrom[chrom_norm].append(info)

    # 3. Processar VCF para Haplótipos e FASTA WT
    if vcf_bytes.startswith(b"\x1f\x8b"):
        vcf_str = gzip.decompress(vcf_bytes).decode('utf-8', errors="ignore")
    else:
        vcf_str = vcf_bytes.decode('utf-8', errors="ignore")

    vcf_variants_by_gene = defaultdict(list)
    vcf_lines = []
    amostras = []
    for line in vcf_str.splitlines():
        if line.startswith("##"): continue
        if line.startswith("#CHROM"):
            parts = line.strip().split("\t")
            if len(parts) > 9: amostras = parts[9:]
            continue
        vcf_lines.append(line)

    for line in vcf_lines:
        parts = line.strip().split("\t")
        chrom_raw = parts[0]
        pos_vcf = int(parts[1])
        ref = parts[3]
        alts = parts[4].split(",")
        chrom_norm = normalizar_chrom(chrom_raw)
        
        promotores_candidatos = promoters_by_chrom.get(chrom_norm, [])
        if not promotores_candidatos: continue
        
        alt_isolates = defaultdict(list)
        if len(parts) > 9 and amostras:
            for i, am in enumerate(amostras):
                gt = parts[9+i].split(":")[0].replace("|", "/")
                for alt_idx in range(len(alts)):
                    if str(alt_idx+1) in gt.split("/"):
                        alt_isolates[alt_idx].append(am)
                        
        for p_info in promotores_candidatos:
            if p_info["Promotor_Start"] <= pos_vcf <= p_info["Promotor_End"]:
                for alt_idx, alt in enumerate(alts):
                    if alt == "*": continue
                    isos = alt_isolates[alt_idx]
                    if not isos: continue
                    vcf_variants_by_gene[p_info["Gene_ID"]].append({
                        "pos_vcf": pos_vcf, "ref": ref, "alt": alt, 
                        "strand": p_info["Strand"], "prom_start": p_info["Promotor_Start"], 
                        "prom_end": p_info["Promotor_End"], "gene_start": p_info["Start"],
                        "gene_end": p_info["End"], "chrom_raw": chrom_raw, "isolados": isos
                    })

    promoters_ref = {}
    wt_promoter_totals = defaultdict(lambda: defaultdict(int))
    mut_promoter_totals = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    
    def open_file_smart_local(filepath):
        with open(filepath, "rb") as test_f: magic_bytes = test_f.read(2)
        if magic_bytes == b"\x1f\x8b": return gzip.open(filepath, "rt", encoding="utf-8", errors="ignore")
        else: return open(filepath, "r", encoding="utf-8", errors="ignore")

    with open_file_smart_local(caminho_fasta) as f_handle:
        for record in SeqIO.parse(f_handle, "fasta"):
            gene_id = record.id
            prom_seq = str(record.seq[:prom_len]).upper()
            promoters_ref[gene_id] = prom_seq
            
            # --- TOTAIS WT REAIS ---
            motivos_wt = []
            for (tf, consensus), regex_obj in tf_patterns.items():
                for match in regex_obj.finditer(prom_seq):
                    motivos_wt.append({"FT": tf, "Pos_5p": match.start(), "Pos_3p": match.end() - 1})
            for s in agrupar_sitios_fisicos(motivos_wt, margem_bp=10, footprint_alvo=12):
                wt_promoter_totals[gene_id][s['FT']] += 1
                
            # --- TOTAIS HAPLÓTIPOS POR PACIENTE ---
            variants = vcf_variants_by_gene.get(gene_id, [])
            isos_in_gene = set()
            for v in variants: isos_in_gene.update(v["isolados"])
            
            for iso in isos_in_gene:
                iso_vars = [v for v in variants if iso in v["isolados"]]
                mut_seq = prom_seq
                vars_with_idx = []
                for v in iso_vars:
                    idx = v["pos_vcf"] - v["prom_start"] if v["strand"] == "+" else v["prom_end"] - v["pos_vcf"]
                    ref_allele = v["ref"] if v["strand"] == "+" else str(Seq.Seq(v["ref"]).reverse_complement())
                    alt_allele = v["alt"] if v["strand"] == "+" else str(Seq.Seq(v["alt"]).reverse_complement())
                    vars_with_idx.append((idx, ref_allele, alt_allele))
                
                vars_with_idx.sort(key=lambda x: x[0], reverse=True)
                for idx, ref_allele, alt_allele in vars_with_idx:
                    if mut_seq[idx:idx+len(ref_allele)] == ref_allele:
                        mut_seq = mut_seq[:idx] + alt_allele + mut_seq[idx+len(ref_allele):]
                        
                motivos_mut = []
                for (tf, consensus), regex_obj in tf_patterns.items():
                    for match in regex_obj.finditer(mut_seq):
                        motivos_mut.append({"FT": tf, "Pos_5p": match.start(), "Pos_3p": match.end() - 1})
                for s in agrupar_sitios_fisicos(motivos_mut, margem_bp=10, footprint_alvo=12):
                    mut_promoter_totals[gene_id][iso][s['FT']] += 1

    # 4. Avaliar Variantes Individuais (Janelas) vs Haplótipo Global
    resultados_variantes = []
    
    for gene_id, variants in vcf_variants_by_gene.items():
        prom_seq = promoters_ref.get(gene_id)
        if not prom_seq: continue
        
        # Agrupar mutações idênticas que afetam múltiplos isolados para evitar repetições na tabela VCF
        unique_vars = {}
        for v in variants:
            k = (v["pos_vcf"], v["ref"], v["alt"])
            if k not in unique_vars:
                unique_vars[k] = v
                unique_vars[k]["all_isos"] = set(v["isolados"])
            else:
                unique_vars[k]["all_isos"].update(v["isolados"])
                
        for (pos_vcf, ref, alt), v in unique_vars.items():
            idx = v["pos_vcf"] - v["prom_start"] if v["strand"] == "+" else v["prom_end"] - v["pos_vcf"]
            pos_relativa = v["pos_vcf"] - v["gene_start"] if v["strand"] == "+" else v["gene_end"] - v["pos_vcf"]
            
            def calc_rel_pos(idx_val):
                # O YEASTRACT utiliza um sistema 1-based, o Python usa 0-based.
                # Adicionamos +1 para calibrar perfeitamente a posição do SNP e os limites do sítio físico.
                if v["strand"] == "+":
                    return (v["prom_start"] + idx_val) - v["gene_start"] + 1
                else:
                    return v["gene_end"] - (v["prom_end"] - idx_val) + 1
                    
            ref_allele = v["ref"] if v["strand"] == "+" else str(Seq.Seq(v["ref"]).reverse_complement())
            alt_allele = v["alt"] if v["strand"] == "+" else str(Seq.Seq(v["alt"]).reverse_complement())
            
            window_start = max(0, idx - 20)
            window_end = min(len(prom_seq), idx + len(ref_allele) + 20)
            
            ref_window = prom_seq[window_start:window_end]
            mut_window = prom_seq[window_start:idx] + alt_allele + prom_seq[idx + len(ref_allele):window_end]
            
            motivos_ref = []
            for (tf, cons), regex_obj in tf_patterns.items():
                for match in regex_obj.finditer(ref_window):
                    motivos_ref.append({"FT": tf, "Pos_5p": window_start + match.start(), "Pos_3p": window_start + match.end() - 1})
            sitios_ref = agrupar_sitios_fisicos(motivos_ref, margem_bp=10, footprint_alvo=12)
            
            motivos_mut = []
            for (tf, cons), regex_obj in tf_patterns.items():
                for match in regex_obj.finditer(mut_window):
                    motivos_mut.append({"FT": tf, "Pos_5p": window_start + match.start(), "Pos_3p": window_start + match.end() - 1})
            sitios_mut = agrupar_sitios_fisicos(motivos_mut, margem_bp=10, footprint_alvo=12)
            
            tfs_in_window = set([s['FT'] for s in sitios_ref] + [s['FT'] for s in sitios_mut])
            
            for tf in tfs_in_window:
                s_ref_tf = [s for s in sitios_ref if s['FT'] == tf]
                s_mut_tf = [s for s in sitios_mut if s['FT'] == tf]
                
                count_ref = len(s_ref_tf)
                count_mut = len(s_mut_tf)
                
                if count_ref == count_mut: continue
                
                impacto = "🔴 PERDA DE FUNÇÃO" if count_ref > count_mut else "🟢 GANHO DE FUNÇÃO"
                dif = abs(count_ref - count_mut)
                
                s_target = s_ref_tf if impacto == "🔴 PERDA DE FUNÇÃO" else s_mut_tf
                s_str_list = []
                for s in s_target:
                    rel_5 = calc_rel_pos(s['inicio'])
                    rel_3 = calc_rel_pos(s['fim'])
                    s_str_list.append(f"[{min(rel_5, rel_3)} a {max(rel_5, rel_3)}]")
                consenso_afetado = "Sítio(s): " + " | ".join(s_str_list) if st.session_state['lang'] == 'PT' else "Site(s): " + " | ".join(s_str_list)
                
                tf_total_ref = wt_promoter_totals[gene_id][tf]
                
                isos_perda_total = []
                isos_ganho_novo = []
                isos_parcial = []
                
                for iso in v["all_isos"]:
                    tf_total_mut = mut_promoter_totals[gene_id][iso][tf]
                    if tf_total_ref > 0 and tf_total_mut == 0:
                        isos_perda_total.append(iso)
                    elif tf_total_ref == 0 and tf_total_mut > 0:
                        isos_ganho_novo.append(iso)
                    else:
                        isos_parcial.append(iso)
                        
                if impacto == "🔴 PERDA DE FUNÇÃO":
                    if isos_perda_total and not isos_parcial: status_exclusividade = "🚨 PERDA TOTAL"
                    elif isos_perda_total and isos_parcial: status_exclusividade = "📉 PERDA MISTA"
                    else: status_exclusividade = "📉 PERDA PARCIAL"
                else:
                    if isos_ganho_novo and not isos_parcial: status_exclusividade = "✨ GANHO DE NOVO"
                    elif isos_ganho_novo and isos_parcial: status_exclusividade = "📈 GANHO MISTO"
                    else: status_exclusividade = "📈 GANHO PARCIAL"
                    
                resultados_variantes.append({
                    "Gene_ID": gene_id,
                    "Fator de Transcrição": tf,
                    "Consenso Afetado": consenso_afetado,
                    "Impacto": impacto,
                    "Qtd_Sítios_Alterados": dif,
                    "TF_Total_Ref_Promotor": tf_total_ref,
                    "Status_Exclusividade": status_exclusividade,
                    "Cromossoma_VCF": v["chrom_raw"],
                    "Posição Absoluta": pos_vcf,
                    "Posição Relativa (TSS)": pos_relativa,
                    "REF": ref,
                    "ALT": alt,
                    "Isolados (Perda Total / De Novo)": ", ".join(sorted(isos_perda_total + isos_ganho_novo)) if (isos_perda_total or isos_ganho_novo) else "-",
                    "Isolados (Impacto Parcial)": ", ".join(sorted(isos_parcial)) if isos_parcial else "-",
                    "Todos Isolados Afetados": ", ".join(sorted(v["all_isos"]))
                })

    # Construir tabela Dashboard Eventos Extremos e Reguloma Dict
    eventos_extremos = []
    clean_regulon = {}
    all_tfs = sorted(list(set([t for t, c in TF_LIST])))
    
    for gene_id in promoters_ref.keys():
        variants = vcf_variants_by_gene.get(gene_id, [])
        isos_in_gene = set()
        for v in variants: isos_in_gene.update(v["all_isos"])
        
        for tf in all_tfs:
            ref_count = wt_promoter_totals[gene_id][tf]
            
            mut_dict_for_tf = {}
            for iso in isos_in_gene:
                mut_count = mut_promoter_totals[gene_id][iso][tf]
                if mut_count != ref_count:
                    mut_dict_for_tf[iso] = mut_count
                    
                if ref_count > 0 and mut_count == 0:
                    eventos_extremos.append({
                        "Gene_ID": gene_id, "Fator de Transcrição": tf, "Isolado": iso,
                        "Tipo": "🚨 PERDA TOTAL", "Sítios_Ref": ref_count
                    })
                elif ref_count == 0 and mut_count > 0:
                    eventos_extremos.append({
                        "Gene_ID": gene_id, "Fator de Transcrição": tf, "Isolado": iso,
                        "Tipo": "✨ GANHO DE NOVO", "Sítios_Mut": mut_count
                    })
                    
            if ref_count > 0 or len(mut_dict_for_tf) > 0:
                if tf not in clean_regulon: clean_regulon[tf] = {}
                clean_regulon[tf][gene_id] = {"WT": ref_count, "Mutados": mut_dict_for_tf}

    df_impacto = pd.DataFrame(resultados_variantes).sort_values(by=["Impacto", "Gene_ID", "Fator de Transcrição"], ascending=[True, True, True]) if resultados_variantes else pd.DataFrame()
    df_extremos = pd.DataFrame(eventos_extremos) if eventos_extremos else pd.DataFrame(columns=["Gene_ID", "Fator de Transcrição", "Isolado", "Tipo", "Sítios_Ref", "Sítios_Mut"])
    
    return df_impacto, df_extremos, clean_regulon, amostras

# ==============================================================================
# FUNÇÕES DE RENDERIZAÇÃO
# ==============================================================================
@st.fragment
def renderizar_painel_clinico(df_clinico_final, df_clinico_externo, df_matriz_atual, df_reg_atual, colunas_mutacoes, mapa_status_paciente, regiao_selecionada, coord_strand="N/A", limiar_mic_atual=32):
    # Definir os nomes exatos das colunas consoante a língua ativa para evitar KeyErrors
    col_diag = tr("Diagnóstico Clínico", "Clinical Diagnosis")
    col_iso = tr("Isolado Clínico", "Clinical Isolate")
    col_nr_mut = tr("Nr. de Mutações", "No. of Mutations")

    tab1, tab2, tab3 = st.tabs([
        tr("🔬 Inspeção por Amostra", "🔬 Sample Inspection"), 
        tr("📈 Predição Multi-FT", "📈 Multi-TF Prediction"), 
        tr("🧮 Matriz: Variações vs Diagnóstico", "🧮 Matrix: Variations vs Diagnosis")
    ])
    
    with tab1:
        st.subheader(tr("📋 Tabela Consolidada (Mutações Originais + Metadados)", "📋 Consolidated Table (Original Mutations + Metadata)"))
        df_tab1_exibicao = df_clinico_final.copy().drop(columns=["Valor_MIC_Numerico", "Idade_Numerica"], errors='ignore')
        st.dataframe(df_tab1_exibicao, use_container_width=True, hide_index=True)
        st.markdown("---")
        st.subheader(tr("📐 Teste de Hipóteses: Resistentes vs. Suscetíveis (MIC50 Fluconazol)", "📐 Hypothesis Testing: Resistant vs. Susceptible (MIC50 Fluconazole)"))
        
        if df_clinico_externo is not None:
            df_teste_mic = df_clinico_final.copy()
            # Usar a col_diag dinâmica em vez de hardcode
            mics_res = df_teste_mic[df_teste_mic[col_diag].str.contains("RESIST")]["Valor_MIC_Numerico"].dropna().values
            mics_sus = df_teste_mic[df_teste_mic[col_diag].str.contains("SUSC|SUSC")]["Valor_MIC_Numerico"].dropna().values
            
            if len(mics_res) > 1 and len(mics_sus) > 1:
                stat_levene, p_levene = levene(mics_res, mics_sus)
                variancias_iguais = p_levene >= 0.05
                if variancias_iguais:
                    nome_teste = tr("Teste t de Student", "Student's t-Test"); nome_stat = tr("Estatística t", "t-Statistic"); stat_teste, p_teste = ttest_ind(mics_res, mics_sus, equal_var=True)
                else:
                    nome_teste = tr("Teste U de Mann-Whitney", "Mann-Whitney U Test"); nome_stat = tr("Estatística U", "U-Statistic"); stat_teste, p_teste = mannwhitneyu(mics_res, mics_sus, alternative='two-sided')
                
                c1, c2, c3 = st.columns(3)
                c1.metric(tr("Variâncias Iguais? (Levene)", "Equal Variances? (Levene)"), tr("Sim", "Yes") if variancias_iguais else tr("Não", "No"), formatar_p_texto(p_levene), delta_color="off")
                c2.metric(f"{nome_stat} ({nome_teste})", f"{stat_teste:.3f}")
                c3.metric("p-value", formatar_p_value(p_teste))
                
                if p_teste < 0.05:
                    st.success(tr(
                        f"✅ **Diferença Significativa ({formatar_p_texto(p_teste)}):** Segundo o {nome_teste}, os valores de MIC50 entre isolados resistentes (média: {np.mean(mics_res):.1f}) e suscetíveis (média: {np.mean(mics_sus):.1f}) são estatisticamente diferentes.",
                        f"✅ **Significant Difference ({formatar_p_texto(p_teste)}):** According to {nome_teste}, the MIC50 values between resistant (mean: {np.mean(mics_res):.1f}) and susceptible isolates (mean: {np.mean(mics_sus):.1f}) are statistically different."
                    ))
                else:
                    st.info(tr(
                        f"ℹ️ **Sem Diferença Significativa (p >= 0.05):** Segundo o {nome_teste}, não foi possível comprovar diferença estatística entre os dois grupos.",
                        f"ℹ️ **No Significant Difference (p >= 0.05):** According to {nome_teste}, statistical difference between the two groups could not be proven."
                    ))
            else:
                st.warning(tr("⚠️ Não existem amostras suficientes em ambos os grupos para realizar os testes estatísticos.", "⚠️ Not enough samples in both groups to perform statistical tests."))
        else:
            st.info(tr("💡 Carrega o 'Ficheiro de Metadados' na barra lateral para realizar o teste estatístico.", "💡 Upload the 'Metadata File' in the sidebar to perform the statistical test."))

        st.markdown("---")
        exibir_botoes_exportacao(df_tab1_exibicao, f"tabela_consolidada_{regiao_selecionada.replace(' ', '_')}")
            
    with tab2:
        col_efeito = next((c for c in df_reg_atual.columns if c in ["Efeito", "Consequencia", "Potencial Regulatório / Efeito Predito"]), None)
        if col_efeito is not None and not df_reg_atual.empty and not df_matriz_atual.empty:
            mapa_nr_mutacoes = dict(zip(df_clinico_final[col_iso], df_clinico_final[col_nr_mut])) if col_nr_mut in df_clinico_final.columns else {}
            mapa_multi_ft = {}
            for _, row_reg in df_reg_atual.iterrows():
                col_nome = row_reg.get('Nome_Coluna')
                if not col_nome: continue
                mut_ex = str(row_reg.get('Mutacao_Exata', ''))
                ref_al = mut_ex.split("->")[0] if "->" in mut_ex else None
                efeito = str(row_reg.get(col_efeito, 'Sem informação'))
                if col_nome not in mapa_multi_ft: mapa_multi_ft[col_nome] = {'ref': ref_al, 'impactos': []}
                mapa_multi_ft[col_nome]['impactos'].append({'efeito': efeito, 'mut_exata': mut_ex})
            
            dados_matriz_tfs = []
            
            lista_tfs_alvo = sorted(list(set(tf for tf, _ in TF_LIST)))
            colunas_mut = [c for c in df_matriz_atual.columns if c != 'Pacientes']
            
            mapa_consensos_colunas = {}
            for tf in lista_tfs_alvo:
                primeiro_consenso = next((c for t, c in TF_LIST if t == tf), "")
                mapa_consensos_colunas[tf] = f"{tf} [{primeiro_consenso}...]"
            
            for _, row_iso in df_matriz_atual.iterrows():
                isolado = str(row_iso['Pacientes']).strip()
                status = mapa_status_paciente.get(isolado, "N/A")
                status_f = tr("🔴 RESISTENTE", "🔴 RESISTANT") if status == "RESISTENTE" else (tr("🟢 SUSCETÍVEL", "🟢 SUSCEPTIBLE") if status == "SUSCETÍVEL" else status)
                nr_mut = mapa_nr_mutacoes.get(isolado, 0)
                linha_matriz_tf = {"Isolado": isolado, col_nr_mut: nr_mut, tr("Diagnóstico", "Diagnosis"): status_f}
                impactos_acumulados = {tf_k: set() for tf_k in lista_tfs_alvo}
                
                for col in colunas_mut:
                    if col in mapa_multi_ft:
                        gt = str(row_iso[col]).strip()
                        if gt not in ["./.", ".", "N/A", "", "WT"] and not gt.startswith("0/0"):
                            for imp in mapa_multi_ft[col]['impactos']:
                                ef_str = str(imp['efeito'])
                                for frag in [f.strip() for f in ef_str.split("||")]:
                                    for tf_k in lista_tfs_alvo:
                                        tf_base = re.sub(r'[pP]$', '', tf_k)
                                        if f"[{tf_k}]" in frag or f"[{tf_base}]" in frag or frag.startswith(f"{tf_k}:") or frag.startswith(f"{tf_base}:") or f" {tf_k} " in frag or f" {tf_base} " in frag:
                                            if any(w in frag for w in ["Perda", "Abolição", "Abolicao", "Loss"]): impacts_add = tr("🔴 Abolição", "🔴 Abolition")
                                            elif any(w in frag for w in ["Ganho", "De Novo", "Gain"]): impacts_add = tr("🟢 Ganho de Sítio", "🟢 Gain of Site")
                                            elif any(w in frag for w in ["Aumento", "Indução", "Increase"]): impacts_add = tr("🟡 Aumento", "🟡 Increase")
                                            elif any(w in frag for w in ["Diminuição", "Redução", "Reduction"]): impacts_add = tr("🟠 Redução", "🟠 Reduction")
                                            elif any(w in frag for w in ["Modulação", "Modulacao", "Modulation"]): impacts_add = tr("🔵 Modulação", "🔵 Modulation")
                                            else: impacts_add = tr("🟣 Alterado", "🟣 Altered")
                                            impactos_acumulados[tf_k].add(impacts_add)
                
                for tf_k in lista_tfs_alvo:
                    col_nome_consenso = mapa_consensos_colunas[tf_k]
                    imps = sorted(list(impactos_acumulados[tf_k]))
                    if not imps: linha_matriz_tf[col_nome_consenso] = "⚪ WT"
                    elif len(imps) == 1: linha_matriz_tf[col_nome_consenso] = imps[0]
                    else: linha_matriz_tf[col_nome_consenso] = " + ".join(imps)
                dados_matriz_tfs.append(linha_matriz_tf)
            
            st.markdown(tr("#### 🏁 Matriz Compacta de Impacto Regulatório", "#### 🏁 Compact Matrix of Regulatory Impact"))
            df_matriz_tfs_exibicao = pd.DataFrame(dados_matriz_tfs)
            
            col_iso_view = "Isolado"
            col_diag_view = tr("Diagnóstico", "Diagnosis")
            
            cfg_colunas = {
                col_iso_view: st.column_config.TextColumn(col_iso_view, pinned=True, width="small"),
                col_nr_mut: st.column_config.NumberColumn(col_nr_mut, width="small"),
                col_diag_view: st.column_config.TextColumn(col_diag_view, width="small")
            }
            for tf_k, col_consenso in mapa_consensos_colunas.items():
                cfg_colunas[col_consenso] = st.column_config.TextColumn(label=col_consenso, width="medium")
            st.dataframe(df_matriz_tfs_exibicao, use_container_width=True, hide_index=True, column_config=cfg_colunas)
            exibir_botoes_exportacao(df_matriz_tfs_exibicao, f"matriz_tfs_{regiao_selecionada.replace(' ', '_')}")
        else:
            st.warning(tr("Não existem dados regulatórios disponíveis para esta região.", "No regulatory data available for this region."))

    with tab3:
        st.subheader(tr("🧮 Matriz: Análise estatística de associação (Mutação vs Resistência)", "🧮 Matrix: Statistical Association Analysis (Mutation vs Resistance)"))
        st.markdown(tr(
            f"Esta tabela avalia individualmente cada variante identificada na região promotora, quantificando o número exacto de isolados resistentes e suscetíveis portadores da alteração de acordo com o **limiar de resistência configurado na barra lateral ({limiar_mic_atual} µg/mL)**, e determina a significância estatística e a força da associação fenotípica.",
            f"This table evaluates each identified variant in the promoter region individually, quantifying the exact number of resistant and susceptible isolates carrying the alteration according to the **resistance threshold configured in the sidebar ({limiar_mic_atual} µg/mL)**, and determines the statistical significance and strength of the phenotypic association."
        ))
        st.markdown(tr(
            "ℹ️ **Nota Metodológica:** A significância estatística da associação entre o genótipo (Presença/Ausência da Variante) e o fenótipo (Resistente/Suscetível) é calculada através do **Teste Exato de Fisher** (quando frequências esperadas < 5 ou n < 20) ou do **Teste Qui-Quadrado de Pearson**.",
            "ℹ️ **Methodological Note:** The statistical significance of the association between genotype (Presence/Absence of Variant) and phenotype (Resistant/Susceptible) is calculated via the **Fisher's Exact Test** (when expected frequencies < 5 or n < 20) or **Pearson's Chi-Square Test**."
        ))

        if not df_clinico_final.empty and not df_matriz_atual.empty:
            mapa_info_var = {}
            if not df_reg_atual.empty:
                for _, r_reg in df_reg_atual.iterrows():
                    col_k = r_reg.get('Nome_Coluna')
                    if col_k:
                        mut_ex = str(r_reg.get('Mutacao_Exata', col_k)).strip()
                        pos_abs = r_reg.get('Posição Absoluta', '')
                        pos_rel_val = r_reg.get('Posição Relativa (TSS)', '')
                        try:
                            pos_num = int(float(pos_rel_val))
                            pos_str = f"Rel {pos_num} bp"
                        except (ValueError, TypeError):
                            pos_num = 0
                            pos_str = f"Pos {pos_abs}" if pos_abs != '' else str(col_k)
                        
                        if "->" in mut_ex:
                            var_ref_alt = mut_ex
                        else:
                            var_ref_alt = str(col_k)
                        
                        var_id = f"Pos {pos_abs} ({pos_str})" if pos_abs != '' else f"{col_k} ({pos_str})"
                        
                        rotulo_limpo = f"{pos_num} ({var_ref_alt})" if "->" in var_ref_alt else f"{pos_num}"
                        
                        mapa_info_var[col_k] = {
                            "identificador": var_id,
                            "rotulo_poster": rotulo_limpo,
                            "pos_num": pos_num,
                            "ref_alt": var_ref_alt
                        }
            
            colunas_mut = [c for c in df_matriz_atual.columns if c != 'Pacientes']
            linhas_matriz_assoc = []
            
            for col_mut in colunas_mut:
                info_var = mapa_info_var.get(col_mut, {"identificador": str(col_mut), "rotulo_poster": str(col_mut), "pos_num": 0, "ref_alt": str(col_mut)})
                
                gts_var = df_clinico_final[col_mut].astype(str).str.strip()
                mascara_mut = ~gts_var.isin(["./.", ".", "N/A", "", "WT"]) & ~gts_var.str.startswith("0/0")
                
                df_com_mut = df_clinico_final[mascara_mut]
                df_sem_mut = df_clinico_final[~mascara_mut]
                
                lista_isolados_mut = df_com_mut[col_iso].tolist()
                isolados_str = ", ".join([str(i) for i in lista_isolados_mut]) if lista_isolados_mut else tr("Nenhum", "None")

                res_mut = len(df_com_mut[df_com_mut[col_diag].str.contains("RESIST")])
                sus_mut = len(df_com_mut[df_com_mut[col_diag].str.contains("SUSC|SUSC")])
                
                res_wt = len(df_sem_mut[df_sem_mut[col_diag].str.contains("RESIST")])
                sus_wt = len(df_sem_mut[df_sem_mut[col_diag].str.contains("SUSC|SUSC")])
                
                total_amostras = len(df_clinico_final)
                total_portadores = res_mut + sus_mut
                freq_conservacao = (total_portadores / total_amostras) * 100 if total_amostras > 0 else 0
                
                if freq_conservacao == 100: status_conservacao = tr("🔒 Fixada / Conservada (100%)", "🔒 Fixed / Conserved (100%)")
                elif freq_conservacao >= 95: status_conservacao = tr(f"🛡️ Altamente Conservada ({freq_conservacao:.1f}%)", f"🛡️ Highly Conserved ({freq_conservacao:.1f}%)")
                else: status_conservacao = tr(f"⚡ Polimórfica ({freq_conservacao:.1f}%)", f"⚡ Polymorphic ({freq_conservacao:.1f}%)")

                if sus_mut > 0 and res_mut == 0: exclusividade = tr("🟢 Exclusivo Suscetíveis", "🟢 Exclusive to Susceptible")
                elif res_mut > 0 and sus_mut == 0: exclusividade = tr("🔴 Exclusivo Resistentes", "🔴 Exclusive to Resistant")
                elif sus_mut > 0 and res_mut > 0: exclusividade = tr("🟡 Misto (Presente em Ambos)", "🟡 Mixed (Present in Both)")
                else: exclusividade = tr("⚪ Sem Isolados Portadores", "⚪ No Carrier Isolates")
                
                tabela_contingencia = [[res_mut, sus_mut], [res_wt, sus_wt]]
                total_res = res_mut + res_wt
                total_sus = sus_mut + sus_wt
                total_mut = res_mut + sus_mut
                total_wt = res_wt + sus_wt
                
                if total_mut == 0 or total_wt == 0 or total_res == 0 or total_sus == 0:
                    teste_nome = "N/A"
                    p_val_str = "N/A"
                    p_val_raw = np.nan
                    or_str = "N/A"
                    phi_str = "0.000"
                else:
                    n_total = total_res + total_sus
                    e11 = (total_mut * total_res) / n_total
                    e12 = (total_mut * total_sus) / n_total
                    e21 = (total_wt * total_res) / n_total
                    e22 = (total_wt * total_sus) / n_total
                    
                    a, b, c, d = res_mut, sus_mut, res_wt, sus_wt
                    if a == 0 or b == 0 or c == 0 or d == 0:
                        or_val = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                        se_log_or = np.sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5))
                    else:
                        or_val = (a * d) / (b * c)
                        se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
                        
                    ic_inf = np.exp(np.log(or_val) - 1.96 * se_log_or)
                    ic_sup = np.exp(np.log(or_val) + 1.96 * se_log_or)
                    or_str = f"{or_val:.2f} (IC 95%: {ic_inf:.2f}-{ic_sup:.2f})"
                    
                    phi_val = (a * d - b * c) / np.sqrt(total_mut * total_wt * total_res * total_sus)
                    phi_str = f"{phi_val:.3f}"
                    
                    if min(e11, e12, e21, e22) < 5 or n_total < 20:
                        teste_nome = tr("Teste Exato de Fisher", "Fisher's Exact Test")
                        try: _, p_val_raw = fisher_exact(tabela_contingencia, alternative='two-sided')
                        except Exception: p_val_raw = np.nan
                    else:
                        teste_nome = tr("Teste Qui-Quadrado", "Chi-Square Test")
                        try: _, p_val_raw, _, _ = chi2_contingency(tabela_contingencia, correction=True)
                        except Exception: p_val_raw = np.nan
                    
                    if pd.isna(p_val_raw): p_val_str = "N/A"
                    elif p_val_raw < 0.001: p_val_str = "< 0.001 (*** Sig.)"
                    elif p_val_raw < 0.01: p_val_str = f"{p_val_raw:.4f} (** Sig.)"
                    elif p_val_raw < 0.05: p_val_str = f"{p_val_raw:.4f} (* Sig.)"
                    else: p_val_str = tr(f"{p_val_raw:.4f} (Não sig.)", f"{p_val_raw:.4f} (Not sig.)")
                
                linhas_matriz_assoc.append({
                    tr("Variante identificador", "Variant Identifier"): info_var["identificador"],
                    tr("Rótulo Poster", "Poster Label"): info_var["rotulo_poster"],
                    tr("Posição Numérica", "Numeric Position"): info_var["pos_num"],
                    tr("Estado de Conservação", "Conservation State"): status_conservacao,
                    tr("Variação REF/ALT", "REF/ALT Variation"): info_var["ref_alt"],
                    tr("Mutados Resistentes (R)", "Mutated Resistant (R)"): int(res_mut),
                    tr("Mutados Suscetíveis (S)", "Mutated Susceptible (S)"): int(sus_mut),
                    tr("Total WT (Frequência)", "Total WT (Frequency)"): int(res_wt + sus_wt),
                    tr("WT Resistentes (R)", "WT Resistant (R)"): int(res_wt),
                    tr("WT Suscetíveis (S)", "WT Susceptible (S)"): int(sus_wt),
                    tr("Exclusividade fenotípica", "Phenotypic Exclusivity"): exclusividade,
                    tr("Odds Ratio (IC 95%)", "Odds Ratio (95% CI)"): or_str,
                    tr("Correlação (Phi Φ)", "Correlation (Phi Φ)"): phi_str,
                    tr("Teste estatístico", "Statistical Test"): teste_nome,
                    tr("P-value e significância", "P-value & Significance"): p_val_str,
                    "_p_raw": p_val_raw if not pd.isna(p_val_raw) else 1.0,
                    tr("Isolados Mutados", "Mutated Isolates"): isolados_str
                })
            
            df_matriz_resumo = pd.DataFrame(linhas_matriz_assoc)
            col_mr = tr("Mutados Resistentes (R)", "Mutated Resistant (R)")
            df_matriz_resumo = df_matriz_resumo.sort_values(by=["_p_raw", col_mr], ascending=[True, False]).drop(columns=["_p_raw"]).reset_index(drop=True)
            
            if not df_matriz_resumo.empty:
                col_exc = tr("Exclusividade fenotípica", "Phenotypic Exclusivity")
                col_pos = tr("Posição Numérica", "Numeric Position")
                col_lbl = tr("Rótulo Poster", "Poster Label")
                col_ms = tr("Mutados Suscetíveis (S)", "Mutated Susceptible (S)")
                col_wtr = tr("WT Resistentes (R)", "WT Resistant (R)")
                col_wts = tr("WT Suscetíveis (S)", "WT Susceptible (S)")
                
                df_graf = df_matriz_resumo[df_matriz_resumo[col_exc] != tr("⚪ Sem Isolados Portadores", "⚪ No Carrier Isolates")].copy()
                df_graf = df_graf.sort_values(by=col_pos, ascending=True).head(20)

                if not df_graf.empty:
                    fig_assoc = go.Figure()
                    
                    fig_assoc.add_trace(go.Scatter(
                        x=df_graf[col_lbl], y=df_graf[col_mr], name="Mutated Resistant",
                        mode="lines+markers", marker=dict(color="#e11d48", size=22, opacity=1.0),
                        line=dict(color="#e11d48", width=1.2), hovertemplate="<b>Relative Coord:</b> %{x}<br><b>Mutated Resistant:</b> %{y}<extra></extra>"
                    ))
                    
                    fig_assoc.add_trace(go.Scatter(
                        x=df_graf[col_lbl], y=df_graf[col_ms], name="Mutated Susceptible",
                        mode="lines+markers", marker=dict(color="#10b981", size=22, opacity=1.0),
                        line=dict(color="#10b981", width=1.2), hovertemplate="<b>Relative Coord:</b> %{x}<br><b>Mutated Susceptible:</b> %{y}<extra></extra>"
                    ))

                    fig_assoc.add_trace(go.Scatter(
                        x=df_graf[col_lbl], y=df_graf[col_wtr], name="WT Resistant",
                        mode="lines+markers", marker=dict(color="#881337", size=22, opacity=1.0),
                        line=dict(color="#881337", width=1.2, dash="dot"), hovertemplate="<b>Relative Coord:</b> %{x}<br><b>WT Resistant:</b> %{y}<extra></extra>"
                    ))

                    fig_assoc.add_trace(go.Scatter(
                        x=df_graf[col_lbl], y=df_graf[col_wts], name="WT Susceptible",
                        mode="lines+markers", marker=dict(color="#064e3b", size=22, opacity=1.0),
                        line=dict(color="#064e3b", width=1.2, dash="dot"), hovertemplate="<b>Relative Coord:</b> %{x}<br><b>WT Susceptible:</b> %{y}<extra></extra>"
                    ))
                    
                    max_y_val = max(df_graf[col_mr].max(), df_graf[col_ms].max(), df_graf[col_wtr].max(), df_graf[col_wts].max())
                    limite_superior_y = max(int(max_y_val * 1.2), 5)

                    fig_assoc.update_layout(
                        plot_bgcolor="white", paper_bgcolor="white", font=dict(color="black", family="Arial, sans-serif"),
                        title=dict(
                            text=tr(
                                f"<b>Distribuição das Variantes Promotoras nos Fenótipos Clínicos (Mutados vs. WT)</b><br><span style='font-size:14px; color:#475569;'>Limiar de Resistência: ≥ {limiar_mic_atual} µg/mL | Ordenado por Posição 5' ➔ 3'</span>",
                                f"<b>Distribution of Promoter Variants Across Clinical Phenotypes (Mutated vs. WT)</b><br><span style='font-size:14px; color:#475569;'>Resistance Threshold: ≥ {limiar_mic_atual} µg/mL | Ordered by 5' ➔ 3' Genomic Position</span>"
                            ),
                            x=0.01, y=0.94, font=dict(color="black", size=19)
                        ),
                        xaxis=dict(
                            title=tr("<b>Mutação Promotora / Coordenada Relativa ao TSS</b>", "<b>Promoter Mutation / Relative Coordinate to TSS</b>"),
                            title_font=dict(color="black", size=15), tickfont=dict(color="black", size=13, family="Arial Black"),
                            showline=True, linewidth=1.5, linecolor='black', mirror=False, showgrid=True, gridcolor="#e2e8f0", gridwidth=1, zeroline=False, tickangle=-45
                        ),
                        yaxis=dict(
                            title=tr("<b>Número de Isolados (Portadores vs WT)</b>", "<b>Number of Isolates (Carriers vs. WT)</b>"),
                            title_font=dict(color="black", size=15), tickfont=dict(color="black", size=13),
                            range=[0, limite_superior_y], showline=True, linewidth=1.5, linecolor='black', mirror=False, showgrid=True, gridcolor="#e2e8f0", gridwidth=1, zeroline=False
                        ),
                        height=640, legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, font=dict(size=13, color="black"), bgcolor="rgba(255,255,255,0.9)", borderwidth=0), margin=dict(t=120, b=130, l=70, r=40)
                    )
                    
                    st.plotly_chart(fig_assoc, use_container_width=True, theme=None, config=config_plotly_export, key=f"plot_assoc_{regiao_selecionada}")

            col_im = tr("Isolados Mutados", "Mutated Isolates")
            cfg_df_resumo = {col_im: st.column_config.TextColumn(col_im, width="large")}
            st.dataframe(df_matriz_resumo.drop(columns=[col_lbl, col_pos], errors='ignore'), use_container_width=True, hide_index=True, column_config=cfg_df_resumo)
            
            st.markdown("---")
            with st.expander(tr("🔍 Painel de Inspeção Avançada e Diagnóstico de Frequências (Amostras vs Phenótipo)", "🔍 Advanced Inspection Panel and Frequency Diagnosis (Samples vs Phenotype)"), expanded=False):
                st.markdown(tr("##### 1. Consulta Interativa por Posição/Variante", "##### 1. Interactive Query by Position/Variant"))
                col_vi = tr("Variante identificador", "Variant Identifier")
                var_selecionada = st.selectbox(
                    tr("Selecione uma variante para inspecionar os isolados portadores em detalhe:", "Select a variant to inspect carrier isolates in detail:"),
                    options=df_matriz_resumo[col_vi].tolist()
                )
                if var_selecionada:
                    linha_var = df_matriz_resumo[df_matriz_resumo[col_vi] == var_selecionada].iloc[0]
                    var_col_vcf = next((k for k, v in mapa_info_var.items() if v["identificador"] == var_selecionada), var_selecionada)
                    
                    c_info1, c_info2, c_info3 = st.columns(3)
                    c_info1.metric(tr("Variante / Posição", "Variant / Position"), var_selecionada, linha_var[tr("Variação REF/ALT", "REF/ALT Variation")])
                    c_info2.metric(tr("Portadores Resistentes", "Resistant Carriers"), f"{linha_var[col_mr]} " + tr("isolados", "isolates"))
                    c_info3.metric(tr("Portadores Suscetíveis", "Susceptible Carriers"), f"{linha_var[col_ms]} " + tr("isolados", "isolates"))
                    
                    df_sub_var = df_clinico_final[~df_clinico_final[var_col_vcf].astype(str).str.strip().isin(["./.", ".", "N/A", "", "WT"]) & ~df_clinico_final[var_col_vcf].astype(str).str.strip().str.startswith("0/0")]
                    if not df_sub_var.empty:
                        st.dataframe(df_sub_var[[col_iso, var_col_vcf, col_diag, "MIC50 Fluconazole Real", "Source", "Gender"]], use_container_width=True, hide_index=True)
                    else:
                        st.info(tr("Nenhum isolado portador detetado com esta variante.", "No carrier isolates detected with this variant."))

                st.markdown("---")
                st.markdown(tr("##### 2. Diagnóstico Exato de Distribuição e Amostras em Falta", "##### 2. Exact Diagnosis of Distribution and Missing Samples"))
                universo_total = len(df_clinico_final)
                st.write(tr(f"**Universo Total da Coorte analisada na aba:** `{universo_total}` isolados.", f"**Total Universe of the Cohort analyzed in the tab:** `{universo_total}` isolates."))
                
                if var_selecionada:
                    var_col_vcf = next((k for k, v in mapa_info_var.items() if v["identificador"] == var_selecionada), var_selecionada)
                    gts = df_clinico_final[var_col_vcf].astype(str).str.strip()
                    tem_mut = ~gts.isin(["./.", ".", "N/A", "", "WT"]) & ~gts.str.startswith("0/0")
                    
                    isolados_com_var = df_clinico_final[tem_mut]
                    isolados_sem_var = df_clinico_final[~tem_mut]
                    
                    cnt_res = len(isolados_com_var[isolados_com_var[col_diag].str.contains("RESIST")])
                    cnt_sus = len(isolados_com_var[isolados_com_var[col_diag].str.contains("SUSC|SUSC")])
                    total_mutados = len(isolados_com_var)
                    em_falta_fenotipo = total_mutados - (cnt_res + cnt_sus)
                    
                    col_d1, col_d2, col_d3 = st.columns(3)
                    col_d1.metric(tr("Com Mutações Registadas", "With Recorded Mutations"), f"{total_mutados}", tr(f"de {universo_total} total", f"out of {universo_total} total"))
                    col_d2.metric(tr("Sem Mutação / Wild-Type (WT)", "No Mutation / Wild-Type (WT)"), f"{len(isolados_sem_var)}", f"{len(isolados_sem_var)/universo_total*100:.1f}%")
                    col_d3.metric(tr("Sem Fenótipo R/S (Em Falta)", "Missing R/S Phenotype"), f"{em_falta_fenotipo}", tr("Não R/S", "Not R/S"), delta_color="inverse")
                    
                    if em_falta_fenotipo > 0:
                        st.warning(tr(
                            f"⚠️ Atenção: Existem **{em_falta_fenotipo}** isolado(s) que têm esta mutação mas cujo diagnóstico não foi classificado formalmente como nem 'RESISTENTE' nem 'SUSCETÍVEL'.",
                            f"⚠️ Warning: There are **{em_falta_fenotipo}** isolate(s) that have this mutation but whose diagnosis was not formally classified as 'RESISTANT' or 'SUSCEPTIBLE'."
                        ))
                        perdidos_df = isolados_com_var[~isolados_com_var[col_diag].str.contains("RESIST|SUSC")]
                        st.dataframe(perdidos_df[[col_iso, var_col_vcf, col_diag, "MIC50 Fluconazole Real"]], use_container_width=True, hide_index=True)
                    else:
                        st.success(tr(
                            f"✅ Todos os **{total_mutados}** isolados portadores desta mutação estão 100% perfeitamente classificados!",
                            f"✅ All **{total_mutados}** carrier isolates for this mutation are 100% perfectly classified!"
                        ))

            st.markdown("---")
            exibir_botoes_exportacao(df_matriz_resumo.drop(columns=[col_lbl, col_pos], errors='ignore'), f"matriz_assoc_mut_res_{regiao_selecionada.replace(' ', '_')}")
        else:
            st.warning(tr("⚠️ Dados clínicos ou de mutações insuficientes para gerar a matriz de associação.", "⚠️ Insufficient clinical or mutational data to generate the association matrix."))

@st.cache_data(show_spinner=False)
def carregar_gff(conteudo_bytes, gene_alvo):
    gff_stream = io.BytesIO(conteudo_bytes)
    df_gff = pd.read_csv(gff_stream, sep='\t', comment='#', header=None, names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    
    # 1. Tenta extrair pela CDS e devolve apenas 1 linha
    df_cds = df_gff[(df_gff['type'] == 'CDS') & (df_gff['attributes'].str.contains(gene_alvo))].copy()
    if not df_cds.empty:
        return df_cds.sort_values(by='start', ascending=True).head(1) if df_cds.iloc[0]['strand'] == '+' else df_cds.sort_values(by='end', ascending=False).head(1)
    
    # 2. Fallback: Extrai gene ou mRNA, mas força a devolução de apenas 1 linha para evitar Regiões duplicadas
    df_fallback = df_gff[(df_gff['type'].isin(['gene', 'mRNA'])) & (df_gff['attributes'].str.contains(gene_alvo))].copy()
    if not df_fallback.empty:
        return df_fallback.head(1)
    
    return pd.DataFrame()

@st.cache_data(show_spinner=False)
def carregar_clinico(conteudo_bytes, nome_ficheiro):
    clinico_stream = io.BytesIO(conteudo_bytes)
    df = pd.read_excel(clinico_stream) if nome_ficheiro.endswith('.xlsx') else pd.read_csv(clinico_stream, sep=None, engine='python')
    df.columns = [c.strip() for c in df.columns]
    
    col_alvo = None
    for cand in ['isolado', 'isolate', 'sample', 'amostra', 'id', 'patient', 'paciente', 'strain', 'cepa']:
        for c in df.columns:
            if cand == c.lower().strip() or cand in c.lower().strip(): col_alvo = c; break
        if col_alvo: break
    if not col_alvo and len(df.columns) > 0: col_alvo = df.columns[0]
    df['id_norm_join'] = df[col_alvo].apply(normalizar_id_amostra)
    return df

@st.cache_data(show_spinner=False)
def processar_vcf_global(conteudo_bytes, lista_coordenadas):
    try: conteudo_texto = conteudo_bytes.decode('utf-8')
    except UnicodeDecodeError: conteudo_texto = gzip.decompress(conteudo_bytes).decode('utf-8')
    vcf_stream = io.StringIO(conteudo_texto)

    skip_rows = 0; colunas_vcf = []
    for line in vcf_stream:
        if line.startswith('#CHROM'): colunas_vcf = line.strip().split('\t'); break
        skip_rows += 1

    vcf_stream.seek(0)
    vcf_reader = pd.read_csv(vcf_stream, sep='\t', skiprows=skip_rows, chunksize=50000, dtype=str)
    
    chunks_uteis = []
    for chunk in vcf_reader:
        chunk = chunk.rename(columns={'#CHROM': '#CHROM'})
        chunk['chrom_normalizado_vcf'] = chunk['#CHROM'].apply(normalizar_chrom)
        chunk['POS'] = chunk['POS'].astype(int)
        mascara_global = pd.Series(False, index=chunk.index)
        for coord in lista_coordenadas:
            mascara_local = (chunk['chrom_normalizado_vcf'] == coord['chrom_normalizado']) & (chunk['POS'] >= coord['start']) & (chunk['POS'] <= coord['end'])
            mascara_global = mascara_global | mascara_local
        sub_chunk = chunk[mascara_global]
        if not sub_chunk.empty: chunks_uteis.append(sub_chunk)
            
    df_filtrado = pd.concat(chunks_uteis, ignore_index=True) if chunks_uteis else pd.DataFrame()
    return df_filtrado, (colunas_vcf[9:] if len(colunas_vcf) > 9 else [])

st.sidebar.header(tr("⚙️ Configurações do Utilizador", "⚙️ User Settings"))

GENE_ALVO = st.sidebar.text_input(
    tr("ID do Gene", "Gene ID"), 
    value="CAGL0M01760g",
    help=tr(
        "Insira o ID do gene. A aplicação irá localizar estas coordenadas no GFF para extrair e analisar automaticamente a região promotora a montante.",
        "Enter the gene ID. The app will locate these coordinates in the GFF to automatically extract and analyze the upstream promoter region."
    ),
    key="input_gene"
)
TAMANHO_PROMOTOR = st.sidebar.number_input(tr("Tamanho do Promotor (bp)", "Promoter Size (bp)"), value=1000, step=100, key="input_prom_size")

st.sidebar.subheader(tr("🎛️ Parâmetros Clínicos de Resistência", "🎛️ Clinical Resistance Parameters"))
limiar_mic = st.sidebar.slider(tr("Limiar de Resistência (MIC50 Fluconazol >= X)", "Resistance Threshold (MIC50 Fluconazole >= X)"), min_value=2, max_value=64, value=32, step=2, key="input_mic")

with st.sidebar.expander(tr("⚙️ Configurações Avançadas (Opcional)", "⚙️ Advanced Settings (Optional)"), expanded=False):
    st.markdown(tr("**Calibração de Coordenadas (YEASTRACT)**", "**Coordinate Calibration (YEASTRACT)**"))
    offset_calibracao = st.number_input(
        tr("Offset / Desvio 5' UTR (bp)", "Offset / 5' UTR Deviation (bp)"), 
        value=0, 
        step=1,
        help=tr(
            "Use este valor para calibrar desfasamentos entre as coordenadas relativas do YEASTRACT e o seu ficheiro GFF (ex: tamanho da 5' UTR).",
            "Use this value to calibrate offsets between YEASTRACT relative coordinates and your GFF file (e.g., 5' UTR size)."
        ),
        key="input_offset"
    )

is_apenas_exatos = True
val_rigor_pwm = 1.0

st.sidebar.subheader(tr("📂 Upload de Ficheiros Genómicos", "📂 Genomic Files Upload"))

gff_file = st.sidebar.file_uploader(tr("1. Ficheiro GFF", "1. GFF File"), type=["gff", "gff3"], key="up_gff")
vcf_original_file = st.sidebar.file_uploader(tr("2. Ficheiro VCF Original (Fonte das Variantes)", "2. Original VCF File (Variant Source)"), type=["vcf", "gz"], key="up_vcf")
dados_clinicos_file = st.sidebar.file_uploader(tr("3. Ficheiro de Metadados (Fluconazol, etc.)", "3. Metadata File (Fluconazole, etc.)"), type=["csv", "txt", "xlsx"], key="up_meta")
fasta_file = st.sidebar.file_uploader(tr("4. Ficheiro FASTA (Promotores/Genoma)", "4. FASTA File (Promoters/Genome)"), type=["fasta", "fa", "gz"], key="up_fasta")


# ==============================================================================
# ABAS DA APLICAÇÃO (UI PRINCIPAL)
# ==============================================================================
aba_analise_variantes, aba_matriz_global, aba_mutagenese, aba_reguloma = st.tabs([
    tr("🧬 Análise de Variantes Clínicas (VCF vs GFF)", "🧬 Clinical Variant Analysis (VCF vs GFF)"), 
    tr("🧩 Matriz Global de TFs (Genoma de Referência)", "🧩 Global TF Matrix (Reference Genome)"),
    tr("⚡ Mutagénese in Silico (VCF Total)", "⚡ In Silico Mutagenesis (Total VCF)"),
    tr("🌍 Pan & Core Regulomas", "🌍 Pan & Core Regulomes")
])

# ----------------- ABA 1: FUNCIONALIDADE ORIGINAL DO VCF -----------------
with aba_analise_variantes:
    if gff_file and vcf_original_file:
        with st.spinner(tr("A processar com regra cega IUPAC e coordenadas separadas...", "Processing with blind IUPAC rule and separated coordinates...")):
            try:
                caminho_fasta_resolvido = resolver_caminho_fasta(fasta_file)
                df_gff_filtrado = carregar_gff(gff_file.getvalue(), GENE_ALVO)

                if df_gff_filtrado.empty:
                    st.error(tr(f"❌ Erro: O gene '{GENE_ALVO}' não foi encontrado no GFF.", f"❌ Error: Gene '{GENE_ALVO}' was not found in the GFF."))
                else:
                    regioes_promotoras = {}
                    for idx, row in df_gff_filtrado.reset_index().iterrows():
                        chrom_gff = str(row['seqid']); chrom_normalizado = normalizar_chrom(chrom_gff)
                        strand = row['strand']; start = int(row['start']); end = int(row['end'])
                        
                        if strand == '+': p_start, p_end = max(1, start - TAMANHO_PROMOTOR), start - 1
                        else: p_start, p_end = end + 1, end + TAMANHO_PROMOTOR
                        
                        nome_regiao = tr(f"Região {idx+1} Chr{chrom_normalizado} ({p_start}-{p_end})", f"Region {idx+1} Chr{chrom_normalizado} ({p_start}-{p_end})")
                        regioes_promotoras[nome_regiao] = {
                            'chrom_original': chrom_gff, 'chrom_normalizado': chrom_normalizado, 
                            'start': p_start, 'end': p_end, 'strand': strand, 'start_gene': start, 'end_gene': end
                        }
                    
                    lista_coordenadas = list(regioes_promotoras.values())
                    df_vcf_global, amostras_originais = processar_vcf_global(vcf_original_file.getvalue(), lista_coordenadas)
                    
                    lista_cg = [str(nome).strip() for nome in amostras_originais]
                    mapa_original_para_cg = {str(nome).strip(): str(nome).strip() for nome in amostras_originais}

                    df_clinico_externo = None
                    if dados_clinicos_file:
                        df_clinico_externo = carregar_clinico(dados_clinicos_file.getvalue(), dados_clinicos_file.name)

                    sumario_por_regiao, matrizes_completas = {}, {}
                    dados_regua_por_regiao, analise_regulatoria_por_regiao = {}, {}
                    motivos_nativos_por_regiao = {}

                    for nome_regiao, coord in regioes_promotoras.items():
                        dados_pacientes = {cg: {} for cg in lista_cg}; total_missing = 0
                        tss = coord['start_gene'] if coord['strand'] == '+' else coord['end_gene']
                        
                        nativos_encontrados = mapear_motivos_nativos_referencia(
                            coord, tss, caminho_fasta=caminho_fasta_resolvido, 
                            gene_alvo=GENE_ALVO, offset_utr=offset_calibracao, tamanho_promotor=TAMANHO_PROMOTOR
                        )
                        motivos_nativos_por_regiao[nome_regiao] = nativos_encontrados
                        
                        df_sub_total = pd.DataFrame()
                        if not df_vcf_global.empty:
                            df_sub_total = df_vcf_global[
                                (df_vcf_global['chrom_normalizado_vcf'] == coord['chrom_normalizado']) & 
                                (df_vcf_global['POS'] >= coord['start']) & (coord['end'] >= df_vcf_global['POS'])
                            ]

                        detalhes_regulatorios, variantes_regua = [], []
                        if not df_sub_total.empty:
                            detalhes_regulatorios, variantes_regua = analisar_impacto_regulatorio(
                                df_sub_total, coord, tss, MUTACOES_FENOTIPO,
                                caminho_fasta=caminho_fasta_resolvido, tamanho_promotor=TAMANHO_PROMOTOR,
                                gene_alvo=GENE_ALVO, rigor_pwm=val_rigor_pwm, apenas_exatos=is_apenas_exatos,
                                offset_utr=offset_calibracao
                            )

                        if detalhes_regulatorios:
                            posicoes_processadas_missing = set()
                            for det in detalhes_regulatorios:
                                pos_m, col_nome = det["Posição Absoluta"], det["Nome_Coluna"]
                                row_vcf_especifica = df_sub_total[df_sub_total['POS'] == pos_m].iloc[0]
                                det["QUAL"] = row_vcf_especifica.get("QUAL", "N/A")
                                det["Mutacao_Exata"] = f"{row_vcf_especifica['REF']}->{row_vcf_especifica['ALT']}"
                                
                                for am_orig, nome_cg in mapa_original_para_cg.items():
                                    gt_traduzido = extrair_e_traduzir_gt(row_vcf_especifica[am_orig], row_vcf_especifica['REF'], row_vcf_especifica['ALT'])
                                    if gt_traduzido == "./." and (pos_m, am_orig) not in posicoes_processadas_missing:
                                        total_missing += 1; posicoes_processadas_missing.add((pos_m, am_orig))
                                    dados_pacientes[nome_cg][col_nome] = gt_traduzido
                            
                            total_snps = sum(1 for r in detalhes_regulatorios if str(r.get("Type")).upper() == "SNP")
                            total_indels = sum(1 for r in detalhes_regulatorios if str(r.get("Type")).upper() == "INDEL")
                            total_variantes = len(detalhes_regulatorios)
                        else:
                            total_snps, total_indels, total_variantes = 0, 0, 0

                        sumario_por_regiao[nome_regiao] = {
                            tr("Total Variantes", "Total Variants"): total_variantes, 
                            "SNPs": total_snps, 
                            "Indels": total_indels, 
                            tr("Valores Omissos", "Missing Values"): total_missing
                        }
                        dados_regua_por_regiao[nome_regiao] = detalhes_regulatorios
                        analise_regulatoria_por_regiao[nome_regiao] = pd.DataFrame(detalhes_regulatorios)
                        
                        matrizes_completas[nome_regiao] = pd.DataFrame.from_dict(dados_pacientes, orient='index').reset_index().rename(columns={'index': 'Pacientes'}) if total_variantes > 0 else pd.DataFrame(columns=['Pacientes'])

                    df_sumario_geral = pd.DataFrame.from_dict(sumario_por_regiao, orient='index').reset_index().rename(columns={'index': tr('Região Alvo Analisada', 'Analyzed Target Region')})
                    st.subheader(tr("📊 Estatísticas: Tabela de Frequências", "📊 Statistics: Frequency Table"))
                    st.dataframe(df_sumario_geral, use_container_width=True, hide_index=True)

                    st.markdown("---")
                    st.subheader(tr("🧠 Análise por Região Promotora", "🧠 Analysis by Promoter Region"))
                    abas_avancadas = st.tabs([f"🧬 {k}" for k in regioes_promotoras.keys()])
                    
                    for idx, (nome_regiao, df_reg) in enumerate(analise_regulatoria_por_regiao.items()):
                        with abas_avancadas[idx]:
                            nativos_lista = motivos_nativos_por_regiao.get(nome_regiao, [])
                            if nativos_lista:
                                sitios_agrupados = agrupar_sitios_fisicos(nativos_lista)
                                total_matches = len(nativos_lista)
                                total_loci = len(sitios_agrupados)
                                
                                exp_title = tr(
                                    f"🛡️ Prova de Base: {total_loci} Sítios de Ligação Reais Encontrados ({total_matches} matches de consensos)",
                                    f"🛡️ Base Proof: {total_loci} Real Binding Sites Found ({total_matches} consensus matches)"
                                )
                                with st.expander(exp_title, expanded=False):
                                    linhas_exibicao = []
                                    for idx_sitio, sitio in enumerate(sitios_agrupados, 1):
                                        consensos = list(set([m.get("Consenso_IUPAC", "") for m in sitio["motivos"]]))
                                        fitas = list(set([str(m.get("Fita", "")) for m in sitio["motivos"]]))
                                        coords_y = list(set([str(m.get("Coordenada_YEASTRACT", "")) for m in sitio["motivos"]]))
                                        
                                        linhas_exibicao.append({
                                            "Locus #": idx_sitio,
                                            tr("Fator de Transcrição", "Transcription Factor"): sitio["FT"],
                                            tr("Coordenada Física", "Physical Coordinate"): f"{sitio['inicio']} a {sitio['fim']} bp",
                                            tr("Posições YEASTRACT", "YEASTRACT Positions"): " / ".join(coords_y),
                                            tr("Consensos Envolvidos", "Involved Consensuses"): " | ".join(consensos),
                                            tr("Fitas DNA", "DNA Strands"): " / ".join(fitas),
                                            tr("Nº de Matches", "No. of Matches"): len(sitio["motivos"])
                                        })
                                        
                                    df_nativos_exib = pd.DataFrame(linhas_exibicao)
                                    st.dataframe(df_nativos_exib, use_container_width=True, hide_index=True)
                            else:
                                st.info(tr("ℹ️ Nenhum consenso nativo do YEASTRACT encontrado no FASTA para esta região promotora.", "ℹ️ No native YEASTRACT consensus found in the FASTA for this promoter region."))
                            
                            if df_reg.empty and not nativos_lista:
                                st.info(tr("ℹ️ Nenhuma variante ou sítio nativo detetado nesta região.", "ℹ️ No variant or native site detected in this region."))
                            else:
                                st.markdown(tr(f"#### 🗺️ Régua Física do Promotor: {nome_regiao}", f"#### 🗺️ Physical Promoter Ruler: {nome_regiao}"))
                                fig_regua = go.Figure()
                                x_inicio_rel = -int(TAMANHO_PROMOTOR)
                                
                                pos_cols_x = []
                                if not df_reg.empty: pos_cols_x.extend(df_reg["Posição Relativa (TSS)"].astype(int).tolist())
                                if nativos_lista: pos_cols_x.extend([m["Pos_5p"] for m in nativos_lista] + [m["Pos_3p"] for m in nativos_lista])
                                    
                                min_pos_rel = min(pos_cols_x) if pos_cols_x else x_inicio_rel
                                max_pos_rel = max(pos_cols_x) if pos_cols_x else 100
                                lim_inf_x = min(x_inicio_rel - 60, min_pos_rel - 50)
                                lim_sup_x = max(int(TAMANHO_PROMOTOR * 0.1), max_pos_rel + 50)

                                tfs_ativos_neste_promotor = set()
                                
                                if nativos_lista:
                                    for m in nativos_lista:
                                        nome_raw = str(m.get("FT", "")).strip()
                                        tf_clean = re.sub(r'[pP]$', '', nome_raw)
                                        tfs_ativos_neste_promotor.add(tf_clean)
                                
                                col_impacto = next((c for c in df_reg.columns if c in ["Efeito", "Consequencia", "Potencial Regulatório / Efeito Predito"]), "Consequencia")
                                if not df_reg.empty:
                                    for _, row_r in df_reg.iterrows():
                                        cons_texto = str(row_r.get(col_impacto, ""))
                                        for tf_item, _ in TF_LIST: 
                                            tf_base = re.sub(r'[pP]$', '', tf_item)
                                            if f"[{tf_base}]" in cons_texto or f"[{tf_base}p]" in cons_texto or f"{tf_base}:" in cons_texto or f" {tf_base} " in cons_texto:
                                                tfs_ativos_neste_promotor.add(tf_base)
                                
                                lista_tfs_ordenados = sorted(list(tfs_ativos_neste_promotor))
                                
                                tracks_tfs = ["Other Variants"] 
                                mapa_tf_track = {}
                                paleta_cores = px.colors.qualitative.Alphabet 
                                
                                for idx_dinamico, tf_chave in enumerate(lista_tfs_ordenados):
                                    track_idx = idx_dinamico + 1
                                    nome_exibicao = f"{tf_chave}p" if not tf_chave.endswith("g") else tf_chave
                                    cor_hex = paleta_cores[idx_dinamico % len(paleta_cores)] 
                                    
                                    tracks_tfs.append(nome_exibicao)
                                    mapa_tf_track[tf_chave] = (nome_exibicao, track_idx, cor_hex)

                                fig_regua.add_shape(type="rect", x0=0, y0=-0.5, x1=lim_sup_x, y1=len(tracks_tfs)-0.5, fillcolor="#e0e7ff", line=dict(color="#6366f1", width=2), layer="below")
                                fig_regua.add_vline(x=0, line_width=2, line_color="#10b981", line_dash="solid")
                                fig_regua.add_vline(x=x_inicio_rel, line_width=2, line_color="#3b82f6", line_dash="dash")

                                fig_regua.add_annotation(
                                    x=x_inicio_rel, y=len(tracks_tfs) - 0.3, text="<b>5' PROMOTER</b>",
                                    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, ax=0, ay=-30,
                                    font=dict(size=12, color="#1e40af")
                                )
                                fig_regua.add_annotation(
                                    x=0, y=len(tracks_tfs) - 0.3, text="<b>TSS (+1) &#10132; GENE (3')</b>",
                                    showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2, ax=0, ay=-30,
                                    font=dict(size=12, color="#065f46")
                                )
                                
                                nativos_por_tf = {k: [] for k in mapa_tf_track.keys()}
                                for m_nat in nativos_lista:
                                    nome_raw = str(m_nat.get("FT", "")).strip()
                                    tf_k = re.sub(r'[pP]$', '', nome_raw).strip()
                                    
                                    matched_tf = None
                                    if tf_k in mapa_tf_track:
                                        matched_tf = tf_k
                                    else:
                                        for k in mapa_tf_track.keys():
                                            if k.lower() == tf_k.lower() or k.lower() == nome_raw.lower():
                                                matched_tf = k; break
                                        if not matched_tf:
                                            for k in mapa_tf_track.keys():
                                                if len(k) >= 3 and (k.lower() in nome_raw.lower() or nome_raw.lower() in k.lower()):
                                                    matched_tf = k; break
                                    
                                    if matched_tf:
                                        nativos_por_tf[matched_tf].append(m_nat)

                                coordenadas_sem_barra = {k: [] for k in mapa_tf_track.keys()}
                                col_impacto = next((c for c in df_reg.columns if c in ["Efeito", "Consequencia", "Potencial Regulatório / Efeito Predito"]), "Consequencia")
                                
                                if not df_reg.empty:
                                    for _, row_r in df_reg.iterrows():
                                        cons_texto = str(row_r.get(col_impacto, ""))
                                        pos_rel_m = int(row_r.get("Posição Relativa (TSS)", 0))
                                        if any(w in cons_texto for w in ["Perda", "Abolição", "Abolicao", "Ganho", "De Novo"]):
                                            for k_tf in mapa_tf_track.keys():
                                                if f"[{k_tf}]" in cons_texto or f"[{k_tf}p]" in cons_texto or cons_texto.startswith(f"{k_tf}:") or cons_texto.startswith(f"{k_tf}p:"):
                                                    coordenadas_sem_barra[k_tf].append(pos_rel_m)
                                
                                x_limite_5p = x_inicio_rel 
                                x_limite_3p = max(0, max_pos_rel) if pos_cols_x else 0 

                                for i, track_nome in enumerate(tracks_tfs):
                                    tf_chave_track = "Other"
                                    for k, (lbl, idx_tr, c_hex) in mapa_tf_track.items():
                                        if idx_tr == i:
                                            tf_chave_track = k
                                            break
                                            
                                    if i == 0:
                                        fig_regua.add_trace(go.Scatter(
                                            x=[x_limite_5p, x_limite_3p], y=[i, i],
                                            mode="lines", line=dict(color="black", width=1.5, dash="solid"),
                                            name="Linha Base SNPs" if st.session_state['lang'] == 'PT' else "SNP Baseline", legendgroup=f"TF_{tf_chave_track}",
                                            showlegend=False, hoverinfo="skip"
                                        ))
                                    else:
                                        fig_regua.add_trace(go.Scatter(
                                            x=[x_limite_5p, x_limite_3p], y=[i + 0.08, i + 0.08],
                                            mode="lines", line=dict(color="black", width=1.5, dash="solid"),
                                            name="Fita Direta (+)" if st.session_state['lang'] == 'PT' else "Forward Strand (+)", legendgroup=f"TF_{tf_chave_track}",
                                            showlegend=False, hoverinfo="skip"
                                        ))
                                        
                                        fig_regua.add_trace(go.Scatter(
                                            x=[x_limite_5p, x_limite_3p], y=[i - 0.08, i - 0.08],
                                            mode="lines", line=dict(color="black", width=1.5, dash="solid"),
                                            name="Fita Reversa (-)" if st.session_state['lang'] == 'PT' else "Reverse Strand (-)", legendgroup=f"TF_{tf_chave_track}",
                                            showlegend=False, hoverinfo="skip"
                                        ))

                                for tf_chave, (track_label, track_idx, cor_hex) in mapa_tf_track.items():
                                    siti_tf = nativos_por_tf.get(tf_chave, [])
                                    if not siti_tf: continue
                                    
                                    for m in siti_tf:
                                        p_min = min(int(m["Pos_5p"]), int(m["Pos_3p"]))
                                        p_max = max(int(m["Pos_5p"]), int(m["Pos_3p"]))
                                        
                                        padding_x = 4 
                                        seletor_largura = p_max - p_min
                                        if seletor_largura < 28:
                                            centro = (p_min + p_max) // 2
                                            x_min_vis = centro - 18
                                            x_max_vis = centro + 18
                                        else:
                                            x_min_vis = p_min - padding_x
                                            x_max_vis = p_max + padding_x

                                        y0_val = track_idx - 0.40
                                        y1_val = track_idx + 0.40

                                        poly_x = [x_min_vis, x_max_vis, x_max_vis, x_min_vis, x_min_vis]
                                        poly_y = [y0_val, y0_val, y1_val, y1_val, y0_val]

                                        fita_raw = str(m.get('Fita', '')).strip()
                                        is_direta = fita_raw == "+" or "dir" in fita_raw.lower() or "forward" in fita_raw.lower()
                                        fita_simbolo = "➡️ Direta (+)" if is_direta and st.session_state['lang'] == 'PT' else ("➡️ Forward (+)" if is_direta else ("⬅️ Reversa (-)" if st.session_state['lang'] == 'PT' else "⬅️ Reverse (-)"))

                                        sobrepostos = [
                                            out for out in siti_tf 
                                            if max(p_min, min(int(out["Pos_5p"]), int(out["Pos_3p"]))) <= min(p_max, max(int(out["Pos_5p"]), int(out["Pos_3p"])))
                                        ]
                                        
                                        if len(sobrepostos) > 1:
                                            hover_text = tr(f"<b>🟢 {len(sobrepostos)} Sítios Nativos Sobrepostos ({track_label}):</b><br><br>", f"<b>🟢 {len(sobrepostos)} Overlapping Native Sites ({track_label}):</b><br><br>")
                                            for sub_i, s in enumerate(sobrepostos, 1):
                                                s_min = min(int(s["Pos_5p"]), int(s["Pos_3p"]))
                                                s_max = max(int(s["Pos_5p"]), int(s["Pos_3p"]))
                                                s_fita = str(s.get('Fita', '')).strip()
                                                s_dir = s_fita == "+" or "dir" in s_fita.lower() or "forward" in s_fita.lower()
                                                s_simbolo = "➡️ Direta (+)" if s_dir and st.session_state['lang'] == 'PT' else ("➡️ Forward (+)" if s_dir else ("⬅️ Reversa (-)" if st.session_state['lang'] == 'PT' else "⬅️ Reverse (-)"))
                                                hover_text += (f"<b>[{sub_i}] {s['FT']} | {s_simbolo}:</b> {s['Consenso_IUPAC']}<br>"
                                                               f"Seq FASTA: {s['Sequencia_FASTA']} ({s_min} a {s_max} bp)<br><br>")
                                            hover_text += tr(f"<b>Coordenada Genómica:</b> {m['Coordenada_YEASTRACT']}<extra></extra>", f"<b>Genomic Coordinate:</b> {m['Coordenada_YEASTRACT']}<extra></extra>")
                                        else:
                                            hover_text = (tr(f"<b>🟢 Sítio Nativo YEASTRACT: {track_label}</b><br>", f"<b>🟢 Native YEASTRACT Site: {track_label}</b><br>") +
                                                          tr(f"<b>Fita Genómica:</b> {fita_simbolo}<br>", f"<b>Genomic Strand:</b> {fita_simbolo}<br>") +
                                                          f"<b>Consenso IUPAC:</b> {m['Consenso_IUPAC']}<br>"
                                                          f"<b>Seq FASTA:</b> {m['Sequencia_FASTA']}<br>" +
                                                          tr(f"<b>Coordenada Real:</b> {p_min} bp a {p_max} bp<br>", f"<b>Real Coordinate:</b> {p_min} bp to {p_max} bp<br>") +
                                                          tr(f"<b>Coordenada Genómica:</b> {m['Coordenada_YEASTRACT']}<extra></extra>", f"<b>Genomic Coordinate:</b> {m['Coordenada_YEASTRACT']}<extra></extra>"))

                                        fig_regua.add_trace(go.Scatter(
                                            x=poly_x, y=poly_y,
                                            mode="lines", fill="toself", fillcolor=cor_hex,
                                            line=dict(width=1.5, color=cor_hex, dash="solid"),
                                            name=track_label, legendgroup=f"TF_{tf_chave}",
                                            showlegend=False, hoveron="fills+points", hovertemplate=hover_text
                                        ))

                                if not df_reg.empty:
                                    posicoes_por_track = {i: [] for i in range(len(tracks_tfs))}
                                    df_regua_graf = df_reg.copy()
                                    df_regua_graf["Posição Relativa (TSS)"] = df_regua_graf["Posição Relativa (TSS)"].astype(int)

                                    for idx_v, row_v in df_regua_graf.iterrows():
                                        pos_rel = int(row_v.get("Posição Relativa (TSS)", 0))
                                        
                                        val_5p_raw = row_v.get("Consenso_YEASTRACT_5p", pos_rel)
                                        try: pos_5p_consenso = int(val_5p_raw)
                                        except ValueError: pos_5p_consenso = pos_rel
                                            
                                        intervalo_consenso = row_v.get("Intervalo_Consenso", f"[{pos_rel}]")
                                        pos_abs = row_v.get("Posição Absoluta", "N/A")
                                        tipo_var = row_v.get("Type", "N/A")
                                        qualidade = row_v.get("QUAL", "N/A")
                                        mutacao = row_v.get("Mutacao_Exata", "N/A")
                                        consequencia = str(row_v.get(col_impacto, ""))
                                        
                                        tfs_impactados_nesta_var = []
                                        fragmentos_cons = [f.strip() for f in consequencia.split("||")]
                                        
                                        for tf_chave, (track_label, track_idx, cor_hex) in mapa_tf_track.items():
                                            frag_tf = next((f for f in fragmentos_cons if f"[{tf_chave}]" in f or f"[{tf_chave}p]" in f or f.startswith(f"{tf_chave}:") or f.startswith(f"{tf_chave}p:") or f" {tf_chave} " in f or f" {tf_chave}p " in f), None)
                                            if frag_tf is not None: tfs_impactados_nesta_var.append((track_label, track_idx, cor_hex, tf_chave, frag_tf))
                                        
                                        if not tfs_impactados_nesta_var:
                                            tfs_impactados_nesta_var.append(("Other Variants", 0, "#64748b", "Other", consequencia))
                                            
                                        for track_label, track_idx, cor_hex, tf_nome, frag_texto in tfs_impactados_nesta_var:
                                            simbolo_marcador = "circle"
                                            tamanho_marcador = 28
                                            cor_borda = "white"
                                            largura_borda = 1.5
                                            cor_preenchimento = cor_hex
                                            
                                            if any(w in frag_texto for w in ["Perda", "Abolição", "Abolicao", "Loss"]): 
                                                simbolo_marcador = "circle"
                                                tamanho_marcador = 28
                                                cor_preenchimento = "#dc2626"
                                                cor_borda = "#7f1d1d"
                                                largura_borda = 2.0
                                            elif any(w in frag_texto for w in ["Ganho", "De Novo", "Gain"]): 
                                                simbolo_marcador = "circle"
                                                tamanho_marcador = 28
                                                cor_preenchimento = "#10b981"
                                                cor_borda = "#064e3b"
                                                largura_borda = 2.0
                                            
                                            vizinhos_proximos = sum(1 for pos_ant in posicoes_por_track[track_idx] if abs(pos_5p_consenso - pos_ant) < 18)
                                            posicoes_por_track[track_idx].append(pos_5p_consenso)
                                            
                                            y_calc = track_idx
                                            if vizinhos_proximos > 0:
                                                sinal = 1 if vizinhos_proximos % 2 == 1 else -1
                                                passo = ((vizinhos_proximos + 1) // 2) * 0.22
                                                y_calc = max(track_idx - 0.38, min(track_idx + 0.38, track_idx + (sinal * passo)))
                                            
                                            fig_regua.add_trace(go.Scatter(
                                                x=[pos_5p_consenso], y=[y_calc], mode="markers",
                                                marker=dict(symbol=simbolo_marcador, size=tamanho_marcador, color=cor_preenchimento, opacity=1.0, line=dict(width=largura_borda, color=cor_borda)),
                                                name=f"Mutation ({track_label})", showlegend=False,
                                                legendgroup=f"TF_{tf_nome}",
                                                customdata=[[tipo_var, mutacao, qualidade, pos_abs, pos_rel, frag_texto, tf_nome, val_5p_raw, intervalo_consenso]],
                                                hovertemplate=("<b>🧬 %{customdata[0]} (%{customdata[6]}): %{customdata[1]}</b><br><br>"
                                                               "<b>YEASTRACT Consensus Target (5'):</b> %{customdata[7]}<br>"
                                                               "<b>Native Site Interval:</b> %{customdata[8]}<br>"
                                                               "<b>Relative SNP Coordinate:</b> %{customdata[4]} bp<br>"
                                                               "<b>Absolute Position (VCF):</b> %{customdata[3]}<br>"
                                                               "<b>TF Track:</b> %{customdata[6]}<br>"
                                                               "<b>Predicted Regulatory Impact:</b><br>%{customdata[5]}<extra></extra>")
                                            ))

                                fig_regua.add_trace(go.Scatter(
                                    x=[None], y=[None], mode="markers",
                                    marker=dict(symbol="circle", size=32, color="#64748b", line=dict(width=1.5, color="white")),
                                    name="Other Variants", showlegend=True
                                ))
                                fig_regua.add_trace(go.Scatter(
                                    x=[None], y=[None], mode="markers",
                                    marker=dict(symbol="circle", size=32, color="#10b981", line=dict(width=1.5, color="white")),
                                    name="Gain of Binding Site", showlegend=True
                                ))
                                fig_regua.add_trace(go.Scatter(
                                    x=[None], y=[None], mode="markers",
                                    marker=dict(symbol="circle", size=32, color="#dc2626", line=dict(width=1.5, color="white")),
                                    name="Loss of Binding Site", showlegend=True
                                ))

                                titulo_x = tr("Coordenadas Relativas ao TSS (bp)", "Relative Coordinates to TSS (bp)")
                                fig_regua.update_layout(
                                    plot_bgcolor="white", paper_bgcolor="white", font=dict(color="black", family="Arial, sans-serif"),
                                    title=dict(text=tr("Mapeamento Físico do Promotor: Sítios Nativos vs Variantes Clínicas", "Physical Mapping of Promoter: Native TF Sites vs. Clinical Variants"), y=0.98, font=dict(color="black", size=18)),
                                    xaxis=dict(title=titulo_x, title_font=dict(color="black", size=15), title_standoff=25, tickfont=dict(color="black", size=13), showgrid=False, gridcolor="#e2e8f0", zeroline=False, range=[lim_inf_x, lim_sup_x], dtick=100, showline=True, linewidth=1, linecolor='black', mirror=False),
                                    yaxis=dict(tickmode="array", tickvals=list(range(len(tracks_tfs))), ticktext=tracks_tfs, tickfont=dict(color="black", size=12, family="Arial Black"), showgrid=False, range=[-0.8, len(tracks_tfs) + 0.3], showline=True, linewidth=1, linecolor='black'),
                                    hovermode="closest", height=max(500, len(tracks_tfs)*70), margin=dict(l=140, r=40, t=70, b=150),
                                    legend=dict(
                                        font=dict(color="black", size=14), orientation="h", yanchor="top", y=-0.35, xanchor="center", x=0.5, 
                                        itemclick="toggle", itemdoubleclick="toggleothers", borderwidth=1, bordercolor="#cbd5e1", 
                                        bgcolor="rgba(255,255,255,0.9)", traceorder="normal"
                                    )
                                )
                                st.plotly_chart(fig_regua, use_container_width=True, theme=None, config=config_plotly_export, key=f"plot_regua_{idx}")

                    st.markdown("---")
                    regiao_selecionada = st.selectbox(tr("Selecione a Região Alvo Principal para Inspeção Clínica:", "Select the Main Target Region for Clinical Inspection:"), list(matrizes_completas.keys()))

                    if regiao_selecionada:
                        df_matriz_atual = matrizes_completas[regiao_selecionada]
                        df_reg_atual = analise_regulatoria_por_regiao.get(regiao_selecionada, pd.DataFrame())
                        
                        if not df_matriz_atual.empty:
                            classificacao_clinica, mapa_status_paciente = [], {} 
                            mapa_clinico_sanitizado = {}
                            if df_clinico_externo is not None and 'id_norm_join' in df_clinico_externo.columns:
                                for _, idx_row in df_clinico_externo.iterrows():
                                    chave_s = idx_row['id_norm_join']
                                    if chave_s: mapa_clinico_sanitizado[chave_s] = idx_row

                            col_diag = tr("Diagnóstico Clínico", "Clinical Diagnosis")
                            col_iso = tr("Isolado Clínico", "Clinical Isolate")
                            col_nr_mut = tr("Nr. de Mutações", "No. of Mutations")

                            for _, row_p in df_matriz_atual.iterrows():
                                paciente = str(row_p['Pacientes']).strip(); paciente_norm = normalizar_id_amostra(paciente)
                                mic_real, idade_real = 0, np.nan; source_val, gender_val = "N/A", "N/A"
                                
                                if df_clinico_externo is not None:
                                    sub_clinico_row = corresponder_amostra(paciente, mapa_clinico_sanitizado)
                                            
                                    if sub_clinico_row is not None:
                                        mapa_cols_min = {str(c).lower().strip(): c for c in df_clinico_externo.columns}
                                        
                                        col_mic = next((orig for min_c, orig in mapa_cols_min.items() if any(k in min_c for k in ['fluc', 'fcz', 'fluconazol']) and any(m in min_c for m in ['mic', 'val', '50', 'ug', 'mg', 'ml', 'num'])), None)
                                        if not col_mic: col_mic = next((orig for min_c, orig in mapa_cols_min.items() if any(k in min_c for k in ['fluc', 'fcz', 'fluconazol'])), None)
                                        if not col_mic: col_mic = next((orig for min_c, orig in mapa_cols_min.items() if 'mic50' in min_c or 'mic' in min_c), None)
                                        if col_mic and pd.notna(sub_clinico_row[col_mic]): mic_real = sub_clinico_row[col_mic]
                                        
                                        col_idade = next((orig for min_c, orig in mapa_cols_min.items() if 'idade' in min_c or 'age' in min_c), None)
                                        if col_idade and pd.notna(sub_clinico_row[col_idade]):
                                            try: idade_real = float(sub_clinico_row[col_idade])
                                            except ValueError: pass
                                            
                                        col_source = next((orig for min_c, orig in mapa_cols_min.items() if any(k in min_c for k in ['source', 'fonte', 'origem', 'specimen', 'material'])), None)
                                        if col_source and pd.notna(sub_clinico_row[col_source]): source_val = str(sub_clinico_row[col_source]).strip()
                                        
                                        col_gender = next((orig for min_c, orig in mapa_cols_min.items() if any(k in min_c for k in ['gender', 'género', 'genero', 'sexo', 'sex'])), None)
                                        if col_gender and pd.notna(sub_clinico_row[col_gender]): gender_val = str(sub_clinico_row[col_gender]).strip()

                                val_mic_num = extrair_valor_mic(mic_real)
                                status_clinico = "RESISTENTE" if val_mic_num >= limiar_mic else "SUSCETÍVEL"
                                mapa_status_paciente[paciente] = status_clinico
                                mutacoes_dict = {col: row_p[col] for col in df_matriz_atual.columns if col != 'Pacientes'}
                                nr_mutacoes = sum(1 for col_mut, gt_val in mutacoes_dict.items() if str(gt_val).strip() not in ["./.", ".", "N/A", "", "WT"] and not str(gt_val).strip().startswith("0/0"))
                                
                                dict_paciente = {
                                    col_iso: paciente, col_nr_mut: nr_mutacoes, "Source": source_val, "Gender": gender_val,
                                    tr("Idade", "Age"): idade_real if not pd.isna(idade_real) else "N/A", "MIC50 Fluconazole Real": f"{mic_real} µg/mL" if pd.notna(mic_real) else "0 µg/mL",
                                    col_diag: "🔴 RESISTENTE" if status_clinico == "RESISTENTE" else "🟢 SUSCETÍVEL", "Valor_MIC_Numerico": val_mic_num, "Idade_Numerica": idade_real
                                }
                                dict_paciente.update(mutacoes_dict); classificacao_clinica.append(dict_paciente)
                                
                            df_clinico_final = pd.DataFrame(classificacao_clinica)
                            renderizar_painel_clinico(df_clinico_final, df_clinico_externo, df_matriz_atual, df_reg_atual, [c for c in df_matriz_atual.columns if c != 'Pacientes'], mapa_status_paciente, regiao_selecionada, regioes_promotoras[regiao_selecionada]['strand'], limiar_mic_atual=limiar_mic)

            except Exception as e: st.error(tr(f"❌ Ocorreu um erro no processamento: {str(e)}", f"❌ Processing error occurred: {str(e)}"))
    else:
        st.info(tr("⚠️ Carregue os ficheiros GFF e VCF na barra lateral para analisar as variantes clínicas.", "⚠️ Upload GFF and VCF files in the sidebar to analyze clinical variants."))

# ----------------- ABA 2: NOVA FUNCIONALIDADE DA MATRIZ GLOBAL -----------------
with aba_matriz_global:
    st.header(tr("🧩 Geração de Matriz Global de Fatores de Transcrição", "🧩 Generation of Global Transcription Factor Matrix"))
    st.markdown(tr("Esta ferramenta varre todo o genoma de referência para mapear os sítios nativos.", "This tool scans the entire reference genome to map native sites."))
    
    if "matriz_global_processada" not in st.session_state:
        st.session_state.matriz_global_processada = False
    
    if gff_file and fasta_file:
        if st.button(tr("🚀 Gerar Matrizes Globais", "🚀 Generate Global Matrices"), type="primary"):
            st.session_state.matriz_global_processada = True 
            
        if st.session_state.matriz_global_processada: 
            with st.spinner(tr("A analisar todas as sequências e a agrupar sítios físicos. Isto pode demorar alguns segundos...", "Analyzing all sequences and grouping physical sites. This may take a few seconds...")):
                caminho_fasta_local = resolver_caminho_fasta(fasta_file)
                df_resumo, df_matriz_fts, df_matriz_detalhada = extrair_matrizes_globais(
                    gff_file.getvalue(), caminho_fasta_local, TAMANHO_PROMOTOR
                )
                
                st.success(tr("🎉 Matrizes geradas com sucesso!", "🎉 Matrices generated successfully!"))
                
                t1, t2, t3 = st.tabs([tr("Resumo", "Summary"), tr("Matriz Agrupada (Genes x FTs)", "Grouped Matrix (Genes x TFs)"), tr("Matriz Detalhada (Consensos)", "Detailed Matrix (Consensuses)")])
                t1.dataframe(df_resumo, hide_index=True, use_container_width=True)
                t2.dataframe(df_matriz_fts, hide_index=True, use_container_width=True)
                t3.dataframe(df_matriz_detalhada, hide_index=True, use_container_width=True)
                
                st.markdown("---")
                st.subheader(tr("📥 Exportar Resultados", "📥 Export Results"))
                
                formato = st.radio(tr("Escolha o formato de download:", "Choose download format:"), [tr("Excel (XLSX - Todas as folhas num ficheiro)", "Excel (XLSX - All sheets in one file)"), tr("CSV (Tabelas individuais)", "CSV (Individual tables)")], horizontal=True)
                
                if "Excel" in formato:
                    buffer = io.BytesIO()
                    with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
                        df_resumo.to_excel(writer, index=False, sheet_name='Summary' if st.session_state['lang'] == 'EN' else 'Resumo_por_Consenso')
                        df_matriz_fts.to_excel(writer, index=False, sheet_name='Matrix_Genes_x_TFs' if st.session_state['lang'] == 'EN' else 'Matriz_Genes_x_FTs')
                        df_matriz_detalhada.to_excel(writer, index=False, sheet_name='Detailed_Matrix' if st.session_state['lang'] == 'EN' else 'Matriz_por_Consenso_Exato')
                    st.download_button(tr("📥 Descarregar Ficheiro XLSX", "📥 Download XLSX File"), data=buffer.getvalue(), file_name="Resultados_Completos_FTs.xlsx" if st.session_state['lang'] == 'PT' else "Complete_TFs_Results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                else:
                    c1, c2, c3 = st.columns(3)
                    c1.download_button(tr("📥 Descarregar Resumo (CSV)", "📥 Download Summary (CSV)"), data=df_resumo.to_csv(index=False).encode('utf-8'), file_name="Resumo.csv", mime="text/csv")
                    c2.download_button(tr("📥 Descarregar Matriz Agrupada (CSV)", "📥 Download Grouped Matrix (CSV)"), data=df_matriz_fts.to_csv(index=False).encode('utf-8'), file_name="Matriz_Agrupada.csv", mime="text/csv")
                    c3.download_button(tr("📥 Descarregar Matriz Detalhada (CSV)", "📥 Download Detailed Matrix (CSV)"), data=df_matriz_detalhada.to_csv(index=False).encode('utf-8'), file_name="Matriz_Detalhada.csv", mime="text/csv")
    else:
        st.warning(tr("⚠️ Para utilizar o gerador de matrizes, precisa de carregar o **Ficheiro GFF** e o **Ficheiro FASTA** na barra lateral.", "⚠️ To use the matrix generator, you need to upload the **GFF File** and the **FASTA File** in the sidebar."))


# ----------------- ABA 3: MUTAGÉNESE IN SILICO (VCF VS GENOMA) -----------------
with aba_mutagenese:
    st.header(tr("⚡ Mutagénese in Silico (Genoma Completo vs. Pacientes)", "⚡ In Silico Mutagenesis (Full Genome vs. Patients)"))
    st.markdown(tr("Cruze todas as variantes dos seus pacientes com os promotores de todo o genoma para identificar ganhos e perdas de locais de ligação de FTs de forma massiva.", "Cross-reference all variants from your patients with the promoters of the entire genome to massively identify gains and losses of TF binding sites."))
    
    if "mutagenese_processada" not in st.session_state:
        st.session_state.mutagenese_processada = False
    
    if gff_file and vcf_original_file and fasta_file:
        if st.button(tr("🚀 Iniciar Análise de Mutagénese", "🚀 Start Mutagenesis Analysis"), type="primary"):
            st.session_state.mutagenese_processada = True
            
        if st.session_state.get("mutagenese_processada", False):
            with st.spinner(tr("A cruzar milhares de variantes com o genoma inteiro... (Footprints biológicos em processamento)", "Cross-referencing thousands of variants with the entire genome... (Biological footprints processing)")):
                caminho_fasta_local = resolver_caminho_fasta(fasta_file)
                
                df_mutagenese, df_extremos, dict_reguloma, amostras_vcf = executar_mutagenese_in_silico(
                    vcf_original_file.getvalue(), 
                    gff_file.getvalue(), 
                    caminho_fasta_local, 
                    TAMANHO_PROMOTOR
                )
                
                st.session_state.df_mutagenese_cache = df_mutagenese
                st.session_state.df_extremos_cache = df_extremos
                st.session_state.dict_reguloma_cache = dict_reguloma
                st.session_state.amostras_vcf_cache = amostras_vcf
                
                if not df_mutagenese.empty:
                    st.success(tr(f"🎉 Análise concluída com sucesso! Foram detetados {len(df_mutagenese)} eventos de disrupção funcional nos promotores.", f"🎉 Analysis completed successfully! {len(df_mutagenese)} functional disruption events in the promoters were detected."))
                    
                    total_perdas = int(df_mutagenese[df_mutagenese["Impacto"] == "🔴 PERDA DE FUNÇÃO"]["Qtd_Sítios_Alterados"].sum())
                    total_ganhos = int(df_mutagenese[df_mutagenese["Impacto"] == "🟢 GANHO DE FUNÇÃO"]["Qtd_Sítios_Alterados"].sum())
                    
                    perdas_totais = df_extremos[df_extremos["Tipo"] == "🚨 PERDA TOTAL"]
                    ganhos_novo = df_extremos[df_extremos["Tipo"] == "✨ GANHO DE NOVO"]
                    
                    # 1. RESUMO METRICAS PRINCIPAIS
                    st.markdown(tr("### 📊 Resumo de Impactos Regulatórios", "### 📊 Summary of Regulatory Impacts"))
                    col_m1, col_m2 = st.columns(2)
                    col_m1.metric(tr("🔴 Sítios Físicos Perdidos (Global)", "🔴 Lost Physical Sites (Global)"), f"{total_perdas}", tr(f"-{len(perdas_totais)} Perdas Totais (Gene-TF-Isolado)", f"-{len(perdas_totais)} Total Losses (Gene-TF-Isolate)"), delta_color="inverse")
                    col_m2.metric(tr("🟢 Sítios Físicos Ganhos (Global)", "🟢 Gained Physical Sites (Global)"), f"{total_ganhos}", tr(f"+{len(ganhos_novo)} Ganhos De Novo (Gene-TF-Isolado)", f"+{len(ganhos_novo)} De Novo Gains (Gene-TF-Isolate)"), delta_color="normal")
                    
                    st.markdown("---")
                    
                    # 2. DASHBOARD DE TOP 10 (GRÁFICOS DE BARRAS)
                    st.markdown(tr("#### 🏆 Top 10 Fatores de Transcrição (Eventos Extremos)", "#### 🏆 Top 10 Transcription Factors (Extreme Events)"))
                    col_chart1, col_chart2 = st.columns(2)
                    
                    with col_chart1:
                        if not perdas_totais.empty:
                            top_perdas = perdas_totais.groupby('Fator de Transcrição')['Gene_ID'].nunique().reset_index(name='Contagem')
                            top_perdas = top_perdas.sort_values(by='Contagem', ascending=False).head(10)
                            fig_p = px.bar(top_perdas, x='Contagem', y='Fator de Transcrição', orientation='h', 
                                           title=tr('🚨 Top 10: Perdas Totais (Por Gene Afetado)', '🚨 Top 10: Total Losses (Per Affected Gene)'), color_discrete_sequence=['#ef4444'])
                            fig_p.update_layout(yaxis={'categoryorder':'total ascending'}, height=400, plot_bgcolor="rgba(0,0,0,0)")
                            st.plotly_chart(fig_p, use_container_width=True, config=config_plotly_export)
                            
                    with col_chart2:
                        if not ganhos_novo.empty:
                            top_ganhos = ganhos_novo.groupby('Fator de Transcrição')['Gene_ID'].nunique().reset_index(name='Contagem')
                            top_ganhos = top_ganhos.sort_values(by='Contagem', ascending=False).head(10)
                            fig_g = px.bar(top_ganhos, x='Contagem', y='Fator de Transcrição', orientation='h', 
                                           title=tr('✨ Top 10: Ganhos De Novo (Por Gene Afetado)', '✨ Top 10: De Novo Gains (Per Affected Gene)'), color_discrete_sequence=['#10b981'])
                            fig_g.update_layout(yaxis={'categoryorder':'total ascending'}, height=400, plot_bgcolor="rgba(0,0,0,0)")
                            st.plotly_chart(fig_g, use_container_width=True, config=config_plotly_export)

                    st.markdown("---")
                    
                    # 3. TABELAS DETALHADAS DE EXCLUSIVIDADE (AGRUPADAS POR GENE-TF)
                    st.markdown(tr("### 🔎 Análise de Exclusividade Detalhada (Gene vs Fator de Transcrição)", "### 🔎 Detailed Exclusivity Analysis (Gene vs Transcription Factor)"))
                    st.markdown(tr(
                        "Identificação de promotores e as respetivas populações de isolados que sofreram **Perda Total** (tinham sítios físicos e perderam todos) ou **Ganho De Novo** (não tinham nenhum sítio físico e ganharam).",
                        "Identification of promoters and their respective isolate populations that suffered **Total Loss** (had physical sites and lost all) or **De Novo Gain** (had no physical site and gained)."
                    ))

                    col_ex1, col_ex2 = st.columns(2)
                    with col_ex1:
                        st.error(tr(f"🚨 **Perdas Totais Consolidadas**", f"🚨 **Consolidated Total Losses**"))
                        if not perdas_totais.empty:
                            df_pt_show = perdas_totais.groupby(["Gene_ID", "Fator de Transcrição", "Sítios_Ref"]).agg(
                                Isolados_Afetados=('Isolado', lambda x: ', '.join(sorted(x))),
                                Qtd_Isolados=('Isolado', 'count')
                            ).reset_index().sort_values(by="Qtd_Isolados", ascending=False)
                            # Translating column names before display
                            df_pt_show = df_pt_show.rename(columns={"Fator de Transcrição": tr("Fator de Transcrição", "Transcription Factor"), "Sítios_Ref": tr("Sítios (Ref)", "Sites (Ref)"), "Isolados_Afetados": tr("Isolados Afetados", "Affected Isolates"), "Qtd_Isolados": tr("Qtd Isolados", "Qty Isolates")})
                            st.dataframe(df_pt_show, use_container_width=True, hide_index=True)
                        else:
                            st.write(tr("Nenhuma perda total detetada.", "No total loss detected."))

                    with col_ex2:
                        st.success(tr(f"✨ **Ganhos De Novo Consolidados**", f"✨ **Consolidated De Novo Gains**"))
                        if not ganhos_novo.empty:
                            df_gn_show = ganhos_novo.groupby(["Gene_ID", "Fator de Transcrição", "Sítios_Mut"]).agg(
                                Isolados_Afetados=('Isolado', lambda x: ', '.join(sorted(x))),
                                Qtd_Isolados=('Isolado', 'count')
                            ).reset_index().sort_values(by="Qtd_Isolados", ascending=False)
                            df_gn_show = df_gn_show.rename(columns={"Fator de Transcrição": tr("Fator de Transcrição", "Transcription Factor"), "Sítios_Mut": tr("Sítios (Mut)", "Sites (Mut)"), "Isolados_Afetados": tr("Isolados Afetados", "Affected Isolates"), "Qtd_Isolados": tr("Qtd Isolados", "Qty Isolates")})
                            st.dataframe(df_gn_show, use_container_width=True, hide_index=True)
                        else:
                            st.write(tr("Nenhum ganho de novo detetado.", "No de novo gain detected."))
                    
                    st.markdown("---")
                    
                    # 4. TABELA GERAL
                    st.markdown(tr("### 📋 Tabela Geral de Impactos Regulatórios (Footprints Físicos YEASTRACT)", "### 📋 General Table of Regulatory Impacts (YEASTRACT Physical Footprints)"))
                    
                    df_mut_show = df_mutagenese.rename(columns={
                        "Fator de Transcrição": tr("Fator de Transcrição", "Transcription Factor"),
                        "Consenso Afetado": tr("Consenso Afetado", "Affected Consensus"),
                        "Impacto": tr("Impacto", "Impact"),
                        "Qtd_Sítios_Alterados": tr("Qtd Sítios Alterados", "Qty Altered Sites"),
                        "Status_Exclusividade": tr("Status Exclusividade", "Exclusivity Status"),
                        "Posição Absoluta": tr("Posição Absoluta", "Absolute Position"),
                        "Posição Relativa (TSS)": tr("Posição Relativa (TSS)", "Relative Position (TSS)"),
                        "Isolados (Perda Total / De Novo)": tr("Isolados (Perda Total / De Novo)", "Isolates (Total Loss / De Novo)"),
                        "Isolados (Impacto Parcial)": tr("Isolados (Impacto Parcial)", "Isolates (Partial Impact)"),
                        "Todos Isolados Afetados": tr("Todos Isolados Afetados", "All Affected Isolates")
                    })
                    st.dataframe(df_mut_show, use_container_width=True, hide_index=True)
                    
                    exibir_botoes_exportacao(df_mut_show, "Impacto_Variantes_Pacientes")
                else:
                    st.info(tr("ℹ️ A análise terminou, mas nenhuma variante no VCF fornecido provocou Ganhos ou Perdas de Função nos FTs mapeados.", "ℹ️ Analysis completed, but no variant in the provided VCF caused Gain or Loss of Function in the mapped TFs."))
    else:
        st.warning(tr("⚠️ Para utilizar a Mutagénese in Silico, é obrigatório carregar os ficheiros **GFF**, **VCF** e **FASTA** na barra lateral.", "⚠️ To use In Silico Mutagenesis, uploading **GFF**, **VCF**, and **FASTA** files in the sidebar is required."))

# ----------------- ABA 4: PAN & CORE REGULOMAS -----------------
with aba_reguloma:
    st.header(tr("🌍 Pan & Core Regulomas por Fator de Transcrição", "🌍 Pan & Core Regulomes by Transcription Factor"))
    st.markdown(tr(
        "Nesta aba, avaliamos o perfil de ligação ao nível do **Isolado e Gene**, permitindo identificar se um Fator de Transcrição regula um conjunto fixo de genes em todos os pacientes (Core Regulon) ou se essa regulação é fluida/variável (Accessory Regulon).",
        "In this tab, we evaluate the binding profile at the **Isolate and Gene** level, allowing us to identify if a Transcription Factor regulates a fixed set of genes across all patients (Core Regulon) or if this regulation is fluid/variable (Accessory Regulon)."
    ))
    
    if st.session_state.get("mutagenese_processada", False) and "dict_reguloma_cache" in st.session_state:
        regulon_dict = st.session_state.dict_reguloma_cache
        amostras_vcf = st.session_state.amostras_vcf_cache
        
        tfs_disponiveis = sorted(list(regulon_dict.keys()))
        if not tfs_disponiveis:
            st.info(tr("Não há dados de reguloma disponíveis para os parâmetros selecionados.", "No regulome data available for the selected parameters."))
        else:
            tf_selecionado = st.selectbox(tr("🎯 Selecione o Fator de Transcrição para construir a matriz do Reguloma:", "🎯 Select the Transcription Factor to build the Regulome matrix:"), tfs_disponiveis)
            
            if tf_selecionado:
                dados_tf = regulon_dict[tf_selecionado]
                
                linhas_reguloma = []
                for gene, info in dados_tf.items():
                    wt_count = info["WT"]
                    mut_dict = info["Mutados"]
                    
                    linha = {"Gene_ID": gene, tr("Sítios Físicos (WT)", "Physical Sites (WT)"): wt_count}
                    
                    counts_isolados = []
                    for iso in amostras_vcf:
                        c = mut_dict.get(iso, wt_count)
                        linha[iso] = c
                        counts_isolados.append(c)
                        
                    tem_em_todos = all(c > 0 for c in counts_isolados)
                    tem_em_nenhum = all(c == 0 for c in counts_isolados)
                    
                    if tem_em_todos:
                        linha[tr("Classificação Regulamentar", "Regulatory Classification")] = tr("🟢 Core Regulon (Presente em todos)", "🟢 Core Regulon (Present in all)")
                    elif tem_em_nenhum:
                        linha[tr("Classificação Regulamentar", "Regulatory Classification")] = tr("⚪ Ausente", "⚪ Absent")
                    else:
                        linha[tr("Classificação Regulamentar", "Regulatory Classification")] = tr("🟡 Accessory Regulon (Variável)", "🟡 Accessory Regulon (Variable)")
                        
                    linhas_reguloma.append(linha)
                    
                df_reguloma = pd.DataFrame(linhas_reguloma)
                if not df_reguloma.empty:
                    df_reguloma = df_reguloma[df_reguloma[tr("Classificação Regulamentar", "Regulatory Classification")] != tr("⚪ Ausente", "⚪ Absent")]
                
                if not df_reguloma.empty:
                    cols = ["Gene_ID", tr("Classificação Regulamentar", "Regulatory Classification"), tr("Sítios Físicos (WT)", "Physical Sites (WT)")] + amostras_vcf
                    df_reguloma = df_reguloma[cols].sort_values(by=[tr("Classificação Regulamentar", "Regulatory Classification"), "Gene_ID"])
                    
                    total_pan = len(df_reguloma)
                    total_core = len(df_reguloma[df_reguloma[tr("Classificação Regulamentar", "Regulatory Classification")] == tr("🟢 Core Regulon (Presente em todos)", "🟢 Core Regulon (Present in all)")])
                    total_acc = len(df_reguloma[df_reguloma[tr("Classificação Regulamentar", "Regulatory Classification")] == tr("🟡 Accessory Regulon (Variável)", "🟡 Accessory Regulon (Variable)")])
                    
                    c1, c2, c3 = st.columns(3)
                    c1.metric(tr("🌍 Pan Regulon (Total de Genes)", "🌍 Pan Regulon (Total Genes)"), total_pan)
                    c2.metric(tr("🟢 Core Regulon (Conservado)", "🟢 Core Regulon (Conserved)"), total_core)
                    c3.metric(tr("🟡 Accessory Regulon (Variável)", "🟡 Accessory Regulon (Variable)"), total_acc)
                    
                    st.markdown(tr("#### 🧬 Matriz de Presença/Quantidade de Sítios por Isolado", "#### 🧬 Matrix of Presence/Quantity of Sites per Isolate"))
                    
                    col_reg = tr("Classificação Regulamentar", "Regulatory Classification")
                    col_sitios = tr("Sítios Físicos (WT)", "Physical Sites (WT)")
                    
                    cfg_colunas = {
                        "Gene_ID": st.column_config.TextColumn("Gene_ID", pinned=True, width="small"),
                        col_reg: st.column_config.TextColumn(col_reg, width="medium"),
                        col_sitios: st.column_config.NumberColumn(col_sitios, width="small")
                    }
                    for iso in amostras_vcf:
                        cfg_colunas[iso] = st.column_config.NumberColumn(iso, width="small")
                        
                    st.dataframe(df_reguloma, use_container_width=True, hide_index=True, column_config=cfg_colunas)
                    
                    exibir_botoes_exportacao(df_reguloma, f"Reguloma_{tf_selecionado}")
                else:
                    st.info(tr(f"O fator de transcrição **{tf_selecionado}** não regula nenhum promotor em nenhuma amostra ou o seu sinal foi completamente extinto pelas mutações em toda a coorte.",
                               f"The transcription factor **{tf_selecionado}** does not regulate any promoter in any sample or its signal was completely extinguished by mutations across the entire cohort."))
    else:
        st.info(tr("Execute a Mutagénese in Silico na aba anterior para processar todo o genoma e gerar os dados do Reguloma.", "Run the In Silico Mutagenesis in the previous tab to process the entire genome and generate Regulome data."))