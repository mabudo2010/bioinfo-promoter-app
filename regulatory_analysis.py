import pandas as pd
import numpy as np
import os
import gzip

# =====================================================================
# 1. MOTOR DE LEITURA FASTA / GZIP NATIVO E À PROVA DE ERROS
# =====================================================================

def carregar_fasta_para_dict(caminho):
    genoma = {}
    if not os.path.exists(caminho): 
        return genoma
    
    try:
        if caminho.endswith('.gz'):
            handle = gzip.open(caminho, 'rt', encoding='utf-8', errors='ignore')
        else:
            handle = open(caminho, 'r', encoding='utf-8', errors='ignore')
            
        chave_atual = None
        seq_acumu = []
        for linha in handle:
            linha = linha.strip()
            if not linha: 
                continue
            if linha.startswith('>'):
                if chave_atual: 
                    genoma[chave_atual] = "".join(seq_acumu).upper()
                chave_atual = linha[1:].strip()
                seq_acumu = []
            else:
                seq_acumu.append(linha)
                
        if chave_atual: 
            genoma[chave_atual] = "".join(seq_acumu).upper()
        handle.close()
    except Exception as e:
        print(f"Erro ao ler FASTA/GZ ({caminho}): {e}")
        
    return genoma

def encontrar_chave_fasta(genoma, alvo):
    if not genoma: 
        return None
    if not alvo: 
        return list(genoma.keys())[0] if len(genoma) == 1 else None
        
    alvo_str = str(alvo).strip()
    if alvo_str in genoma.keys(): 
        return alvo_str
    
    for k in genoma.keys():
        if k.lower() == alvo_str.lower() or alvo_str.lower() in k.lower() or k.lower() in alvo_str.lower(): 
            return k
            
    for k in genoma.keys():
        if k.replace("_", " ").split()[0].lower() == alvo_str.lower(): 
            return k
            
    return list(genoma.keys())[0] if len(genoma) == 1 else None

def procurar_base_com_tolerancia(seq_cromossoma, idx_teorico, ref_vcf, alt_vcf, tolerancia=5):
    tamanho_crom = len(seq_cromossoma)
    if not (0 <= idx_teorico < tamanho_crom): 
        return None, ref_vcf, f"Fora do limite ({idx_teorico}/{tamanho_crom} bp)"
        
    base_atual = str(seq_cromossoma[idx_teorico]).upper()
    if base_atual == ref_vcf.upper(): 
        return idx_teorico, ref_vcf.upper(), alt_vcf.upper()
    elif base_atual == complemento(ref_vcf).upper(): 
        return idx_teorico, complemento(ref_vcf).upper(), complemento(alt_vcf).upper()
            
    for desvio in range(1, tolerancia + 1):
        for sinal in [1, -1]:
            idx_teste = idx_teorico + (desvio * sinal)
            if 0 <= idx_teste < tamanho_crom:
                base_teste = str(seq_cromossoma[idx_teste]).upper()
                if base_teste == ref_vcf.upper(): 
                    return idx_teste, ref_vcf.upper(), alt_vcf.upper()
                elif base_teste == complemento(ref_vcf).upper(): 
                    return idx_teste, complemento(ref_vcf).upper(), complemento(alt_vcf).upper()
                    
    return None, ref_vcf, str(seq_cromossoma[idx_teorico]).upper()

# =====================================================================
# 2. CÓDIGOS DEGENERADOS IUPAC & CONSENSOS EXAUSTIVOS (15 FTs YEASTRACT)
# =====================================================================

IUPAC_TABELA = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
}

CONSENSOS_YEASTRACT_CDR1 = {
    "Amt1":  ["WTWHNGCTGW"],
    "Ap1":   ["TTACAAA", "TTAGTAA", "VTTACWAAB", "VDTASTAA", "GWWSTASWWCT"],
    "Ap5":   ["TTACKAA", "ATTACTA", "CASACAGWAS", "AYAGTGGTTCAT"],
    "Atf2":  ["ACACCC"],
    "CAGL0B03421g": ["CCCCTCC"],
    "Cst6":  ["TGACGTCA"],
    "Hap5":  ["CCAAT"],
    "Msn2":  ["CCYYCCYYYW"],
    "Msn4":  ["CCYYCCYYYW"],
    "Pdr1":  ["BCCRYYRGD", "TCCRYGGA", "TCCACGGA", "HYCCGTGGR", "TCCACGGAA", "RTTCCGTGGAA"],
    "Pho4":  ["CCCGTG", "CACGTG"], 
    "Rpn4":  ["TTGCAAA"],
    "Yap3":  ["ATTACGTAAT"],
    "Yap6":  ["AATKACV", "WATKACGW", "ATGACGT", "ATKACGTMAT", "RTKACGTMAY", "GCACAGA", "ATTACGTW", "YCGWCGS"],
    "Yap7":  ["TTACTAA", "WGMTTACTRR", "TTAGTAA", "ATTACTAA", "MTKASTMA", "DATTACTAAD", "GMTKASTMAK", "TTASTAAK", "AGTCATM", "KTACTAAD"]
}

# =====================================================================
# 3. MATRIZES PWM E FUNÇÕES BÁSICAS
# =====================================================================

def complemento(base):
    mapa = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return mapa.get(base, 'N')

def reverso_complemento(sequencia):
    return "".join([complemento(b) for b in reversed(sequencia)])

def sequencia_combina_com_iupac_flexivel(sequencia_kmer, padrao_iupac, max_mismatches=0):
    if len(sequencia_kmer) != len(padrao_iupac): 
        return False
        
    mismatches = 0
    for base_seq, base_iupac in zip(sequencia_kmer.upper(), padrao_iupac.upper()):
        if base_seq not in IUPAC_TABELA.get(base_iupac, {'N'}):
            mismatches += 1
            if mismatches > max_mismatches: 
                return False
    return True

def calcular_coordenadas_motivo(idx_genomico_start, tamanho_mot, tss, strand, offset=0, is_local_fasta=False, tamanho_promotor=1000, idx_atg_ref=None, nome_ft=""):
    if is_local_fasta:
        ancora = idx_atg_ref if idx_atg_ref is not None else tamanho_promotor
        pos_5prime = (idx_genomico_start - ancora) + offset
    else:
        idx_genomico_end = idx_genomico_start + tamanho_mot - 1
        if strand == '+': 
            pos_5prime = (idx_genomico_start - tss) + offset
        else: 
            pos_5prime = (tss - idx_genomico_end) + offset
            
    largura_visual_padrao = 12
    pos_3prime = pos_5prime + largura_visual_padrao
    p_min = min(pos_5prime, pos_3prime)
    p_max = max(pos_5prime, pos_3prime)
    
    return f"[{p_min} a {p_max}]", pos_5prime

# =====================================================================
# 4. MAPEAMENTO DE SÍTIOS NATIVOS: MODO RIGOROSO COMPUTACIONAL
# =====================================================================

def reordenar_e_limpar_nativos(lista_motivos):
    df_nativos = pd.DataFrame(lista_motivos)
    if not df_nativos.empty:
        df_nativos = df_nativos.drop_duplicates(subset=["FT", "Pos_5p", "Consenso_IUPAC", "Fita"]).sort_values(by="Pos_5p")
        return df_nativos.to_dict('records')
    return []

def mapear_motivos_nativos_referencia(coord_regiao, tss, caminho_fasta="genoma_referencia.fasta", gene_alvo="", offset_utr=0, tamanho_promotor=1000):
    largura_visual_padrao = 12
    
    genoma = carregar_fasta_para_dict(caminho_fasta)
    if not genoma: 
        return []
        
    cromossoma_alvo = encontrar_chave_fasta(genoma, gene_alvo)
    if not cromossoma_alvo: 
        cromossoma_alvo = encontrar_chave_fasta(genoma, str(coord_regiao['chrom_original']))
    if not cromossoma_alvo and len(genoma.keys()) == 1: 
        cromossoma_alvo = list(genoma.keys())[0]
    if not cromossoma_alvo: 
        return []

    start_gen = int(coord_regiao['start'])
    end_gen = int(coord_regiao['end'])
    strand = coord_regiao['strand']
    tamanho_crom = len(genoma[cromossoma_alvo])
    is_local_fasta = tamanho_crom <= 30000
    margem_buffer = 15
    
    if is_local_fasta:
        seq_promotor = str(genoma[cromossoma_alvo][:]).upper()
        idx_atg_local = tamanho_promotor if tamanho_crom > (tamanho_promotor + 200) else tamanho_crom
    else:
        s_buff = max(1, start_gen - margem_buffer)
        e_buff = min(tamanho_crom, end_gen + margem_buffer)
        seq_promotor = str(genoma[cromossoma_alvo][s_buff-1:e_buff]).upper()
        idx_atg_local = None
        
    motivos_nativos_encontrados = []
    
    for ft_clean, consensos in CONSENSOS_YEASTRACT_CDR1.items():
        nome_ft_p = f"{ft_clean}p" if not ft_clean.endswith("g") else ft_clean
        sitios_do_ft = []
        
        for motivo_iupac in consensos:
            tamanho_mot = len(motivo_iupac)
            
            for i in range(len(seq_promotor) - tamanho_mot + 1):
                kmer_dir = seq_promotor[i:i+tamanho_mot]
                kmer_rev = reverso_complemento(kmer_dir)
                
                combina_dir = sequencia_combina_com_iupac_flexivel(kmer_dir, motivo_iupac, max_mismatches=0)
                combina_rev = sequencia_combina_com_iupac_flexivel(kmer_rev, motivo_iupac, max_mismatches=0)
                
                if combina_dir or combina_rev:
                    idx_gen = i if is_local_fasta else (start_gen - 1 - margem_buffer) + i
                    coord_str, pos_5p = calcular_coordenadas_motivo(
                        idx_gen, tamanho_mot, tss, strand, offset_utr, is_local_fasta, tamanho_promotor, idx_atg_local, nome_ft=nome_ft_p
                    )
                    
                    if pos_5p > 0: 
                        continue
                        
                    sitios_do_ft.append({
                        "FT": nome_ft_p, 
                        "FT_Track": nome_ft_p,
                        "Consenso_IUPAC": motivo_iupac, 
                        "Sequencia_FASTA": kmer_dir if combina_dir else kmer_rev,
                        "Fita": 'F' if combina_dir else 'R', 
                        "Pos_5p": pos_5p, 
                        "Pos_3p": pos_5p + largura_visual_padrao,
                        "Coordenada_YEASTRACT": f"{coord_str} (YEASTRACT 5': {pos_5p})"
                    })
                    
        motivos_nativos_encontrados.extend(sitios_do_ft)
                
    return reordenar_e_limpar_nativos(motivos_nativos_encontrados)

# =====================================================================
# 5. MOTOR DE AVALIAÇÃO DE MUTAÇÕES
# =====================================================================

def avaliar_motivo_yeastract_iupac(seq_ref_janela, seq_alt_janela, nome_ft, pos_relativa_snp, idx_alinhado_snp, idx_base_janela, tss, strand, offset=0, is_local_fasta=False, tamanho_promotor=1000, idx_atg_ref=None):
    consensos = CONSENSOS_YEASTRACT_CDR1.get(nome_ft, [])
    if not consensos: 
        return None
        
    for motivo_iupac in consensos:
        tamanho_mot = len(motivo_iupac)
        
        for i in range(len(seq_ref_janela) - tamanho_mot + 1):
            kmer_ref_dir = seq_ref_janela[i:i+tamanho_mot]
            kmer_ref_rev = reverso_complemento(kmer_ref_dir)
            kmer_alt_dir = seq_alt_janela[i:i+tamanho_mot]
            kmer_alt_rev = reverso_complemento(kmer_alt_dir)
            
            ref_e_motivo_dir = sequencia_combina_com_iupac_flexivel(kmer_ref_dir, motivo_iupac, 0)
            ref_e_motivo_rev = sequencia_combina_com_iupac_flexivel(kmer_ref_rev, motivo_iupac, 0)
            alt_e_motivo_dir = sequencia_combina_com_iupac_flexivel(kmer_alt_dir, motivo_iupac, 0)
            alt_e_motivo_rev = sequencia_combina_com_iupac_flexivel(kmer_alt_rev, motivo_iupac, 0)
            
            fita = 'F' if (ref_e_motivo_dir or alt_e_motivo_dir) else 'R'
            idx_real_motivo = idx_base_janela + i
            
            coord_str, pos_5p_nativa = calcular_coordenadas_motivo(
                idx_real_motivo, tamanho_mot, tss, strand, offset, is_local_fasta, tamanho_promotor, idx_atg_ref, nome_ft=nome_ft
            )
            
            coord_completa = f"{coord_str} (Fita {fita}; YEASTRACT 5': {pos_5p_nativa})"
            
            if (ref_e_motivo_dir or ref_e_motivo_rev) and not (alt_e_motivo_dir or alt_e_motivo_rev):
                return f"[{nome_ft}] 🔴 Abolição de Consenso YEASTRACT ({motivo_iupac}) em {coord_completa}", pos_5p_nativa, coord_str
            elif not (ref_e_motivo_dir or ref_e_motivo_rev) and (alt_e_motivo_dir or alt_e_motivo_rev):
                return f"[{nome_ft}] 🟢 Ganho de Consenso YEASTRACT De Novo ({motivo_iupac}) em {coord_completa}", pos_5p_nativa, coord_str
                
    return None

def analisar_impacto_regulatorio(df_vcf, coord_regiao, tss, mutacoes_conhecidas, caminho_fasta="genoma_referencia.fasta", tamanho_promotor=1000, gene_alvo="", rigor_pwm=0.85, apenas_exatos=True, offset_utr=0):
    detalhes = []
    genoma = carregar_fasta_para_dict(caminho_fasta)
    fasta_disponivel = bool(genoma)
            
    cromossoma_alvo = None
    if fasta_disponivel:
        cromossoma_alvo = encontrar_chave_fasta(genoma, gene_alvo)
        if not cromossoma_alvo: 
            cromossoma_alvo = encontrar_chave_fasta(genoma, str(coord_regiao['chrom_original']))
        if not cromossoma_alvo and len(genoma.keys()) == 1: 
            cromossoma_alvo = list(genoma.keys())[0]
    
    OFFSET_UTR = offset_utr
    
    for idx, row in df_vcf.iterrows():
        pos_absoluta = int(row['POS'])
        ref_vcf = str(row['REF']).strip()
        alt_vcf = str(row['ALT']).strip()
        qual = float(row.get('QUAL', 0)) if pd.notna(row.get('QUAL')) else 0.0
        strand = coord_regiao['strand']
        
        tamanho_crom_teste = len(genoma[cromossoma_alvo]) if (fasta_disponivel and cromossoma_alvo) else 0
        is_local_fasta = tamanho_crom_teste <= 30000 if fasta_disponivel else False

        if is_local_fasta:
            pos_relativa = (pos_absoluta - tss) + OFFSET_UTR if strand == '+' else (tss - pos_absoluta) + OFFSET_UTR
            idx_atg_ref = tamanho_promotor if tamanho_crom_teste > (tamanho_promotor + 200) else tamanho_crom_teste
            idx_teorico = int(idx_atg_ref + pos_relativa)
            if not (0 <= idx_teorico < idx_atg_ref): 
                continue
        else:
            pos_relativa = (pos_absoluta - tss) + OFFSET_UTR if strand == '+' else (tss - pos_absoluta) + OFFSET_UTR
            idx_teorico = int(pos_absoluta - 1)
            idx_atg_ref = None
            
        tipo_var = "SNP" if len(ref_vcf) == 1 and len(alt_vcf) == 1 else "INDEL"
        efeitos_totais = []
        consensos_totais = []
        
        consenso_yeastract_5p = None
        coord_consenso_interv = None
        houve_impacto_nativo = False
        
        if tipo_var == "SNP":
            if not fasta_disponivel or not cromossoma_alvo:
                efeitos_totais.append("⚠️ Ficheiro FASTA/Alvo não encontrado.")
            else:
                resultado_busca = procurar_base_com_tolerancia(genoma[cromossoma_alvo], idx_teorico, ref_vcf, alt_vcf, tolerancia=5)
                
                if resultado_busca[0] is not None:
                    idx_alinhado, ref_eval, alt_eval = resultado_busca
                    
                    for nome_ft, consenso_lista in CONSENSOS_YEASTRACT_CDR1.items():
                        if consenso_lista:
                            consensos_totais.append(f"{nome_ft}: {consenso_lista[0]}")
                        
                        try:
                            inicio_janela = max(0, idx_alinhado - 12)
                            fim_janela = min(tamanho_crom_teste, idx_alinhado + 13)
                            seq_ref_real = str(genoma[cromossoma_alvo][inicio_janela:fim_janela]).upper()
                            pos_central_janela = idx_alinhado - inicio_janela
                            
                            if len(seq_ref_real) > pos_central_janela and seq_ref_real[pos_central_janela] == ref_eval:
                                seq_alt_real = seq_ref_real[:pos_central_janela] + alt_eval + seq_ref_real[pos_central_janela + 1:]
                                
                                res_iupac = avaliar_motivo_yeastract_iupac(
                                    seq_ref_real, seq_alt_real, nome_ft, pos_relativa, idx_alinhado, 
                                    inicio_janela, tss, strand, OFFSET_UTR,
                                    is_local_fasta=is_local_fasta, tamanho_promotor=tamanho_promotor, idx_atg_ref=idx_atg_ref
                                )
                                
                                if res_iupac:
                                    efeito_iupac, p5_nat, c_str_nat = res_iupac
                                    efeitos_totais.append(efeito_iupac)
                                    
                                    consenso_yeastract_5p = p5_nat
                                    coord_consenso_interv = c_str_nat
                                    houve_impacto_nativo = True
                        except Exception:
                            pass
                else:
                    efeitos_totais.append(f"⚠️ Erro Alinhamento FASTA (VCF/REF={ref_vcf} vs FASTA={resultado_busca[2]})")

        consenso_final = " | ".join(consensos_totais) if consensos_totais else "N/A"
        efeito_final = " || ".join(efeitos_totais) if efeitos_totais else "Sem impacto termodinâmico direto nos FTs avaliados"
        
        if pos_relativa in mutacoes_conhecidas: 
            efeito_final += f" || 🏥 Anotação Clínica: {mutacoes_conhecidas[pos_relativa]}"

        if not houve_impacto_nativo or consenso_yeastract_5p is None:
            texto_yeastract_5p = "N/A (Sem Sítio Nativo Abalado)"
            intervalo_real = "N/A"
        else:
            texto_yeastract_5p = str(consenso_yeastract_5p)
            intervalo_real = str(coord_consenso_interv)

        detalhes.append({
            "ID_Interno": f"var_{pos_absoluta}_{ref_vcf}_{alt_vcf}",
            "Nome_Coluna": f"(POS={pos_absoluta}; QUAL={qual:.2f})",
            "Posição Absoluta": pos_absoluta,
            "Posição Relativa (TSS)": pos_relativa,
            "Consenso_YEASTRACT_5p": texto_yeastract_5p,
            "Intervalo_Consenso": intervalo_real,
            "Type": tipo_var,
            "REF": ref_vcf,
            "ALT": alt_vcf,
            "Consenso": consenso_final,
            "Consequencia": efeito_final
        })
        
    return detalhes, []