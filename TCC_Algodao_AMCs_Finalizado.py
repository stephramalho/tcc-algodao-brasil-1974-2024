# -*- coding: utf-8 -*-
"""
Da dispersão à concentração: a dinâmica espacial do algodão no Brasil,
1974-2024
================================================================================
TCC USP ESALQ - MBA em Data Science and Analytics
Stéphanni Ramalho de Assis

Script único que reproduz, em sequência, todos os dados, tabelas e figuras
do TCC a partir dos microdados públicos da Pesquisa Agrícola Municipal
[PAM/IBGE]:

    ETAPA 1 - Construção do painel em nível de AMC (1974-2024)
    ETAPA 2 - Quociente Locacional [QL] por corte temporal (Tabela 1)
    ETAPA 3 - Série histórica nacional, quebras estruturais (Figuras 1 e 2)
    ETAPA 4 - Índice de Moran Global e LISA (Tabela 2, Figuras 3 e 4)
    ETAPA 5 - Tipologia territorial (Tabela 3, Figura 5)
    ETAPA 6 - Distribuição por UF das AMCs especializadas (Tabela 4)

ENTRADAS (em PASTA_DADOS):
  - PAM_painel_1974_2024.csv          (painel municipal bruto da PAM/IBGE)
  - br_ibge_amc_municipio_de_para.csv (correspondência município -> AMC)

SAÍDAS (em PASTA_DADOS):
  - PAM_painel_AMC_1974_2024.csv
  - QL_AMCs_cortes.csv, Tabela_1_QL_AMCs.csv
  - Quebras_Estruturais_BIC.csv, Figura_1_Serie_Historica_Nacional.png,
    Figura_2_Quebras_Estruturais.png
  - Tabela_2_Autocorrelacao_Espacial_AMCs.csv, Classificacao_LISA_por_AMC.csv,
    Figura_3_Moran_Scatterplot.png, Figura_4_Mapas_LISA.png
  - Tabela_3_Tipologia_AMCs.csv, Classificacao_Tipologia_por_AMC.csv,
    Figura_5_Tipologia.png
  - Tabela_4_UF_especializadas.csv

Convenção de quadrante LISA usada (esda.Moran_Local.q):
    1 = Alto-Alto (HH) | 2 = Baixo-Alto (LH) | 3 = Baixo-Baixo (LL) | 4 = Alto-Baixo (HL)

Reprodutibilidade: a Etapa 4 usa teste de permutação aleatória (I de Moran e
LISA); a semente SEMENTE_ALEATORIA está fixada para que o resultado seja
idêntico a cada execução.
"""

import os
import sys
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

try:
    import ruptures as rpt
except ImportError:
    print('Instalando pacote "ruptures" (busca de quebras estruturais)...')
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'ruptures'])
    import ruptures as rpt

import esda
import geobr
import libpysal

# =============================================================================
# CONFIGURAÇÃO GERAL
# =============================================================================
PASTA_DADOS = r'D:\USP ESALQ - Data Science & Analytics\TCC\Dados preliminares\DADOS CORRETOS FINAIS ATUAIS'

CORTES = [1985, 1996, 2005, 2024]
LETRAS_PAINEL = ['A', 'B', 'C', 'D']
SEMENTE_ALEATORIA = 12345

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 11

print('=' * 78)
print('DA DISPERSÃO À CONCENTRAÇÃO: A DINÂMICA ESPACIAL DO ALGODÃO NO BRASIL')
print('=' * 78)

# =============================================================================
# ETAPA 1 - CONSTRUÇÃO DO PAINEL EM NÍVEL DE AMC (1974-2024)
# -----------------------------------------------------------------------------
# Agrega o painel municipal da PAM/IBGE para Áreas Mínimas Comparáveis [AMCs],
# usando a harmonização territorial 1970-2010 (IPEA/IBGE), garantindo
# comparabilidade longitudinal da malha territorial ao longo da série.
# =============================================================================
print('\n### ETAPA 1: Construção do painel em nível de AMC ###')

ARQUIVO_PAM_MUNICIPAL = os.path.join(PASTA_DADOS, 'PAM_painel_1974_2024.csv')
ARQUIVO_DE_PARA = os.path.join(PASTA_DADOS, 'br_ibge_amc_municipio_de_para.csv')
ARQUIVO_PAINEL_AMC = os.path.join(PASTA_DADOS, 'PAM_painel_AMC_1974_2024.csv')
ANO_DE_ORIGEM = '1970'
ANO_DE_DESTINO = '2010'

print('1.1. Lendo tabela de correspondência município -> AMC...')
df_de_para = pd.read_csv(ARQUIVO_DE_PARA, dtype=str)
de_para_amc = df_de_para[
    (df_de_para['ano_de'] == ANO_DE_ORIGEM)
    & (df_de_para['ano_para'] == ANO_DE_DESTINO)
][['id_municipio', 'id_amc']].copy()
de_para_amc['id_municipio'] = de_para_amc['id_municipio'].str.zfill(7)
de_para_amc['id_amc'] = de_para_amc['id_amc'].astype(str)
print(
    f'      Mapeamento pronto: {len(de_para_amc)} municípios em '
    f'{de_para_amc["id_amc"].nunique()} AMCs.'
)

print('1.2. Lendo painel municipal da PAM/IBGE...')
df_pam = pd.read_csv(
    ARQUIVO_PAM_MUNICIPAL, dtype={'cod_municipio': str}, low_memory=False
)
df_pam['cod_municipio'] = df_pam['cod_municipio'].str.zfill(7)
df_pam_amc = df_pam.merge(
    de_para_amc, left_on='cod_municipio', right_on='id_municipio', how='left'
)
n_sem_amc = df_pam_amc['id_amc'].isna().sum()
if n_sem_amc > 0:
    print(f'      Aviso: {n_sem_amc} registros municipais sem AMC correspondente.')

print('1.3. Agregando variáveis por AMC e ano...')
df_amc = df_pam_amc.groupby(['id_amc', 'ano'], as_index=False).agg({
    'area_colhida_ha': 'sum',
    'total_area_colhida_ha': 'sum',
    'qtd_produzida_t': 'sum',
    'total_qtd_produzida_t': 'sum',
})
df_amc['rendimento_kg_ha'] = np.where(
    df_amc['area_colhida_ha'] > 0,
    (df_amc['qtd_produzida_t'] * 1000.0) / df_amc['area_colhida_ha'],
    0.0,
)
df_amc.to_csv(ARQUIVO_PAINEL_AMC, index=False, encoding='utf-8-sig')
print(
    f'      {df_amc["id_amc"].nunique()} AMCs x {df_amc["ano"].nunique()} anos. '
    f'Arquivo salvo: {ARQUIVO_PAINEL_AMC}'
)

# =============================================================================
# ETAPA 2 - QUOCIENTE LOCACIONAL [QL] POR CORTE TEMPORAL (TABELA 1)
# -----------------------------------------------------------------------------
# QL_i = (e_i / E_i) / (e / E), conforme a equação (1) do TCC, em que e_i é a
# área colhida de algodão na AMC i, E_i é a área colhida agrícola total na
# AMC i, e e e E são os mesmos totais para o Brasil.
# =============================================================================
print('\n### ETAPA 2: Quociente Locacional [QL] e Tabela 1 ###')

ARQUIVO_QL = os.path.join(PASTA_DADOS, 'QL_AMCs_cortes.csv')
ARQUIVO_TABELA1 = os.path.join(PASTA_DADOS, 'Tabela_1_QL_AMCs.csv')

df_ql = df_amc[df_amc['ano'].isin(CORTES)].copy()
totais_nacionais = df_ql.groupby('ano').agg(
    area_algodao_br=('area_colhida_ha', 'sum'),
    area_total_br=('total_area_colhida_ha', 'sum'),
)
df_ql = df_ql.merge(totais_nacionais, on='ano', how='left')
df_ql['part_local'] = np.where(
    df_ql['total_area_colhida_ha'] > 0,
    df_ql['area_colhida_ha'] / df_ql['total_area_colhida_ha'],
    0.0,
)
df_ql['part_nacional'] = df_ql['area_algodao_br'] / df_ql['area_total_br']
df_ql['QL'] = np.where(
    df_ql['part_nacional'] > 0, df_ql['part_local'] / df_ql['part_nacional'], 0.0
)
df_ql[['id_amc', 'ano', 'area_colhida_ha', 'total_area_colhida_ha', 'QL']].to_csv(
    ARQUIVO_QL, index=False, encoding='utf-8-sig'
)

linhas_tabela1 = []
for ano in CORTES:
    sub = df_ql[df_ql['ano'] == ano]
    especializadas = sub[sub['QL'] > 1]
    linhas_tabela1.append({
        'Corte': ano,
        'AMCs_Especializadas_QL_maior_1': len(especializadas),
        'Total_AMCs': len(sub),
        'Participacao_Area_Algodao_Total_pct': round(
            sub['area_colhida_ha'].sum() / sub['total_area_colhida_ha'].sum() * 100, 2
        ),
        'Mediana_QL_Especializadas': round(especializadas['QL'].median(), 2),
    })
df_tabela1 = pd.DataFrame(linhas_tabela1)
df_tabela1.to_csv(ARQUIVO_TABELA1, index=False, encoding='utf-8-sig')

print('\n--- TABELA 1 ---')
print(df_tabela1.to_string(index=False))
print(f'\nArquivos salvos: {ARQUIVO_QL} | {ARQUIVO_TABELA1}')

# =============================================================================
# ETAPA 3 - SÉRIE HISTÓRICA NACIONAL, QUEBRAS ESTRUTURAIS (FIGURAS 1 E 2)
# -----------------------------------------------------------------------------
# As quebras estruturais são identificadas por busca de segmentação ótima
# (programação dinâmica) com seleção do número de quebras pelo Critério de
# Informação Bayesiano [BIC] e segmento mínimo de 15% da série.
#
# Os quatro cortes temporais da pesquisa (1985, 1996, 2005, 2024) NÃO são
# derivados do Bai-Perron: são a estratégia analítica comparativa da
# pesquisa, já usada na Etapa 2 e usada também nas Etapas 4 e 5.
# =============================================================================
print('\n### ETAPA 3: Quebras estruturais (BIC) e Figuras 1 e 2 ###')

CORTES_TCC = CORTES
NUMERO_MAXIMO_QUEBRAS_TESTADAS = 5
TAMANHO_MINIMO_SEGMENTO_PCT = 0.15

serie_nacional = df_amc.groupby('ano', as_index=False).agg(
    area_colhida_ha=('area_colhida_ha', 'sum'),
    qtd_produzida_t=('qtd_produzida_t', 'sum'),
)
serie_nacional['rendimento_kg_ha'] = (
    serie_nacional['qtd_produzida_t'] * 1000.0 / serie_nacional['area_colhida_ha']
)
serie_nacional['area_mil_ha'] = serie_nacional['area_colhida_ha'] / 1000.0
serie_nacional['producao_mil_t'] = serie_nacional['qtd_produzida_t'] / 1000.0

anos = serie_nacional['ano'].values
n_observacoes = len(anos)
tamanho_minimo_segmento = int(np.ceil(TAMANHO_MINIMO_SEGMENTO_PCT * n_observacoes))
print(
    f'      {n_observacoes} observações (1974-2024) | '
    f'segmento mínimo ({TAMANHO_MINIMO_SEGMENTO_PCT:.0%}) = {tamanho_minimo_segmento} anos'
)


def buscar_quebras_por_bic(serie, tamanho_minimo, numero_maximo_quebras):
    """Testa de 0 a numero_maximo_quebras pontos de quebra (modelo de médias
    por segmento) e seleciona o número de quebras que minimiza o BIC."""
    resultados = []
    for numero_quebras in range(numero_maximo_quebras + 1):
        algoritmo = rpt.Dynp(model='l2', min_size=tamanho_minimo, jump=1).fit(serie)
        pontos_quebra = algoritmo.predict(n_bkps=numero_quebras)

        limites_segmento = [0] + pontos_quebra
        soma_residuos_quadrados = sum(
            ((serie[limites_segmento[i]:limites_segmento[i + 1]]
              - serie[limites_segmento[i]:limites_segmento[i + 1]].mean()) ** 2).sum()
            for i in range(len(limites_segmento) - 1)
        )
        numero_parametros = numero_quebras + (numero_quebras + 1)
        bic = (
            n_observacoes * np.log(soma_residuos_quadrados / n_observacoes)
            + numero_parametros * np.log(n_observacoes)
        )
        anos_quebra = [int(anos[indice - 1]) for indice in pontos_quebra[:-1]]
        resultados.append({
            'numero_quebras': numero_quebras,
            'anos_quebra': anos_quebra,
            'bic': bic,
        })

    melhor = min(resultados, key=lambda r: r['bic'])
    return resultados, melhor


resultados_bic = []
quebras_selecionadas = {}
for nome_variavel, nome_exibicao in [
    ('area_colhida_ha', 'Área colhida'),
    ('rendimento_kg_ha', 'Rendimento médio'),
]:
    serie = serie_nacional[nome_variavel].values.astype(float)
    tabela_bic, melhor = buscar_quebras_por_bic(
        serie, tamanho_minimo_segmento, NUMERO_MAXIMO_QUEBRAS_TESTADAS
    )
    print(f'\n   --- {nome_exibicao} ---')
    for linha in tabela_bic:
        marcador = '  <-- selecionado (menor BIC)' if linha is melhor else ''
        print(
            f"   {linha['numero_quebras']} quebra(s): {linha['anos_quebra']} "
            f"| BIC = {linha['bic']:.2f}{marcador}"
        )
        resultados_bic.append({
            'variavel': nome_exibicao,
            'numero_quebras': linha['numero_quebras'],
            'anos_quebra': linha['anos_quebra'],
            'bic': round(linha['bic'], 2),
            'selecionado': linha is melhor,
        })
    quebras_selecionadas[nome_variavel] = melhor['anos_quebra']

df_bic = pd.DataFrame(resultados_bic)
caminho_bic = os.path.join(PASTA_DADOS, 'Quebras_Estruturais_BIC.csv')
df_bic.to_csv(caminho_bic, index=False, encoding='utf-8-sig')

quebras_area = quebras_selecionadas['area_colhida_ha']
quebras_rendimento = quebras_selecionadas['rendimento_kg_ha']
print(f'\n   Quebras finais - área colhida: {quebras_area}')
print(f'   Quebras finais - rendimento médio: {quebras_rendimento}')

# --- Figura 1: série histórica nacional (3 painéis, sem quebras) ---
fig1, eixos1 = plt.subplots(3, 1, figsize=(6.3, 8.5), sharex=True)
paineis_figura1 = [
    (eixos1[0], serie_nacional['area_mil_ha'], 'Área colhida (mil ha)', '#1f77b4', 'A'),
    (eixos1[1], serie_nacional['producao_mil_t'], 'Quantidade produzida (mil t)', '#d62728', 'B'),
    (eixos1[2], serie_nacional['rendimento_kg_ha'], r'Rendimento médio (kg ha$^{-1}$)', '#2ca02c', 'C'),
]
for eixo, serie, rotulo_y, cor, letra in paineis_figura1:
    eixo.plot(anos, serie, color=cor, linewidth=1.5)
    eixo.set_ylabel(rotulo_y, fontsize=11, color='black')
    for lado in ('top', 'right'):
        eixo.spines[lado].set_visible(False)
    for lado in ('left', 'bottom'):
        eixo.spines[lado].set_color('black')
        eixo.spines[lado].set_linewidth(1.5)
    eixo.tick_params(axis='both', colors='black', labelsize=11)
    for ano_corte in CORTES_TCC:
        eixo.axvline(x=ano_corte, color='black', linestyle=':', linewidth=1.0, alpha=0.6)
    eixo.text(0.01, 0.96, letra, transform=eixo.transAxes, fontsize=13,
               fontweight='bold', color='black', va='top', ha='left')
for ano_corte in CORTES_TCC:
    eixos1[0].text(ano_corte, 1.04, str(ano_corte), transform=eixos1[0].get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=10, fontweight='bold', color='black')
eixos1[-1].set_xlabel('Ano', fontsize=11, color='black')
fig1.tight_layout()
caminho_figura1 = os.path.join(PASTA_DADOS, 'Figura_1_Serie_Historica_Nacional.png')
fig1.savefig(caminho_figura1, dpi=300, bbox_inches='tight')

# --- Figura 2: quebras estruturais (2 painéis) ---
fig2, eixos2 = plt.subplots(2, 1, figsize=(6.3, 6.5), sharex=True)
cores_segmento = ['#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#1b9e77']
paineis_figura2 = [
    (eixos2[0], serie_nacional['area_mil_ha'], quebras_area, 'Área colhida (mil ha)', '#1f77b4', 'A'),
    (eixos2[1], serie_nacional['rendimento_kg_ha'], quebras_rendimento, r'Rendimento médio (kg ha$^{-1}$)', '#2ca02c', 'B'),
]
for eixo, serie, quebras, rotulo_y, cor, letra in paineis_figura2:
    eixo.plot(anos, serie, color=cor, linewidth=1.5)
    eixo.set_ylabel(rotulo_y, fontsize=11, color='black')
    for lado in ('top', 'right'):
        eixo.spines[lado].set_visible(False)
    for lado in ('left', 'bottom'):
        eixo.spines[lado].set_color('black')
        eixo.spines[lado].set_linewidth(1.5)
    eixo.tick_params(axis='both', colors='black', labelsize=11)
    pontos_segmento = [1974] + quebras + [2024]
    for indice in range(len(pontos_segmento) - 1):
        inicio, fim = pontos_segmento[indice], pontos_segmento[indice + 1]
        mascara = (anos >= inicio) & (anos <= fim)
        media_segmento = serie.values[mascara].mean()
        eixo.hlines(y=media_segmento, xmin=inicio, xmax=fim,
                    colors=cores_segmento[indice % len(cores_segmento)],
                    linestyles='--', linewidth=1.8)
    for ano_quebra in quebras:
        eixo.axvline(x=ano_quebra, color='black', linestyle=':', linewidth=1.0, alpha=0.6)
        eixo.text(ano_quebra, 1.04, str(ano_quebra), transform=eixo.get_xaxis_transform(),
                   ha='center', va='bottom', fontsize=9, fontweight='bold', color='black')
    eixo.text(0.01, 0.96, letra, transform=eixo.transAxes, fontsize=13,
               fontweight='bold', color='black', va='top', ha='left')
eixos2[-1].set_xlabel('Ano', fontsize=11, color='black')
fig2.tight_layout()
caminho_figura2 = os.path.join(PASTA_DADOS, 'Figura_2_Quebras_Estruturais.png')
fig2.savefig(caminho_figura2, dpi=300, bbox_inches='tight')

print(f'\nArquivos salvos: {caminho_bic} | {caminho_figura1} | {caminho_figura2}')

# =============================================================================
# ETAPA 4 - ÍNDICE DE MORAN GLOBAL, LISA, TABELA 2 E FIGURAS 3 E 4
# -----------------------------------------------------------------------------
# Restringe, em cada corte, o universo espacial às AMCs com registro de
# produção (QL > 0): AMCs que nunca produziram algodão não entram na matriz
# de vizinhança como observação de "baixa especialização". A matriz de
# contiguidade Rainha [Queen], padronizada por linha, é reconstruída a cada
# corte apenas entre as AMCs produtoras.
#
# Ilhas (AMC produtora sem nenhuma vizinha produtora) são sempre tratadas
# como "Não Significativo", independentemente do p-valor de permutação, pois
# a estatística local de uma ilha é degenerada (sempre igual a zero).
# =============================================================================
print('\n### ETAPA 4: Moran Global, LISA, Tabela 2 e Figuras 3 e 4 ###')

ROTULOS_LISA = {
    0: 'Não Significativo',
    1: 'Alto-Alto (Cluster)',
    2: 'Baixo-Alto',
    3: 'Baixo-Baixo',
    4: 'Alto-Baixo',
}
CORES_LISA = {0: '#bdbdbd', 1: '#d7191c', 2: '#abd9e9', 3: '#2b83ba', 4: '#fdae61'}
COR_SEM_PRODUCAO = '#e8e8e8'

df_ql_lido = pd.read_csv(ARQUIVO_QL, dtype={'id_amc': str})
df_ql_lido['id_amc'] = df_ql_lido['id_amc'].str.split('.').str[0].str.strip()

print('4.1. Baixando malha oficial de AMCs (geobr, 1970-2010)...')
gdf_amc = geobr.read_comparable_areas(start_year=1970, end_year=2010, verbose=False)
gdf_amc['id_amc'] = gdf_amc['code_amc'].astype(str).str.split('.').str[0].str.strip()
gdf_amc['geometry'] = gdf_amc['geometry'].make_valid()
print(f'      Malha carregada: {len(gdf_amc)} AMCs.')

resultados_tabela2 = []
registros_classificacao_lisa = []
legenda_figura3 = []

fig3, eixos3 = plt.subplots(2, 2, figsize=(12, 10))
eixos3 = eixos3.flatten()
fig4, eixos4 = plt.subplots(2, 2, figsize=(14, 12))
eixos4 = eixos4.flatten()

for indice_corte, ano in enumerate(CORTES):
    print(f'\n--- Corte {ano} ---')

    sub_ql = df_ql_lido[df_ql_lido['ano'] == ano][['id_amc', 'QL']].copy()
    gdf_ano = gdf_amc.merge(sub_ql, on='id_amc', how='left')
    gdf_ano['QL'] = gdf_ano['QL'].fillna(0.0)

    produtoras = gdf_ano[gdf_ano['QL'] > 0].copy().reset_index(drop=True)
    n_produtoras = len(produtoras)
    print(f'AMCs produtoras (QL > 0): {n_produtoras} de {len(gdf_ano)}')

    w_prod = libpysal.weights.Queen.from_dataframe(
        produtoras, use_index=False, ids='id_amc', silence_warnings=True
    )
    ilhas = set(w_prod.islands)
    print(f'Ilhas (produtora sem nenhuma vizinha produtora): {len(ilhas)}')
    w_prod.transform = 'R'

    y = produtoras['QL'].values

    np.random.seed(SEMENTE_ALEATORIA)
    moran_global = esda.Moran(y, w_prod, permutations=999)
    np.random.seed(SEMENTE_ALEATORIA)
    moran_local = esda.Moran_Local(y, w_prod, permutations=999, seed=SEMENTE_ALEATORIA)

    significativo = moran_local.p_sim < 0.05
    quadrante = moran_local.q
    ids_produtoras = produtoras['id_amc'].tolist()
    mascara_ilha = np.array([amc_id in ilhas for amc_id in ids_produtoras])
    significativo = significativo & ~mascara_ilha

    alto_alto = int(np.sum(significativo & (quadrante == 1)))
    baixo_alto = int(np.sum(significativo & (quadrante == 2)))
    baixo_baixo = int(np.sum(significativo & (quadrante == 3)))
    alto_baixo = int(np.sum(significativo & (quadrante == 4)))

    print(f'I de Moran = {moran_global.I:.4f} | p-valor = {moran_global.p_sim:.4f}')
    print(
        f'Alto-Alto={alto_alto} | Baixo-Alto={baixo_alto} | '
        f'Baixo-Baixo={baixo_baixo} | Alto-Baixo={alto_baixo}'
    )

    resultados_tabela2.append({
        'Ano': ano,
        'N_AMCs_Produtoras': n_produtoras,
        'N_Ilhas': len(ilhas),
        'I_Moran_Global': round(moran_global.I, 4),
        'p_value': moran_global.p_sim,
        'Clusters_Alto_Alto': alto_alto,
        'Clusters_Baixo_Alto': baixo_alto,
        'Clusters_Baixo_Baixo': baixo_baixo,
        'Clusters_Alto_Baixo': alto_baixo,
    })

    categorias = np.zeros(n_produtoras, dtype=int)
    categorias[significativo] = quadrante[significativo]
    produtoras['lisa_cat'] = categorias
    produtoras['lisa_rotulo'] = produtoras['lisa_cat'].map(ROTULOS_LISA)
    registros_classificacao_lisa.append(
        produtoras[['id_amc', 'QL', 'lisa_cat', 'lisa_rotulo']].assign(ano=ano)
    )

    letra = LETRAS_PAINEL[indice_corte]
    legenda_figura3.append(
        f'{letra}: {ano} (I = {moran_global.I:.3f}, p = {moran_global.p_sim:.3f})'
    )

    # Figura 3: scatterplot de Moran
    y_padronizado = (y - y.mean()) / y.std()
    lag_espacial = libpysal.weights.spatial_lag.lag_spatial(w_prod, y_padronizado)
    cores_pontos = []
    for sig, q in zip(significativo, quadrante):
        if not sig:
            cores_pontos.append('#d9d9d9')
        elif q == 1:
            cores_pontos.append('#de2d26')
        elif q == 3:
            cores_pontos.append('#3182bd')
        else:
            cores_pontos.append('#fd8d3c')
    eixo3 = eixos3[indice_corte]
    eixo3.scatter(y_padronizado, lag_espacial, c=cores_pontos, alpha=0.6, s=25, edgecolor='none')
    coef_angular, intercepto = np.polyfit(y_padronizado, lag_espacial, 1)
    x_reta = np.linspace(y_padronizado.min(), y_padronizado.max(), 100)
    eixo3.plot(x_reta, intercepto + coef_angular * x_reta, color='darkred', linewidth=1.5)
    eixo3.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    eixo3.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    for lado in ('top', 'right'):
        eixo3.spines[lado].set_visible(False)
    for lado in ('left', 'bottom'):
        eixo3.spines[lado].set_color('black')
        eixo3.spines[lado].set_linewidth(1.5)
    eixo3.tick_params(axis='both', colors='black', labelsize=11)
    eixo3.set_xlabel('QL padronizado (AMCs com registro de produção)', fontsize=11, color='black')
    eixo3.set_ylabel('Lag espacial do QL padronizado', fontsize=11, color='black')
    eixo3.text(0.02, 0.97, letra, transform=eixo3.transAxes, fontsize=14,
                fontweight='bold', color='black', va='top', ha='left')

    # Figura 4: mapa LISA
    eixo4 = eixos4[indice_corte]
    nao_produtoras = gdf_ano[gdf_ano['QL'] == 0]
    if not nao_produtoras.empty:
        nao_produtoras.plot(ax=eixo4, color=COR_SEM_PRODUCAO, edgecolor='#aaaaaa', linewidth=0.05)
    for categoria in range(5):
        subconjunto = produtoras[produtoras['lisa_cat'] == categoria]
        if not subconjunto.empty:
            subconjunto.plot(ax=eixo4, color=CORES_LISA[categoria], edgecolor='black', linewidth=0.05)
    eixo4.text(0.02, 0.98, letra, transform=eixo4.transAxes, fontsize=15,
                fontweight='bold', color='black', va='top', ha='left')
    eixo4.axis('off')

df_tabela2 = pd.DataFrame(resultados_tabela2)
caminho_tabela2 = os.path.join(PASTA_DADOS, 'Tabela_2_Autocorrelacao_Espacial_AMCs.csv')
df_tabela2.to_csv(caminho_tabela2, index=False, encoding='utf-8-sig')

df_classificacao_lisa = pd.concat(registros_classificacao_lisa, ignore_index=True)
caminho_classificacao_lisa = os.path.join(PASTA_DADOS, 'Classificacao_LISA_por_AMC.csv')
df_classificacao_lisa.to_csv(caminho_classificacao_lisa, index=False, encoding='utf-8-sig')

fig3.tight_layout()
caminho_figura3 = os.path.join(PASTA_DADOS, 'Figura_3_Moran_Scatterplot.png')
fig3.savefig(caminho_figura3, dpi=300, bbox_inches='tight')

patches_lisa = [mpatches.Patch(color=CORES_LISA[c], label=ROTULOS_LISA[c]) for c in range(5)]
patches_lisa.append(mpatches.Patch(color=COR_SEM_PRODUCAO, label='Sem produção no corte'))
fig4.legend(handles=patches_lisa, loc='lower center', ncol=3, fontsize=11, frameon=True,
            facecolor='white', edgecolor='gray', bbox_to_anchor=(0.5, 0.02))
fig4.tight_layout()
fig4.subplots_adjust(bottom=0.08)
caminho_figura4 = os.path.join(PASTA_DADOS, 'Figura_4_Mapas_LISA.png')
fig4.savefig(caminho_figura4, dpi=300, bbox_inches='tight')

print('\n--- TABELA 2 ---')
print(df_tabela2.to_string(index=False))
print(
    f'\nArquivos salvos: {caminho_tabela2} | {caminho_classificacao_lisa} | '
    f'{caminho_figura3} | {caminho_figura4}'
)
print('\nLegenda da Figura 3 (usar na legenda externa do documento):')
for linha in legenda_figura3:
    print('  ' + linha)

# =============================================================================
# ETAPA 5 - TIPOLOGIA TERRITORIAL (TABELA 3, FIGURA 5)
# -----------------------------------------------------------------------------
# Classifica as AMCs em uma matriz de quadrantes 2x2 cruzando o QL e o
# rendimento médio, entre as AMCs formalmente especializadas (QL > 1). Os
# limiares alto/baixo QL e alto/baixo rendimento são as respectivas medianas
# entre as AMCs especializadas, calculadas em cada corte temporal.
# =============================================================================
print('\n### ETAPA 5: Tipologia territorial, Tabela 3 e Figura 5 ###')

CORES_TIPOLOGIA = {
    'Alta espec. / Alta produt.': '#d73027',
    'Alta espec. / Baixa produt.': '#fee08b',
    'Baixa espec. / Alta produt.': '#91bfdb',
    'Baixa espec. / Baixa produt.': '#4575b4',
    'Não especializado': '#bdbdbd',
}


def classificar_tipologia(linha, mediana_ql, mediana_rendimento):
    if linha['ql'] <= 1:
        return 'Não especializado'
    if linha['ql'] >= mediana_ql and linha['rendimento_kg_ha'] >= mediana_rendimento:
        return 'Alta espec. / Alta produt.'
    if linha['ql'] >= mediana_ql:
        return 'Alta espec. / Baixa produt.'
    if linha['rendimento_kg_ha'] >= mediana_rendimento:
        return 'Baixa espec. / Alta produt.'
    return 'Baixa espec. / Baixa produt.'


linhas_tabela3 = []
registros_classificacao_tipologia = []
legenda_figura5 = []

fig5, eixos5 = plt.subplots(2, 2, figsize=(13, 13))
eixos5 = eixos5.flatten()

for indice_corte, ano in enumerate(CORTES):
    df_ano = df_amc[df_amc['ano'] == ano].copy()

    area_algodao_br = df_ano['area_colhida_ha'].sum()
    area_total_br = df_ano['total_area_colhida_ha'].sum()
    df_ano['part_local'] = np.where(
        df_ano['total_area_colhida_ha'] > 0,
        df_ano['area_colhida_ha'] / df_ano['total_area_colhida_ha'], 0.0,
    )
    part_nacional = area_algodao_br / area_total_br if area_total_br > 0 else 0.0
    df_ano['ql'] = np.where(part_nacional > 0, df_ano['part_local'] / part_nacional, 0.0)

    especializadas = df_ano[df_ano['ql'] > 1]
    mediana_ql = especializadas['ql'].median() if not especializadas.empty else 0.0
    mediana_rendimento = (
        especializadas['rendimento_kg_ha'].median() if not especializadas.empty else 0.0
    )

    df_ano['tipologia'] = df_ano.apply(
        classificar_tipologia, axis=1, mediana_ql=mediana_ql, mediana_rendimento=mediana_rendimento
    )

    contagem = df_ano['tipologia'].value_counts()
    linhas_tabela3.append({
        'Corte': ano,
        'Alta_espec_Alta_produt': contagem.get('Alta espec. / Alta produt.', 0),
        'Alta_espec_Baixa_produt': contagem.get('Alta espec. / Baixa produt.', 0),
        'Baixa_espec_Alta_produt': contagem.get('Baixa espec. / Alta produt.', 0),
        'Baixa_espec_Baixa_produt': contagem.get('Baixa espec. / Baixa produt.', 0),
        'Nao_especializado': contagem.get('Não especializado', 0),
        'Mediana_rendimento_AMCs_especializadas_kg_ha': round(mediana_rendimento, 0),
    })
    registros_classificacao_tipologia.append(
        df_ano[['id_amc', 'ql', 'rendimento_kg_ha', 'tipologia']].assign(ano=ano)
    )

    gdf_ano_tipologia = gdf_amc.merge(df_ano, on='id_amc', how='left')
    gdf_ano_tipologia['tipologia'] = gdf_ano_tipologia['tipologia'].fillna('Não especializado')

    letra = LETRAS_PAINEL[indice_corte]
    legenda_figura5.append(f'{letra}: {ano}')

    eixo5 = eixos5[indice_corte]
    for categoria, cor in CORES_TIPOLOGIA.items():
        subconjunto = gdf_ano_tipologia[gdf_ano_tipologia['tipologia'] == categoria]
        if not subconjunto.empty:
            cor_borda = '#4d4d4d' if categoria == 'Não especializado' else 'darkgrey'
            subconjunto.plot(ax=eixo5, color=cor, edgecolor=cor_borda, linewidth=0.1)
    eixo5.text(0.02, 0.98, letra, transform=eixo5.transAxes, fontsize=15,
                fontweight='bold', color='black', va='top', ha='left')
    eixo5.axis('off')

df_tabela3 = pd.DataFrame(linhas_tabela3)
caminho_tabela3 = os.path.join(PASTA_DADOS, 'Tabela_3_Tipologia_AMCs.csv')
df_tabela3.to_csv(caminho_tabela3, index=False, encoding='utf-8-sig')

df_classificacao_tipologia = pd.concat(registros_classificacao_tipologia, ignore_index=True)
caminho_classificacao_tipologia = os.path.join(PASTA_DADOS, 'Classificacao_Tipologia_por_AMC.csv')
df_classificacao_tipologia.to_csv(caminho_classificacao_tipologia, index=False, encoding='utf-8-sig')

patches_tipologia = [mpatches.Patch(color=cor, label=cat) for cat, cor in CORES_TIPOLOGIA.items()]
fig5.legend(handles=patches_tipologia, loc='lower center', ncol=3, fontsize=10, frameon=True,
            facecolor='white', edgecolor='gray', bbox_to_anchor=(0.5, 0.02))
fig5.tight_layout()
fig5.subplots_adjust(bottom=0.08)
caminho_figura5 = os.path.join(PASTA_DADOS, 'Figura_5_Tipologia.png')
fig5.savefig(caminho_figura5, dpi=300, bbox_inches='tight')

print('\n--- TABELA 3 ---')
print(df_tabela3.to_string(index=False))
print(f'\nArquivos salvos: {caminho_tabela3} | {caminho_classificacao_tipologia} | {caminho_figura5}')
print('\nLegenda da Figura 5 (usar na legenda externa do documento):')
for linha in legenda_figura5:
    print('  ' + linha)


# =============================================================================
# ETAPA 6 - DISTRIBUIÇÃO POR UF DAS AMCs ESPECIALIZADAS (TABELA 4)
# -----------------------------------------------------------------------------
# Como a unidade de análise é a AMC (não o município), qualquer afirmação sobre
# "onde" a especialização se concentra precisa vir de um dado gerado pelo
# próprio código - nunca de conhecimento genérico sobre quais municípios são
# conhecidos por cultivar algodão. Esta etapa identifica a Unidade da Federação
# de cada AMC (a partir do primeiro município que a compõe, listado pelo
# geobr) e conta quantas AMCs especializadas (QL > 1) caem em cada UF, nos
# cortes de 1985 e 2024.
#
# FONTE EXTERNA: a lista de códigos e nomes de município usada para obter a UF
# vem do repositório público "Municípios Brasileiros" (Kelvins, 2024), com
# licença MIT, baixada diretamente do GitHub - não é conhecimento do modelo,
# é um dado versionado e citável.
#
# SAÍDA:
#   - Tabela_4_UF_especializadas.csv
# =============================================================================
print('\n### ETAPA 6: Distribuição por UF das AMCs especializadas (Tabela 4) ###')

URL_MUNICIPIOS = (
    'https://raw.githubusercontent.com/kelvins/municipios-brasileiros/'
    'main/csv/municipios.csv'
)
UF_SIGLA = {
    11: 'RO', 12: 'AC', 13: 'AM', 14: 'RR', 15: 'PA', 16: 'AP', 17: 'TO',
    21: 'MA', 22: 'PI', 23: 'CE', 24: 'RN', 25: 'PB', 26: 'PE', 27: 'AL',
    28: 'SE', 29: 'BA', 31: 'MG', 32: 'ES', 33: 'RJ', 35: 'SP', 41: 'PR',
    42: 'SC', 43: 'RS', 50: 'MS', 51: 'MT', 52: 'GO', 53: 'DF',
}

print('6.1. Baixando a base pública de municípios (Kelvins, 2024)...')
municipios = pd.read_csv(URL_MUNICIPIOS, dtype={'codigo_ibge': str})
mapa_uf_numero = dict(zip(municipios['codigo_ibge'], municipios['codigo_uf']))


def uf_da_amc(lista_codigos_municipio):
    """A UF da AMC é a do primeiro município que a compõe (uma AMC nunca
    atravessa fronteira estadual, então qualquer município da lista serve)."""
    primeiro_codigo = str(lista_codigos_municipio).split(',')[0].strip()
    numero_uf = mapa_uf_numero.get(primeiro_codigo)
    return UF_SIGLA.get(int(numero_uf), '??') if numero_uf else '??'


print('6.2. Cruzando a malha de AMCs com a UF de cada uma...')
gdf_amc_uf = gdf_amc.copy()
gdf_amc_uf['UF'] = gdf_amc_uf['list_code_muni_2010'].apply(uf_da_amc)

df_ql_uf = df_ql.merge(gdf_amc_uf[['id_amc', 'UF']], on='id_amc', how='left')

print('6.3. Contando AMCs especializadas (QL > 1) por UF, em 1985 e 2024...')
linhas_tabela4 = []
for ano in [1985, 2024]:
    especializadas_ano = df_ql_uf[(df_ql_uf['ano'] == ano) & (df_ql_uf['QL'] > 1)]
    contagem_uf = especializadas_ano['UF'].value_counts()
    total_ano = len(especializadas_ano)
    for uf, quantidade in contagem_uf.items():
        linhas_tabela4.append({
            'Ano': ano,
            'UF': uf,
            'AMCs_especializadas': int(quantidade),
            'Total_AMCs_especializadas_no_ano': total_ano,
            'Participacao_pct': round(quantidade / total_ano * 100, 1),
        })

df_tabela4 = pd.DataFrame(linhas_tabela4)
caminho_tabela4 = os.path.join(PASTA_DADOS, 'Tabela_4_UF_especializadas.csv')
df_tabela4.to_csv(caminho_tabela4, index=False, encoding='utf-8-sig')

print('\n--- TABELA 4 (top 5 por ano) ---')
for ano in [1985, 2024]:
    print(f'\n{ano}:')
    sub = df_tabela4[df_tabela4['Ano'] == ano].sort_values(
        'AMCs_especializadas', ascending=False
    ).head(5)
    print(sub[['UF', 'AMCs_especializadas', 'Participacao_pct']].to_string(index=False))

print(f'\nArquivo salvo: {caminho_tabela4}')
print('\n=== ETAPA 6 CONCLUÍDA ===')

print('\n' + '=' * 78)
print('SCRIPT CONCLUÍDO - todas as 6 etapas executadas com sucesso')
print('=' * 78)