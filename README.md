# Da dispersão à concentração: a dinâmica espacial do algodão no Brasil, 1974–2024

Código de suporte ao TCC do MBA em Data Science and Analytics (USP/ESALQ),
de Stéphanni Ramalho de Assis. Reproduz integralmente os dados, tabelas e
figuras do artigo a partir dos microdados públicos da Pesquisa Agrícola
Municipal [PAM/IBGE].

## Como rodar

Um único script, `TCC_Algodao_AMCs_Finalizado.py`, executa as seis etapas
da análise em sequência:

1. Construção do painel em nível de AMC (1974-2024)
2. Quociente Locacional [QL] por corte temporal — Tabela 1
3. Distribuição por UF das AMCs especializadas — Tabela 2
4. Série histórica nacional e quebras estruturais (Bai-Perron/BIC) — Figuras 1 e 2
5. Índice de Moran Global e LISA, restritos a AMCs produtoras — Tabela 3 e Figuras 3 e 4
6. Tipologia territorial (QL x rendimento) — Tabela 4 e Figura 5

Basta ajustar `PASTA_DADOS` no início do script e rodar do início ao fim.

## Estrutura do repositório

- `TCC_Algodao_AMCs_Finalizado.py` — script único com as seis etapas
- `Dados/` — tabelas e classificações geradas pelo script (CSVs)
- `Figuras/` — as cinco figuras do TCC, já geradas pelo script

## Dados de entrada necessários

Devem estar na pasta configurada em `PASTA_DADOS`:

- `PAM_painel_1974_2024.csv` — painel municipal bruto da PAM/IBGE (Tabela 1612)
- `br_ibge_amc_municipio_de_para.csv` — correspondência município → AMC (IPEA/IBGE, 1970-2010)

A Etapa 3 também baixa automaticamente, via URL, a base pública de códigos e
nomes de município (Kelvins, 2024) e a malha oficial de AMCs (geobr) — não
precisa baixar nada manualmente.

## Saídas geradas

| Etapa | Tabelas / classificações | Figuras |
|---|---|---|
| 1 | `PAM_painel_AMC_1974_2024.csv` | — |
| 2 | `QL_AMCs_cortes.csv`, `Tabela_1_QL_AMCs.csv` | — |
| 3 | `Tabela_2_UF_Especializadas.csv`, `Classificacao_UF_por_AMC.csv` | — |
| 4 | `Quebras_Estruturais_BIC.csv` | `Figura_1_Serie_Historica_Nacional.png`, `Figura_2_Quebras_Estruturais.png` |
| 5 | `Tabela_3_Autocorrelacao_Espacial_AMCs.csv`, `Classificacao_LISA_por_AMC.csv` | `Figura_3_Moran_Scatterplot.png`, `Figura_4_Mapas_LISA.png` |
| 6 | `Tabela_4_Tipologia_AMCs.csv`, `Classificacao_Tipologia_por_AMC.csv` | `Figura_5_Tipologia.png` |

## Dependências

```
pandas numpy geopandas libpysal esda geobr matplotlib ruptures
```

O script instala `ruptures` automaticamente se necessário. Os demais pacotes
espaciais (`geobr`, `libpysal`, `esda`, `geopandas`) devem estar instalados
previamente.

## Reprodutibilidade

O cálculo do Índice de Moran Global e do LISA (Etapa 5) usa teste de
permutação aleatória. A semente `SEMENTE_ALEATORIA = 12345` está fixada no
script para garantir que o resultado seja idêntico a cada execução.

## Convenção de quadrante LISA

A biblioteca `esda` numera os quadrantes do `Moran_Local` como
`1 = Alto-Alto`, `2 = Baixo-Alto`, `3 = Baixo-Baixo`, `4 = Alto-Baixo`. Essa
é a convenção usada em todo o script e nas legendas das figuras.
