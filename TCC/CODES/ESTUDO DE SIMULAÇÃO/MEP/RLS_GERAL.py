# %% [markdown]
# # Avaliação Geral dos Resultados do Estudo de Simulação

import os
import numpy as np
import pandas as pd
from scipy.stats import norm

# ---------------------------------
# Configurações globais
# ---------------------------------
ns = 0.05
z = norm.ppf(1 - ns/2)
sizes = [30, 50, 100, 500, 1000]

# ---------------------------------
# Função para avaliar um parâmetro
# ---------------------------------
def analyze_parameter(df, col_est, col_se, true_value, adj_lim=False):
    """Calcula média, erro-padrão e probabilidade de cobertura de um parâmetro."""
    aval = df[[col_est, col_se]].copy()

    # Intervalo de confiança
    aval["LOWER"] = aval[col_est] - z * aval[col_se]
    aval["UPPER"] = aval[col_est] + z * aval[col_se]

    # Ajuste do limite inferior para parâmetros de taxa
    if adj_lim:
        aval["LOWER"] = np.where(aval["LOWER"] <= 0, 0, aval["LOWER"])

    # Indicadora de cobertura
    aval["INDICADORA"] = np.where(
        (true_value >= aval["LOWER"]) & (true_value <= aval["UPPER"]),
        1, 0
    )

    return (
        aval[col_est].mean(),           # Média do estimador
        aval[col_se].mean(),            # Média do erro-padrão
        aval["INDICADORA"].mean() * 100 # Probabilidade de cobertura (%)
    )

# ---------------------------------
# Função geral para processar cenários
# ---------------------------------
def process_scenario(base_dir, params, sizes, file_pattern):
    """
    Lê resultados de simulação, calcula métricas e retorna DataFrame resumo.
    
    base_dir: caminho base da parametrização
    params: lista de dicionários com info dos parâmetros
    sizes: lista de tamanhos amostrais
    file_pattern: padrão do nome de arquivo, deve conter {n} para tamanho
    """
    tbl_summary = pd.DataFrame(index=[p["name"] for p in params])

    for n in sizes:
        # Caminho do arquivo
        file_path = os.path.join(base_dir, file_pattern.format(n=n))
        df_rls = pd.read_excel(file_path)

        # Calcula métricas para cada parâmetro
        for p in params:
            mean_par, se_par, cp_par = analyze_parameter(
                df=df_rls,
                col_est=p["col"],
                col_se=p["se"],
                true_value=p["true"],
                adj_lim=p["adj_lim"]
            )
            tbl_summary.loc[p["name"], f"MEAN_{n}"] = mean_par
            tbl_summary.loc[p["name"], f"SE_{n}"] = se_par
            tbl_summary.loc[p["name"], f"CP_{n}"] = cp_par

    return tbl_summary

# ---------------------------------
# Lista de cenários
# ---------------------------------
scenarios = [
    # MEP - Partição 1
    {
        "sheet_name": "MEP tau=(0.5, 2)",
        "base_dir": r"C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/ESTUDO DE SIMULAÇÃO/MEP/Partição 1 - tau = (0.5, 2)",
        "params": [
            {"name": "λ1", "col": "Lambda1", "se": "SE_L1", "true": 1.1, "adj_lim": True},
            {"name": "λ2", "col": "Lambda2", "se": "SE_L2", "true": 0.8, "adj_lim": True},
            {"name": "λ3", "col": "Lambda3", "se": "SE_L3", "true": 0.5, "adj_lim": True},
            {"name": "β1", "col": "Beta1",   "se": "SE_B1", "true": -0.5, "adj_lim": False},
            {"name": "β2", "col": "Beta2",   "se": "SE_B2", "true":  0.5, "adj_lim": False},
        ],
        "file_pattern": r"TAMANHO {n}/SIMULATION_RESULTS_{n}.xlsx"
    },
    # MEP - Partição 2
    {
        "sheet_name": "MEP tau=(2, 6)",
        "base_dir": r"C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/ESTUDO DE SIMULAÇÃO/MEP/Partição 2 - tau = (2, 6)",
        "params": [
            {"name": "λ1", "col": "Lambda1", "se": "SE_L1", "true": 0.4, "adj_lim": True},
            {"name": "λ2", "col": "Lambda2", "se": "SE_L2", "true": 0.2, "adj_lim": True},
            {"name": "λ3", "col": "Lambda3", "se": "SE_L3", "true": 0.5, "adj_lim": True},
            {"name": "β1", "col": "Beta1",   "se": "SE_B1", "true": -0.5, "adj_lim": False},
            {"name": "β2", "col": "Beta2",   "se": "SE_B2", "true":  0.5, "adj_lim": False},
        ],
        "file_pattern": r"TAMANHO {n}/SIMULATION_RESULTS_P2_{n}.xlsx"
    },
    # MEPP - Partição 1 - α = 0.8
    {
        "sheet_name": "MEPP tau=(0.5, 2) α=0.8",
        "base_dir": r"C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/ESTUDO DE SIMULAÇÃO/MEPP/Partição 1 - tau = (0.5, 2)",
        "params": [
            {"name": "λ1", "col": "Lambda1", "se": "SE_L1", "true": 1.1, "adj_lim": True},
            {"name": "λ2", "col": "Lambda2", "se": "SE_L2", "true": 0.8, "adj_lim": True},
            {"name": "λ3", "col": "Lambda3", "se": "SE_L3", "true": 0.5, "adj_lim": True},
            {"name": "α",  "col": "Power",   "se": "SE_PW", "true": 0.8, "adj_lim": True},
            {"name": "β1", "col": "Beta1",   "se": "SE_B1", "true": -0.5, "adj_lim": False},
            {"name": "β2", "col": "Beta2",   "se": "SE_B2", "true":  0.5, "adj_lim": False},
        ],
        "file_pattern": r"TAMANHO {n}/SIMULATION_RESULTS_{n}_A1.xlsx"
    },
    # MEPP - Partição 1 - α = 1.2
    {
        "sheet_name": "MEPP tau=(0.5, 2) α=1.2",
        "base_dir": r"C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/ESTUDO DE SIMULAÇÃO/MEPP/Partição 1 - tau = (0.5, 2)",
        "params": [
            {"name": "λ1", "col": "Lambda1", "se": "SE_L1", "true": 1.1, "adj_lim": True},
            {"name": "λ2", "col": "Lambda2", "se": "SE_L2", "true": 0.8, "adj_lim": True},
            {"name": "λ3", "col": "Lambda3", "se": "SE_L3", "true": 0.5, "adj_lim": True},
            {"name": "α",  "col": "Power",   "se": "SE_PW", "true": 1.2, "adj_lim": True},
            {"name": "β1", "col": "Beta1",   "se": "SE_B1", "true": -0.5, "adj_lim": False},
            {"name": "β2", "col": "Beta2",   "se": "SE_B2", "true":  0.5, "adj_lim": False},
        ],
        "file_pattern": r"TAMANHO {n}/SIMULATION_RESULTS_{n}_A2.xlsx"
    }
    # MEPP - Partição 2 - α = 0.8
    {
        "sheet_name": "MEPP tau=(0.5, 2) α=0.8",
        "base_dir": r"C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/ESTUDO DE SIMULAÇÃO/MEPP/Partição 2 - tau = (2, 6)",
        "params": [
            {"name": "λ1", "col": "Lambda1", "se": "SE_L1", "true": 0.4, "adj_lim": True},
            {"name": "λ2", "col": "Lambda2", "se": "SE_L2", "true": 0.2, "adj_lim": True},
            {"name": "λ3", "col": "Lambda3", "se": "SE_L3", "true": 0.5, "adj_lim": True},
            {"name": "α",  "col": "Power",   "se": "SE_PW", "true": 0.8, "adj_lim": True},
            {"name": "β1", "col": "Beta1",   "se": "SE_B1", "true": -0.5, "adj_lim": False},
            {"name": "β2", "col": "Beta2",   "se": "SE_B2", "true":  0.5, "adj_lim": False},
        ],
        "file_pattern": r"TAMANHO {n}/SIMULATION_RESULTS_{n}_P2_A1.xlsx"
    },
    # MEPP - Partição 2 - α = 1.2
    {
        "sheet_name": "MEPP tau=(0.5, 2) α=1.2",
        "base_dir": r"C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/ESTUDO DE SIMULAÇÃO/MEPP/Partição 2 - tau = (2, 6)",
        "params": [
            {"name": "λ1", "col": "Lambda1", "se": "SE_L1", "true": 0.4, "adj_lim": True},
            {"name": "λ2", "col": "Lambda2", "se": "SE_L2", "true": 0.2, "adj_lim": True},
            {"name": "λ3", "col": "Lambda3", "se": "SE_L3", "true": 0.5, "adj_lim": True},
            {"name": "α",  "col": "Power",   "se": "SE_PW", "true": 1.2, "adj_lim": True},
            {"name": "β1", "col": "Beta1",   "se": "SE_B1", "true": -0.5, "adj_lim": False},
            {"name": "β2", "col": "Beta2",   "se": "SE_B2", "true":  0.5, "adj_lim": False},
        ],
        "file_pattern": r"TAMANHO {n}/SIMULATION_RESULTS_{n}_P2_A2.xlsx"
    }
]

# ---------------------------------
# Processa todos os cenários e salva
# ---------------------------------
output_file = r"C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/ESTUDO DE SIMULAÇÃO/RLS_GERAL_SIMULATIONS.xlsx"

with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    for sc in scenarios:
        tbl = process_scenario(sc["base_dir"], sc["params"], sizes, sc["file_pattern"])
        tbl.to_excel(writer, sheet_name=sc["sheet_name"])

print(f"Arquivo salvo em: {output_file}")
