# ------------------------------------------------
# [1] CONFIGURAÇÕES INICIAIS - ESTUDO DE SIMULAÇÃO
# ------------------------------------------------

# Setar Diretório de Trabalho
setwd("C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC")

# Pacotes Necessários
library(eha)
library(dplyr)

# Funções Implementadas
source("CODES/MEP & MEPP/ARQUIVOS SUPLEMENTARES/GENERATION_FUNCTIONS.R")
source("CODES/MEP & MEPP/ARQUIVOS SUPLEMENTARES/MODELING_FUNCTIONS.R")

# -----------------------
# [2] VALIDAÇÃO DO MODELO
# -----------------------

# Número de Amostras Retiradas
k <- 1000

# Tamnho de Cada Amostra
n <- 30

# -------------------------------------------------
# [2.1] CONFIGURAÇÕES DE SIMULAÇÃO BASE DOS MODELOS
# -------------------------------------------------

# Setar semente aleatória
set.seed(123456789)

# Parâmetros da Simulação
breaks <- c(0.5, 2)       # Vetor de Partições
rates.true <- c(1.1, 0.8, 0.5) # Vetor Taxas de Falha
cens.rate <- 1            # Taxa de Censura (C ~ Exponencial(1))

# Vetor de Coeficientes
coeffs.true <- c(-0.5, 0.5)

# Matriz Design
x1 <- rbinom(n, size = 1, prob=0.5)
x2 <- rnorm(n)
X <- cbind(x1, x2)

# Criando Lista Vazia
df.samples <- data.frame()

# Criando
df.rls <- tibble()

# -----------------------------------
# [2.2] MODELO EXPONENCIAL POR PARTES
# -----------------------------------

# Setar semente aleatória
set.seed(123456789)

# Iteração de Simulação
for (id in 1:k) {
  # Iniciando o ajuste nulo
  fit.model <- NULL
  
  # Ajuste do Modelo
  while(is.null(fit.model) | any(!is.finite(fit.model$par))) {
    # Gerar Data Frame de Dados
    df <- GEN.IC.MEP.DF(
      n, x.cov = X, coeffs.par = coeffs.true, 
      cuts.points = breaks, rates.par = rates.true, 
      cens.rate = cens.rate
    )
    
    # Variável Resposta; Número de Parâmetros; Chute Inicial
    y <- cbind(df$L, df$R)
    n.par <- length(rates.true) + ncol(X)
    init <- rep(1, n.par)
    
    # Ajuste do Modelo
    fit.model <- tryCatch(
      {
        optim(
          par = init, fn = loglikelihood.MEP, method = "BFGS", hessian = TRUE,
          interval = y, cens = df$cens, X = X, cuts.points = breaks,
          control = list(fnscale = -1)
        )
      },
      error = function(e) return(NULL),
      warning = function(w) invokeRestart("muffleWarning") # ignora warning
    )
  }
  
  # Identificando e Armazenando a k-ésima Amostra
  df.samples <- bind_rows(
    df.samples,
    df %>% mutate(ID_SAMPLE = paste0("AMOSTRA ", id))
  )
  
  # Erros Padrões (Trata o erro de inversão da hessiana com tryCatch)
  SE.hat <- tryCatch(
    sqrt(diag(-solve(fit.model$hessian))),
    # Em caso de erro (matriz singular), retorna NA para os erros padrão
    error = function(e) {rep(NA, n.par)}
  )
  
  # Armazenar Resultados
  df.rls <- bind_rows(
    df.rls,
    tibble(
      Lambda1 = exp(fit.model$par[1]),
      Lambda2 = exp(fit.model$par[2]),
      Lambda3 = exp(fit.model$par[3]),
      Beta1   = fit.model$par[4],
      Beta2   = fit.model$par[5],
      SE_L1   = SE.hat[1],
      SE_L2   = SE.hat[2],
      SE_L3   = SE.hat[3],
      SE_B1   = SE.hat[4],
      SE_B2   = SE.hat[5],
      AIC     = 2 * n.par - 2 * fit.model$value,
      BIC     = n.par * log(nrow(y)) - 2 * fit.model$value,
      ID_SAMPLE = paste0("AMOSTRA ", id)
    )
  )
  
  # Imprimir a iteração
  print(id)
}

# Salvar as Amostra Retiradas
saveRDS(df.samples, file = "CODES/MEP & MEPP/MEP/DADOS SIMULADOS/TAMANHO 30/DATA_SIMULATED_30")

# Salvar Resultados
openxlsx::write.xlsx(df.rls, file = "CODES/MEP & MEPP/MEP/RESULTADOS/TAMANHO 30/SIMULATION_RESULTS_30.xlsx")
