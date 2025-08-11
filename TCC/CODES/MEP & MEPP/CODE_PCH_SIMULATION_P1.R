# ------------------------------------------------
# [1] CONFIGURAÇÕES INICIAIS - ESTUDO DE SIMULAÇÃO
# ------------------------------------------------

# Setar Diretório de Trabalho
setwd("C:/Users/user/Documents/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC")

# Pacotes Necessários
library(eha)

# Funções Implementadas
source("CODES/MEP & MEPP/ARQUIVOS SUPLEMENTARES/GENERATION_FUNCTIONS.R")
source("CODES/MEP & MEPP/ARQUIVOS SUPLEMENTARES/MODELING_FUNCTIONS.R")

# -----------------------
# [2] VALIDAÇÃO DO MODELO
# -----------------------

# Apenas um teste
n <- 1000

# -------------------------------------------------
# [2.1] CONFIGURAÇÕES DE SIMULAÇÃO BASE DOS MODELOS
# -------------------------------------------------

# Setar semente aleatória
set.seed(123456789)

# Parâmetros da Simulação
breaks <- c(0.5, 2)       # Vetor de Partições
rates <- c(1.1, 0.8, 0.5) # Vetor Taxas de Falha
cens.rate <- 1            # Taxa de Censura (C ~ Exponencial(1))

# Vetor de Coeficientes
coeffs <- c(-0.5, 0.5)

# Matriz Design
x1 <- rbinom(n, size = 1, prob=0.5)
x2 <- rnorm(n)
X <- cbind(x1, x2)

# -----------------------------------
# [2.2] MODELO EXPONENCIAL POR PARTES
# -----------------------------------

# Gerar Data Frame de Dados
MEP.df <- GEN.IC.MEP.DF(
  n, x.cov = X, coeffs.par = coeffs, cuts.points = breaks, rates.par = rates, 
  cens.rate = cens.rate
)

# Ajuste do Limite Inferior do Intervalo
MEP.df$L <- ifelse(MEP.df$L <= 0, 1e-04, MEP.df$L)

# Variável Resposta
y.MEP <- cbind(MEP.df$L, MEP.df$R)

# Número de Parâmetros
n.par.MEP <- length(rates) + ncol(X)

# Chute Inicial
init.MEP <- rep(0.1, n.par.MEP)

# Ajustar
fit.MEP <- optim(
  par = init.MEP, fn = loglikelihood.MEP, method = "L-BFGS-B", hessian = TRUE,
  interval = y.MEP, cens = MEP.df$cens, X = X, cuts.points = breaks, 
  control = list(fnscale=-1)
)

# Visualizar o Ajuste
print(fit.MEP)
cat("Taxas Verdadeiras:", rates, "\n")
cat("Taxas Estimadas:", round(exp(fit.MEP$par[1:length(rates)]), 2), "\n")
cat("Coeficientes Verdadeiros:", coeffs, "\n")
cat("Coeficientes Estimados:", round(fit.MEP$par[(length(rates)+1):n.par.MEP], 2))

# --------------------------------------------
# [2.3] MODELO EXPONENCIAL POR PARTES POTÊNCIA
# --------------------------------------------

# -----------------------------------
# [2.3.1] PARÂMETRO DE POTÊNCIA = 0.8 
# -----------------------------------

# Primeiro Parâmetro de Potência
alpha1 <- 0.8


# Setar semente aleatória
set.seed(123456789)

# Gerar Data Frame de Dados
MEPP.df <- GEN.IC.MEPP.DF(
  n, x.cov = X, coeffs.par = coeffs, cuts.points = breaks, rates.par = rates, 
  alpha.par = alpha1, cens.rate = cens.rate
)

# Ajuste do Limite Inferior do Intervalo
MEPP.df$L <- ifelse(MEPP.df$L <= 0, 1e-04, MEPP.df$L)

# Variável Resposta
y.MEPP <- cbind(MEPP.df$L, MEPP.df$R)

# Número de Parâmetros
n.par.MEPP <- length(rates) + ncol(X) + 1

# Chute Inicial
init.MEPP <- rep(1, n.par.MEPP)

# Ajustar
fit.MEPP <- optim(
  par = init.MEPP, fn = loglikelihood.MEPP, method = "L-BFGS-B", hessian = TRUE,
  interval = y.MEPP, cens = MEPP.df$cens, X = X, cuts.points = breaks, 
  control = list(fnscale=-1)
)

# Visualizar o Ajuste
print(fit.MEPP)
cat("Taxas Verdadeiras:", rates, "\n")
cat("Taxas Estimadas:", exp(fit.MEPP$par[1:length(rates)]), "\n")
cat("Parâmetro de Potência Verdadeiro:", alpha1, "\n")
cat("Parâmetro de Potência Estimado:", fit.MEPP$par[length(rates)+1], "\n")
cat("Coeficientes Verdadeiros:", coeffs, "\n")
cat("Coeficientes Estimados:", fit.MEPP$par[(length(rates)+2):n.par.MEPP])

# -----------------------------------
# [2.3.2] PARÂMETRO DE POTÊNCIA = 1.2 
# -----------------------------------

# Segundo Parâmetro de Potência
alpha2 <- 1.2

# Setar semente aleatória
set.seed(123456789)

# Gerar Data Frame de Dados
MEPP.df <- GEN.IC.MEPP.DF(
  n, x.cov = X, coeffs.par = coeffs, cuts.points = breaks, rates.par = rates, 
  alpha.par = alpha1, cens.rate = cens.rate
)

# Ajuste do Limite Inferior do Intervalo
MEPP.df$L <- ifelse(MEPP.df$L <= 0, 1e-04, MEPP.df$L)

# Variável Resposta
y.MEPP <- cbind(MEPP.df$L, MEPP.df$R)

# Número de Parâmetros
n.par.MEPP <- length(rates) + ncol(X) + 1

# Chute Inicial
init.MEPP <- rep(1, n.par.MEPP)

# Ajustar
fit.MEPP <- optim(
  par = init.MEPP, fn = loglikelihood.MEPP, method = "L-BFGS-B", hessian = TRUE,
  interval = y.MEPP, cens = MEPP.df$cens, X = X, cuts.points = breaks, 
  control = list(fnscale=-1)
)

# Visualizar o Ajuste
print(fit.MEPP)
cat("Taxas Verdadeiras:", rates, "\n")
cat("Taxas Estimadas:", exp(fit.MEPP$par[1:length(rates)]), "\n")
cat("Parâmetro de Potência Verdadeiro:", alpha1, "\n")
cat("Parâmetro de Potência Estimado:", fit.MEPP$par[length(rates)+1], "\n")
cat("Coeficientes Verdadeiros:", coeffs, "\n")
cat("Coeficientes Estimados:", fit.MEPP$par[(length(rates)+2):n.par.MEPP])