# --------------------------
# [1] CONFIGURAÇÕES INICIAIS
# --------------------------

# -------------------------
# [1.1] PACOTES NECESSÁRIOS
# -------------------------
library(eha)

# ---------------------------
# [1.2] FUNÇÕES IMPLEMENTADAS
# ---------------------------
source("~/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/MEP & MEPP/ARQUIVOS SUPLEMENTARES/MODELING_FUNCTIONS.R")
source("~/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/MEP & MEPP/ARQUIVOS SUPLEMENTARES/GENERATION_FUNCTIONS.R")

# ----------------------------
# [2] MODELAGEM DE DADOS REAIS 
# ----------------------------

# -----------------------------
# [2.1] DADOS DE CÂNCER DE MAMA
# -----------------------------

# -----------------------------
# [3] MODELAGEM VIA MONTE CARLO
# -----------------------------

# -----------------------------
# [3.1] AMOSTRA DE TAMANHO 1000
# -----------------------------

# Tamanho Amostral
n <- 1000

# ------------------
# [3.1.1] MODELO MEP
# ------------------

# ---------------------
# [3.1.1.1] GERAR DADOS
# ---------------------

# Setar semente aleatória
set.seed(123456789)

# Parâmetros da Simulação
rates <- c(0.25, 0.5, 0.75, 1) # Taxas de Falha
breaks <- c(0.4, 1.2, 1.8)     # Vetor de Partições
cens.rate <- 1                 # Taxa de Censura (C ~ Exponencial(1))

# Matriz Design
x1 <- rbinom(n, size = 1, prob=0.5)
x2 <- rnorm(n)
X <- cbind(x1, x2)

# Vetor de Coeficientes
coeffs <- c(-0.5, 0.5)

# Gerar Data Frame de Dados
MEP.df <- GEN.IC.MEP.DF(
  n, x.cov = X, coeffs.par = coeffs, cuts.points = breaks, rates.par = rates, cens.rate = cens.rate
)

# --------------------
# [3.1.1.2] AJUSTE MEP
# --------------------

# Variável Resposta
y.MEP <- cbind(MEP.df$L, MEP.df$R)

# Número de Parâmetros
n.par.MEP <- length(rates) + ncol(X)

# Chute Inicial
init.MEP <- rep(0.1, n.par.MEP)

# Ajustar
fit.MEP <- optim(
  par = init.MEP, fn = loglikelihood.MEP, method = "BFGS", hessian = TRUE,
  interval = y.MEP, cens = MEP.df$cens, X = X, cuts.points = breaks, 
  control = list(fnscale=-1)
)

# Visualizar o Ajuste
print(fit.MEP)
cat("Taxas Verdadeiras:", rates, "\n")
cat("Taxas Estimadas:", exp(fit.MEP$par[1:length(rates)]), "\n")
cat("Coeficientes Verdadeiros:", coeffs, "\n")
cat("Coeficientes Estimados:", fit.MEP$par[(length(rates)+1):n.par.MEP])

# -------------------
# [3.1.2] MODELO MEPP
# -------------------

# ---------------------
# [3.1.2.1] GERAR DADOS
# ---------------------

# Setar semente aleatória
set.seed(123456789)

# Parâmetros da Simulação
rates <- c(0.25, 0.5, 0.75, 1) # Taxas de Falha
breaks <- c(0.4, 1.2, 1.8)     # Vetor de Partições
alpha <- 1/10                  # Parâmetro de Potência
cens.rate <- 1                 # Taxa de Censura (C ~ Exponencial(1))

# Matriz Design
x1 <- rbinom(n, size = 1, prob=0.5)
x2 <- rnorm(n)
X <- cbind(x1, x2)

# Vetor de Coeficientes
coeffs <- c(-0.5, 0.5)

# Gerar Data Frame de Dados
MEPP.df <- GEN.IC.MEPP.DF(
  n, x.cov = X, coeffs.par = coeffs, cuts.points = breaks, rates.par = rates, 
  alpha.par = alpha, cens.rate = cens.rate
)

# ---------------------
# [3.1.2.2] AJUSTE MEPP
# ---------------------

# Variável Resposta
y.MEPP <- cbind(MEPP.df$L, MEPP.df$R)

# Número de Parâmetros
n.par.MEPP <- length(rates) + ncol(X) + 1

# Chute Inicial
init.MEPP <- rep(1, n.par.MEPP)

# Ajustar
fit.MEPP <- optim(
  par = init.MEPP, fn = loglikelihood.MEPP, method = "BFGS", hessian = TRUE,
  interval = y.MEPP, cens = MEPP.df$cens, X = X, cuts.points = breaks, 
  control = list(fnscale=-1)
)

# Visualizar o Ajuste
print(fit.MEPP)
cat("Taxas Verdadeiras:", rates, "\n")
cat("Taxas Estimadas:", exp(fit.MEPP$par[1:length(rates)]), "\n")
cat("Parâmetro de Potência Verdadeiro:", alpha, "\n")
cat("Parâmetro de Potência Estimado:", fit.MEPP$par[length(rates)+1], "\n")
cat("Coeficientes Verdadeiros:", coeffs, "\n")
cat("Coeficientes Estimados:", fit.MEPP$par[(length(rates)+2):n.par])