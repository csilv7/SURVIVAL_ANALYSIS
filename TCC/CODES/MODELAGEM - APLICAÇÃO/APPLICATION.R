# -----------------------
# [*] *******************
# -----------------------

# Pacotes Necessários
library(dplyr)
library(eha)
library(ggplot2)

# Funções Usadas
time.grid.interval <- function(li, ri, type = "OBS", bmax = 10) {
  # Função para criar uma grade de tempo (partição)
  # li: limite inferior dos intervalos observados
  # ri: limite superior dos intervalos observados
  # type: "OBS" para partições baseadas em tempos observados
  #       "EQUI" para partições equidistantes
  # bmax: número máximo de intervalos (opcional)
  
  # type = "OBS"
  if (type == "OBS") {
    times <- c(0, li, ri[is.finite(ri)], Inf)
    grid.vet <- sort(unique(times))
    size <- length(grid.vet)
    
    if (bmax < size) {
      k <- round((size - 1) / bmax)
      id.grid <- round(seq(k, size - 1, length.out = bmax))
      grid.vet <- c(0, grid.vet[-1][id.grid])
      return(grid.vet)
      
    } else {
      times <- c(0, li, ri[is.finite(ri)], Inf)
      grid.vet <- sort(unique(times))
      return(grid.vet)
    }
  }
  # type = "EQUI"
  if (type == "EQUI") {
    grid.vet <- seq(0, max(ri[ri != Inf]), length.out = bmax)
    grid.vet <- c(grid.vet, Inf)
    
    return(grid.vet)
  }
}

# ------------------------------
# [1] LEITURA E AJUSTE DOS DADOS
# ------------------------------

# Caminho URL do arquivo
url <- "https://docs.ufpr.br/~giolo/asa/dados/breast.txt"

# Leitura do arquivo
breast <- read.table(url, h = TRUE)

# Salvar Conjunto de Dados em CSV
path <- "C:/Users/user/Documents/PROJETOS/R/SURVIVAL_ANALYSIS/TCC/CONJUNTO DE DADOS/breast.csv"
write.csv(breast, file = path, row.names = FALSE)

# Ajute nos dados
breast$left <- ifelse(breast$left == 0, 0.01, breast$left)     # Lim Inf
breast$right <- ifelse(is.na(breast$right), Inf, breast$right) # Lim Sup

# ----------------------
# [2] AJUSTE DOS MODELOS
# ----------------------

# ----------------------------
# [2.1] EXPONENCIAL POR PARTES
# ----------------------------

# ----------------------------
# [2.1.1] DEFINIÇÃO DE FUNÇÕES
# ----------------------------

# Contribuição para Log-verossimilhança
f.surv.MEP <- function(L, R, x.cov, coeffs.par, cuts.points, rates.par, type) {
  # Para o Tempo de Falha não Censurado
  if (type == "interval") {
    # Preditor linear
    effect <- x.cov %*% coeffs.par
    
    # Contribuição para Função Verossimilhança
    surv.L <- 1 - ppch(q = L/exp(effect), cuts = cuts.points, levels = rates.par)
    surv.R <- 1 - ppch(q = R/exp(effect), cuts = cuts.points, levels = rates.par)
    return(surv.L - surv.R)
  }
  
  # Para o Tempo de Falha Censurado
  if (type == "cens") {
    # Preditor Linear
    effect <- x.cov %*% coeffs.par
    
    # Contribuição para Função Verossimilhança
    surv.L <- 1 - ppch(q = L/exp(effect), cuts = cuts.points, levels = rates.par)
    return(surv.L)
  }
}

# Função Log-verossimilhança
loglikelihood.MEP <- function(par, interval, cens, X, cuts.points) {
  # Número de parâmetros e partições
  n.par <- length(par)
  n.cuts <- length(cuts.points)
  
  # Distinção de parâmetros
  rates.par <- exp(par[1:(n.cuts+1)])
  coeffs.par <- par[(n.cuts + 2):n.par]
  
  # Distinção dos Limites
  L <- interval[, 1]
  R <- interval[, 2]
  
  # Função Log-verossimilhança
  S1 <- f.surv.MEP(
    L=L, R=R, x.cov=X, coeffs.par=coeffs.par,
    cuts.points=cuts.points, rates.par=rates.par,
    type = "interval"
  )
  S2 <- f.surv.MEP(
    L=L, x.cov=X, coeffs.par=coeffs.par,
    cuts.points=cuts.points, rates.par=rates.par,
    type = "cens"
  )
  
  flv <- sum(cens * log(S1) + (1 - cens) * log(S2))
  return(-flv)
}

# ---------------------
# [2.1.2] IMPLEMENTAÇÃO
# ---------------------

# Variável Resposta (y)
y <- cbind(breast$left, breast$right)

# Ajuste de escala
y <- y / 12

# Selecionando o Vetor de Partições
cuts.grid <- time.grid.interval(li=y[, 1], ri=y[, 2], type="OBS", bmax = 2)
cuts.grid

# Separação dos Dados
X <- as.matrix(breast$ther, ncol = 1)        # Matrix de Covariáveis
cuts <- cuts.grid[c(-1, -length(cuts.grid))] # Pontos de Corte
n.cuts <- length(cuts)                       # Número de Potnos de Corte

# Número de parâmetros para estimar
n.par.MEP <- length(cuts.grid) + ncol(X) - 1

# Chute Inicial
init.MEP <- rep(1, n.par.MEP)

# Maximização
fit.MEP <- optim(par=init.MEP, fn = loglikelihood.MEP, 
                 gr = NULL, method = "BFGS", hessian = TRUE, 
                 interval=y, cens=breast$cens, X=X, cuts.points=cuts)

# Visualização do Ajuste
fit.MEP

# -------------------------------------------
# [2.1.3] VISUALIZANDO CURVA DE SOBREVIVÊNCIA
# -------------------------------------------

# Captura das estimativas
MEP.r <- list(rates = exp(fit.MEP$par[1:(n.cuts+1)]), coeffs = fit.MEP$par[(n.cuts + 2):n.par.MEP])

# t & S(t)
t <- seq(0, 50, length.out = 94)
st <- 1 - ppch(q = t/exp(X %*% MEP.r$coeffs), cuts = 12 * cuts, levels = MEP.r$rates)

# Curva de Sobrevivência
ggplot(data = NULL, aes(x = t, y = st)) +
  geom_line(color = "blue", lwd = 1.1) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# -------------------------------------
# [2.2] EXPONENCIAL POR PARTES POTÊNCIA
# -------------------------------------

# ----------------------------
# [2.2.1] DEFINIÇÃO DE FUNÇÕES
# ----------------------------

# Contribuição para Log-verossimilhança
f.surv.MEPP <- function(L, R, x.cov, coeffs.par, cuts.points, rates.par, power.par, type) {
  # Para o Tempo de Falha não Censurado
  if (type == "interval") {
    # Preditor linear
    effect <- x.cov %*% coeffs.par
    
    # Contribuição para Função Verossimilhança
    surv.L <- 1 - ppch(q = L/exp(effect), cuts = cuts.points, levels = rates.par)^power.par
    surv.R <- 1 - ppch(q = R/exp(effect), cuts = cuts.points, levels = rates.par)^power.par
    return(surv.L - surv.R)
  }
  
  # Para o Tempo de Falha Censurado
  if (type == "cens") {
    # Preditor Linear
    effect <- x.cov %*% coeffs.par
    
    # Contribuição para Função Verossimilhança
    surv.L <- 1 - ppch(q = L/exp(effect), cuts = cuts.points, levels = rates.par)^power.par
    return(surv.L)
  }
}

# Função Log-verossimilhança
loglikelihood.MEPP <- function(par, interval, cens, X, cuts.points) {
  # Número de parâmetros e partições
  n.par <- length(par)
  n.cuts <- length(cuts.points)
  
  # Distinção de parâmetros
  rates.par <- exp(par[1:(n.cuts+1)])
  power.par <- par[n.cuts+2]
  coeffs.par <- par[(n.cuts + 3):n.par]
  
  # Distinção dos Limites
  L <- interval[, 1]
  R <- interval[, 2]
  
  # Função Log-verossimilhança
  S1 <- f.surv.MEPP(
    L=L, R=R, x.cov=X, coeffs.par=coeffs.par,
    cuts.points=cuts.points, rates.par=rates.par,
    power.par = power.par,
    type = "interval"
  )
  S2 <- f.surv.MEPP(
    L=L, x.cov=X, coeffs.par=coeffs.par,
    cuts.points=cuts.points, rates.par=rates.par,
    power.par = power.par,
    type = "cens"
  )
  
  flv <- sum(cens * log(S1) + (1 - cens) * log(S2))
  return(-flv)
}

# ---------------------
# [2.2.2] IMPLEMENTAÇÃO
# ---------------------

# Variável Resposta (y)
y <- cbind(breast$left, breast$right)

# Ajuste de escala
y <- y / 12

# Selecionando o Vetor de Partições
cuts.grid <- time.grid.interval(li=y[, 1], ri=y[, 2], type="OBS", bmax = 2)
cuts.grid

# Separação dos Dados
X <- as.matrix(breast$ther, ncol = 1)        # Matrix de Covariáveis
cuts <- cuts.grid[c(-1, -length(cuts.grid))] # Pontos de Corte
n.cuts <- length(cuts)                       # Número de Potnos de Corte

# Número de parâmetros para estimar
n.par.MEPP <- length(cuts.grid) + ncol(X)

# Chute Inicial
init.MEPP <- rep(1, n.par.MEPP)

# Maximização
fit.MEPP <- optim(par=init.MEPP, fn = loglikelihood.MEPP, 
                 gr = NULL, method = "BFGS", hessian = TRUE, 
                 interval=y, cens=breast$cens, X=X, cuts.points=cuts)

# Visualização do Ajuste
fit.MEPP

# -------------------------------------------
# [2.2.3] VISUALIZANDO CURVA DE SOBREVIVÊNCIA
# -------------------------------------------

# Captura das estimativas
MEPP.r <- list(rates = exp(fit.MEPP$par[1:(n.cuts+1)]), 
               power = fit.MEPP$par[n.cuts+2], 
               coeffs = fit.MEPP$par[(n.cuts + 3):n.par.MEPP])

# t & S(t)
t <- seq(0, 50, length.out = 94)
st <- 1 - ppch(q = t/exp(X %*% MEP.r$coeffs), cuts = 12 * cuts, levels = MEP.r$rates)

# Curva de Sobrevivência
ggplot(data = NULL, aes(x = t, y = st)) +
  geom_line(color = "blue", lwd = 1.1) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# ---------------------------
# [3] DIAGNÓSTICO DOS MODELOS
# ---------------------------

# -------------------------
# [4] COMPARÇÃO DOS MODELOS
# -------------------------

n <- nrow(y)

k.MEP <- length(fit.MEP$par)
max.loglikelihood.MEP <- fit.MEP$value

AIC.MEP <- 2 * k.MEP - 2 * max.loglikelihood.MEP
BIC.MEP <- k.MEP * log(n) - 2 * max.loglikelihood.MEP

k.MEPP <- length(fit.MEPP$par)
max.loglikelihood.MEPP <- fit.MEPP$value

AIC.MEPP <- 2 * k.MEPP - 2 * max.loglikelihood.MEPP
BIC.MEPP <- k.MEPP * log(n) - 2 * max.loglikelihood.MEPP