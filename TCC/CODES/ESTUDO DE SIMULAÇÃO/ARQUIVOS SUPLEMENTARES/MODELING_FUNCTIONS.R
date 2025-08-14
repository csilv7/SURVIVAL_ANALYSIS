# ----------------------------------------
# [1] FUNÇÕES IMPLEMENTADAS PARA MODELAGEM
# ----------------------------------------

# ----------------------------------------------------------------------
# [1.1] ENCONTRA AS MELHORES PARTIÇÕES DE ACORDO COM OS TEMPOS OBSERVADO
# ----------------------------------------------------------------------
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

# -----------------------------------
# [1.2] MODELO EXPONENCIAL POR PARTES
# -----------------------------------

# Contribuição para Log-verossimilhança
f.surv.MEP <- function(L, R, x.cov, coeffs.par, cuts.points, rates.par, type) {
  # Para o Tempo de Falha não Censurado
  if (type == "interval") {
    L <- as.numeric(L)
    R <- as.numeric(R)
    
    # Preditor linear
    effect <- as.numeric(x.cov %*% coeffs.par)
    
    # Contribuição para Função Verossimilhança
    surv.L <- 1 - ppch(q = L/exp(effect), cuts = cuts.points, levels = rates.par)
    surv.R <- 1 - ppch(q = R/exp(effect), cuts = cuts.points, levels = rates.par)
    return(surv.L - surv.R)
  }
  
  # Para o Tempo de Falha Censurado
  if (type == "cens") {
    L <- as.numeric(L)
    
    # Preditor Linear
    effect <- as.numeric(x.cov %*% coeffs.par)
    
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
  coeffs.par <- par[(n.cuts+2):n.par]
  
  # Distinção dos Limites
  L <- as.numeric(interval[, 1])
  R <- as.numeric(interval[, 2])
  
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
  return(flv)
}

# --------------------------------------------
# [1.3] MODELO EXPONENCIAL POR PARTES POTÊNCIA
# --------------------------------------------

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
  coeffs.par <- par[(n.cuts+3):n.par]
  
  # Distinção dos Limites
  L <- as.numeric(interval[, 1])
  R <- as.numeric(interval[, 2])
  
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
  return(flv)
}