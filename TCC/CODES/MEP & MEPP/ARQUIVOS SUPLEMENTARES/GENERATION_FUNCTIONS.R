# -----------------------------------------------
# [1] FUNÇÕES IMPLEMENTADAS PARA GERAÇÃO DE DADOS
# -----------------------------------------------

library(eha)

# -----------------------------
# [1.1] CENSURA À DIREITA - MEP
# -----------------------------
MEP.TIME.r <- function(t, x.cov, coeffs.par, cuts.points, rates.par, u.unif) {
  # Equação em t: S(t) = u
  # x.cov: Matriz Design (Variáveis Preditoras)
  # coeffs.par: Coeficientes do Modelo
  # cuts.points: Partições no Eixo dos Tempos
  # rates.par: Taxa de Falha em cada intervalo
  # u.unif: Realização aleatória de uma Uniforme(0, 1)
  
  # Preditor Linear
  effect <- x.cov %*% coeffs.par
  
  # S(t)
  surv <- 1 - ppch(q = t/exp(effect), cuts = cuts.points, levels = rates.par)
  
  # Retornar S(t) - U = 0
  return(surv - u.unif)
}

GEN.MEP.r <- function(n, x.cov, coeffs.par, cuts.points, rates.par) {
  # Resolver para t a equação: S(t) = u <=> S^{-1}(u) = t
  # n: Tamanho da Amostra
  # x.cov: Matriz Design (Variáveis Preditoras)
  # coeffs.par: Coeficientes do Modelo
  # cuts.points: Partições no Eixo dos Tempos
  # rates.par: Taxa de Falha em cada intervalo
  
  # Vetor para armazenar os tempos gerados
  MEP.times <- numeric(n)
  
  # Iteração de geração por instância
  for (i in 1:n) {
    # Gerando um único valor de U para cada iteração
    u <- runif(1)
    
    # Encontrando a raiz da equação
    raiz <- uniroot(
      MEP.TIME.r, interval = c(0, 10000), x.cov = x.cov[i,],
      coeffs.par = coeffs.par, cuts.points = cuts.points, 
      rates.par = rates.par,   
      u.unif = u
    )
    
    # Armazenando Resultado
    MEP.times[i] <- raiz$root
  }
  
  # Retornar Tempos Simulados
  return(MEP.times)
}

# ------------------------------
# [1.2] CENSURA INTERVALAR - MEP
# ------------------------------
GEN.IC.MEP.DF <- function(n, x.cov, coeffs.par, cuts.points, rates.par, cens.rate) {
  # Vetores para armazenamento do Tempo Observado
  MEP.IC.times <- numeric(n)
  
  # -------------------------------
  # DADOS NÃO CENSURADOS (cens = 1)
  # -------------------------------
  
  # Tempo de Falha
  t <- GEN.MEP.r(n=n, x.cov=x.cov, coeffs.par=coeffs.par, cuts.points=cuts.points, rates.par=rates.par)
  
  # Tempo de Censura
  c <- rexp(n, rate = cens.rate)
  
  # Tempo Observado
  MEP.IC.times <- pmin(t, c)
  
  # Variável Indicadora
  cens <- ifelse(t <= c, 1, 0)
  
  # ---------------------------
  # DADOS CENSURADOS (cens = 0)
  # ---------------------------
  
  # Limite Inferior (L - Left) e Superior (R - Right)
  L <- R <- MEP.IC.times * NA
  
  # Iteração de Ajustes de (L, R)
  for (i in 1:n) {
    # Verifica se é Tempo de Censura
    if (cens[i] == 0) {
      # Armazenamento de Dados
      L[i] <- MEP.IC.times[i] # Lim Inf
      R[i] <- Inf             # Lim Sup
    }
    # Caso contrário
    else {
      # Armanamento de Dados
      L[i] <- 0                        # Lim Inf         
      add <- stats::runif(1, 0.1, 0.5) # Ajuste do Intervalo
      R[i] <- add                      # Lim Sup
      
      # Checa se o valor pertence ao intervalo gerado
      check <- (L[i] <= MEP.IC.times[i] & MEP.IC.times[i] < R[i])
      
      # Equanto não pertencer refaz os ajustes
      while (!check) {
        # Armanamento de Dados
        L[i] <- L[i] + add               # Lim Inf
        add <- stats::runif(1, 0.1, 0.5) # Ajuste do Intervalo
        R[i] <- R[i] + add               # Lim Sup
        
        # Checa se o valor pertence ao intervalo gerado
        check <- (L[i] <= MEP.IC.times[i] & MEP.IC.times[i] < R[i])
      }
    }
  }
  
  # Ajustando a Formatação de Data Frame
  DF <- data.frame(L, R, MEP.IC.times, cens, x.cov)
  
  # Retornar Dados Simulados
  return(DF)
}

# ------------------------------
# [1.3] CENSURA À DIREITA - MEPP
# ------------------------------
MEPP.TIME.r <- function(t, x.cov, coeffs.par, cuts.points, rates.par, alpha.par, u.unif) {
  # Equação em t: S(t) = u
  # x.cov: Matriz Design (Variáveis Preditoras)
  # coeffs.par: Coeficientes do Modelo
  # cuts.points: Partições no Eixo dos Tempos
  # rates.par: Taxa de Falha em cada intervalo
  # u.unif: Realização aleatória de uma Uniforme(0, 1)
  
  # Preditor Linear
  effect <- x.cov %*% coeffs.par
  
  # S(t)
  surv <- 1 - ppch(q = t/exp(effect), cuts = cuts.points, levels = rates.par)^alpha.par
  
  # Retornar S(t) - U = 0
  return(surv - u.unif)
}

GEN.MEPP.r <- function(n, x.cov, coeffs.par, cuts.points, rates.par, alpha.par) {
  # Resolver para t a equação: S(t) = u <=> S^{-1}(u) = t
  # n: Tamanho da Amostra
  # x.cov: Matriz Design (Variáveis Preditoras)
  # coeffs.par: Coeficientes do Modelo
  # cuts.points: Partições no Eixo dos Tempos
  # rates.par: Taxa de Falha em cada intervalo
  
  # Vetor para armazenar os tempos gerados
  MEPP.times <- numeric(n)
  
  # Iteração de geração por instância
  for (i in 1:n) {
    # Gerando um único valor de U para cada iteração
    u <- runif(1)
    
    # Encontrando a raiz da equação
    raiz <- uniroot(
      MEPP.TIME.r, interval = c(0, 10000), x.cov = x.cov[i,],
      coeffs.par = coeffs.par, cuts.points = cuts.points, 
      rates.par = rates.par, alpha.par = alpha.par,
      u.unif = u
    )
    
    # Armazenando Resultado
    MEPP.times[i] <- raiz$root
  }
  
  # Retornar Tempos Simulados
  return(MEPP.times)
}

# -------------------------------
# [1.4] CENSURA INTERVALAR - MEPP
# -------------------------------
GEN.IC.MEPP.DF <- function(n, x.cov, coeffs.par, cuts.points, rates.par, alpha.par, cens.rate) {
  # Vetores para armazenamento do Tempo Observado
  MEPP.IC.times <- numeric(n)
  
  # -------------------------------
  # DADOS NÃO CENSURADOS (cens = 1)
  # -------------------------------
  
  # Tempo de Falha
  t <- GEN.MEPP.r(n=n, x.cov=x.cov, coeffs.par=coeffs.par, cuts.points=cuts.points, rates.par=rates.par, alpha.par=alpha.par)
  
  # Tempo de Censura
  c <- rexp(n, rate = cens.rate)
  
  # Tempo Observado
  MEPP.IC.times <- pmin(t, c)
  
  # Variável Indicadora
  cens <- ifelse(t <= c, 1, 0)
  
  # ---------------------------
  # DADOS CENSURADOS (cens = 0)
  # ---------------------------
  
  # Limite Inferior (L - Left) e Superior (R - Right)
  L <- R <- MEPP.IC.times * NA
  
  # Iteração de Ajustes de (L, R)
  for (i in 1:n) {
    # Verifica se é Tempo de Censura
    if (cens[i] == 0) {
      # Armazenamento de Dados
      L[i] <- MEPP.IC.times[i] # Lim Inf
      R[i] <- Inf             # Lim Sup
    }
    # Caso contrário
    else {
      # Armanamento de Dados
      L[i] <- 0                        # Lim Inf         
      add <- stats::runif(1, 0.1, 0.5) # Ajuste do Intervalo
      R[i] <- add                      # Lim Sup
      
      # Checa se o valor pertence ao intervalo gerado
      check <- (L[i] <= MEPP.IC.times[i] & MEPP.IC.times[i] < R[i])
      
      # Equanto não pertencer refaz os ajustes
      while (!check) {
        # Armanamento de Dados
        L[i] <- L[i] + add               # Lim Inf
        add <- stats::runif(1, 0.1, 0.5) # Ajuste do Intervalo
        R[i] <- R[i] + add               # Lim Sup
        
        # Checa se o valor pertence ao intervalo gerado
        check <- (L[i] <= MEPP.IC.times[i] & MEPP.IC.times[i] < R[i])
      }
    }
  }
  
  # Ajustando a Formatação de Data Frame
  DF <- data.frame(L, R, MEPP.IC.times, cens, x.cov)
  
  # Retornar Dados Simulados
  return(DF)
}