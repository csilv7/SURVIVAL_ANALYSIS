# -----------------------
# [1] ATIVAÇÃO DE PACOTES
# -----------------------
library(dplyr)
library(ggplot2)

# ---------------------------
# [1] DISTRIBUIÇÃO EXPONENCIAL
# ---------------------------
# -------------
# [1.1] FUNÇÕES
# -------------

density.exp <- function(times, rate.par) dexp(x = times, rate = rate.par)
survival.exp <- function(times, rate.par) 1 - pexp(q = times, rate = rate.par)
hazard.exp <- function(times, rate.par) density.exp(times, rate.par)/survival.exp(times, rate.par)
accumul.hazard.exp <- function(times, rate.par) - log(x = survival.exp(times, rate.par))

# ----------------------------------------
# [1.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------

set.seed(123456789)        # Semente para reprodutibilidade
n <- 1000                  # Tamanho amostral
times <- rexp(n, rate = 1) # Simulando dados de uma exponencial
alphas <- c(1, 1.5, 2, 2.5, 3) # Valores do parâmetro a serem avaliados

# Criando um Data Frame com valores das funções
dados.exp <- do.call(
  rbind, lapply(alphas, function(alpha) {
    data.frame(
      times = sort(times),
      ft = density.exp(sort(times), alpha),
      st = survival.exp(sort(times), alpha),
      ht = hazard.exp(sort(times), alpha),
      Ht = accumul.hazard.exp(sort(times), alpha),
      rate = factor(alpha)
    )
  })
)

# ---------------------------------------------
# [1.3] FUNÇÃO GRÁFICA & VISUALIZAÇÕES GRÁFICAS
# ---------------------------------------------
plot.function.exp <- function(df, f, label) {
  ggplot(data = df, aes_string(x = "times", y = f, color = "rate")) +
    geom_line(size = 1.15) +
    labs(x = "Tempo", y = label, color = expression(alpha)) +
    scale_color_manual(
      values = c("red", "blue", "green", "purple", "orange"), 
      labels = (levels(df$rate))
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "right",
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
}

# Plotando a função densidade de probabilidade
plot.function.exp(dados.exp, "ft", "Função Densidade de Probabilidade")

# Plotando a função de sobrevivência
plot.function.exp(dados.exp, "st", "Função de Sobrevivência")

# Plotando a função de risco
plot.function.exp(dados.exp, "ht", "Função de Risco")

# Plotando a função de risco acumulado
plot.function.exp(dados.exp, "Ht", "Função de Risco Acumulado")

# ------------------------
# [2] DISTRIBUIÇÃO WEIBULL
# ------------------------
# -------------
# [2.1] FUNÇÕES
# -------------

density.weib <- function(times, shape.par, scale.par) dweibull(x=times, shape=shape.par, scale=scale.par)
survival.weib <- function(times, shape.par, scale.par) 1 - pweibull(q=times, shape=shape.par, scale=scale.par)
hazard.weib <- function(times, shape.par, scale.par) density.weib(times,shape.par,scale.par)/survival.weib(times,shape.par,scale.par)
accumul.hazard.weib <- function(times, shape.par, scale.par) - log(x = survival.weib(times,shape.par,scale.par))

# ----------------------------------------
# [2.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------
n <- 1000                                  # Tamanho amostral
times <- rweibull(n, shape = 2, scale = 1) # Simulando dados de uma Weibull
alpha <- 1                                 # Fixo para simplificar
gammas <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)  # # Valores do parâmetro a serem avaliados

# Criando um Data Frame com valores das funções
dados.weib <- do.call(rbind, lapply(gammas, function(gamma) {
  data.frame(
    times = sort(times),
    ft = density.weib(sort(times), gamma, alpha),
    st = survival.weib(sort(times), gamma, alpha),
    ht = hazard.weib(sort(times), gamma, alpha),
    Ht = accumul.hazard.weib(sort(times), gamma, alpha),
    gamma = factor(gamma)
  )
}))

# ---------------------------------------------
# [2.3] FUNÇÃO GRÁFICA & VISUALIZAÇÕES GRÁFICAS
# ---------------------------------------------
plot.function.weib <- function(df, f, label) {
  ggplot(data = df, aes_string(x = "times", y = f, color = "gamma")) +
    geom_line(size = 1.15) +
    labs(x = "Tempo", y = label, color = expression(gamma)) +
    scale_color_manual(
      values = c("red", "blue", "green", "purple", "orange", "brown"), 
      labels = (levels(df$gamma))
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "right",
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
}

# Plotando a função densidade de probabilidade
plot.function.weib(dados.weib, "ft", "Função Densidade de Probabilidade")

# Plotando a função de sobrevivência
plot.function.weib(dados.weib, "st", "Função de Sobrevivência")

# Plotando a função de risco
plot.function.weib(dados.weib, "ht", "Função de Risco")

# Plotando a função de risco acumulado
plot.function.weib(dados.weib, "Ht", "Função de Risco Acumulado")

# ---------------------------
# [3] DISTRIBUIÇÃO LOG-NORMAL
# ---------------------------
# -------------
# [3.1] FUNÇÕES
# -------------

density.lnorm <- function(times, loc.par, scale.par) dlnorm(x=times, loc.par, scale.par)
survival.lnorm <- function(times, loc.par, scale.par) 1 - pnorm(q=(log(times)-loc.par)/scale.par)
hazard.lnorm <- function(times, loc.par, scale.par) density.lnorm(times, loc.par, scale.par)/survival.lnorm(times, loc.par, scale.par)
accumul.hazard.lnorm <- function(times, loc.par, scale.par) - log(x=survival.lnorm(times, loc.par, scale.par))

# ----------------------------------------
# [3.2] SIMULAÇÃO E VARIAÇÃO DE PARÂMETROS
# ----------------------------------------
n <- 1000                                  # Tamanho amostral
times <- rlnorm(n, meanlog = 0, sdlog = 1) # Simulando dados de uma Log-normal
loc.pars <- c(0, 0.5, 1, 1.5, 2, 2.5)      # Valores de mu
scale.par <- 1                             # Valor fixo de sigma

# Criando um Data Frame com valores das funções
dados.lnorm <- do.call(
  rbind, lapply(loc.pars, function(loc.par) {
    data.frame(
      times = sort(times),
      ft = density.lnorm(sort(times), loc.par, scale.par),
      st = survival.lnorm(sort(times), loc.par, scale.par),
      ht = hazard.lnorm(sort(times), loc.par, scale.par),
      Ht = accumul.hazard.lnorm(sort(times), loc.par, scale.par),
      mu = factor(loc.par)
    )
  })
)

# ---------------------------------------------
# [3.3] FUNÇÃO GRÁFICA & VISUALIZAÇÕES GRÁFICAS
# ---------------------------------------------

plot.function.lnorm <- function(df, f, label) {
  ggplot(data = df, aes_string(x = "times", y = f, color = "mu")) +
    geom_line(size = 1.15) +
    labs(x = "Tempo", y = label, color = expression(mu)) +
    scale_color_manual(
      values = c("red", "blue", "green", "purple", "orange", "brown"), 
      labels = (levels(df$mu))
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "right",
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )
}

# Plotando a função densidade de probabilidade
plot.function.lnorm(dados.lnorm, "ft", "Função Densidade de Probabilidade")

# Plotando a função de sobrevivência
plot.function.lnorm(dados.lnorm, "st", "Função de Sobrevivência")

# Plotando a função de risco
plot.function.lnorm(dados.lnorm, "ht", "Função de Risco")

# Plotando a função de risco acumulado
plot.function.lnorm(dados.lnorm, "Ht", "Função de Risco Acumulado")