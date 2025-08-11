# --------------------------
# [1] CONFIGURAÇÕES INICIAIS
# --------------------------

# -------------------------
# [1.1] PACOTES NECESSÁRIOS
# -------------------------
library(dplyr)
library(eha)
library(ggplot2)

# ---------------------------
# [1.2] FUNÇÕES IMPLEMENTADAS
# ---------------------------
source("~/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/MEP & MEPP/ARQUIVOS SUPLEMENTARES/MODELING_FUNCTIONS.R")
#source("~/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CODES/MEP & MEPP/ARQUIVOS SUPLEMENTARES/GENERATION_FUNCTIONS.R")

# ----------------------------
# [2] MODELAGEM DE DADOS REAIS 
# ----------------------------

# -----------------------------
# [2.1] DADOS DE CÂNCER DE MAMA
# -----------------------------

# ------------------------------------
# [2.1.1] LEITURA DE AJUSTES DOS DADOS
# ------------------------------------

# Caminho URL do arquivo
url <- "https://docs.ufpr.br/~giolo/asa/dados/breast.txt"

# Leitura do arquivo
breast <- read.table(url, h = TRUE)

# Salvar Conjunto de Dados em CSV
path.for.save <- "~/PROJETOS/R/ANÁLISE DE SOBREVIVÊNCIA/TCC/CONJUNTO DE DADOS/INTERVAL CENSORED/breast.csv"
write.csv(breast, file = path.for.save, row.names = FALSE)

# Ajute nos dados
breast$left <- ifelse(breast$left == 0, 0.01, breast$left)     # Lim Inf
breast$right <- ifelse(is.na(breast$right), Inf, breast$right) # Lim Sup

# Variáveis fixas
y <- cbind(breast$left, breast$right) / 12 # Escala de: MESES -> ANOS
X <- as.matrix(breast$ther, ncol = 1)
cens <- breast$cens

# Partições a serem testadas
bmax.values <- 2:6

# Data frames para armazenar resultados
results <- data.frame(
  bmax = integer(),
  model = character(),
  AIC = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

# ---------------------------------
# [2.1.2] LOOP PARA AJUSTAR MODELOS
# ---------------------------------

for (bmax in bmax.values) {
  # Criação e Ajustes das Partições
  cuts.grid <- time.grid.interval(li = y[, 1], ri = y[, 2], type = "OBS", bmax = bmax)
  cuts <- cuts.grid[c(-1, -length(cuts.grid))]
  
  # Ajustar Modelo Exponencial por Partes (MEP)
  n.par.MEP <- length(cuts.grid) + ncol(X) - 1
  init.MEP <- rep(1, n.par.MEP)
  fit.MEP <- optim(
    par = init.MEP, fn = loglikelihood.MEP, method = "BFGS", hessian = TRUE,
    interval = y, cens = cens, X = X, cuts.points = cuts, control = list(fnscale=-1)
  )
  
  # Calcular AIC e BIC do MEP
  AIC.MEP <- 2 * n.par.MEP - 2 * fit.MEP$value
  BIC.MEP <- n.par.MEP * log(nrow(y)) - 2 * fit.MEP$value
  
  # Ajustar Modelo Exponencial por Partes Potência (MEPP)
  n.par.MEPP <- length(cuts.grid) + ncol(X)
  init.MEPP <- rep(1, n.par.MEPP)
  fit.MEPP <- optim(
    par = init.MEPP, fn = loglikelihood.MEPP, method = "BFGS", hessian = TRUE,
    interval = y, cens = cens, X = X, cuts.points = cuts, control = list(fnscale=-1)
  )
  
  # Calcular AIC e BIC do MEPP
  AIC.MEPP <- 2 * n.par.MEPP - 2 * fit.MEPP$value
  BIC.MEPP <- n.par.MEPP * log(nrow(y)) - 2 * fit.MEPP$value
  
  # Armazanar resultados no Data Frame
  results <- rbind(
    results,
    data.frame(bmax = bmax, model = "MEP", AIC = AIC.MEP, BIC = BIC.MEP),
    data.frame(bmax = bmax, model = "MEPP", AIC = AIC.MEPP, BIC = BIC.MEPP)
  )
}

# -----------------------------------
# [2.1.3] VISUALIZAÇÃO DOS RESULTADOS
# -----------------------------------

# Resultados Tabelados
print(results)

# AIC versus Nº Partições
ggplot(results, aes(x = bmax, y = AIC, color = model)) +
  geom_line(lwd = 1.1) + geom_point(size = 3) +
  labs(x = "Número de Partições", y = "AIC", color = "") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# BIC versus Nº Partições
ggplot(results, aes(x = bmax, y = BIC, color = model)) +
  geom_line(lwd = 1.1) + geom_point(size = 3) +
  labs(x = "Número de Partições", y = "BIC", color = "") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )