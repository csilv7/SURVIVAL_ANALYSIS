# --------------------------
# [*] CONFIGURAÇÕES INICIAIS
# --------------------------

# Pacotes Necessários
library(dplyr)
library(survival)
library(ggplot2)
library(patchwork)
#library(gridExtra)

# Caminho URL para os dados
url <- "https://docs.ufpr.br/~giolo/asa/dados/leucemia.txt"

# Leitura dos dados
dados <- read.table(url, header = TRUE)

# ----------------
# [1] KAPLAN-MEIER
# ----------------

# Estimador de Kaplan-Meier (Pacote survival)
ekm <- survfit(Surv(tempos, cens)~1, data = dados)

# Visualizar Curva de Sobrevivência (Pacote ggplot2)
ggplot(data = NULL, aes(x = ekm$time, y = ekm$surv)) +
  geom_line(color = "blue", lwd = 1.2) +
  geom_ribbon(aes(ymin = ekm$lower, ymax = ekm$upper), fill = "blue", alpha = 0.2) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# ---------------------
# [2] TESTE DE LOG-RANK (Igual de Curvas de Sobrevivência)
# ---------------------

# Ajustando a coluna r_6
dados.adj <- dados %>%
  mutate(grupo = ifelse(r6 == 0, "Category Zero", "Category One"))

# Cria um objeto de Surv
dados.adj_surv <- Surv(dados.adj$tempos, dados.adj$cens)

# Aplica o Teste Logrank
logrank.test <- survdiff(dados.adj_surv ~ grupo, data = dados.adj)

# Imprimir Resultados do Teste Logrank
print(logrank.test)

# ----------------
# [3] NELSON-AALEN
# ----------------

# Nível de Significância
alpha <- 0.05

# Quantil da Distribuição Normal
z <- qnorm(1 - alpha/2)

# Estatísticas
fit <- survfit(Surv(tempos, cens)~1, data = dados)

# h(t) e Var[h(t)]
ena <- cumsum(fit$n.event/fit$n.risk)
var.ena <- cumsum(fit$n.event/fit$n.risk^2)

# S(t) e Var[S(t)]
ena.st <- exp(-ena)
var.ena.st <- ena.st^2 * var.ena

# Intervalo de Confiança de S(t)
conf.lower <- ena.st - z * sqrt(var.ena.st)
conf.upper <- ena.st + z * sqrt(var.ena.st)

# Risco Acumulado
cumhaz <- ggplot(data = NULL, aes(x = fit$time, y = ena)) +
  geom_line(color = "red", lwd = 1.2) +
  labs(x = "Tempo", y = "Risco Acumulado") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Probabilidade de Sobrevivência
surv <- ggplot(data = NULL, aes(x = fit$time, y = ena.st)) +
  geom_line(color = "blue", lwd = 1.2) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  geom_ribbon(aes(ymin = conf.lower, ymax = conf.upper), fill = "blue", alpha = 0.2) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Visualização Emparelhada
cumhaz | surv
#grid.arrange(cumhaz, surv, ncol = 2)

# ------------------------------
# [4] COMPARAÇÃO DOS ESTIMADORES
# ------------------------------

# Cria um objeto `ggplot`
ggplot(data = NULL, aes(x = ekm$time)) +
  geom_line(aes(y = ekm$surv, color = "Kaplan-Meier"), lwd = 1.2) +
  geom_ribbon(aes(ymin = ekm$lower, ymax = ekm$upper), fill = "blue", alpha = 0.2) +
  
  geom_line(aes(y = ena.st, color = "Nelson-Aalen"), lwd = 1.2) +
  geom_ribbon(aes(ymin = conf.lower, ymax = conf.upper), fill = "red", alpha = 0.2) +
  
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência", color = "Estimador") +
  
  scale_color_manual(values = c("Kaplan-Meier" = "blue", "Nelson-Aalen" = "red")) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )