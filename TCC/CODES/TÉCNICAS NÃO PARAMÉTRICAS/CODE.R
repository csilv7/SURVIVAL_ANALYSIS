library(survival)
library(ggplot2)
library(gridExtra)

# -----------------------------
# [1] TÉCNICAS NÃO PARAMÉTRICAS
# -----------------------------

# Caminho URL para os dados
url <- "https://docs.ufpr.br/~giolo/asa/dados/leucemia.txt"

# Leitura dos dados
dados <- read.table(url, header = TRUE)

# ------------------
# [1.1] KAPLAN-MEIER
# ------------------

# ---------------------------------
# [1.1.2] Estimador de Kaplan-Meier
# ---------------------------------
ekm <- survfit(Surv(tempos, cens)~1, data = dados)

# --------------------
# [1.1.3] Visualização
# -------------------

# Gráfico com ggplot2
ggplot(data = NULL, aes(x = ekm$time, y = ekm$surv)) +
  geom_line(color = "blue", lwd = 1.2) +
  geom_ribbon(aes(ymin = ekm$lower, ymax = ekm$upper), fill = "blue", alpha = 0.2) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# ------------------
# [1.2] NELSON-AALEN
# ------------------

# Nível de Significância
alpha <- 0.05

# Quantil da Distribuição Normal
z <- qnorm(1 - alpha/2)

# ---------------------------------
# [1.2.2] Estimador de Nelson-Aalen
# ---------------------------------
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

# --------------------
# [1.2.3] Visualização
# --------------------
cumhaz <- ggplot(data = NULL, aes(x = fit$time, y = ena)) +
  geom_line(color = "red", lwd = 1.2) +
  labs(x = "Tempo", y = "Risco Acumulado") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

surv <- ggplot(data = NULL, aes(x = fit$time, y = ena.st)) +
  geom_line(color = "blue", lwd = 1.2) +
  labs(x = "Tempo", y = "Probabilidade de Sobrevivência") +
  geom_ribbon(aes(ymin = conf.lower, ymax = conf.upper), fill = "blue", alpha = 0.2) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

grid.arrange(cumhaz, surv, ncol = 2)

# ----------------
# [1.3] COMPARAÇÃO
# ----------------
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
