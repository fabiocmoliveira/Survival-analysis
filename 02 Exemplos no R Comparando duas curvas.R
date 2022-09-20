## Dados de Carciogeneis
## Estudo com ratos
## tempo em dias da indução da doença até a morte por cancer vaginal
## dois grupos diferenciados por pre tratamento
## Kalbfleish and Prentice 2002, seção 1.1.1

t.g1 <- c(143, 164, 188, 188, 190, 192, 206, 209, 213, 216, 220,
          227, 230, 234, 246, 265, 304, 216, 244)
e.g1 <- c(rep(1,17),0,0)  
t.g2 <- c(142, 156, 163, 198, 205, 232, 232, 233, 233, 233, 233,
          239, 240, 261, 280, 280, 296, 296, 323, 204, 344)  
e.g2 <- c(rep(1,19),0,0) 

tempo <- c(t.g1,t.g2)

evento <- c(e.g1,e.g2)

grupo <- c(rep(1,19),rep(2,21))

cbind(tempo, evento, grupo)

require(survival)

surv.data <- Surv(tempo, evento)

km <- survfit(surv.data ~ grupo)

km

summary(km)

plot(km, col = 2:3, lwd = 2)

plot(km, col = 2:3, lwd = 2, conf.int = TRUE)

survdiff(surv.data ~ grupo)

plot(km, col = 2:3, lwd = 2, fun = "cloglog")

## modelo exponencial
fit1 <- survreg(surv.data ~ grupo , dist = "exp")
summary(fit1)


## Gráfico das curvas de sobrevivência
plot(km, lty = 1:2, lwd = 2)
(lmb.exp1 <- exp(-5.391))
(lmb.exp2 <- exp(-5.391 -0.093))
curve(exp(-lmb.exp1*x), 0, 350, lty = 1, lwd = 2, col = "red", add = T)
curve(exp(-lmb.exp2*x), 0, 350, lty = 2, lwd = 2, col = "red", add = T)


# Ajustando modelo Weibull
m3 <- survreg(surv.data ~ grupo , dist = "weibull")
summary(m3)

(lmb.wei1 <- exp(-5.319))
(lmb.wei2 <- exp(-5.391 -0.132))
(gamma <- 1/0.183)
curve(exp(-(lmb.wei1*x)^gamma), 0, 350, lty = 1, lwd = 2, col = "blue", add = T)
curve(exp(-(lmb.wei2*x)^gamma), 0, 350, lty = 2, lwd = 2, col = "blue", add = T)

# Ajustando modelo log-normal
m5 <- survreg(surv.data ~ grupo , dist = "lognormal")

summary(m5)

(lmb.ln1 <- exp(-5.2823))
(lmb.ln2 <- exp(-5.2823 -0.0931))
(gamma <- 1/0.21)

curve(1 - pnorm(gamma*log(lmb.ln1*x)), 0, 350, lty = 1, lwd = 2, col = "orange", add = T)

curve(exp(-(lmb.wei2*x)^gamma), 0, 350, lty = 2, lwd = 2, col = "orange", add = T)

# Ajustando modelo log-logistico
m6 <- survreg(surv.data ~ grupo , dist ="loglogistic")

summary(m6)

(lmb.ll1 <- exp(-5.271))
(lmb.ll2 <- exp(-5.271 -0.105))
(gamma <- 1/0.116)
curve(1/(1+(lmb.ll1*x)^gamma), 0, 350, lty = 1, lwd = 2, col = "purple", add = T)
curve(1/(1+(lmb.ll2*x)^gamma), 0, 350, lty = 2, lwd = 2, col = "purple", add = T)

m6 <- coxph(surv.data ~ grupo)
summary(m6)
