## Exemplo 1
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
grupo <- c(rep(1,19),rep(0,21))
dados <- data.frame(tempo, evento, grupo)
(dados <- dados[order(tempo),])
table(tempo)

require(survival)

# Estimador de Kaplan-Meier
(surv.km <- survfit(Surv(tempo, evento) ~ grupo, dados))
summary(surv.km)

# Gráfico da curva de sobrevivência
plot(surv.km, col = c("red","blue"), lwd = 2)
plot(surv.km, col = c("red","blue"), lwd = 2, fun = "cloglog")

# Teste de log-rank
survdiff(Surv(tempo, evento) ~ grupo, dados)

# Modelo de Cox
# Descrever o modelo e a a função de verossimilhança parcial
m1 <- coxph(Surv(tempo, evento) ~ grupo, dados, ties = "exact")
summary(m1)

m2 <- coxph(Surv(tempo, evento) ~ grupo, dados, ties = "efron")
summary(m2)

m3 <- coxph(Surv(tempo, evento) ~ grupo, dados, ties = "breslow")
summary(m3)

cox.zph(m1)

cox.zph(m2)

cox.zph(m3)

# Gráficos das curvas de sobrevivência

Lambda0.1 <- basehaz(m1, centered = FALSE) 
S0.1 <- exp(-Lambda0$hazard)

Lambda0.2 <- basehaz(m2, centered = FALSE) 
S0.2 <- exp(-Lambda0$hazard)

Lambda0.3 <- basehaz(m2, centered = FALSE) 
S0.3 <- exp(-Lambda0$hazard)

## Curvas de sobreviv?ncia estimada por K-M
plot(surv.km, col = c("blue","red"), lwd = 2)

## Curvas de sobreviv?ncia estimada modelo 1
points(Lambda0$time, S0.1, type = "s", col = "blue", lty = 2, lwd = 2)
points(Lambda0$time, S0.1^exp(coef(m1)), type = "s", col = "red", lty = 2, lwd = 2)

## Curvas de sobreviv?ncia estimada modelo 2
points(Lambda0$time, S0.2, type = "s", col = "black", lty = 2, lwd = 2)
points(Lambda0$time, S0.2^exp(coef(m2)), type = "s", col = "gray", lty = 2, lwd = 2)


# Análise dos resíduos
res.m <- residuals(m1, type = "martingale")
hist(res.m)
plot(res.m)

res.d <- residuals(m1, type = "deviance")
plot(res.d)

res.ce <- evento - res.m
plot(tempo, res.ce)

km.ce <- survfit(Surv(res.ce, evento) ~ 1)

plot(km.ce)
curve(exp(-x), 0 , 4, col = 2, lwd = 2, add = T)

## Exemplo 2
# Sobrevida de pacientes com Leucemi Aguda
# tempo de sobrevivência, em semanas, de 17 pacientes com lecemia aguda
# WBC : contagem de glóbulos brancos
# grupos: 1) apresenta antígeno Calla na superfício dos bastos (Ag+)
#         2) ainda não expressaram o antígeno Calla (Ag-)    

tempo <- c(65, 156, 100, 134, 16, 108, 121, 4, 
           39, 143, 56, 26, 22, 1, 1, 5, 65, 
           56, 65, 17, 7, 16, 22, 3, 4, 2, 3, 8, 4, 3, 30,
           4, 43) 
evento <- rep(1,17+16)
wbc <- c(2300, 750, 4300, 2600, 6000, 10000, 10000, 17000, 5400, 
         7000, 9400, 32000, 35000, 100000, 100000, 52000, 100000,
         4400, 3000, 4000, 1500, 9000, 5300, 10000, 19000,  
         27000, 28000, 31000, 26000, 21000, 79000, 100000, 100000)
lwbc <- log(wbc, b = 10)
grupo <- c(rep(0,17),rep(1,16))
dados <- data.frame(id = 1:33, tempo, evento, wbc, lwbc, grupo)
dados
table(tempo)

## estimador de kaplan-meier para grupos
(km.g <- survfit(Surv(tempo, evento) ~ grupo))
plot(km.g, col = c("blue","red"), lwd = 2)
plot(km.g, col = c("blue","red"), lwd = 2, fun = "cloglog")

## modelo de riscos proporcionais completo
m5 <- coxph(Surv(tempo, evento) ~ lwbc*grupo, dados, ties = "exact")
summary(m5)

m5b <- coxph(Surv(tempo, evento) ~ lwbc*grupo, dados, ties = "efron")
summary(m5b)

m5c <- coxph(Surv(tempo, evento) ~ lwbc*grupo, dados, ties = "breslow")
summary(m5c)

res.m <- residuals(m5)
plot(res.m)
plot(lwbc, res.m)

res.cs <- evento - res.m
plot(res.cs)
plot(tempo,res.cs)

km.cs <- survfit(Surv(res.cs, evento) ~ 1)
plot(km.cs)
curve(exp(-x), 0, 4, lwd = 2, col = "red", add = T)

cox.zph(m5b, transform = "identity")

cox.zph(m5c, transform = "identity")

## testar efeito da interação
m6 <- coxph(Surv(tempo, evento) ~ lwbc + grupo, dados, ties = "exact")
summary(m6)
anova(m5, m6)

## Como interpretar a intera??o?



## Obtendo as estimativas do risco acumulado e sobrevivência
## Qual o grupo de refer?ncia? lwbc = 0 e grupo = 0
Lambda0 <- basehaz(m5, centered = F)
S0 <- exp(-Lambda0$hazard)
summary(lwbc)

## Contruindo as curvas de sobrevivência para comparação dos grupos
## para lwbc = média de lwbc

plot(S0^exp(mean(lwbc)*1.689), type = "s", col = "blue", lwd = 2,
     xlab = "tempo", ylab = "S(t)", font.lab = 2)
points(S0^exp(mean(lwbc)*1.689 + 6.763 + mean(lwbc)*(-1.3357)), 
       type = "s", col = "red", lwd = 2)
legend("topright", c("Grupo 1", "Grupo 2"), col = c("blue","red"), 
       lwd = 2)

## Contruindo as curvas de sobrevivência para comparação de dois 
## niveis de lwbc (4, 5) em cada grupo

plot(S0^exp(4*1.689), type = "s", col = "blue", lwd = 2,
     xlab = "tempo", ylab = "S(t)", font.lab = 2,
     main = "Grupo 1")
points(S0^exp(5*1.589), 
       type = "s", col = "red", lwd = 2)
legend("topright", c("lwbc = 4", "lwbc = 5"), col = c("blue","red"), 
       lwd = 2)

plot(S0^exp(4*1.689 + 6.763 + 4*(-1.335)), type = "s", col = "blue", lwd = 2,
     xlab = "tempo", ylab = "S(t)", font.lab = 2,
     main = "Grupo 2")
points(S0^exp(5*1.689 + 6.763 + 5*(-1.335)), 
       type = "s", col = "red", lwd = 2)
legend("topright", c("lwbc = 4", "lwbc = 5"), col = c("blue","red"), 
       lwd = 2)




