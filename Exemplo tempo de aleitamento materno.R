## Análise de dados de aleitamento materno
## 150 mães de crianças menores de 2 anos de idade
## tempo até o desmame
## indicadora de desmame
## Covariaveis
## V1 : Experiencia anterior (0 sim 1 não) - expe
## v2 : Número de filhos (0 <= 2 1 > 2) - nfil
## V3 : Conceito materno sobre o tempo ideal de amamentação
##      (0 > 6m e 1 <= 6m) - conc
## V4 : Dificuldade para amamentar nos primeiros dias
##      (0 não 1 sim ) - difi
## V5 : tipo servico pre natal (0 publico 1 privada) - serv
## V6 : recebeu exclus leite materno na maternidade 
##     (0 sim 1 não) exclus
## V7 : a criança teve contato com o pai
##      (0 sim 1 não) contpai
## V8 : renda per capita (0 >= 1SM, 1 < 1 SM ) renda
## V9 : peso nascimento (0 >= 2.5kg 1 < 2.5kg) peso
## V10 : tempo de serparacao mae-filho pós parto
##    (0 <= 6h 1 > 6 hs) separ
## V11 : Permanencia no bercario (0 não 1 sim) berc

var.nomes <- c("expe", "nfil", "conc","dific","serv",
               "exclus","contpai","renda","peso","separ",
               "berc")

setwd("G:/Meu Drive/AtividadesAcademicas/CursoEncPos")

dados <- read.table("desmame.txt", h = T)

str(dados)

dados <- dados[,c("id","tempo","cens", paste("V",1:11,sep=""))]

str(dados)

names(dados)[4:14] <- var.nomes

str(dados)

require(survival)

## análise descritiva

table(dados$cens)

table(dados$tempo)

surv.d <- with(dados, Surv(tempo, cens))

km.total <- survfit(surv.d ~ 1)

km.total

plot(km.total, lwd = 2)

for(i in var.nomes){
  print(i)
  print(table(dados[,i]))
  km <- survfit(surv.d ~ dados[,i])
  print(km)
  logrank <- survdiff(surv.d ~ dados[,i], rho = 0)
  print(logrank)
  peto <- survdiff(surv.d ~ dados[,i], rho = 1)
  print(peto)
  par(mfrow = c(1,2))
  plot(km, col = 1:2, lwd = 2, main = i)
  plot(km, col = 1:2, lwd = 2, main = i, fun = "cloglog")
  print("###########################")
}  

## Alguns modelos paramétricos
## Weibull

m.wei <- survreg(surv.d ~ expe + nfil + conc + dific + serv + exclus + 
                          contpai + renda + peso + separ + berc, 
                 dist = "weibul", dados)

l.wei <- exp(-1*predict(m.wei, type = "linear"))
(g.wei <- 1 / m.wei$scale)
res.wei <- (l.wei*dados$tempo)^{g.wei}

km.res.wei <- survfit(Surv(res.wei, dados$cens) ~ 1)
par(mfrow = c(1,2))
plot(km.res.wei, main = "Modelo Weibull")
curve(exp(-x), 0, 3, add = T, col = "red", lwd = 2)
plot(km.res.wei$time, -log(km.res.wei$surv))
abline(a = 0, b = 1, col = "red", lwd = 2)
lm(-log(km.res.wei$surv) ~ km.res.wei$time - 1)

## Exponencial
m.exp <- survreg(surv.d ~ expe + nfil + conc + dific + serv + exclus + 
                          contpai + renda + peso + separ + berc, dist = "exp", 
                   dados)
anova(m.exp, m.wei)

l.exp <- exp(-1*predict(m.exp, type = "linear"))
res.exp <- l.exp*dados$tempo

km.res.exp <- survfit(Surv(res.exp, dados$cens) ~ 1)
par(mfrow = c(1,2))
plot(km.res.exp, main = "Modelo Exponencial")
curve(exp(-x), 0, 3, add = T, col = "red", lwd = 2)
plot(km.res.exp$time, -log(km.res.exp$surv))
abline(a = 0, b = 1, col = "red", lwd = 2)
lm(-log(km.res.exp$surv) ~ km.res.exp$time - 1)

## Log-logistico
m.ll <- survreg(surv.d ~ expe + nfil + conc + dific + serv + exclus + 
                          contpai + renda + peso + separ + berc, 
                      dist = "loglogistic", dados) 

l.ll <- exp(-1*predict(m.ll, type = "linear"))
(g.ll <- 1 / m.ll$scale)
res.ll <- -log(1/(1+(l.ll*dados$tempo)^g.ll))

km.res.ll <- survfit(Surv(res.ll, dados$cens) ~ 1)
par(mfrow = c(1,2))
plot(km.res.ll, main = "Modelo Log-Logístico")
curve(exp(-x), 0, 3, add = T, col = "red", lwd = 2)
plot(km.res.ll$time, -log(km.res.ll$surv))
abline(a = 0, b = 1, col = "red", lwd = 2)
lm(-log(km.res.ll$surv) ~ km.res.ll$time - 1)

## Log-normal
m.ln <- survreg(surv.d ~ expe + nfil + conc + dific + serv + exclus + 
                          contpai + renda + peso + separ + berc, 
                      dist = "lognormal", dados) 

l.ln <- exp(-1*predict(m.ln, type = "linear"))
(g.ln <- 1 / m.ln$scale)
res.ln <- -log(1-pnorm(g.ln*log(l.ln*dados$tempo)))

km.res.ln <- survfit(Surv(res.ln, dados$cens) ~ 1)

par(mfrow = c(1,2))
plot(km.res.ln, main = "Modelo Log-Normal")
curve(exp(-x), 0, 3, add = T, col = "red", lwd = 2)
plot(km.res.ln$time, -log(km.res.ln$surv))
abline(a = 0, b = 1, col = "red", lwd = 2)
lm(-log(km.res.ln$surv[km.res.ln$surv>0]) ~ km.res.ln$time[km.res.ln$surv>0] - 1)

## Modelo de Cox
m.cox <- coxph(surv.d ~ expe + nfil + conc + dific + serv + exclus + 
                          contpai + renda + peso + separ + berc, dados)

res.cox <- dados$cens - residuals(m.cox, type = "martingal")

km.res.cox <- survfit(Surv(res.cox, dados$cens) ~ 1)

par(mfrow = c(1,2))
plot(km.res.cox, main = "Modelo de Cox")
curve(exp(-x), 0, 3, add = T, col = "red", lwd = 2)
plot(km.res.cox$time, -log(km.res.cox$surv))
abline(a = 0, b = 1, col = "red", lwd = 2)
lm(-log(km.res.cox$surv) ~ km.res.cox$time - 1)

## schoenfeld residuals
(schoen <- cox.zph(m.cox, "identity"))
plot(schoen)

#install.packages("timereg")
require(timereg)
tm.cox <- timecox(surv.d ~ expe + nfil + conc + dific + serv + exclus + 
                    contpai + renda + peso + separ + berc, dados, 
                    max.time = 18)
summary(tm.cox)
plot(tm.cox)

## Selecao de variáveis
summary(m.cox)
#summary(m.ln)

#
m.expe <- coxph(surv.d ~ nfil + conc + dific + serv + exclus + 
                  contpai + renda + peso + separ + berc, dados)
anova(m.expe, m.cox)

m.nfil <- coxph(surv.d ~ expe + conc + dific + serv + exclus + 
                  contpai + renda + peso + separ + berc, dados)
anova(m.nfil, m.cox)

m.conc <- coxph(surv.d ~ expe + nfil + dific + serv + exclus + 
                  contpai + renda + peso + separ + berc, dados)
anova(m.conc, m.cox)

m.difi <- coxph(surv.d ~ expe + nfil + conc + serv + exclus + 
                  contpai + renda + peso + separ + berc, dados)
anova(m.difi, m.cox)

m.serv <- coxph(surv.d ~ expe + nfil + conc + dific + exclus + 
                  contpai + renda + peso + separ + berc, dados)
anova(m.serv, m.cox)

m.exclus <- coxph(surv.d ~ expe + nfil + conc + dific + serv + 
                  contpai + renda + peso + separ + berc, dados)
anova(m.exclus, m.cox)

m.contpai <- coxph(surv.d ~ expe + nfil + conc + dific + serv + 
                   exclus + renda + peso + separ + berc, dados)
logLik(m.contpai, m.cox)

m.renda <- coxph(surv.d ~ expe + nfil + conc + dific + serv + 
                  exclus + contpai + peso + separ + berc, dados)
anova(m.renda, m.cox)

m.peso <- coxph(surv.d ~ expe + nfil + conc + dific + serv + 
                   exclus + contpai + renda + separ + berc, dados)
anova(m.peso, m.cox)

m.separ <- coxph(surv.d ~ expe + nfil + conc + dific + serv + 
                  exclus + contpai + renda + peso + berc, dados)
anova(m.separ, m.cox)

m.berc <- coxph(surv.d ~ expe + nfil + conc + dific + serv + 
                   exclus + contpai + renda + peso + separ, dados)
anova(m.berc, m.cox)

## modelo cox reduzido 
m.cox2 <- coxph(surv.d ~ expe + conc + dific + 
                  exclus + contpai + renda, dados)
anova(m.cox2, m.cox)
summary(m.cox2)

m.exp2 <- coxph(surv.d ~ conc + dific + 
                  exclus + contpai + renda, dados)
anova(m.exp2, m.cox2)

m.conc2 <- coxph(surv.d ~ expe + dific + 
                  exclus + contpai + renda, dados)
anova(m.conc2, m.cox2)

m.difi2 <- coxph(surv.d ~ expe + conc + 
                   exclus + contpai + renda, dados)
anova(m.difi2, m.cox2)

m.exclus2 <- coxph(surv.d ~ expe + conc + 
                   dific + contpai + renda, dados)
anova(m.exclus2, m.cox2)

m.contpai2 <- coxph(surv.d ~ expe + conc + 
                     dific + exclus + renda, dados)
anova(m.contpai2, m.cox2)

m.renda2 <- coxph(surv.d ~ expe + conc + 
                      dific + exclus + contpai, dados)
anova(m.renda2, m.cox2)

m.cox3 <- coxph(surv.d ~ conc + dific + exclus, dados)
anova(m.cox3, m.cox2)

m.cox4 <- coxph(surv.d ~ conc + dific + exclus + renda, dados)
anova(m.cox4, m.cox2)
anova(m.cox4, m.cox)

summary(m.cox4)

cox.zph(m.cox4, "identity")

## interpretação do modelo de cor

# variável conceito
# as mães que cujo o conceito materno sobre o tempo ideal
# de amamentação é menor do 6 meses apresentaram risco de 
# interromper a amanetação 1,78 vezes maior que as outras mães

# variável dificuldade
# as mães que disseram ter dificuldade em amamentar nos primeiros
# dias após o parto apresentararam riso de interromper a 
# amementação 2,57 maior que as mães que não tiveram dificuldades.

# variável exclusividade
# as mães que dosseram que não criança recebeu exclusimente leite
# materno na maternidade apresentaram risco 1,75 vezes maior de
# interromper o aleitamento comparadas as demais mães

# variável renda
# as mães com renda per capita inferior a 1 SM, apresentaram risco
# de interromper o aleitamento materno 1,79 vezes maior que as demais
# mães com renda percapita igual ou supeior a 1 SM

## modelo log-linear
m.ln3 <- survreg(surv.d ~ conc + dific + exclus + renda, 
                 dist = "lognormal", dados)
summary(m.ln3)
anova(m.ln3, m.ln)
g.ln3 <- 1/m.ln3$scale
summary(m.ln3)
coef.ln <- m.ln3$coefficients

#Contruindo as curvas de sobrevivência para comparação dos grupos

# obtendo a função de rico de referência / modelo de cox

L0 <- basehaz(m.cox4)
S0 <- exp(-L0$hazard)
coef.cox <- m.cox4$coefficients


S.ln <- function(x,l,g) pnorm(-g*log(l*x))

plot(L0$time, S0, type = "s", lwd = 2, ylim = c(0,1),
     main = "Conceito sobre tempo ideal de aleitamento")
curve(S.ln(x,exp(-coef.ln[1]),g.ln3), 0, 24, lty = 2, lwd = 2, add = T)
points(L0$time, S0^exp(coef.cox[1]), type = "s", lwd = 2, col = 2)
curve(S.ln(x,exp(-sum(coef.ln[c(1,2)])),g.ln3), 0, 24, 
           col = 2, lty = 2, lwd = 2, add = T)
legend("topright", c("> 6 meses"," <= 6 meses"), lty = 1, lwd = 2, 
       col = 1:2 )
legend("bottomleft", c("Modelo de Cox","Modelo log-normal"), 
       lty = 1, lwd = 2, col = 1:2 )

plot(L0$time, S0, type = "s", lwd = 2, ylim = c(0,1),
     main = "Dificuldade em amamentar nos primeiros dias")
curve(S.ln(x,exp(-coef.ln[1]),g.ln3), 0, 24, lty = 2, lwd = 2, add = T)
points(L0$time, S0^exp(coef.cox[2]), type = "s", lwd = 2, col = 2)
curve(S.ln(x,exp(-sum(coef.ln[c(1,3)])),g.ln3), 0, 24, 
      col = 2, lty = 2, lwd = 2, add = T)
legend("topright", c("não","sim"), lty = 1, lwd = 2, 
       col = 1:2 )
legend("bottomleft", c("Modelo de Cox","Modelo log-normal"), 
       lty = 1, lwd = 2, col = 1:2 )

