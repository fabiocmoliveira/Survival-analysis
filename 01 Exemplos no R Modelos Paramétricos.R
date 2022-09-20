#Exemplo : Pacientes com cancer de bexiga
# tempos de reincidência em meses de pacientes com cancer de bexiga

tempo  <- c(3,5,6,7,8,9,10,10,12,15,15,18,19,20,22,25,28,30,40,45)

evento <- c(1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0) 

cbind(tempo, evento)

require(survival)

surv.data <- Surv(tempo, evento)

km <- survfit(surv.data ~ 1)

summary(km)

plot(km)

plot(km, fun = "cloglog")

m.exp <- survreg(surv.data ~ 1, dist = 'exponential') 
summary(m.exp)
#ajusta o modelo para y = log(T)
#fdp = exp(y - \mu - exp(y - \mu))
#\mu = -log(\lambda)
#\lambda = exp(-\mu)
lambda.exp <- exp(-3.02)
lambda.exp

m.wei <- survreg(surv.data ~ 1, dist = 'weibull') 
summary(m.wei)
#ajusta o modelo para y = log(T)
#fdp = (1\sigma)exp( (y - \mu)/\sigma - exp((y - \mu)/ \sigma)
#\sigma = 1\gamma
#\mu = -log(\lambda)
#\lambda = exp(-\mu)
lambda.wei <- exp(-3.061)
lambda.wei
gamma.wei <- 1/0.648
gamma.wei

m.ln <- survreg(surv.data ~ 1, dist = 'lognorm') 
summary(m.ln)
#ajusta o modelo para y = log(T)
#Y ~ normal (\mu, \sigma^2)
#\mu = -log(\lambda)
#\sigma = 1/gamma
lambda.ln = exp(-2.717)
lambda.ln
gamma.ln <- 1/0.765
gamma.ln

## comparando as curvas de sobrevivencia

tempos <- km$time
ekm <- km$surv
eexp <- exp(-lambda.exp*tempos)
ewei <- exp(-(lambda.wei*tempos)^gamma.wei)
eln <- 1 - pnorm(gamma.ln*log(lambda.ln*tempos))

a <- data.frame(tempo = tempos, "Kaplan-Meier" = ekm, 
                Exponencial = eexp, Weibull = ewei, "Log-normal" = eln)
round(a, 3)

plot(km)
curve(exp(-lambda.exp*x), 0, 50, col = "red", lwd = 2, add = T)
curve(exp(-(lambda.wei*x)^gamma.wei), 0, 50, col = "blue", lwd = 2, add = T)
curve(exp(-(lambda.wei*x)^gamma.wei), 0, 50, col = "blue", lwd = 2, add = T)
curve(1 - pnorm(gamma.ln*log(lambda.ln*x)), 0, 50, col = "orange", lwd = 2, add = T)
legend("topright", c("Kaplan-Meier","Exponencial","Weibull","Log-normal"),
       lty = 1, lwd = 2, col = c("black","red","blue","orange"))

