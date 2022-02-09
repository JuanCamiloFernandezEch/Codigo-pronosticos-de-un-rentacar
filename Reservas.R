library(readxl)
library(tseries)
library(forecast)
library(DescTools)
library(TSA)
library(lmtest)
library(car)



library(graphics)
library(aTSA)#Raices Unitarias
library(uroot)#Raices Unitarias

A <- read_excel("C:/Users/jfernandeze/OneDrive - Renting Colombia S.A/JuanCamilo - Personal/Especializacion/Seminario/HistReservas.xlsx")
#A <- read.table("HistReservas", header = TRUE, stringsAsFactors = False)
attach(A)

n <- length(CantReservas)
m <- 30

Fechas <- as.Date(Fecha, format="%y-%m-%d")
ejex.mes <- seq(Fechas[1], Fechas[n], "months")
ejex.año <- seq(Fechas[1], Fechas[n], "years")

plot (Fechas, CantReservas, xaxt="n", panel.first = grid(), type = 'l', ylab = 'ReservasDiarias')
axis.Date(1, at =ejex.mes, format="%m/%y")
axis.Date(1, at=ejex.año, labels = FALSE, tcl = -0.2)

r <- CantReservas
yi <- ts(r[1:(n-m)], frequency = 7)
yf <- ts(r[(n-m+1):n], frequency = 7)
ni <- length(yi)

boxplot(r, ylab="Cantidad Reservas", outpch=30)
summary(r)

#Descomposicion stl
y.stl <- stl(yi,"per")
plot(y.stl)
St.stl <- y.stl$time.series[,1]
Tt.stl <- y.stl$time.series[,2]
Et.stl <- y.stl$time.series[,3]
y1 <- St.stl+Et.stl
ts.plot(y1)

#FAC y FAC-Parcial
par(mfrow=c(1,2))
TSA::acf(yi, lag.max = 60, ci.type = "ma", drop.lag.0 = TRUE)
pacf(yi,60,main="")

#ljung-box
Box.test(x = yi, lag = 60, type="Ljung-Box")

#Raices Unitarias: ordinarias y estaiconales
adf.test(yi)
TestCH <- ch.test(y1,type = "trigonometric", pvalue='raw') #detecta raiz estaiconal
TestHEGY <- hegy.test(yi,deterministic = c(1,1,0), lag.method = "fixed", maxlag = 1)
TestHEGY <- hegy.test(y1,deterministic = c(1,1,1), lag.method = "fixed", maxlag = 1)
TestOCSB <- ocsb.test(yi)

#Modelos lineales
t <- seq(1,length(yi))
t2 <- t*t
t3 <- t*t*t
lyi <- log(yi)
it = seasonaldummy(yi)

mod.lin <- lm(yi ~ t + it)
mod.cuad <- lm(yi ~ t + t2 + it)
mod.cub <- lm(yi ~ t + t2 + t3 + it)

summary(mod.lin)
summary(mod.cuad)
summary(mod.cub)

yhat.lin = mod.lin$fitted.values
yhat.cuad = mod.cuad$fitted.values
yhat.cub = mod.cub$fitted.values

#Modelo BSM
mod.bsm <- StructTS(yi, type = "BSM")
print(mod.bsm$coef)
yhat.bsm = mod.bsm$fitted[,1]+mod.bsm$fitted[,2]+mod.bsm$fitted[,3]


#Modelo ETS(AAA)
modets <- ets(yi, model = "AAA", damped = FALSE)
summary(modets)
yhat.ets <- modets$fitted

#Modelo Neuronal autoregresivo (NNAR)
modnnar <- nnetar(yi,lambda=0)
print(modnnar)
str(modnnar)
yhat.nnar= fitted(modnnar)

#Modelo ARMA-ARIMA-SARIMA
auto.arima(yi)
par(mfrow=c(1,1))
parma<-armasubsets (y=yi, nar=21, nma=21, y.name = "Y", ar.method = "ols")
plot(parma)


modSARIMA <- arima(yi,order=c(0,0,2),seasonal = list(order=c(0,1,2),period=7))
coeftest(modSARIMA)
yhat.SARIMA <- modSARIMA$fitted
yhat.SARIMA1 <- fitted(modSARIMA)


#Gráfico comparativo ajuste entrenamiento 
plot(t,yi,type='l',col='darkgray', ylab = "Cantidad Reservas", xlab='Registro', ylim=c(0,1600))
lines(t,yhat.lin,col='firebrick3')
lines(t,yhat.cuad,col='darkcyan')
lines(t,yhat.cub,col='orange2')
lines(t,yhat.bsm,col='Green4')
lines(t,yhat.ets, col='magenta3')
lines(t,yhat.nnar, col='blue3')
lines(t,yhat.SARIMA1, col='aquamarine4')

legend("topleft", 
       c("Observaciones","lineal","cuadrático","cúbico","BSM","ETS(AAA)","NNAR", "SARIMA"), 
       lty = c(1,1),lwd=c(2,2),
       col=c('gray20','firebrick3','darkcyan','orange2', 'Green4','magenta3','blue3','aquamarine4'))



medidas <- function(m,y,k){
  # m <- objeto producido con lm()
  # y <- variable dependiente
  # k <- numero de coeficiente beta
  T <- length(y)
  yest <- fitted(m)
  sse <- sum((yest - y)^2)
  ssr <- sum((y - mean(y))^2)
  mse <- sse/(T-k)
  R2 <- 1 - sse/ssr
  Ra2 <- 1- (T-1)*(1-R2)/(T-k)
  aic <- log((T-k)*exp(2*k/T)*mse/T)
  bic <- log(T^(k/T)*(T-k)*mse/T)
  M <- c(sqrt(mse), Ra2, aic, bic)
  names(M) <- c("rmse", "R2-adj", "log.aic", "log.bic")
  return(M)
}
(mlin <- medidas(mod.lin,yi,2))
(mcuad <- medidas(mod.cuad,yi,3))
(mcub <- medidas(mod.cub,yi,4))
(mets <- medidas(modets,yi,3))
(msarima <- medidas(modSARIMA,yi,4))


medidas.struc = function(y,yest,k){
  # y = serie, m = modelo, k = numero parametros
  T = length(y)
  sse = sum((yest-y)^2)
  ssr = sum((y-mean(y))^2) 
  mse = sse/(T-k)
  R2 = 1 - sse/ssr
  Ra2 = 1 - (T-1)*(1-R2)/(T-k)
  aic = log((T-k)*exp(2*k/T)*mse/T)
  bic = log(T^(k/T)*(T-k)*mse/T)
  M = c(sqrt(mse),Ra2,  aic, bic)
  names(M) = c("rmse","R2-ad","log.aic","log.bic")
  return(M)
}
(mnnar <- rbind(medidas.struc(yi[!is.na(yhat.nnar)],na.omit(yhat.nnar),10)))
(mbsm = medidas.struc(yi,yhat.bsm,3))

#--------Residuos de los modelos
#Reg.Lineal
rlin <- mod.lin$residuals
plot(rlin)
rlin <- ts(rlin,frequency=7)

tm = seq(1,length(rlin),1)

par(mfrow=c(3,2))
plot(tm,rlin,type='l',ylab='residuo')
abline(h=0,lty=2)
plot(density(rlin),xlab='x',main= '')
acf(rlin,60,ci.type="ma",drop.lag.0=TRUE)
pacf(rlin,60,main="")
qqnorm(rlin)
qqline(rlin,col=2)
cpgram(rlin)

Box.test(x = rlin, lag = 60, type="Ljung-Box")

#Reg.cuadrática
rcuad <- mod.cuad$residuals
plot(rcuad)
rcuad <- ts(rcuad,frequency=7)

tm = seq(1,length(rcuad),1)

par(mfrow=c(3,2))
plot(tm,rcuad,type='l',ylab='residuo')
abline(h=0,lty=2)
plot(density(rcuad),xlab='x',main= '')
acf(rcuad,60,ci.type="ma",drop.lag.0=TRUE)
pacf(rcuad,60,main="")
qqnorm(rcuad)
qqline(rcuad,col=2)
cpgram(rcuad)

Box.test(x = rcuad, lag = 60, type="Ljung-Box")

#Reg.cub
rcub <- mod.cub$residuals
plot(rcub)
rcub <- ts(rcub,frequency=7)

tm = seq(1,length(rcub),1)

par(mfrow=c(3,2))
plot(tm,rcub,type='l',ylab='residuo')
abline(h=0,lty=2)
plot(density(rcub),xlab='x',main= '')
acf(rcub,60,ci.type="ma",drop.lag.0=TRUE)
pacf(rcub,60,main="")
qqnorm(rcub)
qqline(rcub,col=2)
cpgram(rcub)

Box.test(x = rcub, lag = 60, type="Ljung-Box")

#BSM
rbsm <- mod.bsm$residuals
plot(rbsm)
tm = seq(1,length(rbsm),1)

par(mfrow=c(3,2))
plot(tm,rbsm,type='l',ylab='residuo')
abline(h=0,lty=2)
plot(density(rbsm),xlab='x',main= '')
acf(rbsm,60,ci.type="ma",drop.lag.0=TRUE)
pacf(rbsm,60,main="")
qqnorm(rbsm)
qqline(rbsm,col=2)
cpgram(rbsm)

Box.test(x = rbsm, lag = 60, type="Ljung-Box")

#ETS-AAA
rets <- modets$residuals
plot(rets)
rets <- ts(rets,frequency=7)

tm = seq(1,length(rets),1)

par(mfrow=c(3,2))
plot(tm,rets,type='l',ylab='residuo')
abline(h=0,lty=2)
plot(density(rets),xlab='x',main= '')
acf(rets,60,ci.type="ma",drop.lag.0=TRUE)
pacf(rets,60,main="")
qqnorm(rets)
qqline(rets,col=2)
cpgram(rets)

Box.test(x = rets, lag = 60, type="Ljung-Box")

#Neuronal Autorregresivo
rnnal <- modnnar$residuals
plot(rnnal)
rnnal <- ts(rnnal,frequency=7)
tm = seq(1,length(rnnal),1)

par(mfrow=c(3,2))
plot(tm,rnnal,type='l',ylab='residuo')
abline(h=0,lty=2)
plot(density(rnnal),xlab='x',main= '')
acf(rnnal,60,ci.type="ma",drop.lag.0=TRUE)
pacf(rnnal,60,main="")
qqnorm(rnnal)
qqline(rnnal,col=2)
cpgram(rnnal)

Box.test(x = rnnal, lag = 60, type="Ljung-Box")

#SARIMA
rSARIMA <- modSARIMA$residuals
plot(rSARIMA)
rSARIMA <- ts(rSARIMA,frequency=7)
tm = seq(1,length(rSARIMA),1)

par(mfrow=c(3,2))
plot(tm,rSARIMA,type='l',ylab='residuo')
abline(h=0,lty=2)
plot(density(rSARIMA),xlab='x',main= '')
acf(rSARIMA,60,ci.type="ma",drop.lag.0=TRUE)
pacf(rSARIMA,60,main="")
qqnorm(rSARIMA)
qqline(rSARIMA,col=2)
cpgram(rSARIMA)

Box.test(x = rSARIMA, lag = 60, type="Ljung-Box")

#Pronósticos
tt <- seq((n-m+1),n,1)
tt2 <- tt*tt
tt3 <- tt*tt*tt
itp <- seasonaldummy(yi,m)

pr.lin <- predict(mod.lin,data.frame(t=tt, it=I(itp)))
prf.lin <- ts(pr.lin, frequency=7, start=c(2021,11))

pr.cuad <- predict(mod.cuad,data.frame(t=tt, t2=tt2, it=I(itp)))
prf.cuad <- ts(pr.cuad, frequency=7, start=c(2021,11))

pr.cub <- predict(mod.cub,data.frame(t=tt, t2=tt2, t3=tt3, it=I(itp)))
prf.cub <- ts(pr.cub, frequency=7, start=c(2021,11))

pr.bsm <- forecast(mod.bsm,m)
prf.bsm <- ts(pr.bsm$mean,frequency=7)

pr.ets <- forecast(modets,m)$mean
prf.ets <- ts(pr.ets,frequency=7)

pr.nnal = forecast(modnnar,m)$mean
prf.nnal <- ts(pr.nnal,frequency=7)

pr.SARIMA <- predict(modSARIMA,n.ahead=m)$pred
prf.SARIMA <- ts(pr.SARIMA, frequency=7)

ff <- seq(as.Date("2021/11/01"), length.out=m, by="day")

par(mfrow=c(1,1))

plot(ff, yf, type = 'o', ylab='Reservas',lwd = 3,col='gray20',ylim=c(300,1200))
lines (ff,prf.lin, col='firebrick3',lwd=2,lty=1 )
lines(ff,prf.cuad, col='darkcyan',lwd=2,lty=1 )
lines(ff,prf.cub, col='orange2',lwd=2,lty=1 )
legend("topright", 
       c("Observaciones","lineal","cuadrático","cúbico"), 
       lty = c(1,1),lwd=c(2,2),
       col=c('gray20','firebrick3','darkcyan','orange2'))

plot(ff, yf, type = 'o', ylab='Reservas',lwd = 3,col='gray20',ylim=c(300,1100))
lines(ff,prf.bsm, col='Green4',lwd=2,lty=1)
legend("topright", c("Observaciones","BSM"),  lty = c(1,1),lwd=c(2,2),col=c('gray20','Green4'))

plot(ff, yf, type = 'o', ylab='Reservas',lwd = 3,col='gray20',ylim=c(300,1100))
lines(ff,prf.ets, col='magenta3', lwd=2, lty=1)
legend("topright", c("Observaciones","ETS(AAA)"),  lty = c(1,1),lwd=c(2,2),col=c('gray20','magenta3'))

plot(ff, yf, type = 'o', ylab='Reservas',lwd = 3,col='gray20',ylim=c(300,1100))
lines(ff,prf.nnal, col='blue3', lwd=2, lty=1)
legend("topright", c("Observaciones","NNAR"), lty = c(1,1),lwd=c(2,2), col=c('gray20','blue3'))

plot(ff, yf, type = 'o', ylab='Reservas',lwd = 3,col='gray20',ylim=c(300,1100))
lines(ff,pr.SARIMA, col='aquamarine4', lwd=2, lty=1)
legend("topright", c("Observaciones","SARIMA"), lty = c(1,1),lwd=c(2,2), col=c('gray20','aquamarine4'))


#Estadisticos del pronóstico
A <- rbind(accuracy(pr.lin, yf), accuracy(pr.cuad, yf), accuracy(pr.cub, yf),
           accuracy(pr.bsm, yf),accuracy(pr.ets, yf),accuracy(pr.SARIMA, yf),accuracy(pr.nnal, yf) )
rownames(A) <- c("lineal", "cuadrático", "cúbico", "BSM","ETS(AAA)","SARIMA","NNAR")
Utheil <- c(Utheil <- TheilU(yf,prf.lin),Utheil <- TheilU(yf,prf.cuad),Utheil <- TheilU(yf,prf.cub),
            Utheil <- TheilU(yf,prf.bsm),Utheil <- TheilU(yf,prf.ets),Utheil <- TheilU(yf,prf.SARIMA),
            Utheil <- TheilU(yf,prf.nnal))
W <- cbind(A,Utheil)
(W)


#Estadisticos del pronóstico
accuracy(pr.lin, yf)
accuracy(pr.cuad, yf)
accuracy(pr.cub, yf)
accuracy(prf.bsm, yf)
accuracy(prf.ets, yf)
accuracy(prf.SARIMA, yf)
accuracy(prf.nnal, yf)

