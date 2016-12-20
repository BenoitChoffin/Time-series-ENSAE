library(urca)
library(apt)
library(tseries)
library(foreign)
library(forecast)
library(fUnitRoots)

#PARTIE 1
#Extraction d’hydrocarbures
data <- read.csv("time_series.csv")
#transformation en objet "time series"
datats <- ts(data, frequency=12, start=c(1990,1))
#permutation de la serie, sinon elle est dans le mauvais sens
hycarb <- rev(datats)
hycarbts <- ts(hycarb, frequency=12, start=c(1990,1))

#Représentation graphique de la série brute
plot.ts(hycarbts, type="l",lwd=1, col="red", xlab="Time", ylab="Indice")
#on prend le logarithme de la série pour écraser les pics de saisonnalité
loghycarb <- log(hycarbts)
loghycarbts <- ts(loghycarb, frequency = 12, start = c(1990,1))
#Représentation graphique de la série log-transformée
plot.ts(loghycarbts, type="l",lwd=1, col="red", xlab="Time", ylab="Indice")

#désaisonnalisation de la série loghycarbts
df_desais <- decompose(loghycarbts)
plot(df_desais)
hycarb_desais <- loghycarbts - df_desais$seasonal

#Représentation graphique de la série désaisonnalisée
plot.ts(hycarb_desais,type="l",lwd=1, col="red", xlab="Time", ylab="Indice")

#Test ADF sur la série désaisonnalisée
urdfTest(hycarb_desais)

#création de la série différenciée
hycarbdiff1 <- diff(hycarb_desais, differences=1)
#Représentation graphique de la série différenciée
plot.ts(hycarbdiff1,type="l",lwd=1, col="red", xlab="Time", ylab="Indice")

#Test ADF pour déterminer la stationnarité de la série différenciée
urdfTest(hycarbdiff1)
#on obtient bien une stationnarité de la série différenciée (seuil de 1%)


#PARTIE 2
#Détermination des pmax et qmax à partir des ACF/PACF
acf(hycarbdiff1,lag.max=20) #a priori qmax = 2
pacf(hycarbdiff1,lag.max=20) #a priori pmax = 6

#Mise en place d’une stratégie ascendante :
#Estimation des modèles (1,1), (2,2), (3,3), (4,4), (5,5) et (6,6) : minimisation du critère AIC
mod1 <- Arima(hycarbdiff1, order=c(1,0,1)) #AIC = -827.55
mod2 <- Arima(hycarbdiff1, order = c(2,0,2)) #AIC = -828.84
mod3 <- Arima(hycarbdiff1, order = c(3,0,3)) #AIC = -824.87
mod4 <- Arima(hycarbdiff1, order = c(4,0,4)) #AIC = -826.5
mod5 <- Arima(hycarbdiff1, order = c(5,0,5)) #AIC = -821.15
mod6 <- Arima(hycarbdiff1, order = c(6,0,6)) #AIC = -821.6
#on choisit un ARMA(2,2) pour commencer, car il minimise l’AIC entre les 6 modèles
#comme l’AIC privilégie les gros modèles, on va essayer d’enlever des
#paramètres grâce aux tests sur les paramètres
abs(mod2$coef)/sqrt(diag(mod2$var.coef))
#on peut supprimer deux coefficients (pas en même temps bien entendu)
#on teste donc les modèles ARMA(1,2) et ARMA(2,1)

mod12 <- Arima(hycarbdiff1, order=c(1,0,2)) #AIC = -830.12
mod21 <- Arima(hycarbdiff1, order=c(2,0,1)) #AIC = -829.02
#on sélectionne donc un ARMA(1,2) car il minimise l’AIC
#Testons maintenant la blancheur des résidus de ce modèle
tsdiag(mod12, gof.lag=30)
#PARTIE 3
#Prévision
hycarbforecasts <- forecast.Arima(mod12,h=2,level=c(95))
plot.forecast(hycarbforecasts,shadecols = "orange",fcol="red",
              main = "Intervalles de confiance à 95% pour T+1 et T+2")
