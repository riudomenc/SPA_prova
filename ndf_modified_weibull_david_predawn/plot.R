# script to read and plot SPA output

figdir  = commandArgs()[5]
sde     = commandArgs()[6]
correls = commandArgs()[7]
stream  = commandArgs()[8]
library("hydroGOF")

plot.parameters = FALSE
plot.uncertainties= FALSE

gravelvol = 0.34 # value taken from Saxton equations in SPA, with gravelpc = 0.5

source("R_lib/sourceDir.txt")
sourceDir("R_lib/")

if (getwd()!="./output") setwd("./output")
exp = "Prades"
year = 2011
nyears = 1
steps = 96 # nr of timesteps per day
ndim  = 29 # nr of state variables (stocks,fluxes,parameters) in state vector
wi = 1500
he = 750

res = "daily" # plot daily or hourly values
err = 1    # plot error bars if err = 1
png = TRUE
SWC = TRUE
LWP = TRUE

# 15min time variable
timeh = ISOdate(year,1,1,hour=0)+(1:(365*24*4*nyears))*60*15
timed = ISOdate(year,1,1)+(1:(365*nyears))*60*60*24

filterd = matrix(-999., nrow=(365*nyears), ncol=(ndim+2))
filterd[,1] = 1:(365*nyears)
sdd = filterd
sapod = 1:(365*nyears)

type = factor(c("n",rep("f",7),rep("s",3),rep("f",3),"s","s","f","s",rep("n",4),"f","f","p","p","p","p","p","p","f"))
npar = sum(type=="p")

filter = read.csv("filter.csv"   , header=T)
isthereiota = (sum(match(dimnames(filter)[[2]],"iota"),na.rm=TRUE)) == 1
if (isthereiota) {
  iotacolumn = which(dimnames(filter)[[2]]=="iota")
  filter$iota = log10(filter$iota-1.)
}
sd     = read.csv("filter_sd.csv", header=T)
foo = read.csv(paste("../input/Prades_driver",substr(year,3,4),".csv", sep=""))

bar = read.csv(paste("../input/EnKF_obs_", stream, ".csv", sep=""))
sapo = bar$trans
sapo[sapo == -9999] = NA
Bsd    = read.csv("B_sd.csv"     , header=F)
B = 1:35040
spurious = read.csv("spurious.csv", header=T)
spurious[spurious==1.] = NA
setup = read.csv(paste("../input/EnKF_setup_ndf.csv", sep=""), header=F, as.is=T)
dummy = which(setup$V1 == "lo")
fr=dummy; to=dummy+1
bounds = as.numeric(c(setup$V2[fr:to],setup$V3[fr:to],setup$V4[fr:to],setup$V5[fr:to],setup$V6[fr:to],setup$V7[fr:to]))
bounds[1:2] = log10(bounds[1:2]-1.)

varnames = dimnames(filter)[2][[1]]
dimnames(filterd)[2][[1]] = varnames
dimnames(sdd)[2][[1]] = varnames

# read veg file
veg = read.csv(paste("../input/Prades_", stream, "_veg.csv", sep=""), header=F, as.is=T)
LMA = veg[which(veg$V1 == "LMA"),2]

#read MODIS EVI/NDVI time series data for simulation year
modis = read.csv(paste("../input/MODIS", year, ".csv", sep=""))
timem = ISOdate(year,1,1)+modis$doy*60*60*24

soilR = read.csv("soilR.csv", header=T)

if (correls){
if (png) png(paste(figdir, "parcor.png", sep=""), wi=wi, he=he, units="px", poi=20)
parcor = read.csv("parcor.csv", header=T)
k = 4*24*30 # filter variable, determines filter window size
par(mfrow=c(4,4), mar=c(2, 4, 2, 1), cex=1)
ylim = c(-1,1)
for (i in 2:16) {
  plot(timeh,parcor[,i],ty="l",lwd=0.4,ylab=dimnames(parcor)[[2]][i],ylim=ylim)
  abline(h=0.,lwd=3,lty=2)
  lines(timeh,filter(parcor[,i],rep(1/k,k)),col="red",lwd=4)
  }
dev.off()

if (png) png(paste(figdir, "transparcor.png", sep=""), wi=wi, he=he, units="px", poi=20)
parcor = read.csv("transparcor.csv", header=T)
k = 4*24*30 # filter variable, determines filter window size
par(mfrow=c(2,3), mar=c(2, 4, 2, 1), cex=1)
ylim = c(-1,1)
for (i in 2:7) {
  plot(timeh,parcor[,i],ty="l",lwd=0.4,ylab=dimnames(parcor)[[2]][i],ylim=ylim)
  abline(h=0.,lwd=3,lty=2)
  lines(timeh,filter(parcor[,i],rep(1/k,k)),col="red",lwd=4)
  }
dev.off()

if (png) png(paste(figdir, "capacitanceparcor.png", sep=""), wi=wi, he=he, units="px", poi=20)
parcor = read.csv("capacitanceparcor.csv", header=T)
k = 4*24*30 # filter variable, determines filter window size
par(mfrow=c(2,3), mar=c(2, 4, 2, 1), cex=1)
ylim = c(-1,1)
for (i in 2:7) {
  plot(timeh,parcor[,i],ty="l",lwd=0.4,ylab=dimnames(parcor)[[2]][i],ylim=ylim)
  abline(h=0.,lwd=3,lty=2)
  lines(timeh,filter(parcor[,i],rep(1/k,k)),col="red",lwd=4)
  }
dev.off()
}

# create daily values
sapod = rebin(sapo,steps,"sum",NArm=TRUE) # do this only once
for (d in 2:dim(filter)[2]){ # cycle over all variables in filter
                                        # ... for fluxes (sum of all quarter-hourly values per day)

  if (type[d]=="f"){
    filterd[,d] = rebin(filter[,d],steps,"sum",NArm=TRUE)
    sdd[,d]     = rebin(sd[,d]    ,steps,"sum",NArm=TRUE)
  }
  if (type[d]=="s" | type[d]=="p"){
    filterd[,d] = rebin(filter[,d],steps,"mean",NArm=TRUE)
    sdd[,d]     = rebin(sd[,d]    ,steps,"mean",NArm=TRUE)
  }
}

filterd = as.data.frame(filterd)
sdd = as.data.frame(sdd)

if (res == "daily"){
  time = timed
  vars = filterd
  sdev = sdd
} else{
  time = timeh
  vars = filter
  sdev = sd
}

# plot flux variables
n = sum(type == "f")
foo = pnum(n+1, n)
if (png) png(paste(figdir, "fluxes.png", sep=""), wi=wi, he=he, units="px", poi=20)
par(mfrow=c(4,4), mar=c(2, 4, 2, 1), cex=.5)
for (i in 2:dim(vars)[2]){
  ylim = range(c(vars[,i] - sdev[,i] * err, vars[,i] + sdev[,i] * err))
  if (type[i]=="f" ) {
    plot(time, vars[,i], pch=20, cex=0.5, xlab="", ylab=dimnames(vars)[[2]][i], ty="n", ylim=ylim)
    if (err == 1) arrows(time, vars[,i] - sdev[,i], time, vars[,i] + sdev[,i], col="grey", angle=90., length = 0., code=3, lwd=2.)
    points(time, vars[,i], pch=20)
  }
}
if (png) dev.off()

# plot stock variables
if (png) png(paste(figdir, "stocks.png", sep=""), wi=wi, he=he, units="px", poi=20)
par(mfrow=c(2,3), mar=c(2, 4, 2, 1), cex=1)
for (i in 2:dim(vars)[2]){
  ylim = range(c(vars[,i] - sdev[,i] * err, vars[,i] + sdev[,i] * err))
  if (type[i]=="s" ) {
    plot(time, vars[,i], pch=20, cex=0.5, xlab="", ylab=dimnames(vars)[[2]][i], ty="n", ylim=ylim)
    if (err == 1) arrows(time, vars[,i] - sdev[,i], time, vars[,i] + sdev[,i], col="grey", angle=90., length = 0., code=3, lwd=2.)
    lines(time, vars[,i])
  }
}
if (png) dev.off()

# plot parameters
pp = -1
ppp = 0
varnames[25] = expression(paste(italic(iota), " (stom. effic.)",sep=""))
varnames[26] = expression(italic(G)[plant])
varnames[27] = expression(capacitance)
varnames[28] = expression(italic(ratio)[R/L])
varnames[29:30] = ""

if (plot.parameters){
if (png) png(paste(figdir, "parameters.png", sep=""), wi=wi, he=he, units="px", poi=20)
par(mfrow=c(2,3), mar=c(2, 4, 2, 1), cex=1.5)
for (i in 2:dim(vars)[2]){
  #ylim = range(c(vars[,i] - sdev[,i] * err, vars[,i] + sdev[,i] * err))
  if (type[i]=="p" ) {
    pp = pp + 2
    ppp = ppp + 1
    ylim = bounds[pp:(pp+1)]
    #ylim = c(min(vars[,i] - sdev[,i]), max(vars[,i] + sdev[,i]) )
    plot(time, vars[,i], pch=20, cex=1, xlab="", ylab=varnames[i], ty="n", ylim=ylim)
    if (err == 1) arrows(time, vars[,i] - sdev[,i], time, vars[,i] + sdev[,i], col="grey", angle=90., length = 0., code=3, lwd=2.)
    lines(time, vars[,i], lwd=3)
    if (sde) abline( h = priors[ppp] , lwd=2 , lty = 2 , col="red" )
  }
}
if (png) dev.off()
}

# plot NEE and error bars
if (png) png(paste(figdir, "NEE.png", sep=""), wi=wi, he=he, units="px", poi=20)
par(cex=1)
ylim = range(c(filterd$NEE - sdd$NEE, filterd$NEE + sdd$NEE))
plot(timed, filterd$NEE, ty="n", ylim = ylim)
arrows(timed, filterd$NEE - sdd$NEE, timed, filterd$NEE + sdd$NEE, col="grey", angle=90., length = 0.05, code=3, lwd=2.)
points(timed, filterd$NEE, pch=20, cex=1., ty="b", lwd=1.5)
if (png) dev.off()

# stand transpiration (from sapflow data, in mm, daily)

sapmh   = filter$trans + filter$capacitance
transmh = filter$trans
sapmhsd = sd$trans     + sd$capacitance
sapmd   = filterd$trans + filterd$capacitance
sapmdsd = sdd$trans     + sdd$capacitance

foo = sapod
sapod = rebin(sapo, steps, "sum", NArm=FALSE)
if (png) png(paste(figdir, "trans.png", sep=""), wi=wi, he=he, units="px", poi=20)
par(cex=1)
ylim = range(c(sapod, sapmd-sapmdsd, sapmd+sapmdsd), na.rm = T); ylim[1] = 0
plot(timed, sapod, ylab="transpiration [mm]", xlab="", ylim=ylim, ty="n")
#arrows(timed, sapmd-sapmdsd, timed, sapmd+sapmdsd, col="black", angle=90., length=0.05, code=3, lwd=2.)
#lines(timed, sapmd, lty=1, lwd=2)
bfill(timed, sapmd+sapmdsd,col=grey(0.3))
bfill(timed, sapmd)
points(timed, sapod, pch=21, cex=1, bg="lightgrey", lwd=2)
items = c("observed","modelled")
if (sde) items = c("best guess","modelled")
legend("topright",items,pch=c(21,22), col=c("black","black"),
       pt.lwd=c(2,2), pt.bg=c("lightgrey","darkgrey"), pt.cex=c(1.5,1.5),
       box.lwd=1,cex=1.5, inset=0.05)
cor(sapod, sapmd, use="compl")**2
slope = round(lm(sapmd~sapod)$coefficients[2],3)
#par(new=T)
#plot(timem, modis$evi,axes=F,xlab="",ylab="",pch=2,ty="p",col="darkgreen",lwd=3)
#axis(4,col="darkgreen")
par(fig=c(0.05,0.35,0.45,0.95), new=T, cex=1)
plot(sapod, sapmd, pch=20, axes=T, ylab="", xlab="", cex=0.5, bty="o")
abline(0, 1, lty=2, lwd=2)
fit = lm(sapmd ~ sapod)
abline(fit$coefficients, lwd=2)
inset = 0.05
x = range(sapod, na.rm=T); y = range(sapmd)
x = (x[2] - x[1]) * inset; y = (y[2] - y[1]) * (1 - inset)
rsquare = round(cor(sapod, sapmd, use="compl")**2,3)
text(x,y, substitute(R^2 == rsquare, list(rsquare=rsquare)), adj=c(0,0.5),cex=.8)
inset = 0.2
y = range(sapmd); y = (y[2] - y[1]) * (1 - inset)
text(x,y, paste("slope = ", slope, sep=""), adj=c(0,0.5),cex=.8)
if (png) dev.off()
sapod = foo

# plot cumulative transpiration
if (png) png(paste(figdir, "cumulative_trans.png", sep=""), wi=wi, he=he, units="px", poi=20)
foo = sapo #; foo[is.na(foo)] = 0.
bar = sapmh
barsd = sapmhsd
bar[is.na(foo)] = 0.
barsd[is.na(foo)] = 0.
foo[is.na(foo)] = 0.
if (sde) bar[foo < 0.01] = 0.
if (sde) foo[foo < 0.01] = 0.
#for (i in 1:nyears){
#  start = 1+365*(i-1)
#  end = start+364
#  bar[start:end] = cumsum(bar[start:end])
#}
ylim = range(c(0,cumsum(foo),cumsum(bar),cumsum(bar)+cumsum(barsd),cumsum(bar)-cumsum(barsd)))
plot(timeh, cumsum(foo), ty="l", lwd=3, lty=2, ylab="transpiration [mm]",ylim=ylim,xlab="")
lines(timeh, cumsum(bar), lwd=3)
lines(timeh, cumsum(bar)+cumsum(barsd), lwd=3, col="grey")
lines(timeh, cumsum(bar)-cumsum(barsd), lwd=3, col="grey")
axis(2,at=seq(0,1000,by=50),tck=1,labels=F,lwd=0.5,lty=2)
items = c("observed","modelled")
if (sde) items = c("best guess","modelled")
legend("topleft",items,lty=c(2,1),lwd=3,cex=1.5,inset=0.05)
if (png) dev.off()

# uncertainties
if (plot.uncertainties){
ylab="residual uncertainties [obs. - mod.]"
if (sde) ylab="residual uncertainties [RDE - mod.]"
if (png) png(paste(figdir, "uncertainties.png", sep=""), wi=wi, he=he, units="px", poi=20)
ind = which(setup$V1 == "obssde")
sapodsd = as.numeric(setup$V2[ind])*steps
plot(timed, sapodsd - sapmdsd, pch = 20, cex=0.75, xlab="", ylab=ylab)
abline(0,0, lty=2)
if (png) dev.off()
}

# residuals
foo = sapod
sapod = rebin(sapo, steps, "sum", NArm=FALSE)
ylab="obs. - mod."
if (sde) ylab="RDE - mod."
if (png) png(paste(figdir, "residuals.png", sep=""), wi=wi, he=he, units="px", poi=20)
ylim = c(-max(c(sapod,sapmd),na.rm=T),max(c(sapod,sapmd),na.rm=T))
par(cex=1)
plot(timed, sapod - sapmd, xlab = "", ylab=ylab, pch=20, ylim=ylim)
rmse = sqrt(mean((sapod-sapmd)**2, na.rm=T))
inset = 0.05
x = range(timed, na.rm=T); y = ylim
x = (x[2] - x[1]) * inset; y = y[2] - ((y[2] - y[1]) * inset)
meanError = round(mean(sapod-sapmd, na.rm = T),2)
text(timed[18], y, paste("RMSE = ", round(rmse,3), ", mean Error = ", meanError, sep=""), pos=4)
abline(0,0, lty=2, lwd=2)
abline(meanError, 0, lty=2)
if (png) dev.off()
sapod = foo

# stand transpiration (hourly), selected days

if (png) png(paste(figdir, "trans_hourly.png", sep=""), wi=wi, he=he, units="px", poi=20)
par(cex=1)
#start = 96*148; end = 96*155
start = 96*90; end = 96*97
foo = sapo[start:end]; bar = sapmh[start:end]; foobar = sapmhsd[start:end]
foot = timeh[start:end]; fooboo = transmh[start:end]
ylim = range(c(foo, bar+foobar, bar-foobar), na.rm = T); ylim[1] = 0
plot(foot, foo, ylab="transpiration [mm]", xlab="", ylim=ylim, ty="n")
items = c("observed","modelled")
if (sde) items = c("best guess","modelled")
legend("top",items,pch=c(21,22), col=c("black","black"),
       pt.lwd=c(2,2), pt.bg=c("lightgrey","darkgrey"), pt.cex=c(1.5,1.5),
       box.lwd=1, cex=1.5, inset=0.025)
bfill(foot,bar+foobar, col=grey(0.3), lwd=2)
bfill(foot, bar)
points(foot, foo, pch=21, cex=1, bg="lightgrey",lwd=2)
cor(foo, bar, use="compl")**2
if (png) dev.off()

# stand transpiration (hourly), whole year

if (png) png(paste(figdir, "trans_hourly_all.png", sep=""), wi=wi, he=he, units="px", poi=20)
par(cex=1)
foo = sapo; bar = sapmh; foobar = sapmhsd
foot = timeh; fooboo = transmh
ylim = range(c(foo, bar+foobar, bar-foobar), na.rm = T); ylim[1] = 0
plot(foot, foo, ylab="transpiration [mm]", xlab="", ylim=ylim, ty="n")
items = c("observed","modelled")
if (sde) items = c("best guess","modelled")
legend("topleft",items,pch=20, col=c("black","red"),
       pt.cex=c(1.5,1.5),
       box.lwd=1, cex=1.5, inset=0.025)
#bfill(foot,bar+foobar, col=grey(0.3), lwd=2)
#bfill(foot, bar)
points(foot, foo, pch=20, cex=0.5)
points(foot, bar, pch=20, cex=0.5, col="red")
cor(foo, bar, use="compl")**2
if (png) dev.off()

# GPP
if (png) png(paste(figdir, "GPP.png", sep=""), wi=wi, he=he, units="px", poi=20)
plot(timeh, filter$G, pch=20, cex=.5)
if (png) dev.off()

# LAI
if (png) png(paste(figdir, "LAI.png", sep=""), wi=wi, he=he, units="px", poi=20)
plot(timed, filterd$Cf / LMA, ty="l", lwd=2, xlab="", ylab="LAI")
par(new=T)
seasonality = 0.75
foo = modis$evi
foo = foo / max(foo) * 1.31
foo = foo + ( ( max(foo) * seasonality ) - min(foo) ) * ( max(foo) - foo ) / ( max(foo) - min(foo) )
plot(timem, foo,axes=F,xlab="",ylab="",pch=20,ty="p",col="blue")
axis(4,col="blue")
if (png) dev.off()

# Vallcebre SWC plots

if (SWC==TRUE){
  
foo = read.csv("../input/Prades_raw_data_2011.csv")
SWCo30 = foo$Swcd_av *2.6 #/ 0.68

#bar = read.csv("~/SPA/input/Vallcebre/VallcebreDailyData2003_2005.csv", sep=";")
#SWCo60 = bar$SWC_30_60[bar$Year == 2004]
#time60 = bar$DOY[bar$Year == 2004]; time60 = time60[complete.cases(SWCo60)]; sub60 = time60; time60 = ISOdate(2004,1,1)+time60*60*60*24
#SWCo60 = SWCo60[complete.cases(SWCo60)]

foobar = read.csv("SWC.csv",header=F)
SWCm30 = apply(foobar[,2:4], 1, mean) #* (1.-gravelvol)    # each layer has 10cm thickness: 0-30cm
#SWCm30[8640:length(SWCm30)] = SWCm30[11520:length(SWCm30)]-0.1
SWCm60 = apply(foobar[,5:7], 1, mean) #* (1.-gravelvol)   # 30-60 cm

if (png) png(paste(figdir, "SWC.png", sep=""), wi=wi, he=he, units="px", poi=20)
# plot in two panels, top panel 0-30 cm, bottom panel 30-60 cm
par(mfrow=c(2,1), mai=c(1.,0.8,0.1,0.1), oma=c(0,1,0,0))
ylim=c(0,0.5)#range(c(SWCo30, SWCm30), na.rm=T)
plot(timeh, SWCo30, ty="l", ylim=ylim, ylab="SWC 0-30 cm", lwd=2)
lines(timeh, SWCm30, col="red", lwd=2)
axis(4, tcl=.5)
items = c("observed","modelled")
legend("topright", items, lty=c(1,1), col=c("black","red"))
axis(2,at=seq(0,1.,by=0.1),tck=1,labels=F,lwd=0.5,lty=2)

SWC30rmse = sqrt(mean((SWCo30-SWCm30)**2, na.rm=T))
#SWC60rmse = sqrt(mean((SWCm60[sub60*96]-SWCo60)**2))

plot(timeh, SWCm60, ty="l", ylim=ylim, ylab="SWC 30-60 cm", col="red", lwd=2)
#points(timeh[sub60*96], SWCo60, pch=20, col="green")
#points(time60, SWCo60, pch=20)
axis(4, tcl=.5)
items = c("observed","modelled")
legend("topright", items, lty=c(-1,1), col=c("black","red"), pch=c(20,-1))
if (png) dev.off()

}

if (LWP==TRUE){
  
foo = read.table("../input/LWP_LSC.csv", sep=",", as.is=TRUE, header=TRUE)

files = list.files("./",pattern="layer")
filenames = substr(files,1,7)
nfiles = length(files)
bar = read.csv("layer_1.csv")

time = as.numeric(foo$DOY[foo$Year==year])

if (stream == "ndf"){
PDo = foo$wppnd.pd; PDo = PDo[foo$Year == year]
tPDo = time[complete.cases(PDo)]; tPDo = ISOdate(year,1,1)+tPDo*60*60*24
PDo = PDo[complete.cases(PDo)]

MDo = foo$wppnd.md; MDo = MDo[foo$Year == year]
tMDo = time[complete.cases(MDo)];
tMDo = ISOdate(year,1,1)+tMDo*60*60*24
MDo = MDo[complete.cases(MDo)]
}
if (png) png(paste(figdir, "LWP.png", sep=""), wi=wi, he=he, units="px", poi=20)
par(mfrow=c(2,3), mar=c(2, 4, 2, 1), cex=1)

for (i in 1:nfiles){
bar = read.csv(paste("./",files[i], sep=""))
MDm = bar$LWP[complete.cases(match(timeh, tMDo+(2*60*60)))]
PDm = bar$LWP[complete.cases(match(timeh, tMDo-(6*60*60)))]
if (i == 1) {
  gs  = bar$stom_conduct
  MDmAV = MDm
  PDmAV = PDm
}
if (i > 1){
  gs  = gs + bar$stom_conduct #summing gs over layers
  MDmAV = MDmAV + MDm
  PDmAV = PDmAV + PDm
}
if (i == nfiles){
  gs = gs / nfiles #averaging gs over layers
  MDmAV = MDmAV / nfiles
  PDmAV = PDmAV / nfiles
  MDrmse = sqrt(mean((MDo - MDmAV)**2))
  PDrmse = sqrt(mean((PDo - PDmAV)**2))
}

cex = 1.5
xlim = c(ISOdate(year,1,1),ISOdate(year,12,31))
plot(tPDo, PDo, pch=16, ylim=c(-4.,0.), main = filenames[i], xlim = xlim,
     xlab="", ylab="leaf water potential [MPa]", cex=cex)
points(tMDo, MDo, pch=17, cex=cex)
points(tPDo, PDm, pch=21, bg="red", cex=cex)
points(tMDo, MDm, pch=24, bg="red", cex=cex)
}
plot(tPDo, PDo, ty="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
legend("top", c("PD obs.", "PD mod.", "MD obs.", "MD mod."), pch=c(16,21,17,24),cex=1.25,pt.cex=2,pt.bg="red")

if (png) dev.off()
}

gperc = quantile(gs,seq(0,1,0.01))
#print("max. stomatal conductance:")
#print(round(max(gperc,na.rm=T),0))
if (png) png(paste(figdir, "stom_cond.png", sep=""), wi=wi, he=he, units="px", poi=20)
plot(gs, pch=20)
if (png) dev.off()
