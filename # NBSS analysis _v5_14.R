# NBSS analysis
# Size Spectra analysis 
#
# for the manuscript:
#" Size niche interactions (functional replacement and competitive exclusion)
# between mero- and holoplankton shape the size spectrum 
#of tropical estuarine and marine ecosystems "
#
# total Biovolume vs individual biovolume  (NBSS)

# v. 05.01 - final plots 
# Ralf Schwamborn
# August 16, 2025


# Summary of the rationale -----------

# z  -> totalZoo Total zooplankton, including holo- and meroplankton. Includes all multicellular heterotrophic organisms found in the plankton samples (i.e., excludes algae and non-organismic particles, such as aggregates, biogenic detritus and microplastics).  
#
# md -> md = "NBSSnmd". wo_meroDec without meroplanktonic Decapoda. As above, but without meroplanktonic Decapoda (without brachyuran zoea, etc., but with holoplanktonic luciferid and sergestid adults and larvae. 
#                                                                                                  
# m  ->  m = "NBSSho. Holo_only Holoplankton only. As above, but  without meroplankton, without parasites, without  benthopelagic organisms  (e.g., without polychaete and mollusc larvae, without cirripedian nauplii, without parasitic copepods, without benthic  Cumarea and Isopods andAmpipiuds,  without meroplanktonic decapoda larvae, but with holoplanktonic luciferid and sergestid adults and larvae). 
#                                                                                                  
# These three types of NBSS were built for the three sampling areas (estuary, bay,  shelf), 
# totaling 3 x 3 = 9 NBSS built in this study.
                                                                                              

# 0. load libraries

library(scales)
library(robust)
library(permuco)

# package "permuco" ---------

library(permuco)

# ?aovperm

# example
data("emergencycost")
emergencycost$LOSc <- scale(emergencycost$LOS, scale = FALSE)
mod_cost_0 <- aovperm(cost ~ LOSc*sex*insurance, data = emergencycost, np = 2000)
mod_cost_0


# package "robust" 
# example  "lmRob"-------------------

data(stack.dat)
stack.small <- lmRob(Loss ~ Water.Temp + Acid.Conc., data = stack.dat)
anova(stack.small)

summary(stack.small)
stack.full <- lmRob(Loss ~ ., data = stack.dat)
anova(stack.full)
anova(stack.full, stack.small, test = "RWald")#  < 2.2e-16 ***,Wald-test
anova(stack.full, stack.small, test = "RF")#  9.839e-08 ***, Robustified F-test (default)

# Wald-test or robustified F-test
# Test used: Robustified F-test (default)
# Hampel,F.R.,Ronchetti,E.M.,Rousseeuw,P.J.,andStahel,W.A.(1986).Robust statistics:the approachbasedoninfluencefunctions.JohnWiley&Sons.


#library(WRS2)

# package WRS2 
# A collection of robust statistical methods basedon Wilcox' WRS functions.
# With package WRS2 (includes bootstrap)

# head(invisibility) 
#  ancova(mischief2~cloak*mischief1,data=invisibility) 
# summary(otto)
# ##specifyingcovariateevaluationpoints 
# ancova(mischief2~cloak+mischief1,data=invisibility,pts=c(3,4,8,1)) 
# ##bootstrapversion 
# ancboot(mischief2~cloak+mischief1,data=invisibility)


# 1. import data -----------------------------


# set working directory -----------------------

setwd("C:/Users/RALF/Documents/Papers/0 - Paper Denise Size Spectra")

#text1

NBSSz <- read.csv("~/Papers/0 - Paper Denise Size Spectra/NBSSz_v02.csv")
#  View(NBSSz_v01)
NBSSz <- replace (NBSSz, "NaN", NA)
NBSSz[, (5:38)] <- as.numeric(unlist(NBSSz[, (5:38)]))

NBSSm <- read.csv("~/Papers/0 - Paper Denise Size Spectra/NBSSm_v01.csv")
#  View(NBSSm)
NBSSm <- replace (NBSSm, "NaN", NA)
NBSSm[, (5:38)] <- as.numeric(unlist(NBSSm[, (5:38)]))

NBSSmd <- read.csv("~/Papers/0 - Paper Denise Size Spectra/NBSSmd_v01.csv")
#  View(NBSSmd)
NBSSmd <- replace (NBSSmd, "NaN", NA)
NBSSmd[, (5:38)] <- as.numeric(unlist(NBSSmd[, (5:38)]))


summary(NBSSmd)
head(NBSSmd)
names(NBSSmd)

# 2. subsets by region, for NBSSz, NBSSm, NBSSmd ------------------------------------

NBSSz$station <- as.character(NBSSz$station)

e_NBSSz <- NBSSz[NBSSz$station == "e",]
b_NBSSz <-  NBSSz[NBSSz$station == "b",]
s_NBSSz <-  NBSSz[NBSSz$station == "p",]

e_NBSSm <- NBSSm[NBSSm$station == "e",]
b_NBSSm <-  NBSSm[NBSSm$station == "b",]
s_NBSSm <-  NBSSm[NBSSm$station == "p",]

e_NBSSmd <- NBSSmd[NBSSmd$station == "e",]
b_NBSSmd <-  NBSSmd[NBSSmd$station == "b",]
s_NBSSmd <-  NBSSmd[NBSSmd$station == "p",]

#  View(e_NBSSmd) # OK

# 3.Analyze the slopes provided by Catarina -------------------
slope_z1 <- NBSSz$slopes
slope_m1 <- NBSSm$slopes
slope_md1 <- NBSSmd$slopes
t.test(slope_z1, slope_m1) # p-value = 0.02056
t.test(slope_z1, slope_md1) # p-value = 0.03
summary(slope_z1, slope_m1, slope_md1)
boxplot(slope_z1, slope_m1, slope_md1) # all Catarina slopes are around -2 !
# 4. calculate new, robust slopes ---------------------
# 4.1 X-axis size  vector  -----------------------
#x-axis Log biovolume bins (mm3) ------------------------------
# log10 biovolume --------------

Xvec <-   c(-2.274	,-2.098,	-1.9219,	-1.7458,	-1.5697	,-1.3937	,-1.2176	,-1.0415,	-0.86543,	-0.68935	,
            -0.51328,	-0.3372,	-0.16113,	0.014947,	0.19102,	0.3671	,0.54317	,0.71925	,0.89532,	1.0714,
            1.2475,	1.4235,	1.5996,	1.7757,	1.9518	,2.1278,	2.3039	,2.48,	2.6561,	2.8321,	3.0082,	3.1843,
            3.3604,	3.536)

diff(Xvec) # difference is 0.176

length(Xvec)

# ESD in mm -----------------------------------------------

#  d = (6V/π)^(1/3), where d is the diameter of the sphere and V is its volume.

vol.vec <-  10^Xvec

PI <- 3.141592

esd.vec.ok = ((6*vol.vec)/PI)^(1/3)

# ESDvec <-   (10^ Xvec)^(1/3)
# 
# ESDvec <-   (10^ Xvec)^(1/3)
# 
# ESDvec <-   10^ (Xvec^1/3)
# 
# ESDvec

# useful size 
diff(esd.vec.ok)
#
NBSSz <- replace (NBSSz, "NaN", NA)
NBSSz[, (5:38)] <- as.numeric(unlist(NBSSz[, (5:38)]))
Yvec <- unlist(NBSSz[(13), (5:38)])

# plot station by station
plot(Yvec ~ Xvec, main = "Bay, example sample")
lines(Yvec ~ Xvec)
points(Yvec[6:8] ~ Xvec[6:8], col = "navyblue", pch = 16)
#rect(-1.5,1.5, -0.9, 1.6, col= "grey")

# 5. Define useful size range ---------------------------

# 5.1 plot all stations with rectangles --------------------------



length(NBSSz$station) # 121 stations (all 3 areas)
Xlongvec121 <- rep(Xvec, 121)
Ylongvec121 <- as.vector (t((NBSSz[(1:121), (5:38)])))
#as.vector
length(Ylongvec121)

plot(Ylongvec121 ~ Xlongvec121) # plot ALL samples, NBSSz
#rect(-1.5,1.5, -0.95, 1.6, col= "grey") # 3 size classes
#rect(-1.5,1.7, -0.6, 1.8, col= "lightblue") # 5 classes
rect(-1.5,1.9, -0.2, 2, col= "darkorange") # 7 classes

#rect(-1.5,2.1, 0.1, 2.2, col= "purple") # 9 size classes

rect(-1.5,-1.2, 0.1, -1.4, col= gray(0.7, alpha = 0.6)) # 9 small size classes
# from -1.5 to 0.1: 9 medium size classes

rect(0.1, 0.1, 1.7, -1.4, col= alpha("darkgreen",  0.6)) # 9 large size classes
# from 0.1 to 1.7: 9 large size classes


rect(-1.5,-1.5, 0.1, -1.4, col= alpha("darkgreen", 0.6)) # 9 small size classes
# from -1.5 to 0.1: 9 medium size classes

rect(0.1, 0.1, 1.7, -1.4, col= alpha("darkgreen",  0.6)) # 9 large size classes
# from 0.1 to 1.7: 9 large size classes


# 9 size classes, from - 1.5 to 0.1 log10 (mm3)  -------
ini.mm3  <-   10^(-1.5) # 0.032 mm3  (0.032 mm3 biovol)
fin.mm3  <-   10^(0.1) # 1.07 mm3 (1.25 mm2  biovol)
ini.mm3 # 0.03162278
ini.mm3^(1/3)
ini.mm3*10
0.316227^3
fin.mm3^(1/3)
1.079^3

(ini.esd  <-   ((6*( 10^(-1.5) )/PI)^(1/3)  )) # 0.392 mm esd (392 micron ESD)
(fin.esd  <-   ((6*( 10^( .1) )/PI)^(1/3)  ))# 1.34 mm esd (1.1 mm ESD)

#  d = (6V/π)^(1/3), where d is the diameter of the sphere and V is its volume.

vol.vec <-  10^Xvec

PI <- 3.141592

esd.vec.ok = ((6*vol.vec)/PI)^(1/3)


# size range used for NBSS: -1.5 log10 mm3 to 0.1 log10 mm3 (9 size classes)
# size range a priori defined for NBSS:0.32 mm3 biovol to 1.07 mm3 biovol (9 size classes)
# size range a priori for NBSS: 316 micron ESD to 1.08 mm ESD (9 size classes)

# in practice, the size range effectively used was from -1.4 to 0.015 log10 (mm3)  -------
# in practice, the size range effectively used was  from 0.341 mm esd  to 1.07 mm esd (9 size classes)  -------

(ini.esd_IN_P  <-   ((6*( 10^(-1.5) )/PI)^(1/3)  )) # 0.341 mm esd (341 micron ESD)
(fin.esd_IN_P  <-   ((6*( 10^(0.1) )/PI)^(1/3)  )) # 1.07 mm esd (1.01 mm ESD)


(ini.esd_Large  <-   ((6*( 10^(0.1) )/PI)^(1/3)  )) # 1.339 mm esd (1.339 mm ESD)
(fin.esd_Large  <-   ((6*( 10^(1.7) )/PI)^(1/3)  ))# 4.574 mm esd (4.574 mm ESD)



Yvec <- unlist(NBSSz[(13), (5:38)])

# plot station by station
plot(Yvec ~ Xvec, main = "Bay, example sample")
lines(Yvec ~ Xvec)
points(Yvec[6:8] ~ Xvec[6:8], col = "navyblue", pch = 16)

#rect(-1.5,1.5, -0.9, 1.6, col= "grey")


# plot station by station

NBSSz[44,]

Yvec <- unlist(NBSSz[(44), (5:38)])

plot(Yvec ~ Xvec, main = "Estuary, example sample")
lines(Yvec ~ Xvec)
points(Yvec[6:8] ~ Xvec[6:8], col = "navyblue", pch = 16)

#rect(-1.5,1.5, -0.9, 1.6, col= "grey")

# plot random stations

(stat.randm <- sample( (1: length(NBSSz) ), 1) )

Yvec <- unlist(NBSSz[(stat.randm), (5:38)])

plot(Yvec ~ Xvec, main = paste( "Sation", stat.randm ))
lines(Yvec ~ Xvec)
points(Yvec[6:8] ~ Xvec[6:8], col = "navyblue", pch = 16)

rect(-1.5,1.5, -0.9, 1.6, col= "grey")




# 5.2 Define useful size range --------------------------- 

# plot all stations --------------------------
# n stations 
length(NBSSz$station) # 121 stations (all 3 areas)

Xlongvec121 <- rep(Xvec, 121)
Ylongvec121 <- as.vector (t((NBSSz[(1:121), (5:38)])))

#as.vector

length(Ylongvec121)


plot(Ylongvec121 ~ Xlongvec121) # plot ALL samples, NBSSz

#rect(-1.5,1.5, -0.95, 1.6, col= "grey") # 3 size classes

#rect(-1.5,1.7, -0.6, 1.8, col= "lightblue") # 5 classes

rect(-1.5,1.9, -0.2, 2, col= "darkorange") # 7 classes

#rect(-1.5,2.1, 0.1, 2.2, col= "purple") # 9 size classes

rect(-1.5,-1.2, 0.1, -1.4, col= gray(0.7, alpha = 0.6)) # 9 size classes

rect(-1.5,-1.5, 0.1, -1.4, col= alpha("darkgreen", 0.6)) # 9 small size classes
# from -1.5 to 0.1: 9 medium size classes

rect(0.1, 0.1, 1.7, -1.4, col= alpha("darkgreen",  0.6)) # 9 large size classes
# from 0.1 to 1.7: 9 large size classes


# BEST SIZE CLASSRS = 9, from - 1.5 to 0.1 log10 (mm3)
# 7 size classes, from - 1.5 to -0.2 log10 (mm3)  -------
# or
# 9 size classes, from - 1.5 to 0.1 log10 (mm3) BEST!!! -------

ini.mm3  <-   10^(-1.5) # 0.0316 mm3  (0.032 mm3 biovol)# 9 size classes
fin.mm3  <-   10^(0.1) # 1.25 mm3 (1.25 mm3  biovol)# 9 size classes

ini.mm3 # 0.03162278
ini.mm3^(1/3)
ini.mm3*10

0.316227^3


fin.mm3^(1/3)
1.079^3

ini.esd  <-   10^ (-1.5^1/3) # 0.31 mm esd (310 micron ESD)
fin.esd  <-   10^ (0.1 ^1/3) # 1.07 mm esd (1.1 mm ESD)

# size range used for NBSS: -1.5 log10 mm3 to 0.1 log10 mm3 (9 size classes)
# size range used for NBSS:0.32 mm3 biovol to 1.07 mm3 biovol (9 size classes)
# size range used for NBSS: 310 micron ESD to 1.1 mm ESD (9 size classes)

#points(Yvec[14:20] ~ Xvec[14:20], col = "red", pch = 16, cex = 4)



length(e_NBSSz$station) # 38 stations ( area = estuary)

# 6. LONG VECTORS Datasets for analysis ----------- 
# 6.1 LONG VECTORS by area and z,m, md -------

Xlongvec38 <- rep(Xvec, 38)

# Estuary,  LONG VECTORS for z,m, md---------------
e_Ylongvec38 <- as.vector (t((e_NBSSz[(1:38), (5:38)])))
e_Ylongvec38m <- as.vector (t((e_NBSSm[(1:38), (5:38)])))
e_Ylongvec38md <- as.vector (t((e_NBSSmd[(1:38), (5:38)])))

# Bay,  LONG VECTORS for z,m, md---------------

b_Ylongvec38 <- as.vector (t((b_NBSSz[(1:38), (5:38)])))
b_Ylongvec38m <- as.vector (t((b_NBSSm[(1:38), (5:38)])))
b_Ylongvec38md <- as.vector (t((b_NBSSmd[(1:38), (5:38)])))

# Shelf,  LONG VECTORS for z,m, md---------------

s_Ylongvec38 <- as.vector (t((s_NBSSz[(1:38), (5:38)])))
s_Ylongvec38m <- as.vector (t((s_NBSSm[(1:38), (5:38)])))
s_Ylongvec38md <- as.vector (t((s_NBSSmd[(1:38), (5:38)])))


length(e_Ylongvec38)
length(e_Ylongvec38m)
length(e_Ylongvec38md)



# 6.2 Data for analysis ----------- 
# Areas and by z,m, md -------

Xlongvec38 <- rep(Xvec, 38)

# Estuary,  LONG VECTORS for z,m, md---------------
e_Ylongvec38 <- as.vector (t((e_NBSSz[(1:38), (5:38)])))
e_Ylongvec38m <- as.vector (t((e_NBSSm[(1:38), (5:38)])))
e_Ylongvec38md <- as.vector (t((e_NBSSmd[(1:38), (5:38)])))

# Bay,  LONG VECTORS for z,m, md---------------

b_Ylongvec38 <- as.vector (t((b_NBSSz[(1:38), (5:38)])))
b_Ylongvec38m <- as.vector (t((b_NBSSm[(1:38), (5:38)])))
b_Ylongvec38md <- as.vector (t((b_NBSSmd[(1:38), (5:38)])))

# Shelf,  LONG VECTORS for z,m, md---------------

s_Ylongvec38 <- as.vector (t((s_NBSSz[(1:38), (5:38)])))
s_Ylongvec38m <- as.vector (t((s_NBSSm[(1:38), (5:38)])))
s_Ylongvec38md <- as.vector (t((s_NBSSmd[(1:38), (5:38)])))


length(e_Ylongvec38)
length(e_Ylongvec38m)
length(e_Ylongvec38md)




# 7. Many PLOTS with rectangles ---------------

# NBSSz, estuary------
plot(e_Ylongvec38 ~ Xlongvec38,  main = "NBSSz, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6),
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

rect(-1.5,1.5, -0.95, 1.6, col= "grey")

rect(-1.5,1.7, -0.6, 1.8, col= "lightblue")

rect(-1.5,1.9, -0.2, 2, col= "darkorange")

rect(-1.5,-1.2, 0.1, -1.4, col= gray(0.7, alpha = 0.6)) # 9 size classes

# NBSSm, estuary------
plot(e_Ylongvec38m ~ Xlongvec38, main = "NBSSm, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6),
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

rect(-1.5,1.5, -0.95, 1.6, col= "grey")

rect(-1.5,1.7, -0.6, 1.8, col= "lightblue")

rect(-1.5,1.9, -0.2, 2, col= "darkorange")

rect(-1.5,-1.2, 0.1, -1.4, col= gray(0.7, alpha = 0.6)) # 9 size classes

# NBSSmd, estuary------
plot(e_Ylongvec38md ~ Xlongvec38, main = "NBSSmd, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6),
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

rect(-1.5,1.5, -0.95, 1.6, col= "grey")

rect(-1.5,1.7, -0.6, 1.8, col= "lightblue")

rect(-1.5,1.9, -0.2, 2, col= "darkorange")

rect(-1.5,-1.2, 0.1, -1.4, col= gray(0.7, alpha = 0.6)) # 9 size classes


# 8. Boxplots -------------

# 8.1 simple Boxplots -------------


boxplot(e_Ylongvec38 ~ Xlongvec38, main = "NBSSz, estuary",
        cex = 1.2, col = gray(0.7, alpha = 0.6),
        ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
        xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,ylim = c(-2, 5))

boxplot(e_Ylongvec38m ~ Xlongvec38, main = "NBSSm, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6),
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)",ylim = c(-2, 5) )

boxplot(e_Ylongvec38md ~ Xlongvec38, main = "NBSSmd, estuary",
        cex = 1.2, col = gray(0.7, alpha = 0.6),
        ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
        xlab = "Indiv. Biovolume, log10(mm3 ind.-1)",ylim = c(-2, 5) )

# 8.2 two-in-one Boxplots -------------

# two-in-one boxplots
boxplot(e_Ylongvec38 ~ Xlongvec38, main = "NBSSz vs NBSSm, estuary",
        cex = 1, 
        ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
        xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,ylim = c(-2, 5),
        border = alpha("purple",  0.6), col = gray(0.8, alpha = 0.1))

boxplot(e_Ylongvec38m ~ Xlongvec38, 
        cex = 1, 
        ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
        xlab = "Indiv. Biovolume, log10(mm3 ind.-1)",ylim = c(-2, 5),
        border = alpha("black",  0.6), col = gray(0.8, alpha = 0.1),
        add = T)

# Conclusion:
# Estuary
# NMBBSm: less big ones whe we take away the meropl! Meroplankton are the big ones!


boxplot(e_Ylongvec38 ~ Xlongvec38, main = "NBSSz vs NBSSmd, estuary",
        cex = 1, 
        ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
        xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,ylim = c(-2, 5),
        border = alpha("purple",  0.6), col = gray(0.8, alpha = 0.1))

boxplot(e_Ylongvec38md ~ Xlongvec38, 
        cex = 1, 
        ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
        xlab = "Indiv. Biovolume, log10(mm3 ind.-1)",ylim = c(-2, 5),
        border = alpha("darkorange",  0.6), col = gray(0.8, alpha = 0.1),
        add = T)

# Conclusion:
# Estuary
# NMBBSmd: less big ones whe we take away the meropl decapods! Meropl. decapods are the big ones!


# 9. regression models, 9 classes ----------------------------

# 9.1 ordinary least squares regression, ALL size classes, Estuary ------------------

lm1 <- lm(e_Ylongvec38 ~ Xlongvec38)

summary(lm1) # ALL data, OLSR  Slope : -1.216

# choose data with few empty bins

# three size classes only --------------------
e_all<- data.frame (e_Ylongvec38, Xlongvec38)
e_3sizes <- subset (e_all,  Xlongvec38 < -0.95)

e_3sizes <- subset (e_3sizes,  Xlongvec38 > -1.5)
                                      
# NBSSz, estuary------
plot(e_Ylongvec38 ~ Xlongvec38,  main = "NBSSz, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6),
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

points(e_3sizes$e_Ylongvec38 ~e_3sizes$Xlongvec38, main = "NBSSz, estuary", cex = 1.5, col = "navy", pch = 16)
lm3 <- lm(e_3sizes$e_Ylongvec38 ~e_3sizes$Xlongvec38)
abline (lm3, col = "red")
summary(lm3) # slope: -3.5 !!
lm3ROB <- lmRob(e_3sizes$e_Ylongvec38 ~e_3sizes$Xlongvec38)
abline (lm3ROB, col = "darkgreen")
summary(lm3ROB) # slope: -3.2 !!

# 5 size classes  ---------------
e_all<- data.frame (e_Ylongvec38, Xlongvec38)
e_3sizes <- subset (e_all,  Xlongvec38 < -0.6)

e_3sizes <- subset (e_3sizes,  Xlongvec38 > -1.5)

plot(e_Ylongvec38 ~ Xlongvec38,  main = "NBSSz, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6),
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

points(e_3sizes$e_Ylongvec38 ~e_3sizes$Xlongvec38, main = "NBSSz, estuary", cex = 1.5, col = "navy", pch = 16)
lm3 <- lm(e_3sizes$e_Ylongvec38 ~e_3sizes$Xlongvec38)
abline (lm3)
summary(lm3) # slope: -3.2 !!
lm3ROB <- lmRob(e_3sizes$e_Ylongvec38 ~e_3sizes$Xlongvec38)
abline (lm3ROB, col = "darkgreen")
summary(lm3ROB) # slope: -3.15 !!



# 6.  NBSSz, 9 size classes - FOR PAPER ------------------------

# 9 size classes, BEST! ---------------
# 9 size classes, from - 1.5 to 0.1 log10 (mm3)  --

e_all<- data.frame (e_Ylongvec38, Xlongvec38)

e_9sizes <- subset (e_all,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )

summary(e_9sizes)

# In-Practice size range for regression  ---------------
# in practice, the data range from -1.4 to 0.015 log10 (mm3)  -------

plot(e_Ylongvec38 ~ Xlongvec38,  main = "NBSSz, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6),
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))


points(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38, main = "NBSSz, estuary", cex = 1.5, 
       col = alpha ("navy", 0.4) , pch = 16)
lm9 <- lm(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
abline (lm9, col = "red")
summary(lm9) 
# slope: 7 classes,  lm slope = -2.7 !!9 classes,lmRob  slope = -2.41, estuary, TotZOO
lm9ROB <- lmRob(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
abline (lm9ROB, col = "darkgreen")
summary(lm9ROB)
# 7 classes,  lm slope: -2.7 !!  9 classes,lmRob  slope  =  -2.45, , estuary, TotZOO

# FINAL  NBSS plot, ESTUARY, 9 size classes - FOR PAPER ----------
# biovol for classes between -1.5 and -0.25 log10(mm3) - for paper

plot(e_Ylongvec38 ~ Xlongvec38, main = "NBSSz, estuary", 
     cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
# from -1.5 to 0.1: 9 medium size classes
rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
# from 0.1 to 1.7: 9 large size classes

points(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38, main = "NBSSz, estuary", cex = 1.5, 
       col = alpha ("navy",  0.4), pch = 16)
 lm9 <- lm(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
 abline (lm9, col =  "red", lwd =2 )
# summary(lm9) # slope: -2.7 !!
lm9ROB <- lmRob(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
abline (lm9ROB, col = "green", lwd =2)
summary(lm9ROB) # 7 lASSES: slope: -2.7 !! 9 classes: e_9sizes$Xlongvec38  -2.4546

###
# 7. NBSSm, 9 size classes - FOR PAPER ------------------------
# 9 size classes  - FOR PAPER ---------------


e_allm <- data.frame (e_Ylongvec38m, Xlongvec38)

e_9sizesm <- subset (e_allm,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )

#b_Ylongvec38m


plot(e_Ylongvec38m ~ Xlongvec38, main = "NBSSm, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
# from -1.5 to 0.1: 9 medium size classes
rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
# from 0.1 to 1.7: 9 large size classes

points(e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38, main = "NBSSz, estuary",
       cex = 1.5,
       col = alpha ("navy", 0.4), pch = 16)
lm9m <- lm(e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38)
abline (lm9m, col = "red", lwd =1.5)
summary(lm9m) 
# All size classes, estuary
# slope of lm: -2.7 for  NBSSz!!
# slope of lm: -2.97 for NBSSm !! 
# NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)
# All size classes
lm9ROBm <- lmRob(e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38)
abline (lm9ROBm, col =   "green", lwd =1.5)
summary(lm9ROBm) 
# All size classes, estuary
# slope for lmRob: -2.7 for NBSSz in estuary !!
#! 7 clases , slope for lmRob:-3.201 for NBSSm in estuary! e_9sizes$Xlongvec38 , slope =  -2.7158 
# NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)


###
# 7b. NBSSmd, 9 size classes - FOR PAPER ------------------------
# 9 size classes  - FOR PAPER ---------------


e_allmd <- data.frame (e_Ylongvec38md, Xlongvec38)
e_9sizesmd <- subset (e_allmd,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )


b_allmd <- data.frame (b_Ylongvec38md, Xlongvec38)
b_9sizesmd <- subset (b_allmd,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )

s_allmd <- data.frame (s_Ylongvec38md, Xlongvec38)
s_9sizesmd <- subset (s_allmd,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )



#b_Ylongvec38md


plot(e_Ylongvec38md ~ Xlongvec38, main = "NBSSmd, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
# from -1.5 to 0.1: 9 medium size classes
rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
# from 0.1 to 1.7: 9 large size classes

points(e_9sizesmd$e_Ylongvec38m ~e_9sizesmd$Xlongvec38, main = "NBSSmd, estuary",
       cex = 1.5,
       col = alpha ("navy", 0.4), pch = 16)
lm9m <- lm(e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38)
abline (lm9m, col = "red", lwd =1.5)
summary(lm9m) 
# All size classes, estuary
# slope of lm: -2.7 for  NBSSz!!
# slope of lm: -2.97 for NBSSm !! 
# NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)
# All size classes
lm9ROBmd <- lmRob(e_9sizesmd$e_Ylongvec38m ~e_9sizesmd$Xlongvec38)
abline (lm9ROBmd, col =   "green", lwd =1.5)
summary(lm9ROBmd) 
# All size classes, estuary
# slope for lmRob: -2.7 for NBSSz in estuary !!
#! 7 clases , slope for lmRob:-3.201 for NBSSm in estuary! e_9sizes$Xlongvec38 , slope =  -2.7158 
# NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)







# 8. Mero_only - FOR PAPER! -------------------------------

# Plot Mero_only - FOR PAPER! estuary -------------------------------

# Mero_only dataset  (9 classes only)- FOR PAPER! -----------------

# ALL size classes, Meropl_only estuary

e_all.mero_only <- data.frame (e_Ylongvec38m, e_Ylongvec38, Xlongvec38)
e_all.mero_only2 <- e_all.mero_only
head(e_all.mero_only2)
e_all.mero_only2 <- e_all.mero_only2[e_all.mero_only2$e_Ylongvec38m == "NaN",]



# 9 sizeclasses, Meropl_only -------------
e_9sizes.mero_only <- cbind(e_9sizesm, e_9sizes)
e_9sizes.mero_only2 <- e_9sizes.mero_only
e_9sizes.mero_only2 <- e_9sizes.mero_only2[e_9sizes.mero_only2$e_Ylongvec38m == "NaN",]

plot(e_all.mero_only2$e_Ylongvec38 ~ e_all.mero_only2$Xlongvec38, main = "Meropl only, estuary",
     cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5))

points(e_9sizes.mero_only2$e_Ylongvec38 ~ e_9sizes.mero_only2$Xlongvec38, main = "Meropl only, estuary",
     ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
     xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
     ylim = c(-2, 5), xlim = c(-2.5, 2.5), cex = 1.5,
  col = alpha ("darkgreen", 0.4), pch = 16)


rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
# from -1.5 to 0.1: 9 medium size classes
rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
# from 0.1 to 1.7: 9 large size classes


#9.  ANCOVA NBSSm vs NBSSz, robust --------------

# with package "robust"
summary(lm9ROB) #  9 classes: e_9sizes$Xlongvec388  -2.4546
summary(lm9ROBm) # 9 classes: e_9sizes$Xlongvec38 slope   -2.7158 +-0.24

# anova(lm9ROBm, lm9ROB, test = "RF")
# Error in anova.lmRoblist(list(object, ...), cst, ipsi, yc, test = test) : 
#   models were not all fitted to the same size dataset

# with package "permuco"

eallz_m <-  data.frame (e_Ylongvec38m, e_Ylongvec38)
head(eallz_m)
eallz_m.stack <- stack(eallz_m)
head(eallz_m.stack)

eallz_m.stack.X <- cbind( eallz_m.stack,Xlongvec38 )
#View (eallz_m.stack.X)

attach(eallz_m.stack.X)


plot(eallz_m.stack.X$values[eallz_m.stack.X$ind == "e_Ylongvec38m"] ~ eallz_m.stack.X$Xlongvec38[ind == "e_Ylongvec38m"])

plot(eallz_m.stack.X$values[eallz_m.stack.X$ind == "e_Ylongvec38"] ~ eallz_m.stack.X$Xlongvec38[ind == "e_Ylongvec38"])

# select 9 classes for ancova with lmpern -------------

e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38

e_9allz_mX <-  data.frame (NBSSm =  e_9sizesm$e_Ylongvec38m,NBSS =  e_9sizes$e_Ylongvec38 , 
                           Xvec = e_9sizes$Xlongvec38)
e_9allz_m <-  data.frame (NBSSm =  e_9sizesm$e_Ylongvec38m,NBSS =  e_9sizes$e_Ylongvec38 )

head(e_9allz_m)

e_9allz_m.stack <- stack(e_9allz_m)

head(e_9allz_m.stack)
length(e_9allz_m.stack$values)

e_9allz_m.stack.X <- cbind( e_9allz_m.stack, Xvec = e_9sizes$Xlongvec38 )
#View (e_9allz_m.stack.X)

attach(e_9allz_m.stack.X)

plot(e_9allz_m.stack.X$values[e_9allz_m.stack.X$ind == "NBSSm"] ~
       e_9allz_m.stack.X$Xvec[e_9allz_m.stack.X$ind == "NBSSm"])

plot(e_9allz_m.stack.X$values[e_9allz_m.stack.X$ind == "NBSS"] ~
       e_9allz_m.stack.X$Xvec[e_9allz_m.stack.X$ind == "NBSS"])

 lmperm(e_9allz_m.stack.X$values[e_9allz_m.stack.X$ind == "NBSSm"] ~
         e_9allz_m.stack.X$Xvec[e_9allz_m.stack.X$ind == "NBSSm"])

 summary (lm(e_9allz_m.stack.X$values ~
          e_9allz_m.stack.X$Xvec + e_9allz_m.stack.X$ind) )
 
 summary(lmperm(e_9allz_m.stack.X$values ~
          e_9allz_m.stack.X$Xvec + e_9allz_m.stack.X$ind) )
  

 summary (lm(e_9allz_m.stack.X$values ~
               e_9allz_m.stack.X$Xvec * e_9allz_m.stack.X$ind) )
 
 summary(lmperm(e_9allz_m.stack.X$values ~
                  e_9allz_m.stack.X$Xvec * e_9allz_m.stack.X$ind) )
 
 summary(lm(e_9allz_m.stack.X$values ~
                  e_9allz_m.stack.X$Xvec * e_9allz_m.stack.X$ind) )
 
 # ANCOVA RESULT:
 
 # comparison of NBSSz vs NBSSm: -----------------------
 # slopes are NOT significantly different! -------------------------
 # for "lm" and for "lmperm"!  
 # 9 size classes, estuary
 
 
 
 # 10. Compare NAs,   ----------------------------
 # Wilcoxon test  ----------------------------
 
 e_9allz_m <-  data.frame (NBSSm =  e_9sizesm$e_Ylongvec38m,NBSS =  e_9sizes$e_Ylongvec38 )
 
 summary(e_9allz_m)
 
 
# NBSSm: NA's   :186      
 # NA's   :120  
 
 
 wilcox.test(e_9allz_m$NBSSm,e_9allz_m$NBSS, paired = T) # p-value < 2.2e-16
 # Mean ranks are DIFFERENT! p-value < 2.2e-16
 
 wilcox.test(e_9allz_m$NBSSm,e_9allz_m$NBSS, paired = T, na.action = na.exclude) # p-value < 2.2e-16
 # Mean ranks are DIFFERENT! p-value < 2.2e-16
 
 
 wilcox.test(e_9allz_m$NBSSm,e_9allz_m$NBSS, paired = T, na.action = na.pass) # p-value < 2.2e-16
 # Mean ranks are DIFFERENT! p-value < 2.2e-16
 
 
 library(coin)
 e_9allz_m.stack <- stack(e_9allz_m)

 boxplot(e_9allz_m.stack$value ~ e_9allz_m.stack$ind)
  
#  11. Permutation small ones, estuary  -----------------------
 
independence_test(e_9allz_m.stack$value ~ e_9allz_m.stack$ind , paired = T) # n.s.!!
 ## testing zeros and ones
 ## ESTUARY -----------
 # presence-absence transform: zero and one  ------------------------
 
 e_9allz_m <-  data.frame (NBSSm =  e_9sizesm$e_Ylongvec38m,NBSS =  e_9sizes$e_Ylongvec38 )
 
 e_9allz_m.p_a <- e_9allz_m
 
 e_9allz_m.p_a[!is.na(e_9allz_m.p_a)] <- 1
 
   e_9allz_m.p_a[is.na(e_9allz_m.p_a)] <- 0
  
 #  View(e_9allz_m.p_a)
   
   # estuary
  # Wilcoxon test with presence-absence - highly signifficnt difference!!!
   
   wilcox.test(e_9allz_m.p_a$NBSSm,  e_9allz_m.p_a$NBSS,
               paired = T, na.action = na.pass) # p-value =  4.648e-16
   
   wilcox.test(e_9allz_m.p_a$NBSSm,  e_9allz_m.p_a$NBSS,
               paired = T, na.action = na.omit) # p-value =  4.648e-16
   
 
   e_9allz_m.p_astack <- stack(e_9allz_m.p_a)
   
   #boxplot(e_9allz_m.p_astack$value ~ e_9allz_m.p_astack$ind)
   
   independence_test(e_9allz_m.p_astack$value ~ e_9allz_m.p_astack$ind , paired = T) # 
  
   # p-value = 3.944e-07 Highly sgnifficant for permutation test!!! 
   
  summary(e_9allz_m.p_astack)

    summary(as.factor(e_9allz_m.p_a$NBSSm))
    summary(as.factor(e_9allz_m.p_a$NBSS))
   
    df <- data.frame(NBSSz = c(120, 222), NBSSm =c(186, 156))

    df

    
#12. Permutation small ones, Bay  -----------------------
    
   # BAY 
    ## testing zeros and ones
    ## Bay -----------
    # presence-absence transform: zero and one  ------------------------
    
  
    b_all <- data.frame (b_Ylongvec38, Xlongvec38)
    
    b_9sizes <- subset (b_all,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    
    b_allm <- data.frame (b_Ylongvec38m, Xlongvec38)
    
    b_9sizesm <- subset (b_allm,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    #b_Ylongvec38m
    
    
    
    b_9allz_m <-  data.frame (NBSSm =  b_9sizesm$b_Ylongvec38m,NBSS =  b_9sizes$b_Ylongvec38 )
    
    b_9allz_m.p_a <- b_9allz_m
    
    b_9allz_m.p_a[!is.na(b_9allz_m.p_a)] <- 1
    
    b_9allz_m.p_a[is.na(b_9allz_m.p_a)] <- 0
    
    #  View(b_9allz_m.p_a)
    
    # Bay
    # Wilcoxon test with presence-absence - highly signifficnt difference!!!
    
    wilcox.test(b_9allz_m.p_a$NBSSm,  b_9allz_m.p_a$NBSS,
                paired = T, na.action = na.pass) # p-value = 0.0001
    
    wilcox.test(b_9allz_m.p_a$NBSSm,  b_9allz_m.p_a$NBSS,
                paired = T, na.action = na.omit) # p-value = 0.0001
    
    
    b_9allz_m.p_astack <- stack(b_9allz_m.p_a)
    
    head(b_9allz_m.p_astack)
 
    independence_test(b_9allz_m.p_astack$values ~ b_9allz_m.p_astack$ind , paired = T) # 
    # p-value = 0.167
    # p-value = NOT significant for permutation test!!! 
    
    summary(b_9allz_m.p_astack)
    
    summary(as.factor(b_9allz_m.p_a$NBSSm))
    summary(as.factor(b_9allz_m.p_a$NBSS))
    
     
    
    
    
    
    
   #13.Permutation small ones, Shelf  -----------------------
    
    ## SHELF -----------------------------------
    ##
    # presence-absence transform: zero and one  ------------------------
    
    
    s_all <- data.frame (s_Ylongvec38, Xlongvec38)
    
    s_9sizes <- subset (s_all,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    
    s_allm <- data.frame (s_Ylongvec38m, Xlongvec38)
    
    s_9sizesm <- subset (s_allm,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    #s_Ylongvec38m
    
    
    
    s_9allz_m <-  data.frame (NBSSm =  s_9sizesm$s_Ylongvec38m,NBSS =  s_9sizes$s_Ylongvec38 )
    
    s_9allz_m.p_a <- s_9allz_m
    
    s_9allz_m.p_a[!is.na(s_9allz_m.p_a)] <- 1
    
    s_9allz_m.p_a[is.na(s_9allz_m.p_a)] <- 0
    
    #  View(s_9allz_m.p_a)
    
    # SHELF
    # Wilcoxon test with presence-absence - highly significant difference!!!
    
    wilcox.test(s_9allz_m.p_a$NBSSm,  s_9allz_m.p_a$NBSS,
                paired = T, na.action = na.pass) # p-value = 0.0001
    
    wilcox.test(s_9allz_m.p_a$NBSSm,  s_9allz_m.p_a$NBSS,
                paired = T, na.action = na.omit) # p-value = 0.0001
    
    
    s_9allz_m.p_astack <- stack(s_9allz_m.p_a)
    
    head(s_9allz_m.p_astack)
    
    independence_test(s_9allz_m.p_astack$values ~ s_9allz_m.p_astack$ind , paired = T) # 
    # p-value = 0.4
    # p-value = NOT significant for permutation test!!! 
    
    summary(s_9allz_m.p_astack)
    
    summary(as.factor(s_9allz_m.p_a$NBSSm))
    summary(as.factor(s_9allz_m.p_a$NBSS))
    
   
    ####
     
    ## 14. LARGE ONES ------------------------
    # permutation test fpor large size classes NBSSm vs NBSSz
    
    ## 14.1 define LARGE ONES ------------------------
  
    # from 0.1 to 1.7  log10(mm3): 9 large size classes
  
    
    #estuary  
    e_all <- data.frame (e_Ylongvec38, Xlongvec38)
    
    e_9Lsizes <- subset (e_all,   Xlongvec38 > 0.1  & Xlongvec38 < 1.7 )
    
    
    e_allm <- data.frame (e_Ylongvec38m, Xlongvec38)
    
    e_9Lsizesm <- subset (e_allm,   Xlongvec38 > 0.1 & Xlongvec38 < 1.7 )
    
    e_allmd <- data.frame (e_Ylongvec38md, Xlongvec38)
    
    e_9Lsizesmd <- subset (e_allmd,   Xlongvec38 > 0.1 & Xlongvec38 < 1.7 )
    
    
    
     #bay
    
    b_all <- data.frame (b_Ylongvec38, Xlongvec38)
    
    b_9Lsizes <- subset (b_all,   Xlongvec38 > 0.1  & Xlongvec38 < 1.7 )
    
    
    b_allm <- data.frame (b_Ylongvec38m, Xlongvec38)
    
    b_9Lsizesm <- subset (b_allm,   Xlongvec38 > 0.1 & Xlongvec38 < 1.7 )
    
  
     b_allmd <- data.frame (b_Ylongvec38md, Xlongvec38)
    
    b_9Lsizesmd <- subset (b_allmd,   Xlongvec38 > 0.1 & Xlongvec38 < 1.7 )
    
    
    
     # shelf

    s_all <- data.frame (s_Ylongvec38, Xlongvec38)
    
    s_9Lsizes <- subset (s_all,   Xlongvec38 > 0.1  & Xlongvec38 < 1.7 )
    
    
    s_allm <- data.frame (s_Ylongvec38m, Xlongvec38)
    
    s_9Lsizesm <- subset (s_allm,   Xlongvec38 > 0.1 & Xlongvec38 < 1.7 )
    
    
    s_allmd <- data.frame (s_Ylongvec38md, Xlongvec38)
    
    s_9Lsizesmd <- subset (s_allmd,   Xlongvec38 > 0.1 & Xlongvec38 < 1.7 )
    
        
    
    
    
   
    
    ## testing zeros and ones
    
    ## ESTUARY presence-absence transform:  -----------
    # presence-absence transform: zero and one  ------------------------
    
    #15. Permutation large ones, ESTUARY  -----------------------
    
    
    e_9allz_m <-  data.frame (NBSSm =  e_9Lsizesm$e_Ylongvec38m,NBSS =  e_9Lsizes$e_Ylongvec38 )
    
    e_9allz_m.p_a <- e_9allz_m
    
    e_9allz_m.p_a[!is.na(e_9allz_m.p_a)] <- 1
    
    e_9allz_m.p_a[is.na(e_9allz_m.p_a)] <- 0
    
    #  View(e_9allz_m.p_a)
    
    # estuary
    # Wilcoxon test with presence-absence - highly signifficnt difference!!!
    
    wilcox.test(e_9allz_m.p_a$NBSSm,  e_9allz_m.p_a$NBSS,
                paired = T, na.action = na.pass) # p-value = 1.451e-05
    
    wilcox.test(e_9allz_m.p_a$NBSSm,  e_9allz_m.p_a$NBSS,
                paired = T, na.action = na.omit) # p-value = 1.451e-05
    
    
    e_9allz_m.p_astack <- stack(e_9allz_m.p_a)
    
    #boxplot(e_9allz_m.p_astack$value ~ e_9allz_m.p_astack$ind)
    
    independence_test(e_9allz_m.p_astack$value ~ e_9allz_m.p_astack$ind , paired = T) # 
    
    # p-value = p-value = 0.01086  ,  sgnifficant for permutation test!!! 
    
    summary(e_9allz_m.p_astack)
    
    summary(as.factor(e_9allz_m.p_a$NBSSm))
    summary(as.factor(e_9allz_m.p_a$NBSS))
    
    df <- data.frame(NBSSz = c(120, 222), NBSSm =c(186, 156))
    
    df
    
    
    #16. Permutation large ones, Bay  -----------------------
    
    # BAY 
    ## testing zeros and ones
    ## Bay -----------
    # presence-absence transform: zero and one  ------------------------
    
    
    
    
    b_9allz_m <-  data.frame (NBSSm =  b_9Lsizesm$b_Ylongvec38m,NBSS =  b_9Lsizes$b_Ylongvec38 )
    
    b_9allz_m.p_a <- b_9allz_m
    
    b_9allz_m.p_a[!is.na(b_9allz_m.p_a)] <- 1
    
    b_9allz_m.p_a[is.na(b_9allz_m.p_a)] <- 0
    
    #  View(b_9allz_m.p_a)
    
    # Bay
    # Wilcoxon test with presence-absence - highly signifficnt difference!!!
    
    wilcox.test(b_9allz_m.p_a$NBSSm,  b_9allz_m.p_a$NBSS,
                paired = T, na.action = na.pass) # p-value < 2.2e-16
    
    wilcox.test(b_9allz_m.p_a$NBSSm,  b_9allz_m.p_a$NBSS,
                paired = T, na.action = na.omit) # p-value < 2.2e-16
    
    
    b_9allz_m.p_astack <- stack(b_9allz_m.p_a)
    
    head(b_9allz_m.p_astack)
    
    independence_test(b_9allz_m.p_astack$values ~ b_9allz_m.p_astack$ind , paired = T) # 
    # p-value =p-value = 4.658e-13
    # p-value = highly  significant for permutation test!!! 
    
    summary(b_9allz_m.p_astack)
    
    summary(as.factor(b_9allz_m.p_a$NBSSm))
    summary(as.factor(b_9allz_m.p_a$NBSS))
    
  
    
    #17.Permutation large ones, Shelf  -----------------------
    
    ## SHELF -----------------------------------
    ##
    # presence-absence transform: zero and one  ------------------------
    
    
    
    s_9allz_m <-  data.frame (NBSSm =  s_9Lsizesm$s_Ylongvec38m,NBSS =  s_9Lsizes$s_Ylongvec38 )
    
    s_9allz_m.p_a <- s_9allz_m
    
    s_9allz_m.p_a[!is.na(s_9allz_m.p_a)] <- 1
    
    s_9allz_m.p_a[is.na(s_9allz_m.p_a)] <- 0
    
    #  View(s_9allz_m.p_a)
    
    # SHELF
    # Wilcoxon test with presence-absence - highly significant difference!!!
    
    wilcox.test(s_9allz_m.p_a$NBSSm,  s_9allz_m.p_a$NBSS,
                paired = T, na.action = na.pass) # p-value = 1.303e-07
    
    wilcox.test(s_9allz_m.p_a$NBSSm,  s_9allz_m.p_a$NBSS,
                paired = T, na.action = na.omit) # p-value = 1.303e-07
    
    
    s_9allz_m.p_astack <- stack(s_9allz_m.p_a)
    
    head(s_9allz_m.p_astack)
    
    independence_test(s_9allz_m.p_astack$values ~ s_9allz_m.p_astack$ind , paired = T) # 
    # p-value = 0.002461
    # p-value = highly  significant for permutation test!!! 
    
    summary(s_9allz_m.p_astack)
    
    summary(as.factor(s_9allz_m.p_a$NBSSm))
    summary(as.factor(s_9allz_m.p_a$NBSS))
    
    # now with vegan::adonis2
    vegan::adonis2(s_9allz_m.p_astack$values ~ s_9allz_m.p_astack$ind,
    method = "euclidean") # p = 0.002 **
    
   # vegan::adonis2 and coin::independence_test give exactly the same result 
   
    #ok!
    
    
    
    
    ########## A. NBSS plots for ESTUARY - for PAPER ---------  
    
    
    # A. ESTUARY NBSSz, 9 size classes - PLOTS & ANCOVA  FOR PAPER ------------------------
    
    # 9 size classes, BEST! ---------------
    # 9 size classes, from - 1.5 to 0.1 log10 (mm3)  --
    
    e_all<- data.frame (e_Ylongvec38, Xlongvec38)
    
    e_9sizes <- subset (e_all,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    summary(e_9sizes)
    
    # In-Practice size range for regression  ---------------
    # in practice, the data range from -1.4 to 0.015 log10 (mm3)  -------
    
    plot(e_Ylongvec38 ~ Xlongvec38,  main = "NBSSz, estuary",
         cex = 1.2, col = gray(0.7, alpha = 0.6),
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    
    points(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38, main = "NBSSz, estuary", cex = 1.5, 
           col = alpha ("navy", 0.4) , pch = 16)
    lm9 <- lm(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
    abline (lm9, col = "red")
    summary(lm9) 
    # slope: 7 classes,  lm slope = -2.7 !!9 classes,lmRob  slope = -2.41, estuary, TotZOO
    lm9ROB <- lmRob(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
    abline (lm9ROB, col = "darkgreen")
    summary(lm9ROB)
    # 7 classes,  lm slope: -2.7 !!  9 classes,lmRob  slope  =  -2.45, , estuary, TotZOO
    
    # FINAL  NBSS plot, ESTUARY, 9 size classes - FOR PAPER ----------
    # biovol for classes between -1.5 and -0.25 log10(mm3) - for paper
    
    plot(e_Ylongvec38 ~ Xlongvec38, main = "NBSSz, estuary", 
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    points(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38, main = "NBSSz, estuary", cex = 1.5, 
           col = alpha ("navy",  0.4), pch = 16)
    lm9 <- lm(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
    abline (lm9, col =  "red", lwd =2 )
    # summary(lm9) # slope: -2.7 !!
    lm9ROB <- lmRob(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
    abline (lm9ROB, col = "green", lwd =2)
    summary(lm9ROB) # 7 lASSES: slope: -2.7 !! 9 classes: e_9sizes$Xlongvec38  -2.4546
    
    ###
    # 7. NBSSm, 9 size classes - FOR PAPER ------------------------
    # 9 size classes  - FOR PAPER ---------------
    
    
    e_allm <- data.frame (e_Ylongvec38m, Xlongvec38)
    
    e_9sizesm <- subset (e_allm,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    #b_Ylongvec38m
    
    
    plot(e_Ylongvec38m ~ Xlongvec38, main = "NBSSm, estuary",
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    points(e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38, main = "NBSSz, estuary",
           cex = 1.5,
           col = alpha ("navy", 0.4), pch = 16)
    lm9m <- lm(e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38)
    abline (lm9m, col = "red", lwd =1.5)
    summary(lm9m) 
    # All size classes, estuary
    # slope of lm: -2.7 for  NBSSz!!
    # slope of lm: -2.97 for NBSSm !! 
    # NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)
    # All size classes
    lm9ROBm <- lmRob(e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38)
    abline (lm9ROBm, col =   "green", lwd =1.5)
    summary(lm9ROBm) 
    # All size classes, estuary
    # slope for lmRob: -2.7 for NBSSz in estuary !!
    #! 7 clases , slope for lmRob:-3.201 for NBSSm in estuary! e_9sizes$Xlongvec38 , slope =  -2.7158 
    # NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)
    
    
    # 8. Mero_only - FOR PAPER! -------------------------------
    
    # Plot Mero_only - FOR PAPER! estuary -------------------------------
    
    # Mero_only dataset  (9 classes only)- FOR PAPER! -----------------
    
    # ALL size classes, Meropl_only estuary
    
    e_all.mero_only <- data.frame (e_Ylongvec38m, e_Ylongvec38, Xlongvec38)
    e_all.mero_only2 <- e_all.mero_only
    head(e_all.mero_only2)
    e_all.mero_only2 <- e_all.mero_only2[e_all.mero_only2$e_Ylongvec38m == "NaN",]
    
    
    
    # 9 sizeclasses, Meropl_only -------------
    e_9sizes.mero_only <- cbind(e_9sizesm, e_9sizes)
    e_9sizes.mero_only2 <- e_9sizes.mero_only
    e_9sizes.mero_only2 <- e_9sizes.mero_only2[e_9sizes.mero_only2$e_Ylongvec38m == "NaN",]
    
    plot(e_all.mero_only2$e_Ylongvec38 ~ e_all.mero_only2$Xlongvec38, main = "Meropl only, estuary",
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    points(e_9sizes.mero_only2$e_Ylongvec38 ~ e_9sizes.mero_only2$Xlongvec38, main = "Meropl only, estuary",
           ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
           xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
           ylim = c(-2, 5), xlim = c(-2.5, 2.5), cex = 1.5,
           col = alpha ("darkgreen", 0.4), pch = 16)
    
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    
    #9.  ANCOVA NBSSm vs NBSSz, robust --------------
    
    # with package "robust"
    summary(lm9ROB) #  9 classes: e_9sizes$Xlongvec388  -2.4546
    summary(lm9ROBm) # 9 classes: e_9sizes$Xlongvec38 slope   -2.7158 +-0.24
    
    # anova(lm9ROBm, lm9ROB, test = "RF")
    # Error in anova.lmRoblist(list(object, ...), cst, ipsi, yc, test = test) : 
    #   models were not all fitted to the same size dataset
    
    # with package "permuco"
    
    eallz_m <-  data.frame (e_Ylongvec38m, e_Ylongvec38)
    head(eallz_m)
    eallz_m.stack <- stack(eallz_m)
    head(eallz_m.stack)
    
    eallz_m.stack.X <- cbind( eallz_m.stack,Xlongvec38 )
    #View (eallz_m.stack.X)
    
    attach(eallz_m.stack.X)
    
    
    plot(eallz_m.stack.X$values[eallz_m.stack.X$ind == "e_Ylongvec38m"] ~ eallz_m.stack.X$Xlongvec38[ind == "e_Ylongvec38m"])
    
    plot(eallz_m.stack.X$values[eallz_m.stack.X$ind == "e_Ylongvec38"] ~ eallz_m.stack.X$Xlongvec38[ind == "e_Ylongvec38"])
    
    # select 9 classes for ancova with lmpern -------------
    
    e_9sizesm$e_Ylongvec38m ~e_9sizes$Xlongvec38
    
    e_9allz_mX <-  data.frame (NBSSm =  e_9sizesm$e_Ylongvec38m,NBSS =  e_9sizes$e_Ylongvec38 , 
                               Xvec = e_9sizes$Xlongvec38)
    e_9allz_m <-  data.frame (NBSSm =  e_9sizesm$e_Ylongvec38m,NBSS =  e_9sizes$e_Ylongvec38 )
    
    head(e_9allz_m)
    
    e_9allz_m.stack <- stack(e_9allz_m)
    
    head(e_9allz_m.stack)
    length(e_9allz_m.stack$values)
    
    e_9allz_m.stack.X <- cbind( e_9allz_m.stack, Xvec = e_9sizes$Xlongvec38 )
    #View (e_9allz_m.stack.X)
    
    attach(e_9allz_m.stack.X)
    
    plot(e_9allz_m.stack.X$values[e_9allz_m.stack.X$ind == "NBSSm"] ~
           e_9allz_m.stack.X$Xvec[e_9allz_m.stack.X$ind == "NBSSm"])
    
    plot(e_9allz_m.stack.X$values[e_9allz_m.stack.X$ind == "NBSS"] ~
           e_9allz_m.stack.X$Xvec[e_9allz_m.stack.X$ind == "NBSS"])
    
    lmperm(e_9allz_m.stack.X$values[e_9allz_m.stack.X$ind == "NBSSm"] ~
             e_9allz_m.stack.X$Xvec[e_9allz_m.stack.X$ind == "NBSSm"])
    
    summary (lm(e_9allz_m.stack.X$values ~
                  e_9allz_m.stack.X$Xvec + e_9allz_m.stack.X$ind) )
    
    summary(lmperm(e_9allz_m.stack.X$values ~
                     e_9allz_m.stack.X$Xvec + e_9allz_m.stack.X$ind) )
    
    
    summary (lm(e_9allz_m.stack.X$values ~
                  e_9allz_m.stack.X$Xvec * e_9allz_m.stack.X$ind) )
    
    summary(lmperm(e_9allz_m.stack.X$values ~
                     e_9allz_m.stack.X$Xvec * e_9allz_m.stack.X$ind) )
    
    summary(lm(e_9allz_m.stack.X$values ~
                 e_9allz_m.stack.X$Xvec * e_9allz_m.stack.X$ind) )
    
    # ANCOVA RESULT:
    
    # comparison of NBSSz vs NBSSm: -----------------------
    # slopes are NOT significantly different! -------------------------
    # for "lm" and for "lmperm"!  
    # 9 size classes, estuary

    
    ### Bay------------------
    ########## B. NBSS plots for BAY - for PAPER ---------  
    
    # B. BAY NBSSz, 9 size classes - PLOTS & ANCOVA  FOR PAPER ------------------------
    
    # 9 size classes, BEST! ---------------
    # 9 size classes, from - 1.5 to 0.1 log10 (mm3)  --
    
    b_all<- data.frame (b_Ylongvec38, Xlongvec38)
    
    b_9sizes <- subset (b_all,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    summary(b_9sizes)
    
    # In-Practice size range for regression  ---------------
    # in practice, the data range from -1.4 to 0.015 log10 (mm3)  -------
    
    plot(b_Ylongvec38 ~ Xlongvec38,  main = "NBSSz, bay",
         cex = 1.2, col = gray(0.7, alpha = 0.6),
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    
    points(b_9sizes$b_Ylongvec38 ~b_9sizes$Xlongvec38, main = "NBSSz, bay", cex = 1.5, 
           col = alpha ("navy", 0.4) , pch = 16)
    lm9 <- lm(b_9sizes$b_Ylongvec38 ~b_9sizes$Xlongvec38)
    abline (lm9, col = "red")
    summary(lm9) 
    # slope: 7 classes,  lm slope = xx !!9 classes,lmRob  slope = -1.85, bay, TotZOO
    lm9ROB <- lmRob(b_9sizes$b_Ylongvec38 ~b_9sizes$Xlongvec38)
    abline (lm9ROB, col = "darkgreen")
    summary(lm9ROB)
    #9  classes,   lm slope -1.85 9 classes,lmRob  slope  = -1.80, , bay, TotZOO
    
    # FINAL  NBSS plot, Bay, 9 size classes - FOR PAPER ----------
    # biovol for classes between -1.5 and -0.25 log10(mm3) - for paper
    
    plot(b_Ylongvec38 ~ Xlongvec38, main = "NBSSz, bay", 
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    points(b_9sizes$b_Ylongvec38 ~b_9sizes$Xlongvec38, main = "NBSSz, bay", cex = 1.5, 
           col = alpha ("navy",  0.4), pch = 16)
    lm9 <- lm(b_9sizes$b_Ylongvec38 ~b_9sizes$Xlongvec38)
    abline (lm9, col =  "red", lwd =2 )
    # summary(lm9) # slope: -2.7 !!
    lm9ROB <- lmRob(b_9sizes$b_Ylongvec38 ~b_9sizes$Xlongvec38)
    abline (lm9ROB, col = "green", lwd =2)
    summary(lm9ROB) #7 lASSES: slope: 
    
    ###
    # NBSSm, 9 size classes - FOR PAPER ------------------------
    # 9 size classes  - FOR PAPER ---------------
    
    
    b_allm <- data.frame (b_Ylongvec38m, Xlongvec38)
    
    b_9sizesm <- subset (b_allm,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    #b_Ylongvec38m
    
    
    plot(b_Ylongvec38m ~ Xlongvec38, main = "NBSSm, bay",
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    points(b_9sizesm$b_Ylongvec38m ~b_9sizes$Xlongvec38, main = "NBSSz, bay",
           cex = 1.5,
           col = alpha ("navy", 0.4), pch = 16)
    lm9m <- lm(b_9sizesm$b_Ylongvec38m ~b_9sizes$Xlongvec38) # slope lm -1.93 in Bay, NBSSm
    abline (lm9m, col = "red", lwd =1.5)
    summary(lm9m) 
    # All size classes, bay
    # slope of lm: xxx for  NBSSz!!
    # slope of lm: -xxx for NBSSm !! 
    # NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)
    # All size classes
    lm9ROBm <- lmRob(b_9sizesm$b_Ylongvec38m ~b_9sizes$Xlongvec38)
    abline (lm9ROBm, col =   "green", lwd =1.5)
    summary(lm9ROBm)  # -1.899 sllpe lmROB for NBSSm in Bay
    # All size classes, bay
    # slope for lmRob: xx for NBSSz in bay !!
    #! 9 clases , slope for lmRob:-1.8992  for NBSSm in bay! 
    # NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)
    
    
    # 8. Mero_only - FOR PAPER! -------------------------------
    
    # Plot Mero_only - FOR PAPER! bay -------------------------------
    
    # Mero_only dataset  (9 classes only)- FOR PAPER! -----------------
    
    # ALL size classes, Meropl_only bay
    
    b_all.mero_only <- data.frame (b_Ylongvec38m, b_Ylongvec38, Xlongvec38)
    b_all.mero_only2 <- b_all.mero_only
    head(b_all.mero_only2)
    b_all.mero_only2 <- b_all.mero_only2[b_all.mero_only2$b_Ylongvec38m == "NaN",]
    
    
    
    # 9 sizeclasses, Meropl_only -------------
    b_9sizes.mero_only <- cbind(b_9sizesm, b_9sizes)
    b_9sizes.mero_only2 <- b_9sizes.mero_only
    b_9sizes.mero_only2 <- b_9sizes.mero_only2[b_9sizes.mero_only2$b_Ylongvec38m == "NaN",]
    
    plot(b_all.mero_only2$b_Ylongvec38 ~ b_all.mero_only2$Xlongvec38, main = "Meropl only, bay",
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    points(b_9sizes.mero_only2$b_Ylongvec38 ~ b_9sizes.mero_only2$Xlongvec38, main = "Meropl only, bay",
           ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
           xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
           ylim = c(-2, 5), xlim = c(-2.5, 2.5), cex = 1.5,
           col = alpha ("darkgreen", 0.4), pch = 16)
    
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    
    #9.  ANCOVA NBSSm vs NBSSz, robust --------------
    
    # with package "robust"
    summary(lm9ROB) #  9 classes: b_9sizes$Xlongvec388  -1.80606  +-  0.09835
    summary(lm9ROBm) # 9 classes: b_9sizes$Xlongvec38 slope   -1.89927 +-   0.09655
    
    # anova(lm9ROBm, lm9ROB, test = "RF")
    # Error in anova.lmRoblist(list(object, ...), cst, ipsi, yc, test = test) : 
    #   models were not all fitted to the same size dataset
    
    # with package "permuco"
    
    eallz_m <-  data.frame (b_Ylongvec38m, b_Ylongvec38)
    head(eallz_m)
    eallz_m.stack <- stack(eallz_m)
    head(eallz_m.stack)
    
    eallz_m.stack.X <- cbind( eallz_m.stack,Xlongvec38 )
    #View (eallz_m.stack.X)
    
    attach(eallz_m.stack.X)
    
    
    plot(eallz_m.stack.X$values[eallz_m.stack.X$ind == "b_Ylongvec38m"] ~ eallz_m.stack.X$Xlongvec38[ind == "b_Ylongvec38m"])
    
    plot(eallz_m.stack.X$values[eallz_m.stack.X$ind == "b_Ylongvec38"] ~ eallz_m.stack.X$Xlongvec38[ind == "b_Ylongvec38"])
    
    # select 9 classes for ancova with lmpern -------------
    
    b_9sizesm$b_Ylongvec38m ~b_9sizes$Xlongvec38
    
    b_9allz_mX <-  data.frame (NBSSm =  b_9sizesm$b_Ylongvec38m,NBSS =  b_9sizes$b_Ylongvec38 , 
                               Xvec = b_9sizes$Xlongvec38)
    b_9allz_m <-  data.frame (NBSSm =  b_9sizesm$b_Ylongvec38m,NBSS =  b_9sizes$b_Ylongvec38 )
    
    head(b_9allz_m)
    
    b_9allz_m.stack <- stack(b_9allz_m)
    
    head(b_9allz_m.stack)
    length(b_9allz_m.stack$values)
    
    b_9allz_m.stack.X <- cbind( b_9allz_m.stack, Xvec = b_9sizes$Xlongvec38 )
    #View (b_9allz_m.stack.X)
    
    attach(b_9allz_m.stack.X)
    
    plot(b_9allz_m.stack.X$values[b_9allz_m.stack.X$ind == "NBSSm"] ~
           b_9allz_m.stack.X$Xvec[b_9allz_m.stack.X$ind == "NBSSm"])
    
    plot(b_9allz_m.stack.X$values[b_9allz_m.stack.X$ind == "NBSS"] ~
           b_9allz_m.stack.X$Xvec[b_9allz_m.stack.X$ind == "NBSS"])
    
    lmperm(b_9allz_m.stack.X$values[b_9allz_m.stack.X$ind == "NBSSm"] ~
             b_9allz_m.stack.X$Xvec[b_9allz_m.stack.X$ind == "NBSSm"])
    
    summary (lm(b_9allz_m.stack.X$values ~
                  b_9allz_m.stack.X$Xvec + b_9allz_m.stack.X$ind) )
    
    summary(lmperm(b_9allz_m.stack.X$values ~
                     b_9allz_m.stack.X$Xvec + b_9allz_m.stack.X$ind) )
    
    
    summary (lm(b_9allz_m.stack.X$values ~
                  b_9allz_m.stack.X$Xvec * b_9allz_m.stack.X$ind) )
    
    summary(lmperm(b_9allz_m.stack.X$values ~
                     b_9allz_m.stack.X$Xvec * b_9allz_m.stack.X$ind) )
    
    summary(lm(b_9allz_m.stack.X$values ~
                 b_9allz_m.stack.X$Xvec * b_9allz_m.stack.X$ind) )
    
    
    # 
    # ANCOVA RESULT BAY:
    # slopes are NOT significantly different! , p = 0.585, N.S.
    
    # comparison of NBSSz vs NBSSm: -----------------------
    # slopes are NOT significantly different! -------------------------
    # for "lm" and for "lmperm"!  
    # 9 size classes, bay 
    
  #?lmperm  
  
    ########## C. NBSS plots for SHELF - for PAPER ---------  
    
    # SHELF NBSSz, 9 size classes - PLOTS & ANCOVA  FOR PAPER ------------------------
    
    # 9 size classes, BEST! ---------------
    # 9 size classes, from - 1.5 to 0.1 log10 (mm3)  --
    
    s_all<- data.frame (s_Ylongvec38, Xlongvec38)
    
    s_9sizes <- subset (s_all,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    summary(s_9sizes)
    
    # In-Practice size range for regression  ---------------
    # in practice, the data range from -1.4 to 0.015 log10 (mm3)  -------
    
    plot(s_Ylongvec38 ~ Xlongvec38,  main = "NBSSz, shelf",
         cex = 1.2, col = gray(0.7, alpha = 0.6),
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    
    points(s_9sizes$s_Ylongvec38 ~s_9sizes$Xlongvec38, main = "NBSSz, shelf", cex = 1.5, 
           col = alpha ("navy", 0.4) , pch = 16)
    lm9 <- lm(s_9sizes$s_Ylongvec38 ~s_9sizes$Xlongvec38)
    abline (lm9, col = "red")
    summary(lm9) 
    # slope: 7 classes,  lm slope = xx !!9 classes,lmRob  slope = -1.85, shelf, TotZOO
    lm9ROB <- lmRob(s_9sizes$s_Ylongvec38 ~s_9sizes$Xlongvec38)
    abline (lm9ROB, col = "darkgreen")
    summary(lm9ROB)
    #9  classes,   lm slope -1.85 9 classes,lmRob  slope  = -1.80, , shelf, TotZOO
    
    # FINAL  NBSS plot, shelf, 9 size classes - FOR PAPER ----------
    # biovol for classes between -1.5 and -0.25 log10(mm3) - for paper
    
    plot(s_Ylongvec38 ~ Xlongvec38, main = "NBSSz, shelf", 
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    points(s_9sizes$s_Ylongvec38 ~s_9sizes$Xlongvec38, main = "NBSSz, shelf", cex = 1.5, 
           col = alpha ("navy",  0.4), pch = 16)
    lm9 <- lm(s_9sizes$s_Ylongvec38 ~s_9sizes$Xlongvec38)
    abline (lm9, col =  "red", lwd =2 )
    # summary(lm9) # slope: -2.7 !!
    lm9ROB <- lmRob(s_9sizes$s_Ylongvec38 ~s_9sizes$Xlongvec38)
    abline (lm9ROB, col = "green", lwd =2)
    summary(lm9ROB) #7 lASSES: slope: 
    
    ###
    #  NBSSm, 9 size classes - FOR PAPER ------------------------
    # 9 size classes  - FOR PAPER ---------------
    
    
    s_allm <- data.frame (s_Ylongvec38m, Xlongvec38)
    
    s_9sizesm <- subset (s_allm,   Xlongvec38 > -1.5 & Xlongvec38 < 0.1 )
    
    #s_Ylongvec38m
    
    
    plot(s_Ylongvec38m ~ Xlongvec38, main = "NBSSm, shelf",
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    points(s_9sizesm$s_Ylongvec38m ~s_9sizes$Xlongvec38, main = "NBSSz, shelf",
           cex = 1.5,
           col = alpha ("navy", 0.4), pch = 16)
    lm9m <- lm(s_9sizesm$s_Ylongvec38m ~s_9sizes$Xlongvec38) # slope lm -1.93 in shelf, NBSSm
    abline (lm9m, col = "red", lwd =1.5)
    summary(lm9m) 
    # All size classes, shelf
    # slope of lm: xxx for  NBSSz!!
    # slope of lm: -xxx for NBSSm !! 
    # NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)
    # All size classes
    lm9ROBm <- lmRob(s_9sizesm$s_Ylongvec38m ~s_9sizes$Xlongvec38)
    abline (lm9ROBm, col =   "green", lwd =1.5)
    summary(lm9ROBm)  # -1.899 sllpe lmROB for NBSSm in shelf
    # All size classes, shelf
    # slope for lmRob: xx for NBSSz in shelf !!
    #! 9 clases , slope for lmRob:-1.8992  for NBSSm in shelf! 
    # NBSSm é mais ingreme que NBSSz mais ingreme, mais org. pequenos!)
  
    
    
    # 8. Mero_only - FOR PAPER! -------------------------------
    
    # Plot Mero_only - FOR PAPER! shelf -------------------------------
    
    # Mero_only dataset  (9 classes only)- FOR PAPER! -----------------
    
    # ALL size classes, Meropl_only shelf
    
    s_all.mero_only <- data.frame (s_Ylongvec38m, s_Ylongvec38, Xlongvec38)
    s_all.mero_only2 <- s_all.mero_only
    head(s_all.mero_only2)
    s_all.mero_only2 <- s_all.mero_only2[s_all.mero_only2$s_Ylongvec38m == "NaN",]
    
    
    
    # 9 sizeclasses, Meropl_only -------------
    s_9sizes.mero_only <- cbind(s_9sizesm, s_9sizes)
    s_9sizes.mero_only2 <- s_9sizes.mero_only
    s_9sizes.mero_only2 <- s_9sizes.mero_only2[s_9sizes.mero_only2$s_Ylongvec38m == "NaN",]
    
    plot(s_all.mero_only2$s_Ylongvec38 ~ s_all.mero_only2$Xlongvec38, main = "Meropl only, shelf",
         cex = 1.2, col = gray(0.7, alpha = 0.6), pch = 16,
         ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
         xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
         ylim = c(-2, 5), xlim = c(-2.5, 2.5))
    
    points(s_9sizes.mero_only2$s_Ylongvec38 ~ s_9sizes.mero_only2$Xlongvec38, main = "Meropl only, shelf",
           ylab = "Normalized Biovolume, log10(mm3 m-3 mm-3)" ,
           xlab = "Indiv. Biovolume, log10(mm3 ind.-1)" ,
           ylim = c(-2, 5), xlim = c(-2.5, 2.5), cex = 1.5,
           col = alpha ("darkgreen", 0.4), pch = 16)
    
    
    rect(-1.5,-1.9, 0.1, -2.1, col= alpha("navy", 0.6)) # 9 small size classes
    # from -1.5 to 0.1: 9 medium size classes
    rect(0.1, -1.9,  1.7, -2.1, col= alpha("darkorange",  0.6)) # 9 large size classes
    # from 0.1 to 1.7: 9 large size classes
    
    
    #9.  ANCOVA NBSSm vs NBSSz, robust --------------
    
    # with package "robust"
    summary(lm9ROB) #  9 classes: s_9sizes$Xlongvec388  -1.80606  +-  0.09835
    summary(lm9ROBm) # 9 classes: s_9sizes$Xlongvec38 slope   -1.89927 +-   0.09655
    
    # anova(lm9ROBm, lm9ROB, test = "RF")
    # Error in anova.lmRoblist(list(object, ...), cst, ipsi, yc, test = test) : 
    #   models were not all fitted to the same size dataset
    
    # with package "permuco"
    
    eallz_m <-  data.frame (s_Ylongvec38m, s_Ylongvec38)
    head(eallz_m)
    eallz_m.stack <- stack(eallz_m)
    head(eallz_m.stack)
    
    eallz_m.stack.X <- cbind( eallz_m.stack,Xlongvec38 )
    #View (eallz_m.stack.X)
    
    attach(eallz_m.stack.X)
    
    
    plot(eallz_m.stack.X$values[eallz_m.stack.X$ind == "s_Ylongvec38m"] ~ eallz_m.stack.X$Xlongvec38[ind == "s_Ylongvec38m"])
    
    plot(eallz_m.stack.X$values[eallz_m.stack.X$ind == "s_Ylongvec38"] ~ eallz_m.stack.X$Xlongvec38[ind == "s_Ylongvec38"])
    
    # select 9 classes for ancova with lmpern -------------
    
    s_9sizesm$s_Ylongvec38m ~s_9sizes$Xlongvec38
    
    s_9allz_mX <-  data.frame (NBSSm =  s_9sizesm$s_Ylongvec38m,NBSS =  s_9sizes$s_Ylongvec38 , 
                               Xvec = s_9sizes$Xlongvec38)
    s_9allz_m <-  data.frame (NBSSm =  s_9sizesm$s_Ylongvec38m,NBSS =  s_9sizes$s_Ylongvec38 )
    
    head(s_9allz_m)
    
    s_9allz_m.stack <- stack(s_9allz_m)
    
    head(s_9allz_m.stack)
    length(s_9allz_m.stack$values)
    
    s_9allz_m.stack.X <- cbind( s_9allz_m.stack, Xvec = s_9sizes$Xlongvec38 )
    #View (s_9allz_m.stack.X)
    
    attach(s_9allz_m.stack.X)
    
    plot(s_9allz_m.stack.X$values[s_9allz_m.stack.X$ind == "NBSSm"] ~
           s_9allz_m.stack.X$Xvec[s_9allz_m.stack.X$ind == "NBSSm"])
    
    plot(s_9allz_m.stack.X$values[s_9allz_m.stack.X$ind == "NBSS"] ~
           s_9allz_m.stack.X$Xvec[s_9allz_m.stack.X$ind == "NBSS"])
    
    lmperm(s_9allz_m.stack.X$values[s_9allz_m.stack.X$ind == "NBSSm"] ~
             s_9allz_m.stack.X$Xvec[s_9allz_m.stack.X$ind == "NBSSm"])
    
    summary (lm(s_9allz_m.stack.X$values ~
                  s_9allz_m.stack.X$Xvec + s_9allz_m.stack.X$ind) )
    
    summary(lmperm(s_9allz_m.stack.X$values ~
                     s_9allz_m.stack.X$Xvec + s_9allz_m.stack.X$ind) )
    
    
    summary (lm(s_9allz_m.stack.X$values ~
                  s_9allz_m.stack.X$Xvec * s_9allz_m.stack.X$ind) )
    
    summary(lmperm(s_9allz_m.stack.X$values ~
                     s_9allz_m.stack.X$Xvec * s_9allz_m.stack.X$ind) )
    
    summary(lm(s_9allz_m.stack.X$values ~
                 s_9allz_m.stack.X$Xvec * s_9allz_m.stack.X$ind) )
    
    
    # 
    # ANCOVA RESULT shelf:
    # slopes are NOT significantly different! , p = 0.585, N.S.
    
    # comparison of NBSSz vs NBSSm: -----------------------
    # slopes are NOT significantly different! -------------------------
    # for "lm" and for "lmperm"!  
    # 9 size classes, shelf 
    
    
    
    #### ANCOVA  SHELF  vs Bay vs estuary
   
    
    e_lm9ROB <- lmRob(e_9sizes$e_Ylongvec38 ~e_9sizes$Xlongvec38)
    
     b_lm9ROB <- lmRob(b_9sizes$b_Ylongvec38 ~b_9sizes$Xlongvec38)
    
    s_lm9ROB <- lmRob(s_9sizes$s_Ylongvec38 ~s_9sizes$Xlongvec38)
    
    
    
    
    # multiple model
 
    # GLOBAL permmuation ANCOVA, Estuary vs Bay vs Shelf, NBSSz --------------------------
    
    e.vs.b.vs.s.df <-   data.frame( e_data = e_9sizes$e_Ylongvec38, b_data = b_9sizes$b_Ylongvec38,
                                    s_data = s_9sizes$s_Ylongvec38)
    
    e.vs.b.vs.s.df.stack = stack(e.vs.b.vs.s.df)
    plot(e.vs.b.vs.s.df.stack$values ~e.vs.b.vs.s.df.stack$ind)
    
    e.vs.b.vs.s.df.stack.X <- e.vs.b.vs.s.df.stack
    
    e.vs.b.vs.s.df.stack.X$Xlongvec38 <- rep (b_9sizes$Xlongvec38,3) 
    
    summary(lmperm(e.vs.b.vs.s.df.stack.X$values ~
                     e.vs.b.vs.s.df.stack.X$ind + e.vs.b.vs.s.df.stack.X$Xlongvec38) )
    
    summary(aov(lmperm(e.vs.b.vs.s.df.stack.X$values ~
                     e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) ))
    # Global ANCOVA is signifficant ! p  = 0.0001 ***in "lmperm"
    #  1e-04 ***
    
    
    summary(lm(e.vs.b.vs.s.df.stack.X$values ~
                     e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) )
    #p-value: < 2.2e-16, R² = 0.59, n = 796
    
    summary(lm(e.vs.b.vs.s.df.stack.X$values ~
                 e.vs.b.vs.s.df.stack.X$ind + e.vs.b.vs.s.df.stack.X$Xlongvec38) )
    #p-value: < 2.2e-16, R² = 0.58, n = 796
    
    summary(aov(lm(e.vs.b.vs.s.df.stack.X$values ~
                 e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) ))
    # Global ANCOVA is signifficant ! p <2e-16 ***in "lm"
    
  
      ## Now With NBSSho, compare areas with ANCOVA -----------
    
    # GLOBAL permmuation ANCOVA, Estuary vs Bay vs Shelf, NBSSho = NBSSm --------------------------
    
    e.vs.b.vs.s.df <-   data.frame( e_data = e_9sizesm$e_Ylongvec38, b_data = b_9sizesm$b_Ylongvec38,
                                    s_data = s_9sizesm$s_Ylongvec38)
    
    e.vs.b.vs.s.df.stack = stack(e.vs.b.vs.s.df)
    plot(e.vs.b.vs.s.df.stack$values ~e.vs.b.vs.s.df.stack$ind)
    
    e.vs.b.vs.s.df.stack.X <- e.vs.b.vs.s.df.stack
    
    e.vs.b.vs.s.df.stack.X$Xlongvec38 <- rep (b_9sizes$Xlongvec38,3) 
    
    summary(lmperm(e.vs.b.vs.s.df.stack.X$values ~
                     e.vs.b.vs.s.df.stack.X$ind + e.vs.b.vs.s.df.stack.X$Xlongvec38) )
    
    summary(aovperm(lmperm(e.vs.b.vs.s.df.stack.X$values ~
                         e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) ))
    # Global ANCOVA is signifficant ! p  = 0.003 ***in "lmperm"
    #  1e-04 ***
    
    
    summary(lm(e.vs.b.vs.s.df.stack.X$values ~
                 e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) )
    
    summary(lm(e.vs.b.vs.s.df.stack.X$values ~
                 e.vs.b.vs.s.df.stack.X$ind + e.vs.b.vs.s.df.stack.X$Xlongvec38) )
     
    summary(aov(lm(e.vs.b.vs.s.df.stack.X$values ~
                     e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) ))
    # Global ANCOVA is signifficant ! p = 0.0025 ***in "lm", NBSSho = NBSSm
    

    ## Now With NBSSmd, compare areas with ANCOVA -----------
      #md -> md = "NBSSnmd". wo_meroDec without meroplanktonic Decapoda. As above, but without meroplanktonic Decapoda 
    # GLOBAL permmuation ANCOVA, Estuary vs Bay vs Shelf, NBSS = NBSSmd --------------------------
    
    e.vs.b.vs.s.df <-   data.frame( e_data = e_9Lsizesmd $e_Ylongvec38, b_data = b_9sizesmd$b_Ylongvec38,
                                    s_data = s_9sizesmd$s_Ylongvec38)
    
    e.vs.b.vs.s.df.stack = stack(e.vs.b.vs.s.df)
    plot(e.vs.b.vs.s.df.stack$values ~e.vs.b.vs.s.df.stack$ind)
    
    e.vs.b.vs.s.df.stack.X <- e.vs.b.vs.s.df.stack
    
    e.vs.b.vs.s.df.stack.X$Xlongvec38 <- rep (b_9sizes$Xlongvec38,3) 
    
    summary(lmperm(e.vs.b.vs.s.df.stack.X$values ~
                     e.vs.b.vs.s.df.stack.X$ind + e.vs.b.vs.s.df.stack.X$Xlongvec38) )
    
    summary(aovperm(lmperm(e.vs.b.vs.s.df.stack.X$values ~
                             e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) ))
    # Global ANCOVA is signifficant ! p  = 0.0002 in "lmperm", NBSSmd
    
  
    summary(lm(e.vs.b.vs.s.df.stack.X$values ~
                 e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) )
    
    summary(lm(e.vs.b.vs.s.df.stack.X$values ~
                 e.vs.b.vs.s.df.stack.X$ind + e.vs.b.vs.s.df.stack.X$Xlongvec38) )
    
    summary(aov(lm(e.vs.b.vs.s.df.stack.X$values ~
                     e.vs.b.vs.s.df.stack.X$ind * e.vs.b.vs.s.df.stack.X$Xlongvec38) ))
    # Global ANCOVA is signifficant ! p = 8.58e-08 *** ***in "lm", NBSSmd
    
    
    
    
    
    
    # ANCOVA Estuary vs Bay, NBSSz --------------------------
    
 e.vs.b.df <-   data.frame( e_data = e_9sizes$e_Ylongvec38, b_data = b_9sizes$b_Ylongvec38)
   
 e.vs.b.df.stack = stack(e.vs.b.df)
 plot(e.vs.b.df.stack$values ~e.vs.b.df.stack$ind)

 e.vs.b.df.stack.X <- e.vs.b.df.stack
     
 e.vs.b.df.stack.X$Xlongvec38 <- rep (b_9sizes$Xlongvec38,2) 

 summary(lmperm(e.vs.b.df.stack.X$values ~
                  e.vs.b.df.stack.X$ind + e.vs.b.df.stack.X$Xlongvec38) )
 

 summary(lmperm(e.vs.b.df.stack.X$values ~
                  e.vs.b.df.stack.X$ind * e.vs.b.df.stack.X$Xlongvec38) )
 
## Result: BAY vs Estuary  p = 1.246e-03 = 0.00124, ANCOVA with lmperm, NBSSz 

#  BAY vs Estuary is different slope from By, for nbssz!

 # ANCOVA Estuary vs Shelf, NBSSz    --------------------------------------------------------- 
 
 e.vs.s.df <-   data.frame( e_data = e_9sizes$e_Ylongvec38, s_data = s_9sizes$s_Ylongvec38)
 
 e.vs.s.df.stack = stack(e.vs.s.df)
 plot(e.vs.s.df.stack$values ~e.vs.s.df.stack$ind)
 
 e.vs.s.df.stack.X <- e.vs.s.df.stack
 
 e.vs.s.df.stack.X$Xlongvec38 <- rep (s_9sizes$Xlongvec38,2) 
 
 summary(lmperm(e.vs.s.df.stack.X$values ~
                  e.vs.s.df.stack.X$ind + e.vs.s.df.stack.X$Xlongvec38) )
 
 
 summary(lmperm(e.vs.s.df.stack.X$values ~
                  e.vs.s.df.stack.X$ind * e.vs.s.df.stack.X$Xlongvec38) )
 
 # Shelf vs Estuary  p = 7.584 e-05 , ANCOVA with lmperm, NBSSz 
 
 
 
 # ANCOVA Bay vs Shelf, NBSSz    --------------------------------------------------------- 
 
 b.vs.s.df <-   data.frame( b_data = b_9sizes$b_Ylongvec38, s_data = s_9sizes$s_Ylongvec38)
 
 b.vs.s.df.stack = stack(b.vs.s.df)
 plot(b.vs.s.df.stack$values ~b.vs.s.df.stack$ind)
 
 b.vs.s.df.stack.X <- b.vs.s.df.stack
 
 b.vs.s.df.stack.X$Xlongvec38 <- rep (s_9sizes$Xlongvec38,2) 
 
 summary(lmperm(b.vs.s.df.stack.X$values ~
                  b.vs.s.df.stack.X$ind + b.vs.s.df.stack.X$Xlongvec38) )
 
 
 summary(lmperm(b.vs.s.df.stack.X$values ~
                  b.vs.s.df.stack.X$ind * b.vs.s.df.stack.X$Xlongvec38) )
 
 # Bay  vs Estuary  p = n.s., ANCOVA with lmperm, NBSSz 
 # p = 5.569e-01 =  0.5569 = ns
 
 
 
    # TASKS-------------------------------------------
    # (OK)finish plots for  Bay
 
    # (OK)finish plots for  Shelf
    
    # (OK) ANCOVA bay vs shelf vs estuary 
    
    
     # ( ) insert equations in plots - urgently necessary?
   
      
    #### 10. Compare elevations -----------------------------------------
 
 # Source: https://stats.stackexchange.com/questions/435644/is-there-a-method-to-look-for-significant-difference-between-two-linear-regressi
 
 compare.coeff <- function(b1,se1,b2,se2){
   return((b1-b2)/sqrt(se1^2+se2^2))
 }
# where b1, b2 are coefficients of the two lm and se1 and se2 are the standard errors
 
 #We fit two linear models:
   
 #  lm1 = lm(Height ~ Age,data=subset(df,df$Fladen=="A"))
 
 # NBSSz for Estuary, Bay, shelf
 lm1e = lm( e_9sizes$e_Ylongvec38 ~ e_9sizes$Xlongvec38)
 lm1b = lm( b_9sizes$b_Ylongvec38 ~ b_9sizes$Xlongvec38)
 lm1s = lm( s_9sizes$s_Ylongvec38 ~ s_9sizes$Xlongvec38)
 
 # Slope ------------------------
 
# Estuary vs Bay, Slope 
 b1 <- summary(lm1e)$coefficients[2,1]
 se1 <- summary(lm1e)$coefficients[2,2]
 b2 <- summary(lm1b)$coefficients[2,1]
 se2 <- summary(lm1b)$coefficients[2,2]
 #We calculate the p-value using the z statistic from a normal distribution:
   
   p_value = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
 p_value # slope Estuary vs Bay, p = 0.00159


 # Estuary vs Shelf, Slope 
 b1 <- summary(lm1e)$coefficients[2,1]
 se1 <- summary(lm1e)$coefficients[2,2]
 b2 <- summary(lm1s)$coefficients[2,1]
 se2 <- summary(lm1s)$coefficients[2,2]
 #We calculate the p-value using the z statistic from a normal distribution:
 
 p_value = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
 p_value # slope Estuary vs Shelf, p = 0.000166
 
 
 # Bay vs Shelf, Slope 
 b1 <- summary(lm1b)$coefficients[2,1]
 se1 <- summary(lm1b)$coefficients[2,2]
 b2 <- summary(lm1s)$coefficients[2,1]
 se2 <- summary(lm1s)$coefficients[2,2]
 #We calculate the p-value using the z statistic from a normal distribution:
 
 p_value = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
 p_value # slope Estuary vs Shelf, p = 0.55, Not Sign.
 
 # Intercept  ---------------
 # Intecept equals  biovolume at log10Vol = 0,  bio vol = 1 mm3.
 # The intercept is located at the large-sized  end of the useful linear spectrum (approx  1.3 mm esd size, large copepods, large Brachyuran crab larvae, Caridean shrimp larvae, upper end of the useful linear spectrum)
 

# Estuary vs Bay, Intercept 
 b1 <- summary(lm1e)$coefficients[1,1]
 se1 <- summary(lm1e)$coefficients[1,2]
 b2 <- summary(lm1b)$coefficients[1,1]
 se2 <- summary(lm1b)$coefficients[1,2]
 #We calculate the p-value using the z statistic from a normal distribution:
 
 p_value = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
 p_value # Intercept Estuary vs Bay, p = n.s.
 
 
 # Estuary vs Shelf, Intercept 
 b1 <- summary(lm1e)$coefficients[1,1]
 se1 <- summary(lm1e)$coefficients[1,2]
 b2 <- summary(lm1s)$coefficients[1,1]
 se2 <- summary(lm1s)$coefficients[1,2]
 #We calculate the p-value using the z statistic from a normal distribution:
 
 p_value = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
 p_value # Intercept Estuary vs Shelf, p = n.s.
 
 
 # Bay vs Shelf, Intercept 
 b1 <- summary(lm1b)$coefficients[1,1]
 se1 <- summary(lm1b)$coefficients[1,2]
 b2 <- summary(lm1s)$coefficients[1,1]
 se2 <- summary(lm1s)$coefficients[1,2]
 #We calculate the p-value using the z statistic from a normal distribution:
 
 p_value = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
 p_value # Intercept Bay vs Shelf,  0.0085,  Sign.!!
 
 
 
 

 
    # End of Script  -----------------------------------------------------------
    
    