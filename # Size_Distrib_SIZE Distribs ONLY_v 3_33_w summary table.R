# Size Spectra analysis 
#
# for the manuscript:
#" Size niche interactions (functional replacement and competitive exclusion)
# between mero- and holoplankton shape the size spectrum 
#of tropical estuarine and marine ecosystems "
#
# Feret length vs Abundance - distributions
# using normalized Abundance and non-linearly increasing  Size (i.e. Feret length) bins

# version 03.33
# Ralf Schwamborn

#  August 11, 2025


# Summary of the rationale: -----------

# PLOT Size Distribs ONLY
# Size (Feret length, mm) vs  Abundance (ind. per m3)


# 0. load libraries

library(scales)
library(stringr)

# 1. import data -----------------------------


# set working directory -----------------------

setwd("C:/Users/RALF/Documents/Papers/0 - Paper Denise Size Spectra")
 
# space- time factors (as dummies)

spacetime   <- read.csv("~/Papers/0 - Paper Denise Size Spectra/Space_Time_disttam_v2RSb.csv")
# View(spacetime) 


 totZoo <- read.table("~/Papers/0 - Paper Denise Size Spectra/totZoo_disttam_v2RSb.csv",
                     sep = ",", header = T)

 # totZoo <- read.table("~/Papers/0 - Paper Denise Size Spectra/totZoo_disttam_v2RSb_no_head.csv",
 #                      sep = ",", header = F)

#   View(totZoo)
 totZoo <- as.data.frame((totZoo))
 #totZoo <- as.numeric((totZoo))
 
 mean_totZoo <-  colMeans(totZoo)
 
 
 totZoo2 <- cbind(totZoo,spacetime )
 #   View(totZoo2)  
 
 plot (as.numeric(totZoo[1,]))
 diff(as.numeric(totZoo[1,]))
 
 
BraZ <- read.csv("~/Papers/0 - Paper Denise Size Spectra/BraZ_disttam_v2RSb.csv")
 #   View(BraZ)
BraZ2 <- cbind(BraZ,spacetime )
#   View(BraZ2)  
plot (as.numeric(BraZ[43,]))
diff(as.numeric(BraZ[43,]))


mean_BraZ <-  colMeans(BraZ)



Cop  <- read.csv("~/Papers/0 - Paper Denise Size Spectra/Cop_disttam_v2RSb.csv")
#   View(Cop)
Cop2 <- cbind(Cop,spacetime )
#   View(Cop2)

length(Cop$X0.43498)

mean_Cop <-  colMeans(Cop)




# X axis vector -----------------------------------------

char.X.feret.vec.mm <- "0.29,	0.33196,	0.38,	0.43498,	0.49793,	0.56998	0.65245	0.74686	0.85493	0.97864	1.1202	1.2823	1.4679	1.6803	1.9235	2.2018	2.5204	2.8851	3.3025	3.7804	4.3274	4.9536	5.6704	6.4909	7.4302	8.5053	9.736	11.145	12.757	14.603	16.717	19.136	21.904	25.074	28.702"
  
X.feret.vec.mm <-  str_replace_all(char.X.feret.vec.mm, "\t", "," )
  
X.feret.vec.mm <-c(0.29,0.33196,0.38,0.43498,0.49793,0.56998,0.65245,0.74686,0.85493,0.97864,1.1202,1.2823,1.4679,1.6803,1.9235,2.2018,2.5204,2.8851,3.3025,3.7804,4.3274,4.9536,5.6704,6.4909,7.4302,8.5053,9.736,11.145,12.757,14.603,16.717,19.136,21.904,25.074,28.702)


#  ?str_replace

diff(X.feret.vec.mm)

plot(X.feret.vec.mm)


plot((X.feret.vec.mm))

plot(diff(X.feret.vec.mm))

plot(log10(X.feret.vec.mm))




length(totZoo[1,])

# plot example size spectrum ----------------

# linear-liner plot (no log)
plot(as.numeric(totZoo[1,]) ~ (X.feret.vec.mm))
lines(as.numeric(totZoo[1,]) ~ (X.feret.vec.mm))

# linear-liner plot (no log), only small animals < 8 mm
plot(as.numeric(totZoo[1,]) ~ (X.feret.vec.mm), xlim = c(0,8))
lines(as.numeric(totZoo[1,]) ~ (X.feret.vec.mm))


# y is linear, x is log10
plot(as.numeric(totZoo[1,]) ~ log10(X.feret.vec.mm))
lines(as.numeric(totZoo[1,]) ~ log10(X.feret.vec.mm))

# log-log plot
plot(log10(as.numeric(totZoo[1,])) ~ log10(X.feret.vec.mm))
lines(log10(as.numeric(totZoo[1,])) ~ log10(X.feret.vec.mm))

# log-log plot with logx+1
plot(log10(1+as.numeric(totZoo[1,])) ~ log10(X.feret.vec.mm))
lines(log10(1+as.numeric(totZoo[1,])) ~ log10(X.feret.vec.mm))

# 
# linear-liner plot (no log), only small animals < 8 mm
plot(as.numeric(totZoo[1,]) ~ (X.feret.vec.mm), xlim = c(0,8),
     xlab = "Size (Feret length, mm)", ylab = "Abundance (ind. m-3)" ) 
lines(as.numeric(totZoo[1,]) ~ (X.feret.vec.mm))


#  NICE PLOT ------------------
# log-log plot , NORMALIZED Abund
plot(log10(1+as.numeric(totZoo[96,1:30]))/diff(X.feret.vec.mm[1:30]) ~ log10(X.feret.vec.mm[1:30]), 
     xlim = c(-0.7,0.9),
     xlab = "Size, log10(Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ) 
lines(log10(1+as.numeric(totZoo[96,1:30]))/diff(X.feret.vec.mm[1:30]) ~ log10(X.feret.vec.mm[1:30]))



# FINAL NICE PLOTs for paper------------------

# example station sample  43, estuary ----------------

# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot(log10(1+as.numeric(totZoo[43,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Sample 43,estuary") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines(log10(1+as.numeric(BraZ[43,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund
lines(log10(1+as.numeric(Cop[43,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "darkgreen") 



legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)


# FINAL NICE PLOTs for paper------------------
# EXAMPLES
# example station sample  74, estuary ----------------

st.num <- 73
#st.num <- 51
#st.num <- 49
# st.num <- 41


# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot(log10(1+as.numeric(totZoo[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Sample 73, estuary") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines(log10(1+as.numeric(BraZ[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund
lines(log10(1+as.numeric(Cop[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 



legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)


#### Test different bin sizes 
# Scenario 2, bins2 =   bins_0 ^2

# X axis vector -----------------------------------------


X.feret.vec.mm <-c(0.29,0.33196,0.38,0.43498,0.49793,0.56998,0.65245,0.74686,0.85493,0.97864,1.1202,1.2823,1.4679,1.6803,1.9235,2.2018,2.5204,2.8851,3.3025,3.7804,4.3274,4.9536,5.6704,6.4909,7.4302,8.5053,9.736,11.145,12.757,14.603,16.717,19.136,21.904,25.074,28.702)

X2feret.vec.mm <- X.feret.vec.mm^2


# example station sample  74, estuary ----------------

st.num <- 73
#st.num <- 51
#st.num <- 49
# st.num <- 41


# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot(log10(1+as.numeric(totZoo[st.num,1:30]))/diff(X2feret.vec.mm[1:30]) ~ (X2feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Sample 73, estuary") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines(log10(1+as.numeric(BraZ[st.num,1:30]))/diff(X2feret.vec.mm[1:30]) ~ (X2feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund
lines(log10(1+as.numeric(Cop[st.num,1:30]))/diff(X2feret.vec.mm[1:30]) ~ (X2feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 



legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)


#### Test different bin sizes 
# Scenario 3, bins3 =  sqrt(  bins_0 )

# X axis vector -----------------------------------------


X.feret.vec.mm <-c(0.29,0.33196,0.38,0.43498,0.49793,0.56998,0.65245,0.74686,0.85493,0.97864,1.1202,1.2823,1.4679,1.6803,1.9235,2.2018,2.5204,2.8851,3.3025,3.7804,4.3274,4.9536,5.6704,6.4909,7.4302,8.5053,9.736,11.145,12.757,14.603,16.717,19.136,21.904,25.074,28.702)

X2feret.vec.mm <- sqrt(X.feret.vec.mm)


# example station sample  74, estuary ----------------

st.num <- 73
#st.num <- 51
#st.num <- 49
# st.num <- 41


# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot(log10(1+as.numeric(totZoo[st.num,1:30]))/diff(X2feret.vec.mm[1:30]) ~ (X2feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Sample 73, estuary") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines(log10(1+as.numeric(BraZ[st.num,1:30]))/diff(X2feret.vec.mm[1:30]) ~ (X2feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund
lines(log10(1+as.numeric(Cop[st.num,1:30]))/diff(X2feret.vec.mm[1:30]) ~ (X2feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 



legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)





############################
# correlation analysis







# correlation analysis

Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])

A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])


A1contrib <- A1/Total.zoo
B1contrib  <- B1/Total.zoo
                              
cor.test (A1 , B1)
c1 <- cor.test (A1contrib , B1contrib)
p.value_pears <- c1$p.value
cor.value_pears <- c1$estimate
c_sp <- cor.test (A1contrib , B1contrib, method = "spearman")
p.value_spearm <- round( c_sp$p.value, 18)
cor.value_rho_Spear <- round( c_sp$estimate , 3)

text(2.5, 7, paste(p.value_spearm,cor.value_rho_Spear ))





# FINAL NICE PLOTs for paper------------------

### OVERALL MEAN - ALL REGIONS ----------------------------

# FINAL NICE PLOT for paper------------------

# All Areas ----------------

# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot(log10(1+as.numeric(mean_totZoo[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Overall Mean, All Areas") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines(log10(1+as.numeric(mean_BraZ[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund - Copp
lines(log10(1+as.numeric(mean_Cop[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)




  # subsets by region  ------------------------------------

e_totZoo <- totZoo[totZoo2$est == 1,]

mean_e_totZoo <-  colMeans(e_totZoo)

b_totZoo <-  totZoo[totZoo2$bai == 1,]

mean_b_totZoo  <-  colMeans(b_totZoo )

s_totZoo <-  totZoo[totZoo2$pla == 1,]

mean_s_totZoo <-  colMeans(s_totZoo)


e_BraZ2 <- BraZ2[BraZ2$est == 1,]

mean_e_BraZ2 <-  colMeans(e_BraZ2)

b_BraZ2 <-  BraZ2[BraZ2$bai == 1,]

mean_b_BraZ2 <-  colMeans(b_BraZ2)

s_BraZ2 <-  BraZ2[BraZ2$pla == 1,]

mean_s_BraZ2 <-  colMeans(s_BraZ2)


e_Cop2 <- Cop2[Cop2$est == 1,]

mean_e_Cop2 <-  colMeans(e_Cop2)

b_Cop2 <-  Cop2[Cop2$bai == 1,]

mean_b_Cop2 <-  colMeans(b_Cop2)

s_Cop2 <-  Cop2[Cop2$pla == 1,]

mean_s_Cop2 <-  colMeans(s_Cop2)


# plots by region ---------------------------



# Estuary ------------------

# FINAL NICE PLOT for paper------------------

# Log-linear plots -----------

# Estuary ----------------

# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot(log10(1+as.numeric(mean_e_totZoo[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),ylim= c(0,28),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Estuary") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines(log10(1+as.numeric(mean_e_BraZ2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund - Copp
lines(log10(1+as.numeric(mean_e_Cop2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)

# Bay ----------------

# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot(log10(1+as.numeric(mean_b_totZoo[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),ylim= c(0,25),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Bay") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines(log10(1+as.numeric(mean_b_BraZ2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund - Copp
lines(log10(1+as.numeric(mean_b_Cop2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)



# Shelf ----------------

# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot(log10(1+as.numeric(mean_s_totZoo[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5), ylim= c(0,25),
     xlab = "Size, (Feret length, mm)", 
     ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Shelf") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines(log10(1+as.numeric(mean_s_BraZ2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund - Copp
lines(log10(1+as.numeric(mean_s_Cop2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)


###### linear-linear plots (not log) -----------------

# Estuary ------------------

#  NICE PLOTs------------------

# Estuary ----------------

# lin-lin plot , NORMALIZED Abund - TOTAL ZOO
plot((as.numeric(mean_e_totZoo[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size (Feret length, mm)", 
     ylab = "Normalized Abundance (ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Estuary") 

# lin-lin plot , NORMALIZED Abund - BRAZ
lines((as.numeric(mean_e_BraZ2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size (Feret length, mm)", 
      ylab = "Normalized Abundance, (ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# lin-lin plot , NORMALIZED Abund - Cop
lines((as.numeric(mean_e_Cop2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size (Feret length, mm)", 
      ylab = "Normalized Abundance (ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)

# Bay ----------------

# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot((as.numeric(mean_b_totZoo[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size (Feret length, mm)", 
     ylab = "Normalized Abundance (ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Bay") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines((as.numeric(mean_b_BraZ2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund - Copp
lines((as.numeric(mean_b_Cop2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)



# Shelf ----------------

# log-lin plot , NORMALIZED Abund - TOTAL ZOO
plot((as.numeric(mean_s_totZoo[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
     xlim = c(0,5),
     xlab = "Size (Feret length, mm)", 
     ylab = "Normalized Abundance (ind. m-3 mm-1)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Shelf") 

# log-lin plot , NORMALIZED Abund - BRAZ
lines((as.numeric(mean_s_BraZ2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length, mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkorange") 


# log-lin plot , NORMALIZED Abund - Copp
lines((as.numeric(mean_s_Cop2[1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
      xlim = c(0,5),
      xlab = "Size, (Feret length , mm)", 
      ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)




### 3c LOOP - plot All stations ---------------------------------------


# save figs on Disk -----------------------------------

# sise spec not log  ------------------------------------------------

# setwd("C:/Users/RALF/Documents/Papers/0 - Paper Denise Size Spectra/Diverse_figures_fromR/test_Spectra5")

setwd("C:/Users/RALF/Documents/Papers/0 - Paper Denise Size Spectra/Diverse_figures_fromR/test_Spectra4")


st.list <- 1:121

#w = 4

for (w in 1:38)  { 
  
  st.num <- w
  
  #  svg(filename=  paste ("sample_bay",st.num,".svg") ) 
  
  png(filename=  paste ("sample_bay",st.num,".png") ) 
  
#  st.num <- 1
  
  # log-lin plot , NORMALIZED Abund - TOTAL ZOO
  plot(log10(1+as.numeric(totZoo[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
       xlim = c(0,5),
       xlab = "Size, (Feret length, mm)", 
       ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
       type = "l", lty = 1, cex = 15, col = "black",
       main =paste ("Bay, sample no. " , st.num) ) 
  
  # log-lin plot , NORMALIZED Abund - BRAZ
  lines(log10(1+as.numeric(BraZ[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
        xlim = c(0,5),
        xlab = "Size, (Feret length, mm)", 
        ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
        type = "l", lty = 1, cex = 15, col = "darkorange") 
  
  
  # log-lin plot , NORMALIZED Abund
  lines(log10(1+as.numeric(Cop[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
        xlim = c(0,5),
        xlab = "Size, (Feret length, mm)", 
        ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
        type = "l", lty = 1, cex = 15, col = "darkgreen") 
  
  
  
  legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
         text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
         merge = TRUE, bg = "white", trace=F)
  
  
  # correlation analysis
  
  Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  
  A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  
  
  A1contrib <- A1/Total.zoo
  B1contrib  <- B1/Total.zoo
  
  cor.test (A1 , B1)
  c1 <- cor.test (A1contrib , B1contrib)
  p.value_pears <- c1$p.value
  cor.value_pears <- c1$estimate
  c_sp <- cor.test (A1contrib , B1contrib, method = "spearman")
  p.value_spearm <- round( c_sp$p.value, 18)
  cor.value_rho_Spear <- round( c_sp$estimate , 3)
  
  text(2.5, 1, paste(p.value_spearm,cor.value_rho_Spear ))
  
  
  
  dev.off()
}

for (w in 39:76)  { 
  
  st.num <- w
  
  #svg(filename=  paste ("sample_estu",st.num,".svg") ) 
  
  png(filename=  paste ("sample_estu",st.num,".png") ) 
  
  # log-lin plot , NORMALIZED Abund - TOTAL ZOO
  plot(log10(1+as.numeric(totZoo[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
       xlim = c(0,5),
       xlab = "Size, (Feret length, mm)", 
       ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
       type = "l", lty = 1, cex = 15, col = "black",
       main =paste ("Estuary, sample no. ",st.num)) 
  
  # log-lin plot , NORMALIZED Abund - BRAZ
  lines(log10(1+as.numeric(BraZ[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
        xlim = c(0,5),
        xlab = "Size, (Feret length, mm)", 
        ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
        type = "l", lty = 1, cex = 15, col = "darkorange") 
  
  
  # log-lin plot , NORMALIZED Abund
  lines(log10(1+as.numeric(Cop[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
        xlim = c(0,5),
        xlab = "Size, (Feret length, mm)", 
        ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
        type = "l", lty = 1, cex = 15, col = "darkgreen") 
  
  
  
  legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
         text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
         merge = TRUE, bg = "white", trace=F)
  
  
  # correlation analysis
  
  Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  
  A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  
  
  A1contrib <- A1/Total.zoo
  B1contrib  <- B1/Total.zoo
  
  cor.test (A1 , B1)
  c1 <- cor.test (A1contrib , B1contrib)
  p.value_pears <- c1$p.value
  cor.value_pears <- c1$estimate
  c_sp <- cor.test (A1contrib , B1contrib, method = "spearman")
  p.value_spearm <- round( c_sp$p.value, 18)
  cor.value_rho_Spear <- round( c_sp$estimate , 3)
  
  text(2.5, 5, paste(p.value_spearm,cor.value_rho_Spear ))
  
  
  
  dev.off()
}

for (w in 77:121)  { 
  
  st.num <- w
  
  # svg(filename=  paste ("sample_shelf",st.num,".svg") ) 
  
  png(filename=  paste ("sample_shelf",st.num,".png") ) 
  
  
  
  # log-lin plot , NORMALIZED Abund - TOTAL ZOO
  plot(log10(1+as.numeric(totZoo[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
       xlim = c(0,5),
       xlab = "Size, (Feret length, mm)", 
       ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
       type = "l", lty = 1, cex = 15, col = "black",
       main =paste ("Shelf, sample no.",st.num) ) 
  
  # log-lin plot , NORMALIZED Abund - BRAZ
  lines(log10(1+as.numeric(BraZ[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
        xlim = c(0,5),
        xlab = "Size, (Feret length, mm)", 
        ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
        type = "l", lty = 1, cex = 15, col = "darkorange") 
  
  
  # log-lin plot , NORMALIZED Abund
  lines(log10(1+as.numeric(Cop[st.num,1:30]))/diff(X.feret.vec.mm[1:30]) ~ (X.feret.vec.mm[1:30]), 
        xlim = c(0,5),
        xlab = "Size, (Feret length, mm)", 
        ylab = "Normalized Abundance, log10(1+x, ind. m-3 mm-1)" ,
        type = "l", lty = 1, cex = 15, col = "darkgreen") 
  
  
  
  legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
         text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
         merge = TRUE, bg = "white", trace=F)
  
  
  # correlation analysis
  
  Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  
  A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  
  
  A1contrib <- A1/Total.zoo
  B1contrib  <- B1/Total.zoo
  
  cor.test (A1 , B1)
  c1 <- cor.test (A1contrib , B1contrib)
  p.value_pears <- c1$p.value
  cor.value_pears <- c1$estimate
  c_sp <- cor.test (A1contrib , B1contrib, method = "spearman")
  p.value_spearm <- round( c_sp$p.value, 18)
  cor.value_rho_Spear <- round( c_sp$estimate , 3)
  
  text(2.5, 1, paste(p.value_spearm,cor.value_rho_Spear ))
  
  
  
  
  dev.off() 
}



## 4. Table  with  values of "p" and Spearman "rho"--------------------
# identical values in the figures
# LOOP for Table --------


st.list <- 1:121

#w = 4



result.table <- data.frame (  sample.num = rep(NA, 121), 
                              sampling.area = rep(NA, 121),
                              p = NA , rho = rep(NA, 121) )


result.table$sampling.area[1:38] <- "Bay"

result.table$sampling.area[39:76] <- "Estuary"

result.table$sampling.area[77:121] <- "Shelf"


for (w in 1:121)  { 
  
  st.num <- w
  
  #  svg(filename=  paste ("sample_bay",st.num,".svg") ) 
  
  
  
  # correlation analysis
  
  Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  
  A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
  
  
  A1contrib <- A1/Total.zoo
  B1contrib  <- B1/Total.zoo
  
  cor.test (A1 , B1)
  c1 <- cor.test (A1contrib , B1contrib)
  p.value_pears <- c1$p.value
  cor.value_pears <- c1$estimate
  c_sp <- cor.test (A1contrib , B1contrib, method = "spearman")
  p.value_spearm <- round( c_sp$p.value, 18)
  cor.value_rho_Spear <- round( c_sp$estimate , 3)
  
  paste(p.value_spearm,cor.value_rho_Spear )
 
  
  result.table[ st.num, 1] <- st.num
   result.table[st.num, 3] <- p.value_spearm
  result.table[st.num, 4] <- cor.value_rho_Spear
  
    
 
}

# View(result.table)


setwd("C:/Users/RALF/Documents/Papers/0 - Paper Denise Size Spectra")

write.table (result.table, file = "result_table_spearman.csv")

# Analyse the results Table --------------------------------------------
# list significantly negative and significantly	positive results ----------


result.table$sign_negat <- result.table$rho
result.table$sign_posit <- result.table$p


result.table$sign_negat[ result.table$sign_posit > 0.05] <- NA

result.table$sign_rho_values <- result.table$sign_negat

result.table$sign_rho_values

result.table$sign_rho_values[result.table$sign_rho_values > 0] 
# Only one positive correlation

# negative correlations - All samples ------------
negativ_rhos_andNAs <- as.numeric(result.table$sign_rho_values[result.table$sign_rho_values < 0]) 
negativ_rhos_andNAs
length(negativ_rhos_andNAs) # 120 values (negative and NAs)

negativ_rhos_andNAs
sum(!is.na(negativ_rhos_andNAs)) # 34 negative signifficant values

# Bay only - negative correlations -------------------
Bay_negativ_rhos_andNAs <-   negativ_rhos_andNAs[1:38] # 38 samples

length(Bay_negativ_rhos_andNAs) # 38 values (negative and NAs)
sum(!is.na(Bay_negativ_rhos_andNAs)) # 5 negative significant values

# Estuary only - negative correlations -------------------
Estu_negativ_rhos_andNAs <-   negativ_rhos_andNAs[39:76] # 38 samples

length(Estu_negativ_rhos_andNAs) # 38 values (negative and NAs)
sum(!is.na(Estu_negativ_rhos_andNAs)) # 20 negative significant values

+20/38 # 52%


# Shelf only - negative correlations -------------------
Shelf_negativ_rhos_andNAs <-   negativ_rhos_andNAs[77:121] # 45 samples

length(Shelf_negativ_rhos_andNAs) # 45 values (negative and NAs)
sum(!is.na(Shelf_negativ_rhos_andNAs)) # 20 negative significant values

+9/45 # 20%





# positive correlations ------------
negativ_rhos_NAs_only <- is.na (negativ_rhos_NAs[negativ_rhos_NAs =! NA] )

length (negativ_rhos_NAs)
 length (negativ_rhos_NAs  [ negativ_rhos_NAs == NA])
 121-length (negativ_rhos_NAs  [ negativ_rhos_NAs == NA]) # 1
# Only one positive correlation

 
 
 
 
 # 5.correlation with All samples (n = 121) -----------------------------
 
 Total.zoo.vec <-  1:121
 A1.vec <- 1:121
 B1.vec <- 1:121
 A1contrib <- 1:121
 B1contrib <- 1:121
 
 A_contrib_matrix <- matrix ( nrow = 121, ncol = 11 , NA )
 B_contrib_matrix <- matrix( nrow = 121, ncol = 11 , NA )
 
 
# Fill the Table with standardized (% total) size spectra bins
 
  for (w in 1:121)  { 
   
   st.num <- w
  
   # st.num <- 1
   
   Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   
   A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   
   
   A1contrib <- A1/Total.zoo
   B1contrib  <- B1/Total.zoo
   
   
   A_contrib_matrix[st.num,1:11 ] <- A1contrib
   B_contrib_matrix[st.num,1:11 ] <- B1contrib
     
   
 }
      
 A_contrib_matrix
 B_contrib_matrix
 
 
 
 c_sp <- cor.test (A_contrib_matrix , B_contrib_matrix, method = "spearman")
 c_sp
 p.value_spearm <- round( c_sp$p.value, 18)
 cor.value_rho_Spear <- round( c_sp$estimate , 3)
 
 paste(p.value_spearm,cor.value_rho_Spear )
 
 

 # 6. correlation with all samples by Area ( Bay, Estuary, Shelf) -----------
 
 
 result.table$sampling.area[1:38] <- "Bay"
 
 result.table$sampling.area[39:76] <- "Estuary"
 
 result.table$sampling.area[77:121] <- "Shelf"
 
 
 
 # Bay -----------------
 
 Total.zoo.vec <-  1:38
 A1.vec <- 1:38
 B1.vec <- 1:38
 A1contrib <- 1:38
 B1contrib <- 1:38
 
 A_contrib_matrix <- matrix ( nrow = 38, ncol = 11 , NA )
 B_contrib_matrix <- matrix( nrow = 38, ncol = 11 , NA )
 
 
 # Fill the Table with standardized (% total) size spectra bins
 
 for (w in 1:38)  { 
   
   st.num <- w
   
   # st.num <- 1
   
   Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   
   A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   
   
   A1contrib <- A1/Total.zoo
   B1contrib  <- B1/Total.zoo
   
   
   A_contrib_matrix[st.num,1:11 ] <- A1contrib
   B_contrib_matrix[st.num,1:11 ] <- B1contrib
   
   
 }
 
 A_contrib_matrix
 B_contrib_matrix
 
 
 
 c_sp <- cor.test (A_contrib_matrix , B_contrib_matrix, method = "spearman")
 c_sp
 p.value_spearm <- round( c_sp$p.value, 18)
 cor.value_rho_Spear <- round( c_sp$estimate , 3)
 
 paste(p.value_spearm,cor.value_rho_Spear )
 # p-value = 0.0028
 
 # Estuary -----------------
 
 
 length(39:76) # 38 samples 
 
 A_contrib_matrix <- matrix ( nrow = 121, ncol = 11 , NA )
 B_contrib_matrix <- matrix( nrow = 121, ncol = 11 , NA )
 
 
 # Fill the Table with standardized (% total) size spectra bins
 
 for (w in 39:76)  { 
   
   st.num <- w
   
   # st.num <- 39
   
   Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   
   A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   
   
   A1contrib <- A1/Total.zoo
   B1contrib  <- B1/Total.zoo
   
   
   A_contrib_matrix[st.num,1:11 ] <- A1contrib
   B_contrib_matrix[st.num,1:11 ] <- B1contrib
   
   
 }
 
 A_contrib_matrix
 B_contrib_matrix
 
 
 
 c_sp <- cor.test (A_contrib_matrix[39:76,1:11 ] , B_contrib_matrix[39:76,1:11], method = "spearman")
 c_sp
 p.value_spearm <- round( c_sp$p.value, 18)
 cor.value_rho_Spear <- round( c_sp$estimate , 3)
 
 paste(p.value_spearm,cor.value_rho_Spear )
 
 # Shelf -----------------
 
 
 length(77:121) # 45 samples 
 
 A_contrib_matrix <- matrix ( nrow = 121, ncol = 11 , NA )
 B_contrib_matrix <- matrix( nrow = 121, ncol = 11 , NA )
 
 
 # Fill the Table with standardized (% total) size spectra bins
 
 for (w in 77:121)  { 
   
   st.num <- w
   
   # st.num <- 39
   
   Total.zoo  <- log10(1+as.numeric(totZoo[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   
   A1 <- log10(1+as.numeric(Cop[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   B1 <- log10(1+as.numeric(BraZ[st.num,5:15]))/diff(X.feret.vec.mm[5:16])
   
   
   A1contrib <- A1/Total.zoo
   B1contrib  <- B1/Total.zoo
   
   
   A_contrib_matrix[st.num,1:11 ] <- A1contrib
   B_contrib_matrix[st.num,1:11 ] <- B1contrib
   
   
 }
 
 A_contrib_matrix
 B_contrib_matrix
 
 
 
 c_sp <- cor.test (A_contrib_matrix[77:121,1:11 ] , B_contrib_matrix[77:121,1:11], method = "spearman")
 c_sp
 p.value_spearm <- round( c_sp$p.value, 18)
 cor.value_rho_Spear <- round( c_sp$estimate , 3)
 
 paste(p.value_spearm,cor.value_rho_Spear )
 
 
 # Abundance vs rho ---------------------------------
 
 rs <- read.csv("~/Papers/0 - Paper Denise Size Spectra/result_table_spearman_2.csv")
#   View(result_table_spearman_2)
 
 
 attach (rs)
 
 
 plot(rs$rho ~ rs$BrZ_Abundind.m3)
 abline (v = 0.3)
 
 plot ((rs$rho) ~ log(1+rs$BrZ_Abundind.m3))
 abline (v = log(1+0.3))
 lmx <- lm ((rs$rho) ~ log(1+rs$BrZ_Abundind.m3))
 abline (lmx)
 summary (lmx)
 
 # lmperm ------------
 library(lmPerm)
 # 
 # lmp {lmPerm}
 
 lmxp <- lmp ((rs$rho) ~ log(1+rs$BrZ_Abundind.m3))
 abline (lmxp)
 summary (lmxp)
 

 lmxpcop <- lmp ((rs$rho) ~ log10(1+rs$Cop_Abundance_ind.m3))
 
 summary (lmxpcop)# n.s.
 
 ## NICE plot fig 8-------------------------------
 
 library(ggplot2)
 library(cars)
 
 
rs$Coplogabund <- log10(1+rs$Cop_Abundance_ind.m3)
rs$BRZlogabund <- log10(1+rs$BrZ_Abundind.m3)

 


theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

p <- ggplot(rs, aes(BrZ_Abundind.m3, rho)) +
  geom_point()
p


 p <- ggplot(rs, aes(BRZlogabund, rho)) +
   geom_point()
 # Add regression line
 p + geom_smooth(method = lm)
 
 
 
 
 
 theme_set(
   theme_bw() +
     theme(legend.position = "top")
 )
 
 p <- ggplot(rs, aes(Cop_Abundance_ind.m3, rho)) +
   geom_point()
 p
 
 
 
 p <- ggplot(rs, aes(Coplogabund, rho)) +
   geom_point()
 # Add regression line
 p + geom_smooth(method = lm)
 
 
 