# Size Spectra analysis with different bin sizes
#
# for the manuscript:
#" Size niche interactions (functional replacement and competitive exclusion)
# between mero- and holoplankton shape the size spectrum 
#of tropical estuarine and marine ecosystems "
#
#  Feret Size vs Abundance

# v. 3.00 - multiple  plots  - simple tests
# Ralf Schwamborn
# August 16, 2025


#Inputs: ----------------------
  
# 4 files head and log files, large and small fraction, total fractions
# x vec files (bins)
  
# Outputs:  ----------
#graphs of Size Spectra for Brahchyuran zoeae, copepods, Total  
# with differet bn sizs 

# Procedures (workflow):     --------
# 1. read in 4 files: head and log files, large and small fraction  ----------
# 2. calculate different bin vectors -----------------
# 3. sum : total fraction ( large fraction +  small fraction) ----------
# 4. plot ------------------
# 5. repeat and plot with different bins ------------------



### 0. load libraries        

library(stringr)
library(scales)


# set working directory

setwd("~/Papers/0 - Paper Denise Size Spectra")

setwd("~/Papers/0 - Paper Denise Size Spectra/PIDS_OK_NEW_W_TAXA_HEAD_LOG")


### I. Basic Feret Size vs Abundance analysis and one plot

### 1. read in 4 files: head and log files, large and small fraction  ----------

menor.meta.name <- "estuario2_300_07agosto2013_600ml_tot_1_meta.txt"


menor.valid.name <- "Valid_0125_estuario2_300_07agosto2013_600ml.txt"


menor.meta <- readLines(menor.meta.name)

menor.valid <- readLines(menor.valid.name)

head(menor.valid)


df.all <- read.table(menor.valid.name,
                 sep = "\t",                # tab-separated
                 header = TRUE,             # header present
                 skip = 1,                  # skip first line (filename)
                 fill = TRUE,               # fill uneven rows
                 quote = "",                # prevent "!" from being treated as comment
                 check.names = FALSE)       # keep original column names

names(df.all)

dim(df.all) # 988 objects

df.all$Valid <- as.factor(df.all$Valid )



df.all$Fsize2_micron <- df.all$Feret * 10.58


summary(df.all$Valid)
# View (df.all)

df.cop <- subset (df.all, Valid == "03_Calanoida")
df.BracZ <- subset (df.all, Valid == "04_Brachyura_zoe")


dim(df.cop) # 174 calanoida
dim(df.BracZ)  # 345 Br. zoeae


# View (df.all)

# Read metadata -------------------------
# read Vol,  SubPart, and Station ID ------------------


metadata <- menor.meta

lines <- metadata

# split by '=' and trim spaces
kv <- strsplit(lines, "=")
kv <- lapply(kv, trimws)

# convert to named vector
data_list <- setNames(sapply(kv, `[`, 2), sapply(kv, `[`, 1))

data_list <-  as.list (data_list)

names(data_list)

Vol <-  as.numeric(data_list$Vol)

#SubPart <-  as.numeric(data_list$..)
SubPart = 1700/10

SampleId <- data_list$SampleId   #  "estuario1_300_19janeiro2015_1700ml"

 StationId <-  data_list$StationId  # 01


### 2. calculate different bin vectors -----------------


char.X.feret.vec.mm <- "0.29,	0.33196,	0.38,	0.43498,	0.49793,	0.56998	0.65245	0.74686	0.85493	0.97864	1.1202	1.2823	1.4679	1.6803	1.9235	2.2018	2.5204	2.8851	3.3025	3.7804	4.3274	4.9536	5.6704	6.4909	7.4302	8.5053	9.736	11.145	12.757	14.603	16.717	19.136	21.904	25.074	28.702"

X.feret.vec.mm <-  str_replace_all(char.X.feret.vec.mm, "\t", "," )

X.feret.vec.mm <-c(0.29,0.33196,0.38,0.43498,0.49793,0.56998,0.65245,0.74686,0.85493,0.97864,1.1202,1.2823,1.4679,1.6803,1.9235,2.2018,2.5204,2.8851,3.3025,3.7804,4.3274,4.9536,5.6704,6.4909,7.4302,8.5053,9.736,11.145,12.757,14.603,16.717,19.136,21.904,25.074,28.702)


X.feret.vec.micron <-  X.feret.vec.mm * 1000



### 3. sum : histogram , SubPart and Vol  ( large fraction +  small fraction) ----------

Fsize2_micron.all <- df.all$Fsize2_micron
Fsize2_micron.cop <- df.cop$Fsize2_micron
Fsize2_micron.BrZ <- df.BracZ$Fsize2_micron

max_size <- max(Fsize2_micron.all,Fsize2_micron.cop, Fsize2_micron.BrZ )

bins.lin.simple.microns.1 <- seq(0, (100+ max_size))
length(bins.lin.simple.microns.1)
summary(bins.lin.simple.microns.1)
diff(bins.lin.simple.microns.1)


all.hist <- hist(Fsize2_micron.all, breaks = seq(0, (100+ max_size), by = 8))

cop.hist <- hist(Fsize2_micron.cop, breaks = seq(0, (100+ max_size), by = 8))

BrZ.hist <- hist(Fsize2_micron.BrZ, breaks = seq(0, (100+ max_size), by = 8))


##bins.lin.simple.microns.1
mids.lin.simple.microns.1 <- all.hist$mids
diff(mids.lin.simple.microns.1)

#cop.hist <- hist(Fsize2_micron.cop, breaks = seq(0, (100+ max_size), by = 12))

#BrZ.hist <- hist(Fsize2_micron.BrZ, breaks = seq(0, (100+ max_size), by = 12))


# Multiply by SubPart and divide by Vol  = Abundance

all.abund <- (all.hist$counts * SubPart) / Vol

cop.abund <- (cop.hist$counts * SubPart) / Vol

BrZ.abund <- (BrZ.hist$counts * SubPart) / Vol


### 4. plot ------------------



# linear and Log-linear plots -----------

# Estuary ----------------

# linear plot, not NORMALIZED Abund - TOTAL ZOO, Brac.Z, Cop
plot((as.numeric(all.abund)) ~ (mids.lin.simple.microns.1), 
     xlim = c(0,3000), # ylim= c(0,28),
     xlab = "Size, (Feret length, micron)", 
     ylab = "Abundance (ind. m-3)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main = paste ("bin size = 8 micron, Estuary, Station" ,   StationId) )

lines((as.numeric(BrZ.abund)) ~ (mids.lin.simple.microns.1),
        type = "l", lty = 1, cex = 15, col = "darkorange") 

lines((as.numeric(cop.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)



# log-lin plot, Abund - TOTAL ZOO, Brac.Z, Cop
plot((log10(1+all.abund)) ~ (mids.lin.simple.microns.1), 
     xlim = c(0,4000), # ylim= c(0,28),
     xlab = "Size, (Feret length, micron)", 
     ylab = "log10(1+Abundance (ind. m-3 mm-1))" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Estuary") 

lines((log10(1+BrZ.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkorange") 

lines((log10(1+cop.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)


cor.test(cop.abund, BrZ.abund)

cor.test(cop.abund, BrZ.abund,method = "spearman")#  p-value < 2.2e-16

plot(cop.abund ~ BrZ.abund)


# II. repeat and plot with different bins ------------------

# II.1 vec10 ----------

# bin size: 10 micron


char.X.feret.vec.mm <- "0.29,	0.33196,	0.38,	0.43498,	0.49793,	0.56998	0.65245	0.74686	0.85493	0.97864	1.1202	1.2823	1.4679	1.6803	1.9235	2.2018	2.5204	2.8851	3.3025	3.7804	4.3274	4.9536	5.6704	6.4909	7.4302	8.5053	9.736	11.145	12.757	14.603	16.717	19.136	21.904	25.074	28.702"

X.feret.vec.mm <-  str_replace_all(char.X.feret.vec.mm, "\t", "," )

X.feret.vec.mm <-c(0.29,0.33196,0.38,0.43498,0.49793,0.56998,0.65245,0.74686,0.85493,0.97864,1.1202,1.2823,1.4679,1.6803,1.9235,2.2018,2.5204,2.8851,3.3025,3.7804,4.3274,4.9536,5.6704,6.4909,7.4302,8.5053,9.736,11.145,12.757,14.603,16.717,19.136,21.904,25.074,28.702)


X.feret.vec.micron <-  X.feret.vec.mm * 1000



### 3. sum : histogram , SubPart and Vol  ( large fraction +  small fraction) ----------

Fsize2_micron.all <- df.all$Fsize2_micron
Fsize2_micron.cop <- df.cop$Fsize2_micron
Fsize2_micron.BrZ <- df.BracZ$Fsize2_micron

max_size <- max(Fsize2_micron.all,Fsize2_micron.cop, Fsize2_micron.BrZ )

bins.lin.simple.microns.1 <- seq(0, (100+ max_size))
length(bins.lin.simple.microns.1)
summary(bins.lin.simple.microns.1)
diff(bins.lin.simple.microns.1)


all.hist <- hist(Fsize2_micron.all, breaks = seq(0, (100+ max_size), by = 10))

cop.hist <- hist(Fsize2_micron.cop, breaks = seq(0, (100+ max_size), by = 10))

BrZ.hist <- hist(Fsize2_micron.BrZ, breaks = seq(0, (100+ max_size), by = 10))


##bins.lin.simple.microns.1
mids.lin.simple.microns.1 <- all.hist$mids
diff(mids.lin.simple.microns.1)

#cop.hist <- hist(Fsize2_micron.cop, breaks = seq(0, (100+ max_size), by = 12))

#BrZ.hist <- hist(Fsize2_micron.BrZ, breaks = seq(0, (100+ max_size), by = 12))


# Multiply by SubPart and divide by Vol  = Abundance

all.abund <- (all.hist$counts * SubPart) / Vol

cop.abund <- (cop.hist$counts * SubPart) / Vol

BrZ.abund <- (BrZ.hist$counts * SubPart) / Vol


### 4. plot ------------------



# linear and Log-linear plots -----------

# Estuary ----------------

# linear plot, not NORMALIZED Abund - TOTAL ZOO, Brac.Z, Cop
plot((as.numeric(all.abund)) ~ (mids.lin.simple.microns.1), 
     xlim = c(0,3000), # ylim= c(0,28),
     xlab = "Size, (Feret length, micron)", 
     ylab = "Abundance (ind. m-3)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main = paste ("bin size = 10 micron, Estuary, Station" ,   StationId) )

lines((as.numeric(BrZ.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkorange") 

lines((as.numeric(cop.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)



# log-lin plot, Abund - TOTAL ZOO, Brac.Z, Cop
plot((log10(1+all.abund)) ~ (mids.lin.simple.microns.1), 
     xlim = c(0,4000), # ylim= c(0,28),
     xlab = "Size, (Feret length, micron)", 
     ylab = "log10(1+Abundance (ind. m-3 mm-1))" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Estuary") 

lines((log10(1+BrZ.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkorange") 

lines((log10(1+cop.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)


cor.test(cop.abund, BrZ.abund)

cor.test(cop.abund, BrZ.abund,method = "spearman")#  p-value < 2.2e-16

plot(cop.abund ~ BrZ.abund)




# II.2 vec12 ----------

# bin size: 12 micron


char.X.feret.vec.mm <- "0.29,	0.33196,	0.38,	0.43498,	0.49793,	0.56998	0.65245	0.74686	0.85493	0.97864	1.1202	1.2823	1.4679	1.6803	1.9235	2.2018	2.5204	2.8851	3.3025	3.7804	4.3274	4.9536	5.6704	6.4909	7.4302	8.5053	9.736	11.145	12.757	14.603	16.717	19.136	21.904	25.074	28.702"

X.feret.vec.mm <-  str_replace_all(char.X.feret.vec.mm, "\t", "," )

X.feret.vec.mm <-c(0.29,0.33196,0.38,0.43498,0.49793,0.56998,0.65245,0.74686,0.85493,0.97864,1.1202,1.2823,1.4679,1.6803,1.9235,2.2018,2.5204,2.8851,3.3025,3.7804,4.3274,4.9536,5.6704,6.4909,7.4302,8.5053,9.736,11.145,12.757,14.603,16.717,19.136,21.904,25.074,28.702)


X.feret.vec.micron <-  X.feret.vec.mm * 1000



### 3. sum : histogram , SubPart and Vol  ( large fraction +  small fraction) ----------

Fsize2_micron.all <- df.all$Fsize2_micron
Fsize2_micron.cop <- df.cop$Fsize2_micron
Fsize2_micron.BrZ <- df.BracZ$Fsize2_micron

max_size <- max(Fsize2_micron.all,Fsize2_micron.cop, Fsize2_micron.BrZ )

bins.lin.simple.microns.1 <- seq(0, (100+ max_size))
length(bins.lin.simple.microns.1)
summary(bins.lin.simple.microns.1)
diff(bins.lin.simple.microns.1)


all.hist <- hist(Fsize2_micron.all, breaks = seq(0, (100+ max_size), by = 12))

cop.hist <- hist(Fsize2_micron.cop, breaks = seq(0, (100+ max_size), by = 12))

BrZ.hist <- hist(Fsize2_micron.BrZ, breaks = seq(0, (100+ max_size), by = 12))


##bins.lin.simple.microns.1
mids.lin.simple.microns.1 <- all.hist$mids
diff(mids.lin.simple.microns.1)

#cop.hist <- hist(Fsize2_micron.cop, breaks = seq(0, (100+ max_size), by = 12))

#BrZ.hist <- hist(Fsize2_micron.BrZ, breaks = seq(0, (100+ max_size), by = 12))


# Multiply by SubPart and divide by Vol  = Abundance

all.abund <- (all.hist$counts * SubPart) / Vol

cop.abund <- (cop.hist$counts * SubPart) / Vol

BrZ.abund <- (BrZ.hist$counts * SubPart) / Vol


### 4. plot ------------------



# linear and Log-linear plots -----------

# Estuary ----------------

# linear plot, not NORMALIZED Abund - TOTAL ZOO, Brac.Z, Cop
plot((as.numeric(all.abund)) ~ (mids.lin.simple.microns.1), 
     xlim = c(0,3000), # ylim= c(0,28),
     xlab = "Size, (Feret length, micron)", 
     ylab = "Abundance (ind. m-3)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main = paste ("bin size = 12 micron, Estuary, Station" ,   StationId) )

lines((as.numeric(BrZ.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkorange") 

lines((as.numeric(cop.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)



# log-lin plot, Abund - TOTAL ZOO, Brac.Z, Cop
plot((log10(1+all.abund)) ~ (mids.lin.simple.microns.1), 
     xlim = c(0,4000), # ylim= c(0,28),
     xlab = "Size, (Feret length, micron)", 
     ylab = "log10(1+Abundance (ind. m-3 mm-1))" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Estuary") 

lines((log10(1+BrZ.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkorange") 

lines((log10(1+cop.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)


cor.test(cop.abund, BrZ.abund)

cor.test(cop.abund, BrZ.abund,method = "spearman")#  p-value < 2.2e-16

plot(cop.abund ~ BrZ.abund)




# II.3 vec14 ----------

# bin size: 14 micron


char.X.feret.vec.mm <- "0.29,	0.33196,	0.38,	0.43498,	0.49793,	0.56998	0.65245	0.74686	0.85493	0.97864	1.1202	1.2823	1.4679	1.6803	1.9235	2.2018	2.5204	2.8851	3.3025	3.7804	4.3274	4.9536	5.6704	6.4909	7.4302	8.5053	9.736	11.145	12.757	14.603	16.717	19.136	21.904	25.074	28.702"

X.feret.vec.mm <-  str_replace_all(char.X.feret.vec.mm, "\t", "," )

X.feret.vec.mm <-c(0.29,0.33196,0.38,0.43498,0.49793,0.56998,0.65245,0.74686,0.85493,0.97864,1.1202,1.2823,1.4679,1.6803,1.9235,2.2018,2.5204,2.8851,3.3025,3.7804,4.3274,4.9536,5.6704,6.4909,7.4302,8.5053,9.736,11.145,12.757,14.603,16.717,19.136,21.904,25.074,28.702)


X.feret.vec.micron <-  X.feret.vec.mm * 1000



### 3. sum : histogram , SubPart and Vol  ( large fraction +  small fraction) ----------

Fsize2_micron.all <- df.all$Fsize2_micron
Fsize2_micron.cop <- df.cop$Fsize2_micron
Fsize2_micron.BrZ <- df.BracZ$Fsize2_micron

max_size <- max(Fsize2_micron.all,Fsize2_micron.cop, Fsize2_micron.BrZ )

bins.lin.simple.microns.1 <- seq(0, (100+ max_size))
length(bins.lin.simple.microns.1)
summary(bins.lin.simple.microns.1)
diff(bins.lin.simple.microns.1)


all.hist <- hist(Fsize2_micron.all, breaks = seq(0, (100+ max_size), by = 14))

cop.hist <- hist(Fsize2_micron.cop, breaks = seq(0, (100+ max_size), by = 14))

BrZ.hist <- hist(Fsize2_micron.BrZ, breaks = seq(0, (100+ max_size), by = 14))


##bins.lin.simple.microns.1
mids.lin.simple.microns.1 <- all.hist$mids
diff(mids.lin.simple.microns.1)

#cop.hist <- hist(Fsize2_micron.cop, breaks = seq(0, (100+ max_size), by = 12))

#BrZ.hist <- hist(Fsize2_micron.BrZ, breaks = seq(0, (100+ max_size), by = 12))


# Multiply by SubPart and divide by Vol  = Abundance

all.abund <- (all.hist$counts * SubPart) / Vol

cop.abund <- (cop.hist$counts * SubPart) / Vol

BrZ.abund <- (BrZ.hist$counts * SubPart) / Vol


### 4. plot ------------------



# linear and Log-linear plots -----------

# Estuary ----------------

# linear plot, not NORMALIZED Abund - TOTAL ZOO, Brac.Z, Cop
plot((as.numeric(all.abund)) ~ (mids.lin.simple.microns.1), 
     xlim = c(0,3000), # ylim= c(0,28),
     xlab = "Size, (Feret length, micron)", 
     ylab = "Abundance (ind. m-3)" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main = paste ("bin size = 14 micron, Estuary, Station" ,   StationId) )

lines((as.numeric(BrZ.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkorange") 

lines((as.numeric(cop.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)



# log-lin plot, Abund - TOTAL ZOO, Brac.Z, Cop
plot((log10(1+all.abund)) ~ (mids.lin.simple.microns.1), 
     xlim = c(0,4000), # ylim= c(0,28),
     xlab = "Size, (Feret length, micron)", 
     ylab = "log10(1+Abundance (ind. m-3 mm-1))" ,
     type = "l", lty = 1, cex = 15, col = "black",
     main ="Estuary") 

lines((log10(1+BrZ.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkorange") 

lines((log10(1+cop.abund)) ~ (mids.lin.simple.microns.1),
      type = "l", lty = 1, cex = 15, col = "darkgreen") 

legend("topright", c("total Zoopl.", "Copep.", "Brac. Z."), col = c("black", "darkgreen", "darkorange"),
       text.col = "darkgrey", lty = c(1, 1, 1), pch = c(NA, NA, NA),
       merge = TRUE, bg = "white", trace=F)


cor.test(cop.abund, BrZ.abund)

cor.test(cop.abund, BrZ.abund,method = "spearman")#  p-value < 2.2e-16

plot(cop.abund ~ BrZ.abund)


