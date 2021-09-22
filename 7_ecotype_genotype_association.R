#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("psych")
#install.packages("qqman")
#install.packages("raster")
#install.packages("lfmm")
#install.packages("flashpcaR")
#install.packages("gradientForest")
#install.packages("rasterVis")

#install.packages("devtools")
library("devtools")
#devtools::install_github("bcm-uga/lfmm")

#install.packages("remotes")
remotes::install_github("gabraham/flashpca/flashpcaR")
#install.packages("C:/Users/gzl02/Downloads/flashpcaR_1.2.5")




library(tidyverse) # For data processing
library(flashpcaR) # For scaling genetic data
library(lfmm)      # For LFMM
library(vegan)     # For RDA
library(psych)     # For correlation plots
library(qqman)     # For manhattan plots
library(gradientForest) # For gradient forest analysis
library(raster)    # For working with rasters
library(rasterVis) # For visualizing rasters
#ap <- available.packages("lfmm")


setwd("C:/Users/gzl02/Desktop/Yi2019Research/Research materials")

# Use read_tsv function to read in tab-separated files
# We also assign our own column names since the .012.pos file does not include column names
snp_pos <- readr::read_tsv("out.012.pos",
                           col_names = c("chrom", "pos"))

# Provide column names since they are not included in vcftools output
# First column is sequence of numbers for each individual, which we will remove later
col_names <- c("NUM", 
               paste0(snp_pos$chrom, "_", snp_pos$pos))
# paste0 joins together strings with no spaces in between

col_names[1:10] # Take a look at the output

# Read in genotype matrix
genotypes <- readr::read_tsv("out.012",
                             col_names = col_names, # Set our column names
                             na = "-1") # Reads in -1 values as missing data

# Remove first column

genotypes[, 1:5] # See what the data looks like


indv <- readr::read_tsv("out.012.indv", 
                        col_names = "ID") # Provide our own column name
indv # Should be 126 rows

# Join individual names and genotype data- Must be in same row order!!
gen_dat <- dplyr::bind_cols(indv, genotypes)

#read the table with climate data and the locality
#df <- read.csv("All_variables_A_glandulosa1.csv", header = TRUE)
#str(df)
clim_dat <- readr::read_csv("Selected_variables_A_glandulosa1.csv")
clim_dat
str(clim_dat)
# Create locality/population column that we will use later on
#gen_dat$Locality <- gsub('[0-9]+', '', gen_dat$ID)
#gen_dat$Locality <- df$Locality


# Make sure sample IDs between genetic and climate dataframes are in same order
all(gen_dat$ID == clim_dat$ID) # Should return TRUE

ifelse(gen_dat$ID==clim_dat$ID,0,NA)


# Plot correlations among climate variables
psych::pairs.panels(clim_dat[, -c(1)], # Don't include ID column
                    ellipses = FALSE, lm = TRUE)
x <- clim_dat[, -c(1)]
c <- cor(x,x)
str(c)
c <- data.frame(c)
write.csv(cor(x,x), "correlation_between_soil_bioclim_solar.csv")
str(clim_dat)

#define the function of scale2
scale2 <- function(X, type=c("2", "1"), impute=TRUE)
{
  type <- match.arg(type)
  mult <- ifelse(type == "1", 1, 2)
  
  sum2 <- nrow(X) - colSums(apply(X, 2, is.na))
  p <- colSums(X, na.rm=TRUE) / (2 * sum2)
  xsd <- sqrt(mult * p * (1 - p))
  names(p) <- names(xsd) <- colnames(X)
  
  s <- sweep(
    sweep(X, MARGIN=2, STATS=2 * p, FUN="-"),
    MARGIN=2, STATS=xsd, FUN="/"
  )
  if(impute) {
    s[is.na(s)] <- 0
  }
  attr(s, "scaled:center") <- 2 * p
  attr(s, "scaled:scale") <- xsd
  s
}

# Use flashpcaR package function to scale and impute missing data
#gen_dat_scaled <- flashpcaR::scale2(gen_dat[, -c(1,2)],
#                                    impute = TRUE) # Impute missing data
gen_dat_scaled <- scale2(gen_dat[, -c(1,2)],
                                    impute = TRUE) # Impute missing data

colMeans(gen_dat_scaled[, 1:10]) # Make sure column means are 0 for the first few columns

# Scale and center climate data and convert back to data frame
clim_dat_scaled <- as.data.frame(scale(clim_dat[, -1]))

# Check col means are 0
colMeans(clim_dat_scaled)
apply(clim_dat_scaled, MARGIN = 2, sd) # SDs should all equal 1

# K = 1

# Estimate relationships
lfmm_out_k1 <- lfmm::lfmm_ridge(Y = gen_dat_scaled, # Response variable
                                X = clim_dat_scaled, # Climate + geographic variables
                                K = 1) # 
# Run statistical tests
pv_k1 <- lfmm::lfmm_test(Y = gen_dat_scaled,
                         X = clim_dat_scaled,
                         lfmm = lfmm_out_k1)

# Plot histogram of p values for each climate variable
pval_df_k1 <- as.data.frame(pv_k1$calibrated.pvalue) # Extract p vals and convert

tidyr::gather(pval_df_k1,  # Convert from wide to long format
              key = "clim_var", value = "p_val") %>%
  ggplot(. , aes(x = p_val)) +
  facet_wrap(~clim_var) +
  geom_histogram(fill = "skyblue2", color = "black", binwidth = 0.05) + 
  ggtitle("K = 1") +
  theme_bw()

#K = 2
# Estimate relationships
lfmm_out_k2 <- lfmm::lfmm_ridge(Y = gen_dat_scaled, # Response variable
                                X = clim_dat_scaled, # Climate + geographic variables
                                K = 2)
# Run statistical tests
pv_k2 <- lfmm::lfmm_test(Y = gen_dat_scaled,
                         X = clim_dat_scaled,
                         lfmm = lfmm_out_k2)

# Plot histogram of p values for each climate variable
pval_df_k2 <- as.data.frame(pv_k2$calibrated.pvalue)

tidyr::gather(pval_df_k2,  # Convert from wide to long format
              key = "clim_var", value = "p_val") %>%
  ggplot(. , aes(x = p_val)) +
  facet_wrap(~clim_var) +
  geom_histogram(fill = "skyblue2", color = "black", binwidth = 0.05) + 
  ggtitle("K = 2") +
  theme_bw()

# K =3
# Estimate relationships
lfmm_out_k3 <- lfmm::lfmm_ridge(Y = gen_dat_scaled, # Response variable
                                X = clim_dat_scaled, # Climate + geographic variables
                                K = 3) # 
# Run statistical tests
pv_k3 <- lfmm::lfmm_test(Y = gen_dat_scaled,
                         X = clim_dat_scaled,
                         lfmm = lfmm_out_k3)

# Plot histogram of p values for each climate variable
pval_df_k3 <- as.data.frame(pv_k3$calibrated.pvalue)

tidyr::gather(pval_df_k3,  # Convert from wide to long format
              key = "clim_var", value = "p_val") %>%
  ggplot(. , aes(x = p_val)) +
  facet_wrap(~clim_var) +
  geom_histogram(fill = "skyblue2", color = "black", binwidth = 0.05) + 
  ggtitle("K = 3") +
  theme_bw()


# K =4
# Estimate relationships
lfmm_out_k4 <- lfmm::lfmm_ridge(Y = gen_dat_scaled, # Response variable
                                X = clim_dat_scaled, # Climate + geographic variables
                                K = 4) # 
# Run statistical tests
pv_k4 <- lfmm::lfmm_test(Y = gen_dat_scaled,
                         X = clim_dat_scaled,
                         lfmm = lfmm_out_k4)

# Plot histogram of p values for each climate variable
pval_df_k4 <- as.data.frame(pv_k4$calibrated.pvalue)

tidyr::gather(pval_df_k4,  # Convert from wide to long format
              key = "clim_var", value = "p_val") %>%
  ggplot(. , aes(x = p_val)) +
  facet_wrap(~clim_var) +
  geom_histogram(fill = "skyblue2", color = "black", binwidth = 0.05) + 
  ggtitle("K = 4") +
  theme_bw()

# K =5
# Estimate relationships
lfmm_out_k5 <- lfmm::lfmm_ridge(Y = gen_dat_scaled, # Response variable
                                X = clim_dat_scaled, # Climate + geographic variables
                                K = 5) # 
# Run statistical tests
pv_k5 <- lfmm::lfmm_test(Y = gen_dat_scaled,
                         X = clim_dat_scaled,
                         lfmm = lfmm_out_k5)

# Plot histogram of p values for each climate variable
pval_df_k5 <- as.data.frame(pv_k5$calibrated.pvalue)

tidyr::gather(pval_df_k5,  # Convert from wide to long format
              key = "clim_var", value = "p_val") %>%
  ggplot(. , aes(x = p_val)) +
  facet_wrap(~clim_var) +
  geom_histogram(fill = "skyblue2", color = "black", binwidth = 0.05) + 
  ggtitle("K = 5") +
  theme_bw()


# First, select the variable we are interested in
# You can change this to your liking to explore relationships with various climate variables
clim_var = "soilval"

str(clim_dat)
# Create dataframe to use in manhattan plot
# Needs columns BP, CHR, P, SNP
manhattan_df <- data.frame(CHR = snp_pos$chrom, 
                           BP = snp_pos$pos,
                           P = pval_df_k2[, clim_var],
                           SNP = rownames(pval_df_k2))
head(manhattan_df)
str(manhattan_df)
summary(manhattan_df)
which(manhattan_df$P < 1e-5)
soilval_associated_SNPs <- manhattan_df[which(manhattan_df$P < 1e-5), ]
write.csv(soilval_associated_SNPs, "soilph_associated_SNPs.csv")

# Make vector of chromosome names and format for input into manhattan plot function
#chr_names <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6","chr7","chr8","chr9", "chr10", "chr11", "chr12", "Scq3eQI_1866", "Scq3eQI_1901", "Scq3eQI_1999", "Scq3eQI_2008", "Scq3eQI_2022")

#manhattan_df$CHR <- factor(manhattan_df$CHR, levels =  chr_names)
manhattan_df$CHR <- factor(manhattan_df$CHR)
manhattan_df$CHR <- as.numeric(manhattan_df$CHR)


# Create manhattan plot
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = chr_names)
windows()
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "soilpH")
manhattan(manhattan_df, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Soil pH")

#For bio 3 
clim_var = "bio3val"

str(clim_dat)
# Create dataframe to use in manhattan plot
# Needs columns BP, CHR, P, SNP
manhattan_df <- data.frame(CHR = snp_pos$chrom, 
                           BP = snp_pos$pos,
                           P = pval_df_k2[, clim_var],
                           SNP = rownames(pval_df_k2))
head(manhattan_df)
str(manhattan_df)
summary(manhattan_df)
which(manhattan_df$P < 1e-5)
soilval_associated_SNPs <- manhattan_df[which(manhattan_df$P < 1e-5), ]
write.csv(soilval_associated_SNPs, "bio3_SNPs.csv")
str(gen_dat_scaled)
# Make vector of chromosome names and format for input into manhattan plot function
#chr_names <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6","chr7","chr8","chr9", "chr10", "chr11", "chr12", "Scq3eQI_1866", "Scq3eQI_1901", "Scq3eQI_1999", "Scq3eQI_2008", "Scq3eQI_2022")

#manhattan_df$CHR <- factor(manhattan_df$CHR, levels =  chr_names)
manhattan_df$CHR <- factor(manhattan_df$CHR)
manhattan_df$CHR <- as.numeric(manhattan_df$CHR)


# Create manhattan plot
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = chr_names)
windows()
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio3 Isothermality")
manhattan(manhattan_df, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio3 Isothermality")


#For bio 5 
clim_var = "bio5val"
# Create dataframe to use in manhattan plot
# Needs columns BP, CHR, P, SNP
manhattan_df <- data.frame(CHR = snp_pos$chrom, 
                           BP = snp_pos$pos,
                           P = pval_df_k2[, clim_var],
                           SNP = rownames(pval_df_k2))
which(manhattan_df$P < 1e-5)
soilval_associated_SNPs <- manhattan_df[which(manhattan_df$P < 1e-5), ]
write.csv(soilval_associated_SNPs, "bio5_SNPs.csv")
str(gen_dat_scaled)
# Make vector of chromosome names and format for input into manhattan plot function
manhattan_df$CHR <- factor(manhattan_df$CHR)
manhattan_df$CHR <- as.numeric(manhattan_df$CHR)
# Create manhattan plot
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = chr_names)
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio5  Max Temperature of Warmest Month")
manhattan(manhattan_df, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio5  Max Temperature of Warmest Month")



#For bio 9 
clim_var = "bio9val"
# Create dataframe to use in manhattan plot
# Needs columns BP, CHR, P, SNP
manhattan_df <- data.frame(CHR = snp_pos$chrom, 
                           BP = snp_pos$pos,
                           P = pval_df_k2[, clim_var],
                           SNP = rownames(pval_df_k2))
which(manhattan_df$P < 1e-5)
soilval_associated_SNPs <- manhattan_df[which(manhattan_df$P < 1e-5), ]
write.csv(soilval_associated_SNPs, "bio9_SNPs.csv")
# Make vector of chromosome names and format for input into manhattan plot function
manhattan_df$CHR <- factor(manhattan_df$CHR)
manhattan_df$CHR <- as.numeric(manhattan_df$CHR)
# Create manhattan plot
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = chr_names)
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio9 Mean Temperature of Driest Quarter")
manhattan(manhattan_df, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio9 Mean Temperature of Driest Quarter")



#For solar_radiation 
clim_var = "solar_radiation_val"
# Create dataframe to use in manhattan plot
# Needs columns BP, CHR, P, SNP
manhattan_df <- data.frame(CHR = snp_pos$chrom, 
                           BP = snp_pos$pos,
                           P = pval_df_k2[, clim_var],
                           SNP = rownames(pval_df_k2))
which(manhattan_df$P < 1e-5)
soilval_associated_SNPs <- manhattan_df[which(manhattan_df$P < 1e-5), ]
write.csv(soilval_associated_SNPs, "solar_radiation_SNPs.csv")
# Make vector of chromosome names and format for input into manhattan plot function
manhattan_df$CHR <- factor(manhattan_df$CHR)
manhattan_df$CHR <- as.numeric(manhattan_df$CHR)
# Create manhattan plot
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = chr_names)
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Solar_radiation")
manhattan(manhattan_df, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Solar Radiation")

#For bio12 
clim_var = "bio12val"
# Create dataframe to use in manhattan plot
# Needs columns BP, CHR, P, SNP
manhattan_df <- data.frame(CHR = snp_pos$chrom, 
                           BP = snp_pos$pos,
                           P = pval_df_k2[, clim_var],
                           SNP = rownames(pval_df_k2))
which(manhattan_df$P < 1e-5)
soilval_associated_SNPs <- manhattan_df[which(manhattan_df$P < 1e-5), ]
str(soilval_associated_SNPs)
write.csv(soilval_associated_SNPs, "bio12_SNPs.csv")
# Make vector of chromosome names and format for input into manhattan plot function
manhattan_df$CHR <- factor(manhattan_df$CHR)
manhattan_df$CHR <- as.numeric(manhattan_df$CHR)
# Create manhattan plot
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = chr_names)
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio12 Annual Precipitation")
manhattan(manhattan_df, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio12 Annual Precipitation")

#For bio14 
clim_var = "bio14val"
# Create dataframe to use in manhattan plot
# Needs columns BP, CHR, P, SNP
manhattan_df <- data.frame(CHR = snp_pos$chrom, 
                           BP = snp_pos$pos,
                           P = pval_df_k2[, clim_var],
                           SNP = rownames(pval_df_k2))
which(manhattan_df$P < 1e-5)
soilval_associated_SNPs <- manhattan_df[which(manhattan_df$P < 1e-5), ]
str(soilval_associated_SNPs)
write.csv(soilval_associated_SNPs, "bio14_SNPs.csv")
# Make vector of chromosome names and format for input into manhattan plot function
manhattan_df$CHR <- factor(manhattan_df$CHR)
manhattan_df$CHR <- as.numeric(manhattan_df$CHR)
# Create manhattan plot
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = chr_names)
#manhattan(manhattan_df, annotatePval = 1e-5, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio14 Precipitation of Driest Month")
manhattan(manhattan_df, ylim = c(0, 8), chrlabs = c(manhattan_df$CHR), main = "Bio14 Precipitation of Driest Month")

# Choose number of top SNPs to investigate
n_top_snps = 7

# Extract genetic and climate data for top snps  
top_snps <-  pval_df_k2 %>% # Data frame of p values
  tibble::rownames_to_column(var = "snp") %>% # Make snp names a column
  dplyr::select(snp, clim_var) %>% # Select just climate variable of interest
  dplyr::top_n(-n_top_snps) # Selects the lowest valuesps

top_snps_df <- dplyr::bind_cols(clim_dat[, clim_var],
                                dplyr::select(gen_dat_scaled, top_snps$snp)) 
str(top_snps_df)

# Plot relationships between climate and frequencies of top snps
as.data.frame(top_snps_df) %>%
  tidyr::gather(key = "snp", value = "genotype", -clim_var) %>%
  ggplot(aes(x = .[, clim_var], y = genotype, group = factor(snp))) + 
  geom_jitter(height = .1, width = 0) + facet_wrap(~snp) + geom_smooth(method="lm") +
  xlab(clim_var) + ylab("genotype (scaled)") +
  theme_bw() 

#Read the SNPs associated with the environment 
SNPs <- read.csv("Environment_associated_SNPs.csv", header = TRUE)
str(SNPs)
SNPs1 <- SNPs[!duplicated(SNPs$SNP), ]
str(SNPs1)
library(data.table) ## v 1.9.6+ 
setDT(SNPs1)[, paste0("SNP", 1:2) := tstrsplit(SNP, "_")]
SNPs1
write.csv(SNPs1, "Remove_duplicate_Env_SNPs1.csv")
# Get map of California - already loaded in ggplot2 package
library(maps)
windows()
ca_map <- dplyr::filter(ggplot2::map_data("state"), region == "california")

# Plot genotypes at the top SNP on a map

ggplot(ca_map, aes(x = long, y = lat)) +
  geom_polygon(color = "black", fill = "grey80") +
  geom_point(data = clim_dat, aes(x = longitude, y = Latitude, 
                                  fill = factor(dplyr::pull(gen_dat, as.character(SNPs1$SNP[1])))), # Extract just top SNP
             color = "black", pch = 21, size = 2, position = position_jitter(w = 0.1, h = 0)) +
  scale_fill_discrete(name = "Allele copies") +
  theme_bw() + coord_equal(ratio = 1) 

list <- as.character(SNPs1$SNP)
library(dplyr)
gen_dat_selected <- subset(gen_dat, select= list)
dim(gen_dat_selected)
gen_dat_selected <- cbind(gen_dat$ID, gen_dat_selected)
str(gen_dat_selected)
