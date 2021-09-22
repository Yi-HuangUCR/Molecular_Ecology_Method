# NOTE: Set the maximum percent missing data you will allow. Give this in percentage, not proportion.
max_percent_missing <- 20

# Read in simplified VCF file.
RAD <- as.matrix( read.table("output_from_step_one.txt", 
	header = F, stringsAsFactor = F) )

# Clip off the first 10 columns of the matrix. These are the position and quality info, etc.
RAD2 <- RAD[ , 10:dim(RAD)[2] ]

# Save REF and ALT columns as vectors.
REF <- RAD[, 4]
ALT <- RAD[, 5]

# Pull the first 7 characters out of the genotype strings (just the genotypes).

GT_mat <- RAD2[-1, ]
GT_nrow <- nrow(GT_mat)
GT_ncol <- ncol(GT_mat)
GT_mat_as_vector <- as.vector(GT_mat)
GT_mat_as_vector_missing <- which(GT_mat_as_vector == ".")

library(stringr)
GT_mat_as_vector[GT_mat_as_vector_missing] <- ".:.:.:."
GT_mat_as_vector_ONLY_GT <- str_sub(GT_mat_as_vector, 1, 7)
GT_mat_only_GT <- matrix(nrow = GT_nrow, ncol = GT_ncol, GT_mat_as_vector_ONLY_GT)

RAD3 <- rbind( RAD2[1, ], GT_mat_only_GT )

# Make names vector for value matching.
sample_names <- RAD3[1, ]
# Note for Yi: make sure the sample names in your VCF and the sample names ..
# .. in ref_table.csv are the same. I may have altered the names so that ..
# .. dashes (-) were replaced with underscores (_). If the sample names are ..
# .. not the same, these table references won't work.


# Read in table of population/taxa data.
ref_table <- read.csv("ref_table.csv")

# Eliminate columns that are duplicates.
duplicates <- which( ref_table$Is_duplicate[ match(sample_names, as.character(ref_table$Name_in_SNP_spread)) ] )
RAD4 <- RAD3[, -duplicates]
# Note for Yi: This step is to remove samples that are duplicates in the ddRAD ..
# .. data. These were the sample we got from Dylan that had different tube codes ..
# .. but came from the same plant. You can keep the duplicates by commenting out ..
# .. the line above and uncommenting the next line.
#RAD4 <- RAD3

# Recreate sample names vector.
sample_names <- RAD4[1, ]

# Create species vector.
Species <- ref_table$Species[ match(sample_names, as.character(ref_table$Name_in_SNP_spread)) ]

# Eliminate columns that are not glandulosa.
RAD5 <- RAD4[ , Species == "glandulosa" ]

# Recreate samples names vector.
sample_names <- RAD5[1, ]

# Drop the sample names from the matrix.
RAD6 <- RAD5[-1, ]

# Count the cells with missing data for each locus, and record the ..
# .. number of unique values to later test whether the locus is invariant.
# Also calculate the frequency with which the major allele is the only allele.
n_missing <- NA
n_unique <- NA
n_with_minor_allele_present <- NA
prop_with_minor_allele_present_non_missing <- NA
for(i in 1:dim(RAD6)[1])
{
	row_i <- RAD6[i, ]
	n_unique[i] <- length( unique(row_i) )
	n_missing[i] <- ifelse(".:.:.:." %in% row_i, length( which( row_i == ".:.:.:." ) ), 0)
	n_not_missing <- length(row_i) - n_missing[i]
	n_with_minor_allele_present[i] <- length( which( (row_i != "0/0/0/0") & (row_i != ".:.:.:.") ) )
	if( is.na(n_with_minor_allele_present[i]) )
	{
		n_with_minor_allele_present[i] <- 0 
	}
	prop_with_minor_allele_present_non_missing[i] <- n_with_minor_allele_present[i] / n_not_missing
}


# Calculate the proportion of samples missing data for each locus.
prop_missing <- n_missing / dim(RAD6)[2]

# Test whether the locus is invariant.
is_invariant <- (n_unique == 1) | ((n_unique == 2) & (n_missing > 0))

# Filter data down to just SNPs that have less than the specified percent missing data, and are not invariant.
prop_max <- max_percent_missing / 100
df <- RAD6[ which( (prop_missing <= prop_max) & (!is_invariant) ),  ]
REF <- REF[ (prop_missing <= prop_max) & (!is_invariant) ]
ALT <- ALT[ (prop_missing <= prop_max) & (!is_invariant) ]
df <- cbind(REF, ALT, df)
colnames(df) <- c("REF", "ALT", sample_names)

# Save the filtered data as a table.
write.table(file = "output_from_step_two.txt",
	x = df, row.names = F, col.names = T, quote = F)

