
library("tidyverse")
library("edgeR")
library("limma")
library("psych")
library("gridExtra")
library("stringr")
library("scales")

# read in the data
temp <- read_tsv("ave_missing.txt")

# make a basic plot
ggplot(temp, aes(x = Ave, y = FracMissing)) +
  geom_point() +
  ggtitle("Missing Fraction versus Average SpC") + 
  labs(x = "Ave SpC", y = "Missing Fraction")

# expanded x-axis plot
ggplot(temp, aes(x = Ave, y = FracMissing)) +
  geom_line() + 
  coord_cartesian(xlim = c(0, 8)) +
  ggtitle("Missing Fraction versus Average SpC") + 
  labs( x = "Ave SpC", y = "Missing Fraction") + 
  geom_vline(xintercept = 2.5, linetype = "dotted")

# read in the prepped data
paw_spc <- read_tsv("edgeR_input.txt")

# save accessions in vector and remove from data table
accession <- paw_spc$Accession
paw_spc <- select(paw_spc, -Accession)
head(paw_spc)
nrow(paw_spc)

freq <- read_tsv("Fractions_SpC_ID.txt")
ggplot(freq, aes(InHowMany, Fraction, fill = Measure)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_text(aes(label = Fraction), vjust = 1.6, color = "black",
            position = position_dodge(0.9), size = 2.5) +
  scale_x_continuous( breaks = 1:10) +
  ggtitle("Fractions of total spectra counts or of total identifications") +
  xlab("In How Many Samples") + ylab("Fraction (%)")

# function for simple normalization
SL_Norm <- function(df, color = NULL, plot = TRUE) {
    # This makes each channel sum to the average grand total
        # df - data frame of TMT intensities
        # returns a new data frame with normalized values
    
    # compute scaling factors to make colsums match the average sum
    norm_facs <- mean(c(colSums(df))) / colSums(df)
    cat("SL Factors:\n", sprintf("%-5s -> %f\n", colnames(df), norm_facs))
    df_sl  <- sweep(df, 2, norm_facs, FUN = "*")
    
    # visualize results and return data frame
    if(plot == TRUE) {
        boxplot(log10(df_sl + 1), col = color, notch = TRUE, main = "SL Normalized data")
    }
    df_sl
}

# normalize the data before plotting
color <- c(rep('blue', 5), rep('red', 5))
paw_sl <- SL_Norm(paw_spc, color)

# shortcuts for the cell types
C <- 1:5
R <- 6:10

pairs.panels(paw_sl[C], main = "Choroidal")
pairs.panels(log2(paw_sl[C]+1), main = "Choroidal")
pairs.panels(paw_sl[R], main = "Retinal")
pairs.panels(log2(paw_sl[R]+1), main = "Retinal")

# load the data into edgeR data structures
# group labels need to be factors
y <- DGEList(counts = paw_spc, genes = accession)

# run the TMM normalization (and library size corrections)
y <- calcNormFactors(y)

apply_tmm_factors <- function(y, color = NULL, plot = TRUE) {
    # computes the tmm normalized data from the DGEList object
        # y - DGEList object
        # returns a dataframe with normalized intensities
    
    # compute grand total (library size) scalings
    lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size

    # the TMM factors are library adjustment factors (so divide by them)
    norm_facs <- lib_facs / y$samples$norm.factors
    cat("Overall Factors (lib.size+TMM):\n", sprintf("%-5s -> %f\n", 
                                                     colnames(y$counts), norm_facs))

    # compute the normalized data as a new data frame
    df_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN = "*"))
    colnames(df_tmm) <- str_c(colnames(y$counts), "_tmm")
    
    # visualize results and return data frame
    if(plot == TRUE) {
        boxplot(log10(df_tmm + 1), col = color, notch = TRUE, main = "TMM Normalized data")
    }
    df_tmm
}

paw_spc_tmm <- apply_tmm_factors(y, color)

# check clustering (6 different out of 1748 may not do much)
plotMDS(y)

CV <- function(df) {
    # Computes CVs of data frame rows
        # df - data frame, 
        # returns vector of CVs (%)
    ave <- rowMeans(df)    # compute averages
    sd <- apply(df, 1, sd) # compute standard deviations
    cv <- 100 * sd / ave   # compute CVs in percent (last thing gets returned)
}

labeled_boxplot <- function(df, ylim, title) {
    # Makes a box plot with the median value labeled
        # df - data frame with data to compute CVs of
        # ylim - upper limit for y-axis
        # title - plot title
    cv = CV(df)
    boxplot(cv, ylim = c(0, ylim), notch = TRUE, main = title)
    text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
         labels = round(boxplot.stats(cv)$stats[3], 1))
}

# see what effect TMM had on CV distributions
par(mfrow = c(2, 2))
labeled_boxplot(paw_spc[C], 150, "Choroid before")
labeled_boxplot(paw_spc[C], 150, "Retina before")
labeled_boxplot(paw_spc_tmm[C], 150, "Choroid after")
labeled_boxplot(paw_spc_tmm[C], 150, "Retina after")
par(mfrow = c(1, 1))

group <- factor(c(rep("C", 5), rep("R", 5)))
yy <-  DGEList(counts = paw_spc, group = group, genes = accession)
yy <- calcNormFactors(yy)
yy <- estimateDisp(yy)
et <- exactTest(yy)
topTags(et)$table
summary(decideTests(et, p.value = 0.10))

# create the experimental design matrix
donor <- factor(rep(c(189, 191, 194, 195, 199), 2))
cell <- factor(c(rep("C", 5), rep("R", 5)))

# Example 4.1 in edgeR user's guide
design <- model.matrix(~donor+cell)
rownames(design) <- colnames(y)
design

# extimate the dispersion parameters and check
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
plotBCV(y, main = "Variance Trends")

# fit statistical models (design matrix already in y$design)
fit <- glmQLFit(y, design, robust = TRUE)
plotQLDisp(fit)

# if we do not specify a contrast, the default is the last column
# of the design matrix - a C versus R comparison
paired <- glmQLFTest(fit) # default comparison

# check test results
topTags(paired)
tt <- topTags(paired, n = Inf, sort.by = "none")$table
summary(decideTests(paired, p.value = 0.10))

collect_results <- function(df, tt, x, xlab, y, ylab) {
    # Computes new columns and extracts some columns to make results frame
        # df - data in data.frame
        # tt - top tags from edgeR test
        # x - columns for first condition
        # xlab - label for x
        # y - columns for second condition
        # ylab - label for y
        # returns a new dataframe
    
    # condition average vectors
    ave_x <- rowMeans(df[x])
    ave_y <- rowMeans(df[y])
    
    # FC, direction, candidates
    fc <- ifelse(ave_y > ave_x, (ave_y / ave_x), (-1 * ave_x / ave_y))
    direction <- ifelse(ave_y > ave_x, "up", "down")
    candidate <- cut(tt$FDR, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0), 
                     labels = c("high", "med", "low", "no"))
    
    # make data frame
    temp <- cbind(df[c(x, y)], data.frame(logFC = tt$logFC, FC = fc, 
                                          PValue = tt$PValue, FDR = tt$FDR, 
                                          ave_x = ave_x, ave_y = ave_y, 
                                          direction = direction, candidate = candidate, 
                                          Acc = tt$genes)) 
    
    # fix column headers for averages
    names(temp)[names(temp) %in% c("ave_x", "ave_y")]  <- str_c("ave_", c(xlab, ylab))    
    
    temp # return the data frame
}

# get the results summary
results <- collect_results(paw_spc_tmm, tt, C, "choroid", R, "retina")

pvalue_plots <- function(results, ylim, title) {
    # Makes p-value distribution plots
        # results - results data frame
        # ylim - ymax for expanded view
        # title - plot title
    p_plot <- ggplot(results, aes(PValue)) + 
        geom_histogram(bins = 100, fill = "white", color = "black") +
        geom_hline(yintercept = mean(hist(results$PValue, breaks = 100, 
                                     plot = FALSE)$counts[26:100]))

    # we will need an expanded plot
    p1 <- p_plot + ggtitle(str_c(title, " p-value distribution"))
    p2 <- p_plot + coord_cartesian(xlim = c(0, 1.0), ylim = c(0, ylim)) + ggtitle("p-values expanded")
    grid.arrange(p1, p2, nrow = 2) # from gridExtra package
}

# check the p-value distribution
pvalue_plots(results, 100, "Retina vs Choroid - SpC")

log2FC_plots <- function(results, range, title) {
    # Makes faceted log2FC plots by candidate
        # results - results data frame
        # range - plus/minus log2 x-axis limits
        # title - plot title
    ggplot(results, aes(x = logFC, fill = candidate)) +
        geom_histogram(binwidth=0.1, color = "black") +
        facet_wrap(~candidate) +
        ggtitle(title) + 
        coord_cartesian(xlim = c(-range, range))
}

# make log2FC plots
log2FC_plots(results, 4, "Faceted log2FC")

# see how many candidates are in each category
results %>% count(candidate)

transform <- function(results, x, y) {
    # Make data frame with some transformed columns
        # results - results data frame
        # x - columns for x condition
        # y - columns for y condition
        # return new data frame
    df <- data.frame(log10((results[x] + results[y])/2), 
                     log2(results[y] / results[x]), 
                     results$candidate,
                     -log10(results$FDR))
    colnames(df) <- c("A", "M", "candidate", "P")
    
    df # return the data frame
}

MA_plots <- function(results, x, y, title) {
    # makes MA-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # returns a list of plots 
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # 2-fold change lines
    ma_lines <- list(geom_hline(yintercept = 0.0, color = "black"),
                     geom_hline(yintercept = 1.0, color = "black", linetype = "dotted"),
                     geom_hline(yintercept = -1.0, color = "black", linetype = "dotted"))

    # make main MA plot
    ma <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("logFC (", y, "/", x, ")")) +
        scale_x_continuous("Ave_intensity") +
        ggtitle(title) + 
        ma_lines
    
    # make separate MA plots
    ma_facet <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("log2 FC (", y, "/", x, ")")) +
        scale_x_continuous("log10 Ave_intensity") +
        ma_lines +
        facet_wrap(~ candidate) +
        ggtitle(str_c(title, " (separated)"))

    # make the plots visible
    print(ma)
    print(ma_facet)
}    

# make MA plots
MA_plots(results, "ave_choroid", "ave_retina", "Choroid vs Retina")

scatter_plots <- function(results, x, y, title) {
    # makes scatter-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # returns a list of plots
    
    # 2-fold change lines
    scatter_lines <- list(geom_abline(intercept = 0.0, slope = 1.0, color = "black"),
                          geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          scale_y_log10(),
                          scale_x_log10())

    # make main scatter plot
    scatter <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        ggtitle(title) + 
        scatter_lines

    # make separate scatter plots
    scatter_facet <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scatter_lines +
        facet_wrap(~ candidate) +
        ggtitle(str_c(title, " (separated)")) 

    # make the plots visible
    print(scatter)
    print(scatter_facet)
}

# make scatter plots
scatter_plots(results, "ave_choroid", "ave_retina", "Choroid vs Retina")

volcano_plot <- function(results, x, y, title) {
    # makes a volcano plot
        # results - a data frame with edgeR results
        # x - string for the x-axis column
        # y - string for y-axis column
        # title - plot title string
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # build the plot
    ggplot(temp, aes(x = M, y = P)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        xlab("log2 FC") +
        ylab("-log10 FDR") +
        ggtitle(str_c(title, " Volcano Plot"))
}

# make a volcano plot
volcano_plot(results, "ave_choroid", "ave_retina", "Choroid vs Retina")

# write results
write.table(results, "edgeR_results.txt", sep = "\t", row.names = FALSE, na = " ")

# log the session
sessionInfo()


