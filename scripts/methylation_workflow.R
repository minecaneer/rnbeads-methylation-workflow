##################################################
# RnBeads Methylation Workflow
#
# Structured R-based workflow for DNA methylation
# import, quality control, preprocessing,
# cell type estimation, inference, and
# exploratory analysis.
#
# Note:
# This script is a cleaned and documented version
# of a previously completed local analysis workflow.
# Depending on the available input files and package
# versions, some adaptation may be required before
# re-running the full pipeline.
##################################################

suppressPackageStartupMessages({
  library(RnBeads)
  library(RnBeads.hg19)
})

##################################################
# 1. Project setup
##################################################

# Example project paths
project.dir <- getwd()
idat.dir <- file.path(project.dir, "data")
report.dir <- file.path(project.dir, "results")
sample.annotation <- file.path(project.dir, "data", "SampleSheet_SelectedSamples.txt")

data.type <- "idat.dir"

##################################################
# 2. RnBeads options
##################################################

# General RnBeads settings
rnb.options(assembly = "hg19")
rnb.options(identifiers.column = "SampleID")
rnb.options(import.table.separator = "\t")

# Export current options to XML (optional)
rnb.options2xml(pretty = TRUE)

##################################################
# 3. Data import
##################################################

# Import methylation data and sample annotation
rnbSet <- rnb.run.import(
  data.source = c(idat.dir, sample.annotation),
  data.type = data.type,
  dir.reports = report.dir
)

rnb <- rnbSet$rnb.set

##################################################
# 4. Initial data inspection
##################################################

# Inspect phenotype information
head(pheno(rnb))

# Retrieve raw methylation values
mm <- meth(rnb)

# Plot beta value histogram for one example sample
# Change the column index to inspect another sample
hist(mm[, 5], col = "steelblue", breaks = 50)

# Inspect methylation matrix
head(mm)
dim(mm)

# Check summarized regions
summarized.regions(rnb)

##################################################
# 5. Probe and region annotation
##################################################

# Inspect annotation for 450k probes
anno <- rnb.annotation2data.frame(rnb.get.annotation("probes450"))
head(anno)

# Inspect promoter annotation
annot.promoters <- annotation(rnb, type = "promoters")
head(annot.promoters)

# Example promoter methylation values
meth(rnb, type = "promoters", row.names = TRUE, i = 1:5, j = 1:3)

##################################################
# 6. Coverage and detection p-values
##################################################

# Bead coverage per probe
nbead <- covg(rnb, row.names = TRUE)
nbead[1:5, 1:3]

# Detection p-values
pvals <- dpval(rnb, row.names = TRUE)

##################################################
# 7. Control probe diagnostics
##################################################

# General control probe boxplot
rnb.plot.control.boxplot(rnb)

# Example: bisulfite conversion control
rnb.plot.control.boxplot(rnb, "BISULFITE CONVERSION I")

# Negative control probe boxplot
rnb.plot.negative.boxplot(rnb)

# Barplot for one selected control probe
control.meta.data <- rnb.get.annotation("controls450")
ctrl.probe <- paste0(unique(control.meta.data[["Target"]])[2], ".1")
rnb.plot.control.barplot(rnb, ctrl.probe)

##################################################
# 8. QC data access and SNP-based similarity
##################################################

# Access QC information
qc_data <- qc(rnb)

# Use genotyping probes for sample similarity inspection
snp.probes <- anno[grep("rs", rownames(anno)), ]

# SNP heatmap
rnb.plot.snp.heatmap(rnb)

##################################################
# 9. Quality control settings
##################################################

rnb.options(qc.snp.boxplot = TRUE)
rnb.options(import.sex.prediction = TRUE)
rnb.options(qc.cnv = TRUE)

rnb.options(qc.boxplots = TRUE)
rnb.options(qc.barplots = TRUE)
rnb.options(qc.negative.boxplot = TRUE)
rnb.options(qc.snp.distances = TRUE)
rnb.options(qc.snp.boxplot = TRUE)
rnb.options(qc.snp.barplot = TRUE)
rnb.options(qc.sample.batch.size = 50)
rnb.options(qc.coverage.plots = FALSE)
rnb.options(qc.coverage.threshold.plot = 1:10)
rnb.options(qc.coverage.histograms = FALSE)
rnb.options(qc.coverage.violins = FALSE)

##################################################
# 10. Run complete quality control
##################################################

# Generate the full RnBeads QC report
rnb.run.qc(rnb.set = rnb, dir.reports = report.dir)

##################################################
# 11. Store unprocessed object
##################################################

# Keep a backup of the imported unprocessed dataset
rnb.unprocessed <- rnb

# Save unprocessed object
# Note: this path assumes write access to the results folder
save.rnb.set(
  rnb.unprocessed,
  file.path(report.dir, "rnb_unprocessed"),
  archive = TRUE
)

# Number of probes in the unfiltered object
nrow(meth(rnb.unprocessed))

##################################################
# 12. Filtering
##################################################

# Remove probes outside CpG context
rnb.set.filtered <- rnb.execute.context.removal(
  rnb.set = rnb.unprocessed,
  contexts = "Other"
)$dataset

nrow(meth(rnb.set.filtered))

# Remove probes overlapping SNPs within 3 bp
rnb.set.filtered <- rnb.execute.snp.removal(
  rnb.set = rnb.set.filtered,
  snp = "3"
)$dataset

nrow(meth(rnb.set.filtered))

# Remove probes on sex chromosomes
rnb.set.filtered <- rnb.execute.sex.removal(
  rnb.set = rnb.set.filtered
)$dataset

nrow(meth(rnb.set.filtered))

# Greedycut-based probe filtering
greedycut.results <- rnb.execute.greedycut(
  rnb.set = rnb.set.filtered,
  pval.threshold = 0.05
)

to_remove <- rownames(meth(object = rnb.set.filtered, row.names = TRUE))[
  greedycut.results[["sites"]]
]

rnb.set.filtered <- remove.sites(
  object = rnb.set.filtered,
  probelist = to_remove
)

nrow(meth(rnb.set.filtered))

# Remove probes containing NA values
rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
nrow(meth(rnb.set.filtered))

# Remove low-variability probes
rnb.set.filtered <- rnb.execute.variability.removal(
  rnb.set.filtered,
  0.005
)$dataset

nrow(meth(rnb.set.filtered))

# Remove probes with problematic values
mm.filtered <- meth(rnb.set.filtered, row.names = TRUE)
invalid.probes <- apply(mm.filtered, 1, function(x) any(x <= 0 | is.na(x)))
invalid.probe.names <- rownames(mm.filtered)[invalid.probes]

if (length(invalid.probe.names) > 0) {
  rnb.set.filtered <- remove.sites(
    object = rnb.set.filtered,
    probelist = invalid.probe.names
  )
}

nrow(meth(rnb.set.filtered))

# Save filtered object
save.rnb.set(
  rnb.set.filtered,
  file.path(report.dir, "rnb_filtered"),
  archive = TRUE
)

##################################################
# 13. Preprocessing settings
##################################################

rnb.options(filtering.whitelist = NULL)
rnb.options(filtering.blacklist = NULL)
rnb.options(filtering.snp = "3")
rnb.options(filtering.cross.reactive = FALSE)
rnb.options(filtering.greedycut = TRUE)
rnb.options(filtering.greedycut.pvalue.threshold = 0.05)
rnb.options(filtering.greedycut.rc.ties = "row")
rnb.options(filtering.sex.chromosomes.removal = TRUE)
rnb.options(filtering.missing.value.quantile = 0.8)
rnb.options(filtering.coverage.threshold = 3)
rnb.options(filtering.low.coverage.masking = FALSE)
rnb.options(filtering.high.coverage.outliers = FALSE)
rnb.options(filtering.deviation.threshold = 0)
rnb.options(filtering.context.removal = "Other")

# Normalization settings
rnb.options(normalization = NULL)
rnb.options(normalization.method = "wm.dasen")
rnb.options(normalization.background.method = "methylumi.noob")
rnb.options(normalization.plot.shifts = TRUE)

##################################################
# 14. Run preprocessing
##################################################

preprocessed <- rnb.run.preprocessing(
  rnb.set = rnb.set.filtered,
  dir.reports = report.dir
)

rnb <- preprocessed$rnb.set

##################################################
# 15. Cell type estimation
##################################################

# Estimate cell type composition using the Houseman method
ct <- rnb.execute.ct.estimation(
  rnb,
  cell.type.column = "CellType",
  test.max.markers = 10000,
  top.markers = 500
)

# Plot estimated cell type heatmap
rnb.plot.ct.heatmap(ct.obj = ct)

# Estimated cell type contributions
ct$contributions

ctc <- as.data.frame(ct$contributions)
head(ctc)

# Basic heatmap of estimated proportions
heatmap(as.matrix(ctc), scale = "none")

##################################################
# 16. Add estimated neural content
##################################################

# Add estimated neuron content to phenotype data
phenoInfo <- pheno(rnb)
phenoInfo$NeuralContent <- ctc[
  match(rownames(phenoInfo), rownames(ctc)),
  "Neuron"
]

# Note:
# Depending on the RnBeads object structure, phenotype updates
# may require adaptation in different versions of the package.
rnb <- addPheno(object = rnb, phenoInfo$NeuralContent)

##################################################
# 17. Immune cell content and epigenetic age
##################################################

# Estimate immune cell content
immune.content <- rnb.execute.lump(rnb)

# Estimate epigenetic age
rnb.execute.age.prediction(object = rnb)

##################################################
# 18. Inference settings
##################################################

rnb.options(inference.age.prediction = TRUE)
rnb.options(inference.age.column = "Age")
rnb.options(inference.age.prediction.training = FALSE)
rnb.options(inference.age.prediction.cv = FALSE)
rnb.options(inference.immune.cells = TRUE)
rnb.options(inference.genome.methylation = "Genome-wide methylation")
rnb.options(inference.targets.sva = character())
rnb.options(inference.reference.methylome.column = "CellType")
rnb.options(inference.max.cell.type.markers = 10000)
rnb.options(inference.top.cell.type.markers = 500)
rnb.options(inference.sva.num.method = "leek")

##################################################
# 19. Run inference analysis
##################################################

rnb_inference <- rnb.run.inference(
  rnb.set = rnb,
  dir.reports = report.dir
)

rnb <- rnb_inference$rnb.set

##################################################
# 20. Surrogate variable analysis (SVA)
##################################################

sva.obj <- rnb.execute.sva(
  rnb,
  cmp.cols = "Group",
  numSVmethod = "be"
)

rnb <- set.covariates.sva(rnb, sva.obj)

##################################################
# 21. Exploratory analysis
##################################################

rnb.run.exploratory(rnb.set = rnb, dir.reports = report.dir)

##################################################
# 22. Save final object
##################################################

Sys.setenv("R_ZIPCMD" = "/usr/bin/zip")

save.rnb.set(
  rnb,
  file.path(report.dir, "rnb_exploratory"),
  archive = TRUE
)