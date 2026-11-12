# Step 24: Differential expression using limma

library(limma)

# Load normalized expression and metadata
expression <- read.csv(
  "data/expression_normalized.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "data/metadata_normalized.csv",
  stringsAsFactors = FALSE
)

# Make sample IDs consistent
colnames(expression) <- as.character(colnames(expression))
metadata$Sample_ID <- as.character(metadata$Sample_ID)

# Match metadata to expression columns
metadata <- metadata[
  match(colnames(expression), metadata$Sample_ID),
]

stopifnot(
  !anyNA(metadata$Sample_ID),
  identical(metadata$Sample_ID, colnames(expression))
)

# Convert expression to numeric matrix
expression <- as.matrix(expression)
storage.mode(expression) <- "numeric"

stopifnot(
  nrow(expression) == 730,
  ncol(expression) == 32,
  all(is.finite(expression)),
  all(expression >= 0)
)

# Log2 transform normalized NanoString values
log_expression <- log2(expression + 1)

# Set study groups
metadata$Group <- factor(
  metadata$Group,
  levels = c(
    "Healthy",
    "Early-onset PE",
    "Late-onset PE"
  ),
  labels = c("Healthy", "Early", "Late")
)

stopifnot(!anyNA(metadata$Group))

# Design matrix
design <- model.matrix(~ 0 + Group, data = metadata)
colnames(design) <- c("Healthy", "Early", "Late")

# Fit gene-wise linear models
fit <- lmFit(log_expression, design)

# Define comparisons
contrasts <- makeContrasts(
  Early_vs_Late = Early - Late,
  Early_vs_Healthy = Early - Healthy,
  Late_vs_Healthy = Late - Healthy,
  levels = design
)

# Apply empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Save results for all comparisons
dir.create("results", showWarnings = FALSE)

for (comparison in colnames(contrasts)) {

  results <- topTable(
    fit2,
    coef = comparison,
    number = Inf,
    adjust.method = "BH"
  )

  write.csv(
    results,
    file = paste0("results/limma_", comparison, ".csv"),
    row.names = TRUE
  )

  cat("\n", comparison, "\n")
  cat("Genes tested:", nrow(results), "\n")
  cat("Significant genes (FDR < 0.05):",
      sum(results$adj.P.Val < 0.05), "\n")

  print(head(results, 10))
}

# ==========================================
# Gestational Age Sensitivity Analysis
# ==========================================

library(limma)

# Load clinical metadata
clinical <- read.csv(
  "data/clinical_metadata.csv",
  stringsAsFactors = FALSE
)

clinical$Sample_ID <- as.character(clinical$Sample_ID)

# Match clinical information to expression data
clinical <- clinical[
  match(colnames(expression), clinical$Sample_ID),
]

stopifnot(
  identical(clinical$Sample_ID, colnames(expression))
)

# Keep early-onset and late-onset samples
pe_samples <- metadata$Group %in% c("Early", "Late")

pe_expression <- log_expression[, pe_samples]

pe_metadata <- metadata[pe_samples, ]

pe_clinical <- clinical[pe_samples, ]

# Create numeric gestational age variable
pe_clinical$Gestational_Age <- as.numeric(
  pe_clinical$Gestational_Age
)

stopifnot(
  !anyNA(pe_clinical$Gestational_Age)
)

# Create statistical model
design_adjusted <- model.matrix(
  ~ Group + Gestational_Age,
  data = data.frame(
    Group = droplevels(pe_metadata$Group),
    Gestational_Age = pe_clinical$Gestational_Age
  )
)

# Fit adjusted model
fit_adjusted <- lmFit(
  pe_expression,
  design_adjusted
)

fit_adjusted <- eBayes(fit_adjusted)

# Extract early-versus-late group effect
# With Early as reference, GroupLate = Late - Early
adjusted_results <- topTable(
  fit_adjusted,
  coef = "GroupLate",
  number = Inf,
  adjust.method = "BH"
)

# Reverse coefficient for Early - Late interpretation
adjusted_results$logFC <- -adjusted_results$logFC
adjusted_results$t <- -adjusted_results$t

# Save results
write.csv(
  adjusted_results,
  "results/limma_gestational_age_sensitivity.csv"
)

cat("\nGestational-age sensitivity analysis\n")

cat(
  "Significant genes (FDR < 0.05):",
  sum(adjusted_results$adj.P.Val < 0.05),
  "\n"
)

print(head(adjusted_results, 10))

Rscript -e 'x <- read.csv("results/limma_gestational_age_sensitivity.csv"); cat("Genes tested:", nrow(x), "\nSignificant genes (FDR < 0.05):", sum(x$adj.P.Val < 0.05, na.rm=TRUE), "\n")'
Rscript -e 'x <- read.csv("results/limma_gestational_age_sensitivity.csv", row.names=1); print(x[c("ENG","MME","CEACAM1","BATF","IDO1","IL12RB2"), c("logFC","P.Value","adj.P.Val")])'