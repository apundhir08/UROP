suppressPackageStartupMessages(library(edgeR))
counts <- read.csv("sig/pseudobulk_counts.csv", row.names = 1, check.names = FALSE)
donors <- read.csv("sig/donor_table.csv", row.names = 1)
genes  <- read.csv("sig/genes.csv")
donors <- donors[donors$n_cells >= 20, ]
counts <- as.matrix(counts[, rownames(donors)])
group  <- factor(ifelse(donors$disease == "normal", "normal", "RA"), levels = c("normal", "RA"))
sex    <- factor(donors$sex)
y <- DGEList(counts = counts, genes = genes[match(rownames(counts), genes$feature_id), ])
y <- y[filterByExpr(y, group = group, min.count = 1), , keep.lib.sizes = FALSE]
y <- normLibSizes(y)

fit_coef <- function(g, coef) {
  design <- model.matrix(~ sex + g)
  fit <- glmQLFit(estimateDisp(y, design, robust = TRUE), design, robust = TRUE)
  topTags(glmQLFTest(fit, coef = coef), n = Inf, sort.by = "none")$table
}

sx <- fit_coef(group, 2)
cat("POSITIVE CONTROL (sex coefficient), top 5:\n")
print(head(sx[order(sx$PValue), c("feature_name", "logFC", "PValue")], 5), row.names = FALSE)

frac_sig <- function(g) mean(fit_coef(g, 3)$PValue < 0.05)
obs <- frac_sig(group)
set.seed(1)
perm <- replicate(40, {
  g <- group
  for (s in levels(sex)) { i <- which(sex == s); g[i] <- sample(group[i]) }
  frac_sig(g)
})
cat("PERMUTATION NULL: observed", round(obs, 3), "| permuted mean", round(mean(perm), 3),
    "range", round(min(perm), 3), "-", round(max(perm), 3), "| perms <= observed:", sum(perm <= obs), "of 40\n")
write.csv(data.frame(perm_frac_p05 = perm, observed = obs), "sig/permutation_null.csv", row.names = FALSE)

pc <- prcomp(t(cpm(y, log = TRUE, prior.count = 2)))
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
cat("PC variance %:", round(ve[1:4], 1), "\n")
cov <- data.frame(disease = as.numeric(group), sex = as.numeric(sex), log_cells = log10(donors$n_cells))
for (k in 1:3) cat("PC", k, "|r| with:", paste(names(cov), round(abs(sapply(cov, function(v) cor(pc$x[, k], v))), 2), collapse = "  "), "\n")
write.csv(data.frame(donor = rownames(donors), disease = donors$disease, sex = donors$sex,
                     n_cells = donors$n_cells, PC1 = pc$x[, 1], PC2 = pc$x[, 2], ve1 = ve[1], ve2 = ve[2]),
          "sig/pca_scores.csv", row.names = FALSE)