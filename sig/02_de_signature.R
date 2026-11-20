suppressPackageStartupMessages(library(edgeR))

MIN_CELLS <- 20
N_LIST    <- 150

run_de <- function(min_cells) {
  counts <- read.csv("sig/pseudobulk_counts.csv", row.names = 1, check.names = FALSE)
  donors <- read.csv("sig/donor_table.csv", row.names = 1)
  genes  <- read.csv("sig/genes.csv")

  donors <- donors[donors$n_cells >= min_cells, ]
  counts <- as.matrix(counts[, rownames(donors)])
  group  <- factor(ifelse(donors$disease == "normal", "normal", "RA"), levels = c("normal", "RA"))
  sex    <- factor(donors$sex)

  y <- DGEList(counts = counts, genes = genes[match(rownames(counts), genes$feature_id), ])
  y <- y[filterByExpr(y, group = group, min.count = 1), , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)

  design <- model.matrix(~ sex + group)
  y   <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- topTags(glmQLFTest(fit, coef = "groupRA"), n = Inf, sort.by = "none")$table

  res$t_signed <- sign(res$logFC) * sqrt(pmax(res$F, 0))
  res      <- res[order(-res$t_signed), ]
  res$rank <- seq_len(nrow(res))
  attr(res, "n_donors") <- table(group)
  res
}

sig <- run_de(MIN_CELLS)
n_d <- attr(sig, "n_donors")
cat("donors used: normal =", n_d[["normal"]], " RA =", n_d[["RA"]], "| genes tested:", nrow(sig), "\n")
cat("genes with FDR < 0.05:", sum(sig$FDR < 0.05), "| smallest FDR:", signif(min(sig$FDR), 3), "\n")

write.csv(sig, "sig/ra_vs_normal_classical_monocyte_ranked.csv", row.names = FALSE)
up   <- head(sig, N_LIST)
down <- tail(sig, N_LIST)[N_LIST:1, ]
write.csv(rbind(cbind(list = "up", up), cbind(list = "down", down)),
          "sig/query_up_down_lists.csv", row.names = FALSE)

sens <- do.call(rbind, lapply(c(1, 20, 50, 150), function(m) {
  r <- run_de(m); nd <- attr(r, "n_donors"); shared <- intersect(sig$feature_id, r$feature_id)
  data.frame(min_cells = m, donors_normal = nd[["normal"]], donors_RA = nd[["RA"]],
             genes_tested = nrow(r), min_FDR = signif(min(r$FDR), 3),
             spearman_t_vs_main = round(cor(sig$t_signed[match(shared, sig$feature_id)],
                                            r$t_signed[match(shared, r$feature_id)], method = "spearman"), 3),
             up_overlap   = length(intersect(head(r$feature_id, N_LIST), up$feature_id)),
             down_overlap = length(intersect(tail(r$feature_id, N_LIST), down$feature_id)))
}))
print(sens, row.names = FALSE)
write.csv(sens, "sig/min_cells_sensitivity.csv", row.names = FALSE)