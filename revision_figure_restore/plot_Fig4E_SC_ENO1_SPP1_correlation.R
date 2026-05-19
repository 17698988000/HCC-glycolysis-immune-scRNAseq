# Figure 4E: ENO1-SPP1 correlation in tumor-derived hepatocytes
# Output:
#   Fig4E_SC_ENO1_SPP1_correlation_vector_clean2.pdf
#   Fig4E_SC_ENO1_SPP1_correlation_check.txt

source("00_config_paths.R")
ensure_packages(c("Seurat", "ggplot2"))
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})

candidate_rds <- c(
  "seurat_final.rds", "seurat_annotated.rds", "seurat_harmony.rds",
  "seurat_filtered.rds", "hcc_seurat.rds", "hcc_seurat_copykat.rds"
)
rds_file <- candidate_rds[file.exists(candidate_rds)][1]
if (is.na(rds_file)) stop("No candidate Seurat RDS file found.")
seu <- readRDS(rds_file)

meta_cols <- colnames(seu@meta.data)
site_col <- intersect(c("site", "Site", "tissue", "Tissue", "sample_type", "orig.ident"), meta_cols)[1]
celltype_col <- intersect(c("celltype", "CellType", "cell_type", "celltype_major", "major_celltype", "annotation", "cell_annotation"), meta_cols)[1]
if (is.na(site_col) || is.na(celltype_col)) stop("Could not detect site or celltype metadata columns.")

site_value <- as.character(seu@meta.data[[site_col]])
celltype_value <- as.character(seu@meta.data[[celltype_col]])
tumor_flag <- grepl("Tumor|PT|Primary", site_value, ignore.case = TRUE)
hep_flag <- grepl("Hepatocyte|Hepatocytes|Malignant|Tumor", celltype_value, ignore.case = TRUE)
cells_use <- rownames(seu@meta.data)[tumor_flag & hep_flag]
if (length(cells_use) < 1000) stop("Too few tumor hepatocyte cells selected: ", length(cells_use))

tumor_hep <- subset(seu, cells = cells_use)
missing_genes <- setdiff(c("ENO1", "SPP1"), rownames(tumor_hep))
if (length(missing_genes) > 0) stop("Missing genes: ", paste(missing_genes, collapse = ", "))
expr_data <- tryCatch(GetAssayData(tumor_hep, assay = "RNA", layer = "data"),
                      error = function(e) GetAssayData(tumor_hep, assay = "RNA", slot = "data"))
plot_df <- data.frame(ENO1 = as.numeric(expr_data["ENO1", ]),
                      SPP1 = as.numeric(expr_data["SPP1", ]))
plot_df <- plot_df[complete.cases(plot_df), ]
pearson_test <- cor.test(plot_df$ENO1, plot_df$SPP1, method = "pearson")
spearman_test <- cor.test(plot_df$ENO1, plot_df$SPP1, method = "spearman")

check_file <- "Fig4E_SC_ENO1_SPP1_correlation_check.txt"
sink(check_file)
cat("===== Fig4E ENO1-SPP1 correlation check =====\n\n")
cat("Seurat object:", rds_file, "\n")
cat("site column:", site_col, "\n")
cat("celltype column:", celltype_col, "\n")
cat("Candidate tumor hepatocyte cells:", length(cells_use), "\n")
cat("Cells used for correlation:", nrow(plot_df), "\n")
cat("Pearson R:", unname(pearson_test$estimate), "\n")
cat("Pearson p:", pearson_test$p.value, "\n")
cat("Spearman rho:", unname(spearman_test$estimate), "\n")
cat("Spearman p:", spearman_test$p.value, "\n")
sink()

label_text <- paste0("Pearson R = ", round(unname(pearson_test$estimate), 3),
                     "\nSpearman rho = ", round(unname(spearman_test$estimate), 3))
p <- ggplot(plot_df, aes(x = ENO1, y = SPP1)) +
  geom_point(size = 0.35, alpha = 0.35) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  annotate("label", x = Inf, y = Inf, label = label_text, hjust = 1.05,
           vjust = 1.15, size = 4, fill = "white", label.size = 0.25) +
  labs(x = "ENO1 expression", y = "SPP1 expression") +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 12),
        plot.margin = margin(12, 18, 10, 10)) +
  coord_cartesian(clip = "off")

ggsave("Fig4E_SC_ENO1_SPP1_correlation_vector_clean2.pdf", p, device = cairo_pdf,
       width = 5.2, height = 4.6, units = "in")
cat("Saved figure and check file.\n")
