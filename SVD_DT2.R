# SVD (Singular Value Decomposition)
# Data Ekspresi Gen - Diabetes Tipe 2 (GSE25724)

library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(writexl)
library(gridExtra)
library(scales)

# LANGKAH 1: Load Data
DT2 <- "D:/Statistika Part 2/MODERN PREDIKSI DAN MACHINE LEARNING/Data Diabetes.xlsx"
df_raw <- read_excel(DT2, col_names = TRUE)
print(head(df_raw, 5))

# LANGKAH 2: Cleaning Data
# Hapus baris metadata yang diawali "!"
col1    <- colnames(df_raw)[1]
df_clean <- df_raw %>%
  filter(!grepl("^!", .data[[col1]]))

# Jadikan baris pertama sebagai header nama sampel
colnames(df_clean) <- as.character(df_clean[1, ])
df_clean <- df_clean[-1, ]

# Jadikan kolom ID_REF sebagai rownames (nama gen)
df_clean <- as.data.frame(df_clean)

rownames(df_clean) <- df_clean[["ID_REF"]]
df_clean[["ID_REF"]] <- NULL

# Simpan nama gen
gene_names <- rownames(df_clean)

# Konversi ke numerik
df_clean <- as.data.frame(apply(df_clean, 2, as.numeric))

# Kembalikan nama gen
rownames(df_clean) <- gene_names

# Hapus baris yang semua nilainya NA
df_clean <- df_clean[rowSums(is.na(df_clean)) != ncol(df_clean), ]

cat(sprintf("\n→ Data bersih: %d gen × %d sampel\n",
            nrow(df_clean), ncol(df_clean)))

# LANGKAH 3: Transpose
df_T <- t(df_clean)  # baris = sampel, kolom = gen
cat(sprintf("   → Setelah transpose: %d sampel × %d gen\n",
            nrow(df_T), ncol(df_T)))

# LANGKAH 4: Seleksi 500 Gen Paling Variatif
var_gen  <- apply(df_T, 2, var, na.rm = TRUE)
top500   <- names(sort(var_gen, decreasing = TRUE))[1:500]
df_top   <- df_T[, top500]

cat(sprintf("   → Data siap SVD: %d sampel × %d gen\n",
            nrow(df_top), ncol(df_top)))
View(df_top)
head(df_top)

# LANGKAH 5: Label Kelompok
n_sampel  <- nrow(df_top)
labels    <- c(rep(0, 7), rep(1, 6))          # 0=Non-Diabetik, 1=T2D
kelompok  <- ifelse(labels == 0, "Non-Diabetik", "Diabetes T2D")
nama_samp <- rownames(df_top)

data.frame(
  Sampel = nama_samp,
  Label = labels,
  Kelompok = kelompok
)

cat(sprintf("   → Label: Non-Diabetik = %d | Diabetes T2D = %d\n",
            sum(labels == 0), sum(labels == 1)))

# LANGKAH 6: Standarisasi
X_scaled <- scale(df_top)   # zero mean, unit variance
View(X_scaled)
head(X_scaled)
cat("   → Standarisasi selesai (zero mean, unit variance)\n")

# LANGKAH 7: SVD — Singular Value Decomposition
cat("\n→ Menjalankan SVD...\n")

svd_result <- svd(X_scaled)
str(svd_result)

# Komponen SVD:
#   U     : Left singular vectors  (sampel × sampel)  → representasi sampel
#   D (Σ) : Singular values        (vektor diagonal)  → "kepentingan" tiap dimensi
#   V     : Right singular vectors (gen × sampel)     → representasi gen (loading)

U  <- svd_result$u     # 13 × 13 (Posisi/pola setiap sampel)
View(U)
D  <- svd_result$d     # vektor singular values (Besarnya informasi tiap komponen)
D
V  <- svd_result$v     # 500 × 13 (Kontribusi setiap gen)
View(V)

# Ambil 10 komponen pertama
n_comp <- 10
U10    <- U[, 1:n_comp]
U10
D10    <- D[1:n_comp]
D10
V10    <- V[, 1:n_comp]
V10

# Skor SVD = U × Σ  (setara skor PCA)
skor_svd <- U10 %*% diag(D10)
uji_sv1 <- t.test(skor_svd[,1] ~ kelompok)
uji_sv2 <- t.test(skor_svd[,2] ~ kelompok)

print(uji_sv1)
print(uji_sv2)

# Variance explained dari singular values
var_exp     <- (D^2 / sum(D^2)) * 100
var_exp10   <- var_exp[1:n_comp]
var_cum     <- cumsum(var_exp10)

cat("\n→ Variance yang dijelaskan tiap Komponen SVD:\n")
for (i in 1:n_comp) {
  bar <- paste(rep("█", floor(var_exp10[i] / 2)), collapse = "")
  cat(sprintf("   SV%2d: %5.2f%%  %s\n", i, var_exp10[i], bar))
}
cat(sprintf("\n   → Total SV1 + SV2: %.2f%%\n", sum(var_exp10[1:2])))

# LANGKAH 8: Loading Vector & Kandidat Biomarker
# V = right singular vectors = setara eigenvector PCA
loading_df <- as.data.frame(V10)
colnames(loading_df) <- paste0("SV", 1:n_comp)
rownames(loading_df) <- top500

# Ranking gen berdasarkan kontribusi absolut di SV1
loading_df$Abs_SV1 <- abs(loading_df$SV1)
loading_df$Abs_SV2 <- abs(loading_df$SV2)
loading_sorted     <- loading_df[order(loading_df$Abs_SV1, decreasing = TRUE), ]

# Top 20 kandidat biomarker
top20_biomarker <- head(loading_sorted, 20)
top20_biomarker$Gen <- rownames(top20_biomarker)

cat("\n", strrep("=", 55), "\n")
cat("  TOP 20 KANDIDAT BIOMARKER GEN (berdasarkan SV1)\n")
cat(strrep("=", 55), "\n")
print(top20_biomarker[, c("Gen", "SV1", "SV2", "Abs_SV1")],
      row.names = FALSE)

# LANGKAH 9: Visualisasi (4 Panel)
warna_kel <- c("Non-Diabetik" = "#3B8BD4", "Diabetes T2D" = "#D85A30")

df_plot <- data.frame(
  SV1      = skor_svd[, 1],
  SV2      = skor_svd[, 2],
  Kelompok = kelompok,
  Sampel   = nama_samp
)

# Panel 1: Score Plot SV1 vs SV2 
p1 <- ggplot(df_plot, aes(x = SV1, y = SV2, color = Kelompok)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
  scale_color_manual(values = warna_kel) +
  labs(
    title = "Score Plot SV1 vs SV2",
    x     = sprintf("SV1 (%.1f%%)", var_exp10[1]),
    y     = sprintf("SV2 (%.1f%%)", var_exp10[2]),
    color = "Kelompok"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )
print(p1)

# --- Panel 2: Scree Plot ---
df_scree <- data.frame(
  Komponen   = paste0("SV", 1:n_comp),
  Var_Exp    = var_exp10,
  Kumulatif  = var_cum
)
df_scree$Komponen <- factor(df_scree$Komponen, levels = df_scree$Komponen)

p2 <- ggplot(df_scree, aes(x = Komponen)) +
  geom_col(aes(y = Var_Exp), fill = "#5DCAA5", color = "white", width = 0.7) +
  geom_line(aes(y = Var_Exp, group = 1), color = "#0F6E56", linewidth = 1.2) +
  geom_point(aes(y = Var_Exp), color = "#0F6E56", size = 2.5) +
  labs(
    title = "Scree Plot (Singular Values)",
    x     = "Komponen SVD",
    y     = "Variance Explained (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p2)

# Panel 3: Bar Chart Kandidat Biomarker
top20_plot          <- top20_biomarker
top20_plot$Gen      <- factor(top20_plot$Gen,
                              levels = top20_plot$Gen[order(top20_plot$Abs_SV1)])
top20_plot$Warna    <- ifelse(top20_plot$SV1 > 0, "Positif", "Negatif")

p3 <- ggplot(top20_plot, aes(x = Gen, y = SV1, fill = Warna)) +
  geom_col(width = 0.75, color = "white") +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Positif" = "#D85A30", "Negatif" = "#3B8BD4")) +
  labs(
    title = "Top 20 Kandidat Biomarker Gen (SV1)",
    x     = "Gen",
    y     = "Loading SV1",
    fill  = "Arah Loading"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )
print(p3)

# Panel 4: Biplot SVD
# Skala loading agar proporsional dengan skor
skala   <- D10[1] * 0.6
top10g  <- rownames(head(loading_sorted, 10))

df_vec  <- data.frame(
  x   = loading_df[top10g, "SV1"] * skala,
  y   = loading_df[top10g, "SV2"] * skala,
  Gen = top10g
)

p4 <- ggplot() +
  geom_point(data = df_plot,
             aes(x = SV1, y = SV2, color = Kelompok),
             size = 3.5, alpha = 0.85) +
  geom_segment(data = df_vec,
               aes(x = 0, y = 0, xend = x, yend = y),
               arrow     = arrow(length = unit(0.25, "cm"), type = "closed"),
               color     = "#2E7D32", linewidth = 0.8, alpha = 0.8) +
  geom_label_repel(data = df_vec,
                   aes(x = x, y = y, label = Gen),
                   size = 2.8, color = "#1B5E20",
                   box.padding = 0.3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  scale_color_manual(values = warna_kel) +
  labs(
    title = "Biplot SVD + Kandidat Biomarker Gen",
    x     = sprintf("SV1 (%.1f%%)", var_exp10[1]),
    y     = sprintf("SV2 (%.1f%%)", var_exp10[2]),
    color = "Kelompok"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )
print(p4)

# Gabungkan 4 panel
plot_gabung <- grid.arrange(
  p1, p2, p3, p4,
  ncol = 2,
  top  = grid::textGrob(
    "SVD - Data Ekspresi Gen Diabetes Tipe 2 (GSE25724)",
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
)

PLOT_FILE <- "D:/Statistika Part 2/MODERN PREDIKSI DAN MACHINE LEARNING/SVD_DT2.png"
ggsave(PLOT_FILE, plot = plot_gabung,
       width = 16, height = 12, dpi = 150, bg = "white")
cat(sprintf("   → Plot disimpan: %s\n", PLOT_FILE))
