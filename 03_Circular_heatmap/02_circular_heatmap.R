library(circlize)
library(dplyr)
library(ComplexHeatmap)
library(grid)

# Set working directory & load data from previous pipeline step
data_matrix <- read.csv('/path/to/output/sig_table.tsv', sep='\t')

# Format labels & set order
data_matrix <- data_matrix %>%
  mutate(type = case_when(
    type == "metabolites" ~ "Metabolites",
    type == "taxa" ~ "Microbiome",  
    type == "gmm" ~ "GMM",
    type == "metacyc" ~ "MetaCyc",
    TRUE ~ type
  ))

desired_order <- c("Microbiome", "Metabolites", "MetaCyc", "GMM")
data_matrix$type <- factor(data_matrix$type, levels = desired_order)
data_matrix <- data_matrix %>% arrange(type)

track_cols <- c("HFrEF", "ACS", "SA")
vmax_abundance <- 0.7
vmax_coef <- 1.5

# Define colors
zero_col <- "#EAEAEA" 
col_fun_list <- list(
  Metabolites = colorRamp2(c(-vmax_abundance, 0, vmax_abundance), c("#3182BD", zero_col, "#F768A1")),
  Microbiome  = colorRamp2(c(-vmax_abundance, 0, vmax_abundance), c("#54278F", zero_col, "#00A087")),
  GMM         = colorRamp2(c(-vmax_abundance, 0, vmax_abundance), c("#238B45", zero_col, "#CE1256")),
  MetaCyc     = colorRamp2(c(-vmax_abundance, 0, vmax_abundance), c("#3F51B5", zero_col, "#FF7F00"))
)

type_colors <- c(
  Microbiome  = "#4DAF4A", 
  Metabolites = "#E41A1C", 
  MetaCyc     = "#377EB8", 
  GMM         = "#FF7F00"  
)

coef_col_fun <- colorRamp2(c(-vmax_coef, 0, vmax_coef), c("#4575B4", "#FFFFBF", "#D73027"))

# Insert gap sectors
new_rows <- list()
for (i in 1:(nrow(data_matrix) - 1)) {
  new_rows[[length(new_rows) + 1]] <- data_matrix[i, ]
  if (data_matrix$type[i] != data_matrix$type[i + 1]) {
    gap_row <- data_matrix[i, ]
    gap_row$features <- paste0("gap_", i)
    gap_row$type <- "Gap"
    gap_row[track_cols] <- NA
    gap_row$coef <- NA 
    new_rows[[length(new_rows) + 1]] <- gap_row
  }
}
last_row <- data_matrix[nrow(data_matrix), ]
new_rows[[length(new_rows) + 1]] <- last_row
gap_row <- last_row
gap_row$features <- "gap_end"
gap_row$type <- "Gap"
gap_row[track_cols] <- NA
gap_row$coef <- NA
new_rows[[length(new_rows) + 1]] <- gap_row

data_matrix2 <- bind_rows(new_rows)
n2 <- nrow(data_matrix2)
types2 <- data_matrix2$type
all_types_unique <- levels(data_matrix$type) 

gap_after2 <- rep(0.15, n2)
for (i in 1:(n2 - 1)) {
  if (!is.na(types2[i]) && types2[i] == "Gap") {
    gap_after2[i] <- 3
  }
}
gap_after2[n2] <- 60

# Define Legends
legend_gap <- 2

title_legend <- Legend(
  labels = "Feature Categories",
  legend_gp = gpar(fill = NA),
  labels_gp = gpar(fontsize = 16, fontface = "bold")
)

type_legend_items <- list(
  Legend(labels = "Microbiome", legend_gp = gpar(fill = type_colors["Microbiome"]), labels_gp = gpar(fontsize = 14)),
  Legend(labels = "Metabolome", legend_gp = gpar(fill = type_colors["Metabolites"]), labels_gp = gpar(fontsize = 14)),
  Legend(labels = "MetaCyc", legend_gp = gpar(fill = type_colors["MetaCyc"]), labels_gp = gpar(fontsize = 14)),
  Legend(labels = "GMM", legend_gp = gpar(fill = type_colors["GMM"]), labels_gp = gpar(fontsize = 14))
)

type_legend_items_packed <- packLegend(
  list = type_legend_items,
  direction = "vertical",
  gap = unit(2, "mm"),
  align = TRUE
)

color_legends <- lapply(desired_order, function(t) {
  Legend(
    col_fun = col_fun_list[[t]],
    title = t,
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    direction = "horizontal",
    legend_width = unit(4, "cm"),
    at = seq(-vmax_abundance, vmax_abundance, length.out = 3) 
  )
})

coef_legend <- Legend(
  col_fun = coef_col_fun,
  title = "Coef",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  direction = "horizontal",
  legend_width = unit(4, "cm"),
  at = c(-vmax_coef, 0, vmax_coef)
)

color_legends <- c(list(coef_legend), color_legends)

lgd_colors <- packLegend(
  list = color_legends,
  direction = "vertical",
  gap = unit(2, "mm")
)

# ========================
# Wrap into a single function (logic preserved, margins tightened)
# ========================
draw_my_circos <- function() {
  par(mar = c(0.5, 0.5, 0.5, 0.5)) # Forcefully remove default R margins
  
  lim_val <- 0.75
  circos.clear()
  
  circos.par(
    gap.after = gap_after2,
    start.degree = 90,
    clock.wise = FALSE, 
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.001, 0.001),
    canvas.xlim = c(-lim_val, lim_val),
    canvas.ylim = c(-lim_val, lim_val)
  )
  
  circos.initialize(
    factors = data_matrix2$features, 
    xlim = cbind(rep(0, n2), rep(2.5, n2))
  )
  
  # Track 1: Labels
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.5,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      if (!startsWith(sector.index, "gap_")) {
        circos.text(
          x = 1, y = 0.01,
          labels = sector.index,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = 0.7,
          col = "black"
        )
      }
    }
  )
  
  # Track 2: Coef
  circos.par(track.margin = c(0.015, 0.001)) 
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.02,   
    bg.border = NA,
    panel.fun = function(x, y) {
      i <- get.cell.meta.data("sector.numeric.index")
      sector_type <- data_matrix2$type[i]
      
      if (sector_type == "Gap") {
        circos.rect(0, 0, 2.5, 1, col = NA, border = NA)
      } else {
        value <- data_matrix2$coef[i]
        if(!is.na(value)) {
          circos.rect(0, 0, 2.5, 1, col = coef_col_fun(value), border = "white", lwd = 0.2)
        }
      }
    }
  )
  circos.par(track.margin = c(0.001, 0.001)) 
  
  # Track 3: Abundance
  for (colname in track_cols) {
    values <- data_matrix2[[colname]]
    
    circos.trackPlotRegion(
      ylim = c(0, 1),
      track.height = 0.055,   
      bg.border = NA,
      panel.fun = function(x, y) {
        i <- get.cell.meta.data("sector.numeric.index")
        sector_type <- data_matrix2$type[i]
        
        if (sector_type == "Gap") {
          circos.rect(0, 0, 2.5, 1, col = NA, border = NA)
        } else {
          value <- values[i]
          col_fun <- col_fun_list[[sector_type]]
          if (is.null(col_fun)) {
            col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
          }
          if(!is.na(value)) {
            circos.rect(0, 0, 2.5, 1,
                        col = col_fun(value), border = "white", lwd = 0.2)
          }
        }
      }
    )
  }
  
  # Track 4: Category highlight
  circos.par(track.margin = c(0.001, 0.02))
  circos.track(
    track.height = 0.02,
    ylim = c(0, 1),
    bg.border = NA,
    panel.fun = function(x, y) {}
  )
  circos.par(track.margin = c(0.001, 0.001))
  
  for (t in all_types_unique) {
    features_in_type <- data_matrix$features[data_matrix$type == t] 
    highlight.sector(
      sector.index = features_in_type,
      track.index = get.current.track.index(),
      col = type_colors[t],
      text = "",
      border = NA
    )
  }
  
  draw(title_legend, x = unit(0.5, "npc"), y = unit(0.535, "npc"), just = c("center", "bottom"))
  draw(type_legend_items_packed, x = unit(0.5, "npc"), y = unit(0.525, "npc"), just = c("center", "top"))
  draw(lgd_colors, x = unit(0.692, "npc"), y = unit(0.837, "npc"), just = c("right", "top"))
  
  circos.clear()
}

# ========================
# 1. Export to PNG (12x12 inches, 600dpi)
# ========================
png("/path/to/output/plot/circular_heatmap.png", width = 7200, height = 7200, res = 600)
draw_my_circos()
dev.off()

# ========================
# 2. Export to PDF (12x12 inches)
# ========================
pdf("/path/to/output/plot/circular_heatmap.pdf", width = 12, height = 12)
draw_my_circos()
dev.off()