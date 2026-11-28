# Figure fixes after revision D (run in the same R session right after script 07)
# Figure 4: aligned text panel
t4[, hdr := "Adjusted OR (95% CI), q"]
g4t <- ggplot(t4, aes(0, gen, label = txt)) + geom_text(hjust = 0, size = 3) + facet_wrap(~hdr) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0.02)) + labs(x = " ", y = NULL) + thm +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none")
save_d(g4 + g4t + plot_layout(widths = c(4, 1.2)), "Figure4_coalterations", 14, 5)
# Figure 6: P value labels
br[, lab := sapply(c(p_mss, p_msi, p_alt), function(p) if (startsWith(p, "<")) paste("P", p) else paste("P =", p))]
g6 <- ggplot(hd, aes(grpl, HRD_Score, fill = BRCA2_any)) + geom_hline(yintercept = 42, linetype = 2, colour = "grey40") +
  geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
  geom_text(data = br, aes(x = (x + xend) / 2, y = y + 3.5, label = lab), inherit.aes = FALSE, size = 3) +
  annotate("text", x = 4.45, y = 44.5, label = "HRD score 42", size = 2.8, hjust = 1, colour = "grey30") +
  scale_fill_manual(values = c("grey85", pal[["red"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Genomic HRD score", caption = sprintf("BRCA2 x MSI/HM-SNV interaction in a linear model adjusted for site: P = %s", pfm(pint))) +
  theme_bw(base_size = 10)
save_d(g6, "Figure6_TCGA_HRD", 8, 5.5)
# Figure 7: distinct colors per row, legend in row order
cols <- c(MSI = "#C0392B", `HM-SNV` = "#F1948A", CIN = "#5D6D7E", GS = "#76D7C4", EBV = "#8E7CC3", Missing = "#D5D8DC",
          `Deep deletion` = "#17202A", Truncating = "#2874A6", `Non-truncating only` = "#AED6F1",
          `Homozygous deletion` = "#4A235A", `Pathogenic mutation with locus LOH` = "#1E8449",
          `Multiple truncating hits, phase unresolved` = "#D68910", `Single truncating or curated event` = "#F7DC6F",
          `Non-truncating variant(s) only` = "#E5E7E9", Yes = "#7B241C", No = "grey92", `Not available` = "white")
tiles[, val := factor(val, levels = names(cols))]
g7b <- ggplot(tiles, aes(id, row, fill = val)) + geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(~grp, scales = "free_x", space = "free_x") + scale_fill_manual(values = cols, breaks = names(cols), name = NULL, drop = TRUE) +
  labs(x = "BRCA2-altered TCGA tumors (n = 43)", y = NULL) + theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), strip.text = element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, byrow = FALSE))
save_d(g7a / g7b + plot_layout(heights = c(1, 1.6)), "Figure7_TCGA_oncoprint", 11, 7.5)
cat("Figure fixes done\n")
