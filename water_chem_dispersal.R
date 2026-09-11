library(vegan)      
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(lubridate)

well_order <- c("H14", "H32", "H41", "H51", "H43", "H52", "H53",  "S1", "S2")

##PCA for water chemistry

# gw chem data
GW_chem_TOC <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/GW_chem.csv")

## GW physicho-chem data:
GW_pc <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/GW_Physikochemistry.csv")

##prep pc data to merge
GW_pc_forcor <- GW_pc%>%
  rename(Sample_name = sample_name)%>%
  dplyr::select(-date, -PNK)

env_data_all <- GW_pc_forcor%>%
  dplyr::select(-Well)%>%
  inner_join(GW_chem_TOC, by = "Sample_name")

env_data_for_PCA <- env_data_all %>%
  drop_na() %>%
  distinct(Sample_name, .keep_all = TRUE)%>%
  tibble::column_to_rownames(var = "Sample_name")%>%
  scale()%>%# z-scaling
  as.data.frame()

env_data_PCA <- env_data_for_PCA%>%
  prcomp(center = TRUE, scale. = FALSE)

env_data_PCA_sum <- summary(env_data_PCA)

env_data_scree <- data.frame(
  PC = 1:length(env_data_PCA$sdev),
  Variance = env_data_PCA$sdev^2,
  Proportion = env_data_PCA_sum$importance[2, ],
  Cumulative = env_data_PCA_sum$importance[3, ]
)

# PCA biplot
env_pca_scores <- as.data.frame(env_data_PCA$x[, 1:2])%>%
  mutate(Sample = rownames(env_data_PCA$x),
         well = str_split_i(Sample, "_", 2),
         well = factor(well, levels = well_order),
         Site = case_when(grepl("^H\\d{2}$", well) ~ "Hainich CZE",
                          grepl("^S\\d{1}$", well) ~ "SESO"))

env_pca_loadings <- as.data.frame(env_data_PCA$rotation[, 1:2]) %>%
  mutate(
    Variable = rownames(.),
    Variable = gsub("Result_", "", Variable),
    
    label = case_when(
      Variable == "O2"  ~ "O[2]",
      Variable == "NO3" ~ "NO[3]^{\"-\"}",
      Variable == "SO4" ~ "SO[4]^{2*\"-\"}",
      Variable == "K"   ~ "K^{\"+\"}",
      Variable == "Na"  ~ "Na^{\"+\"}",
      Variable == "Cl"  ~ "Cl^{\"-\"}",
      Variable == "Ca"  ~ "Ca^{2*\"+\"}",
      Variable == "Mg"  ~ "Mg^{2*\"+\"}",
      Variable == "TOC" ~ "DOC",
      TRUE ~ Variable)
  )

## FIG. 2 B
PCA_GWchem <- ggplot(env_pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 4, aes(color = well, shape = Site)) +
  geom_segment(data = env_pca_loadings, 
               aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5),
               arrow = arrow(length = unit(0.3, "cm")), 
               alpha = 0.7) +
  geom_text_repel(
    data = env_pca_loadings,
    aes(x = PC1 * 5, y = PC2 * 5, label = label),
    parse = TRUE,
    size = 5
  )+
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(env_data_scree$Proportion[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(env_data_scree$Proportion[2] * 100, 1), "%)")) +
  coord_fixed()+
  scale_color_manual(values = colors)

##### WATER LEVEL DEVIATION PLOT

#HAI water level
HAI_wl <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/HAI_waterlevel.csv")

HAI_wl_mean <- HAI_wl %>%
  pivot_longer(cols = -Date,
               names_to = "Well_Level",
               values_to = "Value") %>%
  separate(Well_Level, into = c("Well", "GW_level"), sep = "_") %>%
  dplyr::select(-GW_level)%>%
  mutate(Date = as.Date(Date,format="%d/%m/%Y"),
         Year = year(Date))%>%
  filter(!is.na(Value)) %>%
  group_by(Well, Year) %>%
  summarise(Avg = mean(Value),
            Max = max(Value),
            Min = min(Value),
            n_meas = n(),
            sd_y = sd(Value),
            Variance = stats::var(Value),
            .groups = 'drop')%>%
  filter(!Year %in% c("2019", "2024"))%>%
  filter(!Well %in% c("H13", "H31", "H42"))%>%
  mutate(CV_y = sd_y/Avg,
         CV_y_n = CV_y/sqrt(n_meas))

SESO_wl <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/SESO_waterlevel.csv")
SESO_wl_mean <- SESO_wl%>%
  pivot_longer(cols = -Date,
               names_to = "Well_Level",
               values_to = "Value") %>%
  separate(Well_Level, into = c("Well", "GW_level"), sep = "_") %>%
  dplyr::select(-GW_level)%>%
  mutate(Date = as.Date(Date,format="%d/%m/%Y"),
         Year = year(Date))%>%
  filter(!is.na(Value)) %>%
  group_by(Well, Year) %>%
  summarise(Avg = mean(Value),
            Max = max(Value),
            Min = min(Value),
            n_meas = n(),
            sd_y = sd(Value),
            Variance = stats::var(Value),
            .groups = 'drop')%>%
  filter(!Year %in% c("2019", "2024"))%>%
  mutate(CV_y = sd_y/Avg,
         CV_y_n = CV_y/sqrt(n_meas))


well_order <- c("H14","H32","H41","H51","H43","H52","H53","S1","S2" )
xpos <- seq_along(well_order)
xmap <- data.frame(Well = well_order, x = xpos)
H_range <- range(xmap$x[grepl("^H", xmap$Well)])
S_range <- range(xmap$x[grepl("^S", xmap$Well)])

wl_all_mean <- rbind(HAI_wl_mean, SESO_wl_mean)%>%
  mutate(Well = factor(Well, levels = well_order))

colors <- c( "H14" = "#871818", "H32" = "#cd5c5c", "H41" = "#ff6833", "H43" = "#1eb879",
             "H51" = "#f5e23b", "H52" = "#2395d6", "H53" = "#233ed6", "S1" = "#d678dc", "S2" = "#da2f94")

## FIG. 2C

rect_df <- data.frame(
  xmin = c(H_range[1] - 0.5, S_range[1]-0.5),
  xmax = c(H_range[2]+0.5, S_range[2]+0.5),
  fill = c("#40E0D0", "violet")
)

label_df <- data.frame(
  Well = c(H_range[1]+3, S_range[1]),
  label = c("Hainich CZE", "SESO"),
  y = max(wl_all_mean$sd_y, na.rm = TRUE)
)

wl_plot <- ggplot(wl_all_mean, aes(x = Well, y = sd_y))+
  geom_rect(
    data = rect_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = rect_df$fill,
    alpha = 0.15
  ) +
  
  # boxplot + points
  geom_boxplot(aes(fill = Well), outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(color = Well), alpha = 0.2, width = 0.2) +
  
  # labels (top-left inside rectangles)
  geom_text(
    data = label_df,
    aes(x = Well, y = y, label = label),
    inherit.aes = FALSE,
    color = "darkgray",
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  geom_boxplot(aes(fill = Well), outliers = FALSE, alpha = 0.9)+
  geom_jitter(aes(color = Well), alpha = 0.2)+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  theme_minimal()+
  labs(y = "Intraannual deviation \nof water level, m", x = "")+
  theme(legend.position = "none",
        strip.text.x = element_text(size = 10),
        axis.text = element_text(size = 9))
wl <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/wl_forcor_witherrors.csv")

wl_forcor <- wl %>%
  mutate(Sample_name = paste(PNK, well, sep = "_"))%>%
  dplyr::select(-date, -well) %>%
  drop_na()%>%
  mutate(WL_dev_y_abs = abs(WL_dev_y))

wl_forsel <- wl_forcor %>%
  dplyr::select(WL, Sample_name, WL_dev_y_abs)

env_data_NO3 <- env_data_all%>%
  dplyr::select(Sample_name, Result_NO3) %>%
  mutate(
    Well =  str_split(string = Sample_name, simplify = TRUE, 
                      pattern = "_")[, 2],
    Well = factor(Well, levels = well_order))

##Suppl. Fig. 10
NO3_plot <- ggplot(env_data_NO3, aes(x = Well, y = Result_NO3))+
  geom_boxplot(aes(fill=Well), alpha=0.9, outlier.shape = NA) + 
  geom_jitter(aes(color = Well), alpha = 0.2, show.legend = FALSE)+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  labs(y = bquote(NO[3]^"-"))+           
  theme_minimal()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        legend.position = "none")

env_data_for_PCA <- env_data_all %>%
  dplyr::select( -WL)%>%
  drop_na() %>%
  distinct(Sample_name, .keep_all = TRUE)%>%
  tibble::column_to_rownames(var = "Sample_name")%>%
  scale()%>%# z-scaling
  as.data.frame()

chem_dist <- dist(
  #env_data_all_sc, 
  env_data_for_PCA,
  method = "euclidean")

well_chem_all <- env_data_for_PCA%>%rownames_to_column(var = "Sample_name")%>%mutate(well= str_split_i(Sample_name, "_", 2))
chem_disp <- betadisper(chem_dist, well_chem_all$well)

chem_sample_dist <- data.frame(Sample = attr(chem_dist, "Labels"),
                               dist_to_centroid = chem_disp$distances)%>%
  mutate(
    Well =  str_split(string = Sample, simplify = TRUE, 
                      pattern = "_")[, 2],
    Well = factor(Well, levels = well_order))

chem_sample_dist_av <- chem_sample_dist %>%
  group_by(Well)%>%
  dplyr::summarise(disp_mean = mean(dist_to_centroid),
                   disp_median = median(dist_to_centroid),
                   disp_sd = sd(dist_to_centroid))

chem_disp_ubiq <- chem_sample_dist %>%
  inner_join(MF_ubiq_RA, by = c("Sample" = "Sample_name", "Well" = "well"))%>%
  mutate(Site = case_when(Well%in%c("S1", "S2") ~ "SESO", TRUE ~ "Hainich CZE"))

ubiq_RA_av <- MF_ubiq_RA %>%
  group_by(well)%>%
  dplyr::summarise(ubiq_RA_mean = mean(ubiq_RA),
                   ubiq_RA_median = median(ubiq_RA),
                   ubiq_RA_sd = sd(ubiq_RA))

ubiq_vs_disp <- chem_sample_dist_av %>%
  inner_join(ubiq_RA_av, by = c("Well" = "well"))%>%
  mutate(Site = case_when(Well%in%c("S1", "S2") ~ "SESO", TRUE ~ "Hainich CZE"))

kt <- cor.test(x = ubiq_vs_disp$disp_median, 
               y = ubiq_vs_disp$ubiq_RA_median, 
               method = "kendall") #for small # of observations

# Format label
cor_label <- sprintf("τ = %.2f, p = %.3f", kt$estimate, kt$p.value)

# Fig. 5B:
ubiq_chem <- ggplot() +
  geom_smooth(data = ubiq_vs_disp,
              mapping = aes(x = disp_median, y = ubiq_RA_median),
              method = "lm", se = FALSE, color = "grey60", linetype = "dashed", linewidth = 0.8) +
  geom_errorbar(data = ubiq_vs_disp, aes(x = disp_median,  ymin = ubiq_RA_median - ubiq_RA_sd, ymax = ubiq_RA_median + ubiq_RA_sd), color = "darkgray") +
  geom_errorbar(data = ubiq_vs_disp, aes(y = ubiq_RA_median, xmin = disp_median - disp_sd,  xmax = disp_median + disp_sd), orientation = "y",  color = "darkgray") +
  # geom_point(data = chem_disp_ubiq,
  #            mapping = aes(x = dist_to_centroid, y = ubiq_RA, color = Well, shape = Site),
  #            size = 3, alpha = 0.4) +
  geom_point(data = ubiq_vs_disp,
             mapping = aes(x = disp_median, y = ubiq_RA_median, color = Well, shape = Site),
             size = 5) +
  annotate("text", 
           x = Inf, y = Inf,                          
           hjust = 1.1, vjust = 1.5,                 
           label = cor_label, size = 5) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  ylim(0.83, 1)+
  labs(x = "Temporal variability of hydrochemical properties", y = "Ubiquitous DOM, RA") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16))