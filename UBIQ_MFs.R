## Scripts used in data processing and visualization for manuscript
## Ivanova, A. et al: "Converging molecular composition of groundwater dissolved organic matter across contrasting aquifer systems"

library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan) 
library(ggrepel)
library(ggpubr)
library(patchwork)
library(tibble)
library(ggbreak)
library(ggExtra)
library(ggpattern)

# load data frame with assigned molecular formulae (filtered to remove singlets)

load("C:/Users/aivanova/Nextcloud/PhD/WP2/averaged_MF_no_singlets.RData")

# well order for plotting and well colors
well_order <- c("H14", "H32",  "H41", "H51", "H43", "H52", "H53", "S1", "S2")
colors <- c( "H14" = "#871818", "H32" = "#cd5c5c", "H41" = "#ff6833", "H43" = "#1eb879",
             "H51" = "#f5e23b", "H52" = "#2395d6", "H53" = "#233ed6", "S1" = "#d678dc", "S2" = "#da2f94")
# for time series: dates of sampling campaigns
dataref_new <- read.csv("C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/dateref_new.csv")
dataref_new$date <- as.Date(dataref_new$date, format = "%d/%m/%Y")

## PCA on DOM composition:
averaged_MF_for_PCA <- filtered_MF %>%
  dplyr::select(formula, Sample_name, Intensity) %>%
  arrange(formula, Sample_name) %>%
  pivot_wider(
    names_from = Sample_name,
    values_from = Intensity,
    values_fill = list(Intensity = 0)
  ) %>%
  tibble::column_to_rownames("formula")

forPCA_MF_hel <- averaged_MF_for_PCA %>% 
  filter(rowSums(.) != 0) %>%
  t() %>%
  decostand(method = "hellinger")

PCA <- prcomp(forPCA_MF_hel, center = TRUE, scale. = TRUE)
PCA_variance <- PCA$sdev^2 / sum(PCA$sdev^2)

PCA_scores <- PCA$x %>%
  as.data.frame() %>%
  mutate(sample_name = rownames(.),
         well = str_split(sample_name, "_", simplify = TRUE)[,2],
         well = factor(well, levels = well_order),
         site = case_when(
           well %in% c("H14","H32","H41","H51","H43","H52","H53") ~ "Hainich CZE",
           well %in% c("S1","S2") ~ "SESO"
         ))

### DOM bulk properties (to fit on PCA plot)

DOM_bulk_PCA <- filtered_MF %>%
  mutate(AI_mod = ((1+ C- 0.5*O- S -0.5*(H +N +P))/(C- 0.5*O- S- N- P)),
         AI_mod = case_when(AI_mod == -Inf ~ 0,
                            AI_mod == Inf  ~ NA_real_,
                            TRUE ~ AI_mod),
         N_C = N/C) %>%
  group_by(Sample_name) %>%
  summarise(
    AI_wm   = weighted.mean(AI_mod, Intensity, na.rm = TRUE),
    H_C_wm  = weighted.mean(H_C, Intensity, na.rm = TRUE),
    O_C_wm  = weighted.mean(O_C, Intensity, na.rm = TRUE),
    N_C_wm  = weighted.mean(N_C, Intensity, na.rm = TRUE),
    DBE_wm  = weighted.mean(DBE, Intensity, na.rm = TRUE),
    mz_wm   = weighted.mean(mz, Intensity, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )%>%
  column_to_rownames("Sample_name")

bulk_all <- filtered_MF %>%
  group_by(Sample_name) %>%
  summarise(
    shannon = vegan::diversity(Intensity, index="shannon"),
    .groups = "drop"
  )%>%
  inner_join(DOM_bulk_PCA%>%rownames_to_column(var = "Sample_name"), by = "Sample_name")%>%
  column_to_rownames("Sample_name")


DOM_bulk_PCA <- DOM_bulk_PCA[ rownames(PCA_scores), ]

envfit_res <- envfit(PCA, DOM_bulk_PCA, permutations = 999,  choices = c(2, 3))

sig_vectors <- as.data.frame(vegan::scores(envfit_res, display = "vectors"))
sig_vectors$pval <- envfit_res$vectors$pvals
sig_vectors <- sig_vectors[sig_vectors$pval <= 0.05, ]  # keep only significant vectors
sig_vectors$var <- rownames(sig_vectors)

new_names <- c(
  "AI_wm"  = "AI[wm]",
  "H_C_wm" = "H/C[wm]",
  "O_C_wm" = "O/C[wm]",
  "N_C_wm" = "N/C[wm]",
  "DBE_wm" = "DBE[wm]",
  "mz_wm"  = "m/z[wm]",
  "n" = "Richness"
)

sig_vectors$var <- new_names[ sig_vectors$var ]

arrow_scaling = 40

## FIG. 2A:
bulk_PCA_DOM_MF <- ggplot(data = PCA_scores, aes(x = PC2,y = PC3)) + #axis 1 corresponds to time shift due to measurements instability, see Supplementary Information
  geom_point(aes(color = well, shape = site), size = 4)+
  scale_colour_manual(values = colors)+
  theme_minimal()+
  labs(x = paste0("PC2: ", signif(PCA_variance[2]*100,2),"%"),
       y = paste0("PC3: ", signif(PCA_variance[3]*100,2),"%")) +
  geom_segment(data = sig_vectors,
               aes(x=0, xend = PC2*arrow_scaling, y=0, yend = PC3*arrow_scaling),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black") +
  geom_text_repel(data = sig_vectors,
                  aes(x=PC2*arrow_scaling, y=PC3*arrow_scaling, label = var),
                  parse = TRUE,
                  size = 5)+
  theme(legend.position = "none")

#PLOTS for FIGURE 2B (PCA_GWchem) and 2C (wl_plot) ARE CREATED IN THE SCRIPT water_chem_dispersal.R found in https://github.com/arina-iva/Groundwater_DOM_convergence/ 

## FIG. 2:

(bulk_PCA_DOM_MF | (PCA_GWchem / wl_plot))+ plot_layout(guides = "collect", width = c(2, 1)) + plot_annotation(tag_levels = 'A')

## Comparison of N-containing DOM and NO3 levels for Supplementary Information:
DOM_bulk_NC <- DOM_bulk_PCA %>%
  rownames_to_column(var = "Sample_name")%>%
  dplyr::select(Sample_name, N_C_wm) %>%
  mutate(
    Well =  str_split(string = Sample_name, simplify = TRUE, 
                      pattern = "_")[, 2],
    Well = factor(Well, levels = well_order))

Bulk_NC_plot <- ggplot(DOM_bulk_NC, aes(x = Well, y = N_C_wm))+
  geom_boxplot(aes(fill=Well), alpha=0.9, outlier.shape = NA) + 
  geom_jitter(aes(color = Well), alpha = 0.2, show.legend = FALSE)+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  labs(y = bquote(N/C[wm]))+           
  theme_minimal()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        legend.position = "none")

## Supplementary fig. 10:
## NO3_plot is created in the script water_chem_dispersal.R found in https://github.com/arina-iva/Groundwater_DOM_convergence/ 
( NO3_plot | Bulk_NC_plot) + plot_annotation(tag_levels = 'A')

##Assigning CRAM formulae and calculating DOM properties

averaged_MF_CRAM <- filtered_MF%>%
  mutate(CRAM = case_when((DBE/C >= 0.3 & DBE/C <= 0.68 & 
                             DBE/H >= 0.2 & DBE/H <= 0.95 &
                             DBE/O >= 0.77 & DBE/O <= 1.75) ~ "CRAM",
                          TRUE~"not_CRAM"),
         Labile = case_when(H/C >= 1.5 ~ "Labile", 
                            TRUE ~ "not_labile"),
         AI_mod = ((1+ C- 0.5*O- S -0.5*(H +N +P))/(C- 0.5*O- S- N- P)),
         AI_mod = case_when(AI_mod == -Inf ~ 0,
                            AI_mod == Inf  ~ NA_real_,
                            TRUE ~ AI_mod),
         NOSC = 4 - (4*C + H - 2*O -3 *N +5*P - 2*S)/C,
         GFE = 60.3 - 28.5*NOSC,
         N_C = N/C)

formula_info <- averaged_MF_CRAM %>%
  dplyr::select(formula, C, H, O, N, P, S, mz, group, CRAM, Labile, AI_mod, NOSC, GFE, O_C, H_C, N_C)%>%
  distinct()
################################################################################
## FINDING UBIQUITOUS FORMULAE

## Finding unique and ubiquitous formulae

unique_MFs <- averaged_MF_CRAM %>%
  mutate(presence = Intensity > 0) %>%
  group_by(formula, wells) %>%
  dplyr::summarize(detections = sum(presence), .groups = "drop") %>%
  filter(detections >= 2)%>% # only keep formulae that are detected at least twice in each well
  group_by(formula) %>%
  filter(sum(detections > 0) == 1) %>%       # detected in only 1 well
  ungroup()%>%
  left_join(formula_info, by = "formula")

samples_per_well <- averaged_MF_CRAM%>%
  group_by(wells)%>%
  dplyr::summarize(n_samples = length(unique(Sample_name)))

ubiq_MFs <- averaged_MF_CRAM %>%
  mutate(presence = Intensity > 0) %>%
  group_by(formula, wells) %>% # how many times is formula detected per well
  dplyr::summarize(detections = sum(presence), .groups = "drop") %>% 
  ungroup()%>%
  left_join(samples_per_well, by = "wells")%>%
  mutate(ratio_present = detections/n_samples)%>%
  filter(ratio_present >= 0.6)%>%  #molecules should be detected in 60% of sample at least
  group_by(formula) %>%
  filter(
    sum(detections > 0) == 9) %>%       # detected in all nine wells
  ungroup()%>%
  dplyr::select(formula)%>%
  distinct()%>%
  left_join(formula_info, by = "formula")%>%
  mutate(type = "ubiquitous")

all_MFs <- formula_info %>% mutate(type = "all")

### How much of ubiq per Sample?
MF_type <- unique_MFs %>%
  dplyr::select(-wells, -detections)%>%
  mutate(type = "unique")%>%
  rbind(ubiq_MFs)%>%
  dplyr::select(formula, type)

MF_ubiq_RA <- averaged_MF_CRAM %>%
  left_join(MF_type, by = "formula")%>%
  filter(!is.na(type))%>%
  group_by(Sample_name)%>%
  dplyr::summarise(
    unique_RA = sum(Intensity[type == "unique"], na.rm = TRUE),
    ubiq_RA = sum(Intensity[type == "ubiquitous"], na.rm = TRUE))%>%
  ungroup()%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))

## Fig. 3A:
ubiqDOM_RA <- ggplot(MF_ubiq_RA, aes(x = well, y = ubiq_RA))+
  annotate(
    "rect",
    xmin = 0.5,
    xmax = 7.5,
    ymin = -Inf, ymax = Inf,
    fill = "#40E0D0", alpha = 0.15   # turquoise
  ) +
  annotate(
    "rect",
    xmin = 7.5,
    xmax =9.5,
    ymin = -Inf, ymax = Inf,
    fill = "violet", alpha = 0.15)+
  geom_text(data = data.frame(well = c("H51", "S1"), y = c(1, 1), label = c("Hainich CZE", "SESO")),
            aes(x = well, y = y, label = label),
            inherit.aes = FALSE,
            color = "darkgray",
            hjust = 0,
            vjust = 1,
            size = 7)+
  geom_boxplot(aes(fill=well), alpha=0.9, outlier.shape = NA) + 
  geom_jitter(aes(color = well), alpha = 0.2, show.legend = FALSE)+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  labs(y = "Ubiquitous DOM, RA", x = "Well")+            
  theme_minimal()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
ubiqDOM_RA + labs( y = "Ubiquitous DOM, RA")

#### COMPARE BULK PARAMETERS OF UBIQUITOUS AND ALL FORMULAE

ubiq_weigh_data <- averaged_MF_CRAM %>%
  inner_join((ubiq_MFs%>%dplyr::select(formula, type)), by = "formula")%>%
  group_by(Sample_name)%>%
  dplyr::summarise(
    AI_wm    = weighted.mean(AI_mod, Intensity, na.rm = TRUE),
    H_C_wm   = weighted.mean(H_C, Intensity, na.rm = TRUE),
    O_C_wm   = weighted.mean(O_C, Intensity, na.rm = TRUE),
    N_C_wm   = weighted.mean(N_C, Intensity, na.rm = TRUE),
    DBE_wm   = weighted.mean(DBE, Intensity, na.rm = TRUE),
    mz_wm    = weighted.mean(mz, Intensity, na.rm = TRUE),
    GFE_wm = weighted.mean(GFE, Intensity, na.rm = TRUE),
    .groups = "drop"
  )%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))%>%
  pivot_longer(cols = ends_with("wm"), names_to = "variable", values_to = "variable_wm" )%>%
  mutate(type = "ubiq")

all_weigh_data <- averaged_MF_CRAM %>%
  group_by(Sample_name)%>%
  dplyr::summarise(
    AI_wm    = weighted.mean(AI_mod, Intensity, na.rm = TRUE),
    H_C_wm   = weighted.mean(H_C, Intensity, na.rm = TRUE),
    O_C_wm   = weighted.mean(O_C, Intensity, na.rm = TRUE),
    N_C_wm   = weighted.mean(N_C, Intensity, na.rm = TRUE),
    DBE_wm   = weighted.mean(DBE, Intensity, na.rm = TRUE),
    mz_wm    = weighted.mean(mz, Intensity, na.rm = TRUE),
    GFE_wm = weighted.mean(GFE, Intensity, na.rm = TRUE),
    .groups = "drop"
  )%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))%>%
  pivot_longer(cols = ends_with("wm"), names_to = "variable", values_to = "variable_wm" )%>%
  mutate(type = "all")%>%
  rbind(ubiq_weigh_data)%>%
  mutate(variable_plot = recode(variable, 
                                "O_C_wm" = "O/C", 
                                "H_C_wm" = "H/C",
                                "N_C_wm" = "N/C",
                                "DBE_wm" = "DBE",
                                "DBEO_wm" = "DBE-O",
                                "AI_wm" = "AI",
                                "O_wm" = "nO",
                                "mz_wm" = "m/z",
  ))

ubiq_vs_all_aver <- all_weigh_data %>%
  group_by(variable) %>%
  dplyr::summarise(
    p_value_wilcox = wilcox.test(variable_wm ~ type)$p.value,
    p_adjusted_BH = p.adjust(p_value_wilcox, method = "BH"),
    significance = case_when(
      is.na(p_adjusted_BH) ~ "N/A",
      p_adjusted_BH < 0.001 ~ "***",
      p_adjusted_BH < 0.01 ~ "**",
      p_adjusted_BH < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    .groups = 'drop'
  )
## the only statistically significant difference: N/C

## Fig 3B:
N_C_wm <- ggplot(all_weigh_data%>%filter(variable == "N_C_wm"), aes(x = type, 
                                                                    y = variable_wm, pattern = type, fill = type)) + 
  geom_boxplot_pattern(outlier.shape = NA,pattern_fill = "white",pattern_density = 0.3,pattern_spacing = 0.025
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 6) +
  scale_pattern_manual(values = c("none", "stripe"))+
  scale_fill_manual(values = c("#2395d6", "#d678dc"))+
  theme_minimal()+
  theme(strip.text = element_text(size = 11, face = "bold"), axis.text = element_text(size = 11),
        legend.position = "none")+
  scale_x_discrete(
    labels = c(ubiq = "Ubiquitous DOM", all = "Full DOM")
  )+
  labs(x = "", y = bquote(N/C[wm]))+
  coord_cartesian(clip = "off") 

##### COMPARE RELATIVE ABUNDANCES OF DIFFERENT DOM CLASSES IN UBIQUITOUS AND ALL FORMULAE
averaged_MF_CRAM <- averaged_MF_CRAM %>%
  mutate(MF_class = case_when(
    H_C >= 1.53 & H_C <= 2.20 & O_C >= 0.56 & O_C <= 1.23 ~ "carbohydrate",
    H_C >= 0.86 & H_C <= 1.34 & O_C >= 0.21 & O_C <= 0.44 ~ "lignins",
    H_C >= 0.70 & H_C <= 1.01 & O_C >= 0.16 & O_C <= 0.84 ~ "tannins",
    H_C >= 1.34 & H_C <= 2.18 & O_C >= 0.01 & O_C <= 0.35 ~ "lipid-like",
    H_C >= 1.33 & H_C <= 1.84 & O_C >= 0.17 & O_C <= 0.48 & N > 0 ~ "peptide-like",
    AI_mod > 0.66 ~ "condensed aromatics",
    TRUE ~ "other"))

ubiq_DOM_classes_RA <- averaged_MF_CRAM %>%
  dplyr::filter(formula%in%ubiq_MFs$formula)%>%
  group_by(Sample_name)%>%
  mutate(Int_norm = Intensity/sum(Intensity)) %>% # re-normalize intensity after subsetting only ubiquitous MFs
  ungroup()%>%
  group_by(Sample_name, MF_class) %>%
  dplyr::summarize(intensity_cl = sum(Int_norm), .groups="drop")%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))%>%
  mutate(type = "ubiquitous")

all_DOM_classes_RA <- averaged_MF_CRAM%>%
  group_by(Sample_name, MF_class) %>%
  dplyr::summarize(intensity_cl = sum(Intensity), .groups="drop")%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))%>%
  mutate(type = "all")%>%
  rbind(ubiq_DOM_classes_RA)

all_DOM_CRAM <- averaged_MF_CRAM %>%
  group_by(Sample_name, CRAM)%>%
  dplyr::summarize(intensity_CRAM = sum(Intensity), .groups="drop")%>%
  mutate(type = "Full DOM")

ubiq_DOM_CRAM <- averaged_MF_CRAM %>%
  dplyr::filter(formula%in%ubiq_MFs$formula)%>%
  group_by(Sample_name)%>%
  mutate(Int_norm = Intensity/sum(Intensity)) %>% # re-normalize intensity after subsetting only ubiquitous MFs
  ungroup()%>%
  group_by(Sample_name, CRAM) %>%
  dplyr::summarize(intensity_CRAM = sum(Int_norm), .groups="drop")%>%
  mutate(type = "Ubiquitous DOM")%>%
  rbind(all_DOM_CRAM)


#Fig. 3E :
Fig_3_E <- ggplot(ubiq_DOM_CRAM%>%filter(CRAM == "CRAM"), aes(x = type, y = intensity_CRAM, pattern = type, fill = type)) + 
  geom_boxplot_pattern(outlier.shape = NA,pattern_fill = "white",pattern_density = 0.3,pattern_spacing = 0.025
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 6) +
  scale_pattern_manual(values = c("none", "stripe"))+
  scale_fill_manual(values = c("#2395d6", "#d678dc"))+
  theme_minimal()+
  theme(strip.text = element_text(size = 11, face = "bold"), axis.text = element_text(size = 11),
        legend.position = "none")+
  labs(x = "", y = "CRAM, RA")+
  coord_cartesian(clip = "off") 


## supplementary figure 5 (comparing RA of all classes in ubiq. and full DOM)
ggplot(all_DOM_classes_RA, aes(x = interaction(well, type, lex.order = TRUE), 
                               y = intensity_cl, fill = well, pattern = type)) + 
  geom_boxplot_pattern(outlier.shape = NA,pattern_fill = "white",pattern_density = 0.3,pattern_spacing = 0.025
  ) +
  scale_fill_manual(values = colors) +
  scale_pattern_manual(values = c("none", "stripe")) +
  facet_wrap(~MF_class, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))+
  labs(y = "RA of a class", x = "")

##is the differences of relative abundances of DOM classes statistically significant?
ubiq_vs_all_classes_RA <- all_DOM_classes_RA %>%
  group_by(MF_class, type) %>%
  dplyr::summarise(
    mean_int = mean(intensity_cl),
    sd_int = sd(intensity_cl))

ubiq_vs_all_classes <- all_DOM_classes_RA %>%
  dplyr:: filter(MF_class != "condensed aromatics")%>% #present only in full DOM pool
  group_by(MF_class) %>%
  dplyr::summarise(
    p_value_wilcox = wilcox.test(intensity_cl ~ type)$p.value,
    p_adjusted_Bon = p.adjust(p_value_wilcox, method = "bonferroni"),
    significance = case_when(
      is.na(p_adjusted_Bon) ~ "N/A",
      p_adjusted_Bon < 0.001 ~ "***",
      p_adjusted_Bon < 0.01 ~ "**",
      p_adjusted_Bon < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    .groups = 'drop'
  )

imp_classes <- ubiq_vs_all_classes%>%
  filter(p_adjusted_Bon < 0.05)%>%
  pull(MF_class)

# Fig. 3C
imp_cl1 <- ggplot(all_DOM_classes_RA%>%filter(MF_class == imp_classes[1]), aes(x = type, 
                                                                               y = intensity_cl, pattern = type, fill = type)) + 
  geom_boxplot_pattern(outlier.shape = NA,pattern_fill = "white",pattern_density = 0.3,pattern_spacing = 0.025
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 6) +
  scale_pattern_manual(values = c("none", "stripe"))+
  scale_fill_manual(values = c("#2395d6", "#d678dc"))+
  theme_minimal()+
  theme(strip.text = element_text(size = 11, face = "bold"), axis.text = element_text(size = 11),
        legend.position = "none")+
  scale_x_discrete(
    labels = c(ubiquitous = "Ubiquitous DOM", all = "Full DOM")
  )+
  labs(x = "", y = paste0(imp_classes[1], ", RA"))+
  coord_cartesian(clip = "off") 
# Fig. 3D
imp_cl2 <- ggplot(all_DOM_classes_RA%>%filter(MF_class == imp_classes[2]), aes(x = type, 
                                                                               y = intensity_cl, pattern = type, fill = type)) + 
  geom_boxplot_pattern(outlier.shape = NA,pattern_fill = "white",pattern_density = 0.3,pattern_spacing = 0.025
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 6) +
  scale_pattern_manual(values = c("none", "stripe"))+
  scale_fill_manual(values = c("#2395d6", "#d678dc"))+
  theme_minimal()+
  theme(strip.text = element_text(size = 11, face = "bold"), axis.text = element_text(size = 11),
        legend.position = "none")+
  scale_x_discrete(
    labels = c(ubiquitous = "Ubiquitous DOM", all = "Full DOM")
  )+
  labs(x = "", y = paste0(imp_classes[2], ", RA"))+
  coord_cartesian(clip = "off") 

# Compare CRAM and Ubiquitous and all DOM on van Krevelin diagram
CRAM_MFs <- all_MFs%>%
  filter(CRAM == "CRAM")%>%
  mutate(type = "CRAM")

MFs_for_VK <- rbind(all_MFs, CRAM_MFs, ubiq_MFs)

colour_scale <- scale_color_manual(
  values = c("all" = "#2395d6", "CRAM" = "#333333", "ubiquitous" = "#d678dc")
)

p_1 <- ggplot(MFs_for_VK, aes(x = O_C, y = H_C)) +
  geom_point(aes(color = type), size = 1, alpha = 0.7) +
  xlim(0, 1.5) +
  theme_minimal() +
  theme(legend.position = "none",  legend.title = element_blank()) + 
  colour_scale +
  labs(x = "O/C", y = "H/C")

p_2 <- ggplot(MFs_for_VK, aes(x = mz, y = H_C)) +
  geom_point(aes(color = type), size = 1, alpha = 0.7) +
  theme_minimal() +
  colour_scale +
  xlim(0, 1050) +
  labs(x = "m/z", y = "H/C")+
  theme(legend.position = "none",  legend.title = element_blank())+
  theme(
    legend.position = c(.99, .0),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.title = element_blank(),
    legend.key.spacing.y = unit(1, "pt")#,
    #legend.margin = margin(t=0,b=0,unit='pt')
  )

p_1_dens <- ggMarginal(p_1, type = "density", margins = "x",    size = 5, groupFill = TRUE)
p_2_dens <- ggMarginal(p_2, type = "density", margins = "both", size = 5, groupFill = TRUE)

#Supplementary Fig. 6:

wrap_elements(wrap_elements(p_1_dens) | wrap_elements(p_2_dens))

# Fig. 3:

((ubiqDOM_RA + theme(legend.position = "none")) | ((N_C_wm| imp_cl1)/(imp_cl2 |Fig_3_E )))+ 
  plot_layout(width = c(2, 3))  + plot_annotation(tag_levels = 'A')

####################################
## Gibbs free energy distribution of different formulae:

gfe_combined <- bind_rows(
  all_MFs %>% dplyr::select(formula, GFE) %>% mutate(Pool = "Full DOM"),
  ubiq_MFs %>% dplyr::select(formula, GFE) %>% mutate(Pool = "Ubiquitous DOM")
)

# Kernel density plot
#Fig. 4A:
GFE <- ggplot(gfe_combined, aes(x = GFE, fill = Pool, color = Pool)) +
  geom_density(alpha = 0.4) +
  geom_vline(data = gfe_combined %>% group_by(Pool) %>% summarise(mean_GFE = mean(GFE)),
             aes(xintercept = mean_GFE, color = Pool),
             linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = c("Full DOM" = "#2395d6", "Ubiquitous DOM" = "#d678dc")) + 
  scale_color_manual(values = c("Full DOM" = "#2395d6", "Ubiquitous DOM" = "#d678dc")) +
  labs(x = "Gibbs free energy (kJ/mol C)", y = "Relative frequency \nof molecular formulae",
       fill = "DOM pool", color = "DOM pool") +
  theme_minimal()+
  theme(legend.position = c(.0, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "left")

## Temporal stability of ubiquitous vs full DOM pools: dispersal of composition between different time points:

## dispersal of full DOM dataset:
dist_mat_all <- dist(PCA_scores[, 2:155 ], method = "euclidean")

samples_wells <- PCA_scores%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample_name")%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))%>%
  dplyr::select(Sample_name, well)

betadisper_full <- betadisper(dist_mat_all,  samples_wells$well)

##dispersal on ubiq. DOM only:
ubiq_for_PCA <- filtered_MF %>%
  dplyr::filter(formula%in%ubiq_MFs$formula)%>% ##keep only ubiquitous formulae
  dplyr::select(formula, Sample_name, Intensity) %>%
  group_by(Sample_name)%>%
  mutate(Int_norm = Intensity/sum(Intensity))%>% # re-normalize intensity after subsetting only ubiquitous formulae
  ungroup()%>%
  dplyr::select(-Intensity)%>%
  arrange(formula, Sample_name) %>%
  pivot_wider(
    names_from = Sample_name,
    values_from = Int_norm,
    values_fill = list(Int_norm = 0)
  ) %>%
  tibble::column_to_rownames("formula") %>% 
  dplyr::filter(rowSums(.) != 0) %>%
  decostand(method = "total", MARGIN = 2) %>%
  t() %>%
  decostand(method = "hellinger")

ubiq_PCA  <- prcomp(ubiq_for_PCA, center = TRUE, scale. = TRUE)
ubiq_PCA_scores <- ubiq_PCA$x

dist_mat_ubiq <- dist(ubiq_PCA_scores[, 2:155 ], method = "euclidean")

betadisper_ubiq <- betadisper(dist_mat_ubiq, samples_wells$well)

dispersion_data <- data.frame(
  Sample_name = c(names(betadisper_ubiq$distances), 
                  names(betadisper_full$distances)),
  Distance = c(betadisper_ubiq$distances, 
               betadisper_full$distances),
  Pool = rep(c("Ubiquitous", "Full"), 
             each = length(betadisper_ubiq$distances))
)

dispersion_data_av <- dispersion_data%>%
  group_by(Pool)%>%
  summarize(Distance_av = mean(Distance))


dispersion_data_p <- dispersion_data %>%
  dplyr::summarise(
    p_value_wilcox = wilcox.test(Distance ~ Pool)$p.value,
    p_adjusted_Bon = p.adjust(p_value_wilcox, method = "bonferroni"),
    significance = case_when(
      is.na(p_adjusted_Bon) ~ "N/A",
      p_adjusted_Bon < 0.001 ~ "***",
      p_adjusted_Bon < 0.01 ~ "**",
      p_adjusted_Bon < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    .groups = 'drop'
  )

## Fig. 4B
disp <- ggplot(dispersion_data, aes(x = Pool, y = Distance, fill = Pool, pattern = Pool)) +
  geom_boxplot_pattern(outlier.shape = NA,pattern_fill = "white",pattern_density = 0.3,pattern_spacing = 0.025
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 8) +
  scale_pattern_manual(values = c("none", "stripe"))+
  scale_fill_manual(values = c("Ubiquitous" = "#d678dc", "Full" = "#2395d6")) +
  labs(y = "Temporal variability of \nDOM composition",
       x = "DOM pool") +
  theme_minimal() +
  theme(legend.position = "none", axis.text=element_text(size=12)) +
  guides(fill = "none")+
  coord_cartesian(clip = "off") 

#Fig 4C is created in a separate data file: ms1_ms2_pipeline.R available at https://github.com/arina-iva/Groundwater_DOM_convergence 

## Fig 4:
((GFE | disp) / wrap_elements(plot_343)) + plot_annotation(tag_levels = "A")+plot_layout(heights = c(2, 3))

##################################################################################################################
## DOM age and presence of ubiquitous formulae

#load 14C data
GW_14C_av <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/GW_14C_av.csv")%>%
  dplyr::select(F14C_aver, Sample_name)%>%
  filter(!is.na(F14C_aver))%>%
  mutate(del_14C = (F14C_aver-1)*1000)

DOM_bulk_14C_av <- DOM_bulk_PCA%>%
  rownames_to_column(var = "Sample_name")%>%
  inner_join(GW_14C_av, by = "Sample_name")%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))%>%
  group_by(well)%>%
  summarize(n_av= mean(n),
            F14C_av = mean(F14C_aver),
            del_14C_av = mean(del_14C))

DOM_bulk_14C_av_H <- DOM_bulk_14C_av %>% filter(!well%in%c("S1", "S2"))

ubiq_MFs_14C <- MF_ubiq_RA%>% 
  inner_join(GW_14C_av, by = "Sample_name")%>%
  mutate(Site = case_when(well%in%c("S1", "S2") ~ "SESO", TRUE ~ "Hainich CZE"))

ubiq_MFs_14C_aver <- ubiq_MFs_14C%>%
  group_by(well)%>%
  summarize(
    ubiq_RA_av = median(ubiq_RA),
    ubiq_RA_sd = sd(ubiq_RA),
    del_14C_av = median(del_14C),
    del_14C_sd = sd(del_14C)
  )%>%
  mutate(Site = case_when(well%in%c("S1", "S2") ~ "SESO", TRUE ~ "Hainich CZE"))

## Fig. 5, A:
RA_14C <- ggplot()+
  geom_errorbar(data = ubiq_MFs_14C_aver, aes(x = del_14C_av,  ymin = ubiq_RA_av - ubiq_RA_sd, ymax = ubiq_RA_av + ubiq_RA_sd), color = "darkgray") +
  geom_errorbar(data = ubiq_MFs_14C_aver, aes(y = ubiq_RA_av, xmin = del_14C_av - del_14C_sd,  xmax = del_14C_av + del_14C_sd), orientation = "y", color = "darkgray") +
  geom_point(data = ubiq_MFs_14C, aes(y = ubiq_RA , x =  del_14C, color = well, shape = Site), size = 3, alpha = 0.2)+
  geom_point(data = ubiq_MFs_14C_aver, aes(y = ubiq_RA_av , x =  del_14C_av, color = well, shape = Site), size = 5)+
  scale_color_manual(values = colors)+
  labs(y = "Ubiquitous DOM, RA", 
       x = bquote(Delta~{}^{14}*C~", ‰"))+
  theme_minimal()+
  theme(axis.text = element_text(size = 12), legend.position = "none")

## CRAM and 14C age
formula_info <- averaged_MF_CRAM %>%
  dplyr::select(formula, C, H, O, N, P, S, mz, group, CRAM, Labile)%>%
  distinct()%>%
  mutate(H_C = H/C,
         O_C = O/C)

## calculate summary per well: CRAM, labile - # of MF and assigned intensity

CRAM_per_sample <- averaged_MF_CRAM%>%
  group_by(Sample_name)%>%
  summarise(
    # CRAM metrics
    CRAM_RA = sum(Intensity[CRAM == "CRAM"], na.rm = TRUE),
    not_CRAM_RA = sum(Intensity[CRAM == "not_CRAM"], na.rm = TRUE),
    CRAM_n = sum(CRAM == "CRAM"),
    not_CRAM_n = sum(CRAM == "not_CRAM")
  ) %>%
  ungroup()%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))%>%
  inner_join(GW_14C_av, by = "Sample_name")%>%
  mutate(site = case_when(well%in%c("S1", "S2") ~ "SESO",
                          TRUE ~ "Hainich CZE"))

CRAM_14C_p <- ggplot(CRAM_per_sample, aes(x = del_14C, y = CRAM_RA))+
  geom_point(aes(color = well, shape = site), size = 3)+
  scale_color_manual(values = colors)+
  labs(y = "CRAM, RA", x = bquote(Delta~{}^{14}*C~", ‰"))+
  theme_minimal()

DOM_divers <- bulk_all%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample_name")%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         well = factor(well, levels = well_order))


### TIME SERIES for well H41 only:
## DOM aromaticity, DOM bulk age, DOM diversity and microbial diversity
# DOM aromaticity plot: Fig.5, D

H41_divers_time <- bulk_all%>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample_name")%>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         PNK = as.numeric(str_split(string = Sample_name, simplify = TRUE, 
                                    pattern = "_")[, 1])
  )%>%
  inner_join(dataref_new, by = "PNK")%>%
  dplyr::filter(well == "H41")%>%
  mutate(date = as.Date(date, format = "%d/%m/%Y"))

H41_14C <- GW_14C_av %>%
  mutate(well = str_split(string = Sample_name, simplify = TRUE, 
                          pattern = "_")[, 2],
         PNK = as.numeric(str_split(string = Sample_name, simplify = TRUE, 
                                    pattern = "_")[, 1])
  )%>%
  dplyr::filter(well == "H41")

H41_divers_time <- H41_divers_time %>%
  inner_join(H41_14C, by = c("PNK", "well"))
# Fig. D/ Supplementary Fig. S12 C: DOM diversity 
H41_shan_p <- ggplot(H41_divers_time, aes(x = date, y = shannon, group = 1))+
  geom_point(color = "#ff6833")+
  geom_line(color = "#ff6833")+
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  labs(y = "DOM diversity\n (Shannon's H')")
H41_AI_p <- ggplot(H41_divers_time, aes(x = date, y = AI_wm, group = 1))+
  geom_point(color = "#ff6833")+
  geom_line(color = "#ff6833")+
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  labs(y = bquote(AI[wm]))
# DOM bulk 14C values plot: Fig.5, D
H41_14C_p <- ggplot(H41_divers_time, aes(x = date, y = del_14C, group = 1))+
  geom_point(color = "#ff6833")+
  geom_line(color = "#ff6833")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y  = bquote(Delta~{}^{14}*C~", ‰"))
# Are DOM diversity, aromaticity and 14C values and/or thei change concordant in time?
x1 <- H41_divers_time$shannon #absolute values
x2 <- H41_divers_time$AI_wm
x3 <- H41_divers_time$del_14C

d1 <- diff(x1) #change
d2 <- diff(x2)
d3 <- diff(x3)
kendall.global(cbind(d1, d2, d3), nperm = 9999)

kw_result <- kendall.global(cbind(x1, x2, x3), nperm = 9999)

kw_W    <- kw_result$Concordance_analysis["W", "Group.1"]
kw_p    <- kw_result$Concordance_analysis["Prob.perm", "Group.1"]

kw_label <- sprintf("Kendall's W = %.2f, p = %.4f", kw_W, kw_p)

##Fig. 5D:
H41_p <- wrap_elements((H41_AI_p / H41_shan_p / H41_14C_p) + plot_layout(guides = "collect") + plot_annotation(title = "DOM temporal dynamics in well H41", subtitle = kw_label))

## Fig. 5:
p5 <- "ABC
DD#"

(RA_14C +(CRAM_14C_p + theme(legend.position = "none"))+ ubiq_chem+ H41_p) + plot_annotation(tag_levels = "A")+plot_layout(design = p5)

## Including microbial diversity

microb_shan <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/microbial_data/Microb_com_shannon.csv")%>%
  mutate(well = str_split(string = sample, simplify = TRUE, 
                          pattern = "_")[, 1],
         well = factor(well, levels = well_order),
         PNK = str_split(string = sample, simplify = TRUE, 
                         pattern = "_")[, 2],
         Sample_name = paste0(PNK, "_", well))

micr_divers <- microb_shan 

#Suppl. Fig. 12A
DOM_div <- ggplot(DOM_divers, aes(x = well, y = shannon))+
  geom_boxplot(aes(fill = well), alpha = 0.9)+
  geom_jitter(aes(color = well), size = 2)+
  theme_minimal()+
  labs(y = "DOM diversity\n (Shannon's H')")+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  theme(legend.position = "none")
#Suppl. Fig. 12B
Mic_div <- ggplot(micr_divers, aes(x = well, y = shannon))+
  geom_boxplot(aes(fill = well), alpha = 0.9)+
  geom_jitter(aes(color = well), size = 2)+
  theme_minimal()+
  labs(y = "Microbial diversity\n (Shannon's H')")+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  theme(legend.position = "none")

microb_shan_H41 <- microb_shan %>%
  dplyr::filter(well == "H41")%>%
  mutate(PNK=as.numeric(PNK))%>%
  rename(shannon_micr = shannon)%>%
  inner_join(H41_divers_time, by = c("PNK", "well"))

# Are DOM and microbial diversities and/or their change concordant in time?
m1 <- microb_shan_H41$shannon #absolute values
m2 <- microb_shan_H41$shannon_micr
dm1 <- diff(m1) #change
dm2 <- diff(m2)

kendall.global(cbind(m1, m2), nperm = 9999) #are absolute values concordant?
kendall.global(cbind(dm1, dm2), nperm = 9999)  # is the change concordant?

kw_result <- kendall.global(cbind(x1, x2, x3), nperm = 9999)

# Supl. Fig. Fig. S12 C: MICROBIAL diversity 
H41_shan_micr_p <- ggplot(microb_shan_H41, aes(x = date, y = shannon_micr, group = 1))+
  geom_point(color = "#ff6833")+
  geom_line(color = "#ff6833")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "Microbial diversity \n(Shannon's H')")

# Supplemetary Fig.12
H41_div_DOM_micr_p <- wrap_elements(H41_shan_p / H41_shan_micr_p)
Fig_12 <- ((DOM_div | Mic_div) / H41_div_DOM_micr_p)+ plot_layout(guides = "collect")+ plot_annotation(tag_levels = "A")

## MICROBIAL DATA 

#microbial data:

microbe_wide  <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/microbial_data/Rel_abund_16S.csv")%>%
  dplyr::select(-partial_16S)

tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus", "ASV")
sample_cols_micr <- setdiff(colnames(microbe_wide), tax_cols)

colnames(microbe_wide)[match(sample_cols_micr, colnames(microbe_wide))] <- 
  sapply(strsplit(sample_cols_micr, "_"), function(x) paste0(x[2], "_", x[1]))

## samples that are both in DOM data and in microbial data:
samples_all <- Reduce(intersect, list(unique(filtered_MF$Sample_name),
                                      colnames(microbe_wide)[!(colnames(microbe_wide) %in% 
                                                                 tax_cols)]))

microb_wide_order <- microbe_wide %>%
  mutate(
    Order_new = case_when(
      !is.na(Order)   ~ Order,
      is.na(Order) & !is.na(Class)   ~ paste0(Class, "_unclassified_order"),
      is.na(Order) & is.na(Class) & !is.na(Phylum) ~ paste0(Phylum, "_unclassified_order"),
      is.na(Order) & is.na(Class) & is.na(Phylum) & !is.na(Kingdom) ~ paste0(Kingdom, "_unclassified_order"),
      TRUE ~ "Unclassified"
    )
  )

# Prevalence  filtering for order-level microbial data

microbe_order_new <- microb_wide_order %>%
  group_by(Order_new) %>%
  dplyr::summarize(across(all_of(samples_all), sum), .groups="drop")

abund_mat <- as.matrix(microbe_order_new[, samples_all])
prev <- apply(abund_mat > 0, 1, mean)
mean_abund <- rowMeans(abund_mat)
keep <- which(prev >= 0.05 & mean_abund >= 0.001 
)

microbe_order_filt_new <- microbe_order_new[keep, ]
# Renormalize to relative abundance after filtering

X <- microbe_order_filt_new %>%
  column_to_rownames(var = "Order_new")%>%
  sweep(2, colSums(.), "/")
X[is.na(X)] <- 0
microbe_mat <- t(X)

# CLR-transformation function
clr_safe <- function(mat){
  repX <- zCompositions::cmultRepl(t(mat), label=0, method="CZM", z.delete = FALSE) # replace 0s for CLR transformations by small pseudo-counts
  clrX <- compositions::clr(t(repX))
  return(clrX)
}
microbe_clr  <- clr_safe(microbe_mat)
microbe_dist <- dist(microbe_clr, method="euclidean") # Aitchison

### Mantel tests on DOM compositional data from PCA: 
### DOM formulas => re-normalize if subsetting => Hellinger => PCA => Axes 2-5 => eucledian distance => mantel, partial mantel

#1) for the total DOM pool, we can take PCA axes from the calculations above
PCA_all_DOM_a25 <- PCA$x[, 2:5] 
summary(PCA)

setdiff(rownames(microbe_mat),rownames(PCA_all_DOM_a25)) # 1 sample difference 

PCA_all_DOM_a25 <-PCA_all_DOM_a25 %>%as.data.frame()%>% dplyr::filter(rownames(.)%in%rownames(microbe_mat))
PCA_all_DOM_a25 <- PCA_all_DOM_a25[match(rownames(microbe_mat), rownames(PCA_all_DOM_a25)), ]
identical(rownames(PCA_all_DOM_a25), rownames(microbe_mat))

DOM_all_PCA_dist <- dist(PCA_all_DOM_a25, method = "euclidean")

#2) ubiquitous DOM pool: full DOM => subset => re-normalize => hellinger => PCA => axes 2-5 => eucledian distance
PCA_ub_DOM <- averaged_MF_CRAM %>% 
  dplyr::filter(formula %in% ubiq_MFs$formula) %>%
  dplyr::filter(Sample_name %in% rownames(microbe_mat)) %>%
  group_by(Sample_name) %>%
  mutate(Int_norm = Intensity / sum(Intensity))%>%
  dplyr::select(formula, Sample_name, Int_norm) %>%
  arrange(formula, Sample_name) %>%
  pivot_wider(
    names_from = Sample_name,
    values_from = Int_norm,
    values_fill = list(Int_norm = 0)
  ) %>%
  tibble::column_to_rownames("formula")%>%
  filter(rowSums(.) != 0) %>%
  t() %>%
  decostand(method = "hellinger")

PCA_ub_DOM <- prcomp(PCA_ub_DOM, center = TRUE, scale. = TRUE)
summary(PCA_ub_DOM)
PCA_ub_DOM_a25 <- PCA_ub_DOM$x[, 2:5]
PCA_ub_DOM_a25 <- PCA_ub_DOM_a25[match(rownames(microbe_mat), rownames(PCA_ub_DOM_a25)), ]
identical(rownames(PCA_ub_DOM_a25), rownames(microbe_mat))

DOM_ub_PCA_dist <- dist(PCA_ub_DOM_a25, method = "euclidean")

#3) same for non-ubiquitous DOM pool

PCA_nu_DOM <- averaged_MF_CRAM %>% 
  dplyr::filter(!formula %in% ubiq_MFs$formula) %>%
  dplyr::filter(Sample_name %in% rownames(microbe_mat)) %>%
  group_by(Sample_name) %>%
  mutate(Int_norm = Intensity / sum(Intensity))%>%
  dplyr::select(formula, Sample_name, Int_norm) %>%
  arrange(formula, Sample_name) %>%
  pivot_wider(
    names_from = Sample_name,
    values_from = Int_norm,
    values_fill = list(Int_norm = 0)
  ) %>%
  tibble::column_to_rownames("formula")%>%
  filter(rowSums(.) != 0) %>%
  t() %>%
  decostand(method = "hellinger")

PCA_nu_DOM <- prcomp(PCA_nu_DOM, center = TRUE, scale. = TRUE)
summary(PCA_nu_DOM)
PCA_nu_DOM_a25 <- PCA_nu_DOM$x[, 2:5]
PCA_nu_DOM_a25 <- PCA_nu_DOM_a25[match(rownames(microbe_mat), rownames(PCA_nu_DOM_a25)), ]
identical(rownames(PCA_nu_DOM_a25), rownames(microbe_mat))

DOM_nu_PCA_dist <- dist(PCA_nu_DOM_a25, method = "euclidean")

#4) prepare sample distance matrix for partial mantel
samples_mic_DOM <- rownames(PCA_nu_DOM_a25)
well_ids <- sub(".*_", "", samples_mic_DOM)  

# Binary same-well matrix: 1 = same well, 0 = different well
well_mat <- outer(well_ids, well_ids, FUN = "==") * 1
well_dist <- as.dist(1 - well_mat)

#5) full mantel test

mantel_PCA_all   <- mantel(DOM_all_PCA_dist,   microbe_dist, 
                           method = "spearman", permutations = 9999)

mantel_PCA_ub   <- mantel(DOM_ub_PCA_dist,   microbe_dist, 
                          method = "spearman", permutations = 9999)

mantel_PCA_nu   <- mantel(DOM_nu_PCA_dist,   microbe_dist, 
                          method = "spearman", permutations = 9999)

#6) partial mantel test - temporal effect (within each well)

mantel_part_PCA_all   <- mantel.partial(DOM_all_PCA_dist,   microbe_dist, well_dist,
                                        method = "spearman", permutations = 9999)

mantel_part_PCA_ub   <- mantel.partial(DOM_ub_PCA_dist,   microbe_dist, well_dist,
                                       method = "spearman", permutations = 9999)

mantel_part_PCA_nu   <- mantel.partial(DOM_nu_PCA_dist,   microbe_dist, well_dist,
                                       method = "spearman", permutations = 9999)

#### PROCRUSTES TEST TO CHECK SPATIAL EFFECTS (alignment on well-level)

DOM_ub_coords <- PCA_ub_DOM_a25
DOM_nu_coords <- PCA_nu_DOM_a25

# microbial side - ordination
microbe_PCA <- prcomp(microbe_clr, center = TRUE, scale. = FALSE)  # CLR data shouldn't be re-scaled
microbe_coords <- microbe_PCA$x[, 1:4] #same number of axes as in DOM data
summary(microbe_PCA)
identical(rownames(DOM_ub_coords), rownames(microbe_coords)) #checking that the order is the same
well_ids <- sub(".*_", "", rownames(DOM_ub_coords))

# average position per well - only testing spatial effect, so remove temporal variability
collapse_to_well <- function(score_matrix, well_ids) {
  df <- as.data.frame(score_matrix)
  df$well <- well_ids
  well_means <- df %>%
    group_by(well) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    tibble::column_to_rownames("well")
  as.matrix(well_means)
}

PCA_ub_DOM_well <- collapse_to_well(PCA_ub_DOM_a25, well_ids)
PCA_nu_DOM_well <- collapse_to_well(PCA_nu_DOM_a25, well_ids)
microbe_well    <- collapse_to_well(microbe_coords[, 1:2 ], well_ids)

protest_well_ub_microbe <- protest(X = PCA_ub_DOM_well, Y = microbe_well, permutations = 9999)
protest_well_nu_microbe <- protest(X = PCA_nu_DOM_well, Y = microbe_well, permutations = 9999)

## chemodiversity vs radiocarbon age

bulk_vs_age <- bulk_all %>%
  as.data.frame()%>%
  rownames_to_column(var = "Sample_name")%>%
  inner_join(GW_14C_av, by = "Sample_name")%>%
  mutate(well = str_split(Sample_name, "_", simplify = TRUE)[,2],
         well = factor(well, levels = well_order),
         site = case_when(
           well %in% c("H14","H32","H41","H51","H43","H52","H53") ~ "Hainich CZE",
           well %in% c("S1","S2") ~ "SESO"
         ))

bulk_vs_age_H_aver <- bulk_vs_age %>%dplyr::filter(site == "Hainich CZE")%>%
  group_by(well)%>%
  summarize(
    shannon_av = mean(shannon),
    shannon_med = median(shannon),
    richness_av = mean(richness),
    richness_med = median(richness),
    del_14C_av = mean (del_14C),
    del_14C_med = median(del_14C)
  )

ct <- cor.test(x = bulk_vs_age_H_aver$shannon_med, 
               y =  bulk_vs_age_H_aver$del_14C_med, 
               method = "kendall") #for small # of observations

# Format label
cortest_label <- sprintf("τ = %.2f, p = %.3f", ct$estimate, ct$p.value)

## Supplementary fig. 8:

ggplot(bulk_vs_age_H_aver, aes(x = del_14C_med, y =shannon_med))+
  geom_smooth(data = bulk_vs_age_H_aver,
              mapping = aes(x = del_14C_med, y = shannon_med),
              method = "lm", se = FALSE, color = "grey60", linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(color = well), size = 4)+
  annotate("text",
           x = -500, y = 6.8,
           hjust = 1.1, vjust = 1.5,
           label = cortest_label, size = 4) +
  theme_minimal()+
  scale_color_manual(values = colors)+
  labs(y = "Median molecular diversity (Shannon's H')", x = "Δ14C, ‰")


#########################################################################################################################################
### Variance partitioning

## microbial data (aggregated to order level, CLR-transformed, same samples as DOM data):
microbe_clr

## water chemistry data (scaled, made in the script cor_wat_stability_vs_ubiq_DOM.R, can be found at ):
env_data_for_PCA

## determine shared sample pool:
samples_varpart <- intersect(rownames(microbe_clr), rownames(env_data_for_PCA))

## PCA for microbial data:
microbe_clr_vp_PCA <- prcomp(microbe_clr_vp, center = TRUE, scale. = FALSE)
summary(microbe_clr_vp_PCA) # how much variability is explained
microbe_clr_vp_PCA_scores <- microbe_clr_vp_PCA$x[ , 1:2]

## PCA for water chemistry data: 
env_data_for_vp_PCA <- env_data_for_vp %>% prcomp(center = TRUE, scale. = FALSE)
env_data_for_vp_PCA_scores <- env_data_for_vp_PCA$x[ , 1:2]
summary(env_data_for_vp_PCA) # how much variability is explained
colnames(env_data_for_vp_PCA_scores) <- paste("WC", colnames(env_data_for_vp_PCA_scores), sep = "_")

################################################################################
## VARIANCE PARTITIONING BASED ON PCA DATA FOR DOM
## New PCAs because the set of samples is different:
#1) ubiquitous DOM
PCA_ub_DOM_f <- averaged_MF_CRAM %>% 
  dplyr::filter(formula %in% ubiq_MFs$formula) %>%
  dplyr::filter(Sample_name %in% samples_varpart) %>%
  group_by(Sample_name) %>%
  mutate(Int_norm = Intensity / sum(Intensity))%>%
  dplyr::select(formula, Sample_name, Int_norm) %>%
  arrange(formula, Sample_name) %>%
  pivot_wider(
    names_from = Sample_name,
    values_from = Int_norm,
    values_fill = list(Int_norm = 0)
  ) %>%
  tibble::column_to_rownames("formula")%>%
  filter(rowSums(.) != 0) %>%
  t() %>%
  decostand(method = "hellinger")

PCA_ub_DOM_f <- prcomp(PCA_ub_DOM_f, center = TRUE, scale. = TRUE)
summary(PCA_ub_DOM_f) # how much variability is explained
PCA_ub_DOM_a25_f <- PCA_ub_DOM_f$x[, 2:5]

# 2) non-ubiquitous DOM
PCA_nu_DOM_f <- averaged_MF_CRAM %>% 
  dplyr::filter(!formula %in% ubiq_MFs$formula) %>%
  dplyr::filter(Sample_name %in% samples_varpart) %>%
  group_by(Sample_name) %>%
  mutate(Int_norm = Intensity / sum(Intensity))%>%
  dplyr::select(formula, Sample_name, Int_norm) %>%
  arrange(formula, Sample_name) %>%
  pivot_wider(
    names_from = Sample_name,
    values_from = Int_norm,
    values_fill = list(Int_norm = 0)
  ) %>%
  tibble::column_to_rownames("formula")%>%
  filter(rowSums(.) != 0) %>%
  t() %>%
  decostand(method = "hellinger")

PCA_nu_DOM_f <- prcomp(PCA_nu_DOM_f, center = TRUE, scale. = TRUE)
summary(PCA_nu_DOM_f) # how much variability is explained
PCA_nu_DOM_a25_f <- PCA_nu_DOM_f$x[, 2:5]
PCA_nu_DOM_a25_f <- PCA_nu_DOM_a25_f[match(rownames(PCA_ub_DOM_a25_f), rownames(PCA_nu_DOM_a25_f)), ]

# 3) Environmental data - match sample order:

env_data_for_vp_PCA_scores <- env_data_for_vp_PCA_scores[match(rownames(PCA_ub_DOM_a25_f), rownames(env_data_for_vp_PCA_scores)), ]
identical(rownames(PCA_ub_DOM_a25_f), rownames(env_data_for_vp_PCA_scores))
identical(rownames(PCA_nu_DOM_a25_f), rownames(env_data_for_vp_PCA_scores))

# 4) microbial data - match sample order
identical(rownames(env_data_for_vp_PCA_scores), rownames(microbe_clr_vp_PCA_scores))
microbe_clr_vp_PCA_scores <- microbe_clr_vp_PCA_scores[match(rownames(PCA_nu_DOM_a25_f), rownames(microbe_clr_vp_PCA_scores)), ]

# 5) site data - a vector for permutations:
well_vp <- sub(".*_", "", rownames(env_data_for_vp_PCA_scores))%>%as.factor()

## variance partitioning:

vp_ub_DOM_PCA <- varpart(PCA_ub_DOM_a25_f, 
                         microbe_clr_vp_PCA_scores, 
                         env_data_for_vp_PCA_scores) 


vp_nu_DOM_PCA <- varpart(PCA_nu_DOM_a25_f, 
                         microbe_clr_vp_PCA_scores, 
                         env_data_for_vp_PCA_scores) 


##visualizing

vp_results <- data.frame(
  pool     = rep(c("Ubiquitous", "Non-ubiquitous"), each = 4),
  fraction = rep(c("Microbial community", 
                   "Water chemistry", 
                   "Shared",
                   "Unexplained"), 2),
  value    = c(
    vp_ub_DOM_PCA$part$indfract[1, 3],   
    vp_ub_DOM_PCA$part$indfract[2, 3],   
    vp_ub_DOM_PCA$part$indfract[3, 3],   
    vp_ub_DOM_PCA$part$indfract[4, 3],  
    vp_nu_DOM_PCA$part$indfract[1, 3],
    vp_nu_DOM_PCA$part$indfract[2, 3],
    vp_nu_DOM_PCA$part$indfract[3, 3],
    vp_nu_DOM_PCA$part$indfract[4, 3]))%>%
  mutate(label = paste0(round(value * 100), "%"),
         fraction = factor(fraction, levels = c("Unexplained", "Water chemistry", "Shared", "Microbial community"))
  )

# Fig. 6:
ggplot(vp_results, aes(x = pool, y = value, fill = fraction)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = c(
    "Microbial community"        = "#2166ac",   # blue
    "Shared"                  = "#762a83",   # purple
    "Water chemistry" = "#da2f94",   # red
    "Unexplained"             = "grey80"
  )) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  labs(x = "DOM pool",
       y = "Proportion of variance",
       fill = NULL)+
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "right",
    axis.text        = element_text(color = "black", size = 12)
  )

h <- how(
  nperm = 999,
  blocks = well_vp
)

## testing significance for ubiquitous and non-ubiquitous pools:
full_nu <-  capscale(PCA_nu_DOM_a25_f ~ ., data =  cbind(env_data_for_vp_PCA_scores%>%as.data.frame(),microbe_clr_vp_PCA_scores%>%as.data.frame()))
full_ub <-  capscale(PCA_ub_DOM_a25_f ~ ., data =  cbind(env_data_for_vp_PCA_scores%>%as.data.frame(),microbe_clr_vp_PCA_scores%>%as.data.frame()))
set.seed(12345)

sign_full_nu <- anova(
  full_nu,
  permutations = 999
)
# restricted permutations (per well) - removing pseudoreplication
sign_full_nu_r <- anova(
  full_nu,
  permutations = h
)

sign_full_ub <- anova(
  full_ub,
  permutations = 999
)
# restricted permutations (per well) - removing pseudoreplication
sign_full_ub_r <- anova(
  full_ub,
  permutations = h
)

## significance of individual fractions - water chemistry and microbes

pure_chem_nu <- capscale(PCA_nu_DOM_a25_f ~ ., data =  env_data_for_vp_PCA_scores%>%as.data.frame(),
                         Condition = microbe_clr_vp_PCA_scores%>%as.data.frame())
pure_chem_ub <- capscale(PCA_ub_DOM_a25_f ~ ., data =  env_data_for_vp_PCA_scores%>%as.data.frame(),
                         Condition = microbe_clr_vp_PCA_scores%>%as.data.frame())

sign_chem_nu <- anova(
  pure_chem_nu,
  permutations = 999
)

sign_chem_nu_r <- anova(
  pure_chem_nu,
  permutations = h
)

sign_chem_ub <- anova(
  pure_chem_ub,
  permutations = 999
)

sign_chem_ub_r <- anova(
  pure_chem_ub,
  permutations = h
)


pure_micr_nu <- capscale(PCA_nu_DOM_a25_f ~ ., data =  microbe_clr_vp_PCA_scores%>%as.data.frame(),
                         Condition = env_data_for_vp_PCA_scores%>%as.data.frame())

pure_micr_ub <- capscale(PCA_ub_DOM_a25_f ~ ., data =  microbe_clr_vp_PCA_scores%>%as.data.frame(),
                         Condition = env_data_for_vp_PCA_scores%>%as.data.frame())

sign_micr_nu <- anova(pure_micr_nu, 
                      permutations = 999
)

sign_micr_nu_r <- anova(pure_micr_nu,
                        permutations = h
)
sign_micr_ub <- anova(pure_micr_ub,
                      permutations = 999
)

sign_micr_ub_r <- anova(pure_micr_ub,
                        permutations = h
)

## visualizing average class composition of microbial community per well (sup. fig. S11)
## Visualizing the most abundant taxa (averaged per well):
sample_cols <- setdiff(names(microb_wide_order%>%dplyr::select(-Order_new)), tax_cols)

# Convert from wide to long format
microbe_vis <- microb_wide_order %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  mutate(
    # Extract everything after the first "_"
    well = sub("^[^_]+_", "", Sample)
  )

paired_12 <- brewer.pal(12, "Paired")

class_vis <- microbe_vis %>%
  mutate(
    Class = replace_na(Class, "Unclassified")
  ) %>%
  group_by(Sample, well, Class) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )
class_well_vis <- class_vis %>%
  group_by(well, Class) %>%
  summarise(
    MeanAbundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )
top_class <- class_well_vis %>%
  group_by(Class) %>%
  summarise(
    OverallMean = mean(MeanAbundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(OverallMean)) %>%
  slice_head(n = 12)

class_plot <- class_well_vis %>%
  mutate(
    Class = if_else(
      Class %in% top_class$Class,
      Class,
      "Other"
    )
  ) %>%
  group_by(well, Class) %>%
  summarise(
    MeanAbundance = sum(MeanAbundance),
    .groups = "drop"
  )

microb_color_cl <- c(
  setNames(paired_12, setdiff(unique(class_plot$Class), "Other")),
  "Other" = "darkgrey"
)
# Suppl. Fig. S11:
ggplot(class_plot,
       aes(x = well,
           y = MeanAbundance,
           fill = Class)) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = "Well",
    y = "Mean relative abundance",
    fill = "Class"
  ) +
  theme_minimal()+
  scale_fill_manual(values = microb_color_cl)
