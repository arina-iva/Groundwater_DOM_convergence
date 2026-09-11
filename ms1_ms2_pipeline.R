
## MS1/MS2 peak uniqueness + spectral similarity pipeline for manuscript by A. Ivanova et al. 
#"Converging molecular composition of groundwater dissolved organic matter across contrasting aquifer systems"

library(dplyr)
library(tidyr)
library(mzR)
library(ggplot2)
library(patchwork)

samples          <- c("H32", "H43", "H53", "S1")
mz_targets       <- c(343.0, 383.0, 467.0)
precursor_tol    <- 0.4     # +/- window for MS1 precursor selection and MS2 precursor/fragment split
mz_tolerance_ppm <- 0.6     # alignment tolerance across samples

# MS1 precursor filtering: keep peaks above this fraction of the local base peak
ms1_rel_intensity_cutoff <- 0.005   # confirmed reasonable on real data

# MS2 fragment filtering: individual_base only, i.e. keep fragments whose intensity
# exceeds this % of the fragment base peak (max fragment intensity in that sample/mz)
fragment_base_peak_pct <- 1

ms1_dir <- "C:/Users/aivanova/Nextcloud/PhD/WP2/MS2/averaged_MS1"
ms2_dir <- "C:/Users/aivanova/Nextcloud/PhD/WP2/MS2/averaged_MS2"

# File prefix differs by site -- extend this if you add more sample groups later
get_prefix <- function(sample) if (sample == "S1") "SESO34" else "PNK156"

## MS1: precursor info loading

load_ms1_precursors <- function(sample, target_mz,
                                tolerance = precursor_tol,
                                rel_cutoff = ms1_rel_intensity_cutoff,
                                ms1_dir) {
  prefix <- get_prefix(sample)
  path <- file.path(ms1_dir, paste0(prefix, "_", sample, "_MS1-qb.mzML"))
  
  msExp <- openMSfile(path)
  peaks_df <- as.data.frame(spectra(msExp)) %>%
    dplyr::filter(mz > target_mz - tolerance, mz < target_mz + tolerance)
  
  peaks_df %>%
    dplyr::filter(intensity > rel_cutoff * max(intensity)) %>%
    mutate(ion_type = "precursor", sample = sample, target_mz = target_mz) %>%
    dplyr::select(mz, intensity, ion_type, sample)
}

## MS2: raw peak loading + ion_type split

load_ms2_peaks <- function(sample, target_mz, tolerance = precursor_tol, ms2_dir) {
  prefix   <- get_prefix(sample)
  mz_label <- sprintf("%d", round(target_mz))
  path <- file.path(ms2_dir, paste0(prefix, "_", sample, "_", mz_label, "-qb.mzML"))
  
  msExp <- openMSfile(path)
  peaks_df <- as.data.frame(spectra(msExp)) %>%
    dplyr::filter(mz < target_mz + tolerance)
  
  precursor_range <- c(target_mz - tolerance, target_mz + tolerance)
  peaks_df$ion_type <- ifelse(peaks_df$mz >= precursor_range[1] & peaks_df$mz <= precursor_range[2],
                              "precursor", "fragment")
  peaks_df
}

## MS2: filtered fragments only 
##  Precursor rows from MS2 are dropped
## precursor uniqueness comes from MS1.

get_filtered_fragments <- function(sample, target_mz,
                                   base_peak_pct = fragment_base_peak_pct,
                                   ms2_dir) {
  peaks_df <- load_ms2_peaks(sample, target_mz, ms2_dir = ms2_dir)
  fragments <- peaks_df %>% dplyr::filter(ion_type == "fragment")
  frag_base_peak <- max(fragments$intensity)
  
  fragments %>%
    dplyr::filter(intensity > (base_peak_pct / 100) * frag_base_peak) %>%
    mutate(sample = sample) %>%
    dplyr::select(mz, intensity, ion_type, sample)
}

## Combine MS1 precursor + MS2 fragments, per sample 
combine_and_align <- function(filter_results_list, 
                              mz_tolerance_ppm = 0.6,
                              precursor_range = NULL) {
  
  # Combine all samples into one data frame
  all_samples <- bind_rows(filter_results_list, .id = "sample_id")
  
  # Ensure we have sample column
  if (!"sample" %in% colnames(all_samples)) {
    all_samples$sample <- all_samples$sample_id
  }
  
  # Add ion_type if not present
  if (!"ion_type" %in% colnames(all_samples) && !is.null(precursor_range)) {
    all_samples$ion_type <- ifelse(all_samples$mz >= precursor_range[1] & 
                                     all_samples$mz <= precursor_range[2], 
                                   "precursor", "fragment")
  }
  
  # Sort by m/z
  all_samples <- all_samples %>%
    arrange(mz) %>%
    dplyr::select(mz, intensity, sample, ion_type, everything())
  
  # Group m/z values within tolerance
  all_samples$mz_group <- 1
  current_group <- 1
  reference_mz <- all_samples$mz[1]
  
  for (i in 2:nrow(all_samples)) {
    ppm_diff <- abs(all_samples$mz[i] - reference_mz) / reference_mz * 1e6
    
    if (ppm_diff <= mz_tolerance_ppm) {
      all_samples$mz_group[i] <- current_group
    } else {
      current_group <- current_group + 1
      all_samples$mz_group[i] <- current_group
      reference_mz <- all_samples$mz[i]
    }
  }
  
  # Assign representative m/z
  aligned_data <- all_samples %>%
    group_by(mz_group) %>%
    mutate(aligned_mz = mean(mz)) %>%
    ungroup()
  
  # Summary of alignment
  cat("=== M/Z Alignment Summary ===\n")
  cat(sprintf("Total peaks across all samples: %d\n", nrow(all_samples)))
  cat(sprintf("Unique m/z groups after alignment: %d\n", 
              n_distinct(aligned_data$mz_group)))
  
  return(aligned_data)
}

build_sample_data <- function(sample, target_mz,
                              base_peak_pct = fragment_base_peak_pct,
                              ms1_dir, ms2_dir) {
  ms1_prec <- load_ms1_precursors(sample, target_mz, ms1_dir = ms1_dir)
  ms2_frag <- get_filtered_fragments(sample, target_mz, base_peak_pct, ms2_dir = ms2_dir)
  bind_rows(ms1_prec, ms2_frag)
}

## Run for one m/z target across all 4 samples

run_pipeline_for_mz <- function(target_mz,
                                base_peak_pct = fragment_base_peak_pct,
                                ms1_dir, ms2_dir,
                                mz_tolerance_ppm = 0.7) {
  precursor_range <- c(target_mz - precursor_tol, target_mz + precursor_tol)
  
  sample_list <- lapply(samples, build_sample_data,
                        target_mz = target_mz,
                        base_peak_pct = base_peak_pct,
                        ms1_dir = ms1_dir, ms2_dir = ms2_dir)
  names(sample_list) <- samples
  
  aligned <- combine_and_align(filter_results_list = sample_list,
                               mz_tolerance_ppm = mz_tolerance_ppm,
                               precursor_range = precursor_range)
  aligned$target_mz <- target_mz
  aligned
}

## Sharing / uniqueness stats 

compute_sharing_stats <- function(aligned_data) {
  sharing_stats <- aligned_data %>%
    group_by(mz_group, aligned_mz, ion_type) %>%
    summarise(
      n_samples = n_distinct(sample),
      samples = paste(unique(sample), collapse = ", "),
      mean_intensity = mean(intensity),
      sd_intensity = sd(intensity),
      cv_intensity = sd(intensity) / mean(intensity) * 100,
      total_intensity_across_samples = sum(intensity),
      .groups = "drop"
    ) %>%
    arrange(desc(n_samples), desc(total_intensity_across_samples))
  
  sample_tic <- aligned_data %>%
    group_by(sample) %>%
    summarise(sample_tic = sum(intensity))
  
  intensity_contribution <- aligned_data %>%
    left_join(sharing_stats %>% dplyr::select(mz_group, n_samples), by = "mz_group") %>%
    group_by(sample, n_samples, ion_type) %>%
    summarise(category_intensity = sum(intensity), n_peaks = n_distinct(mz_group), .groups = "drop") %>%
    left_join(sample_tic, by = "sample") %>%
    mutate(percent_of_sample_tic = category_intensity / sample_tic * 100)
  
  peak_count_contribution <- aligned_data %>%
    distinct(sample, mz_group, ion_type, .keep_all = FALSE) %>%
    left_join(sharing_stats %>% dplyr::select(mz_group, n_samples), by = "mz_group") %>%
    group_by(sample, n_samples, ion_type) %>%
    summarise(n_peaks = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(total_peaks_per_sample = sum(n_peaks),
           percent_of_sample_peaks = n_peaks / total_peaks_per_sample * 100) %>%
    ungroup()
  
  list(sharing_stats = sharing_stats,
       intensity_contribution = intensity_contribution,
       peak_count_contribution = peak_count_contribution)
}

##Plots

ion_type_labels <- as_labeller(c(fragment = "Fragments", precursor = "Precursors"))
grad_colors <- c("#00429E", "#5B86D7", "#A181DA", "#EF7EF7")

make_sharing_plots <- function(stats_list, title = NULL) {
  p_tic <- ggplot(stats_list$intensity_contribution,
                  aes(x = sample, y = percent_of_sample_tic, fill = factor(n_samples))) +
    geom_col(position = "stack") +
    facet_wrap(~ion_type, labeller = ion_type_labels) +
    labs(x = "Sample origin", y = "% of TIC per sample", title = title) +
    theme_minimal() +
    theme(strip.text.x = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = grad_colors, name = "In how many samples\nwas the same peak\ndetected")
  
  p_n <- ggplot(stats_list$peak_count_contribution,
                aes(x = sample, y = n_peaks, fill = factor(n_samples))) +
    geom_col(position = "stack") +
    facet_wrap(~ion_type, scales = "free_y", labeller = ion_type_labels) +
    labs(x = "Sample origin", y = "Number of peaks") +
    theme_minimal() +
    theme(strip.text.x = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = grad_colors, name = "In how many samples\nwas the same peak\ndetected")
  
  (p_tic | p_n) + plot_layout(guides = "collect")
}

## Spectral similarity (cosine)

cosine_similarity <- function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))

compute_cosine_matrix <- function(aligned_data, ion_types) {
  sim_data <- aligned_data %>%
    dplyr::filter(ion_type %in% ion_types) %>%
    group_by(sample, mz_group) %>%
    summarise(intensity = sum(intensity), .groups = "drop")
  
  intensity_matrix <- sim_data %>%
    pivot_wider(names_from = mz_group, values_from = intensity, values_fill = 0)
  
  sample_names <- intensity_matrix$sample
  intensity_data <- as.matrix(intensity_matrix[, -1])
  rownames(intensity_data) <- sample_names
  
  n <- length(sample_names)
  cosine_mat <- matrix(1, n, n, dimnames = list(sample_names, sample_names))
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      cosine_mat[i, j] <- cosine_mat[j, i] <- cosine_similarity(intensity_data[i, ], intensity_data[j, ])
    }
  }
  round(cosine_mat, 3)
}

compute_similarity <- function(aligned_data) {
  list(
    fragments_only      = compute_cosine_matrix(aligned_data, "fragment"),
    fragments_precursor = compute_cosine_matrix(aligned_data, c("fragment", "precursor"))
  )
}


## MAIN EXECUTION

mz_aligned_343 <- run_pipeline_for_mz(343.0, ms1_dir = ms1_dir, ms2_dir = ms2_dir, mz_tolerance_ppm = mz_tolerance_ppm)
mz_aligned_383 <- run_pipeline_for_mz(383.0, ms1_dir = ms1_dir, ms2_dir = ms2_dir, mz_tolerance_ppm = mz_tolerance_ppm)
mz_aligned_467 <- run_pipeline_for_mz(467.0, ms1_dir = ms1_dir, ms2_dir = ms2_dir, mz_tolerance_ppm = mz_tolerance_ppm)

sharing_343 <- compute_sharing_stats(mz_aligned_343)
sharing_383 <- compute_sharing_stats(mz_aligned_383)
sharing_467 <- compute_sharing_stats(mz_aligned_467)

sharing_343_st <- sharing_343[[2]]
sharing_383_st <- sharing_383[[2]]
sharing_467_st <- sharing_467[[2]]



plot_343 <- make_sharing_plots(sharing_343)
plot_383 <- make_sharing_plots(sharing_383, title = "m/z 383")
plot_467 <- make_sharing_plots(sharing_467, title = "m/z 467")
## Supplementary Figure 7:
(wrap_elements(plot_383) / wrap_elements(plot_467) ) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")

sim_343 <- compute_similarity(mz_aligned_343)
sim_383 <- compute_similarity(mz_aligned_383)
sim_467 <- compute_similarity(mz_aligned_467)


# plot_383
# sim_383$fragments_only
# sim_383$fragments_precursor
# sharing_383$sharing_stats