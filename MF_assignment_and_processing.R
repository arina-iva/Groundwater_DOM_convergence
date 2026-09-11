library(devtools)
install_github('skschum/MFAssignR/MFAssignR')

library(lubridate)
library(MALDIquant)
library(MALDIquantForeign)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MFAssignR)
library(stringr)
library(future.apply)
library(parallel)
library(progressr)
library(fs)
library(vegan)
library(labdsv)
library(ggrepel)
library(tibble)
library(ggpubr)
library(purrr)
library(multcompView)
library(patchwork)
library(ggrepel)


#NOTE: This script uses averaged spectra.
#NECESSARY PREPROCESSING:
#1) average ~100 scans in Thermo QualBrowser.
#2) export as a raw file.
#3) do so for all samples.
#4) use Thermo RawFileParser to do peak picking and export the results as mzML files.

# Set your input file paths here
filePaths <- list.files("C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/mzml/SESO_HAI/", full.names = T, pattern = '.mzML$')

#Set your output file path here
outputPath <- "C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/mzml/SESO_HAI/output/"

# PARALLELIZATION SETTINGS
parallel <-TRUE #would you like to use parallel processing?

# GENERAL SETTINGS####
ionizationMode <- "neg"  #'neg' or 'pos'
LowMassLimit <- 100  #What is the lower limit of the mass range in your scans?
HighMassLimit <- 1000  #What is the upper limit of the mass range in your scans?
CHO_assign_error <- 2.5 #allowed error range (in ppm) for CHO sum formula assignment. Needs to be kind of high, because the data is not recalibrated at this point.
assignment_errorrange <- 0.6 #error range (in ppm) for final sum formula assignment
SN_import <- 7 #Signal to noise ratio during import
SN_recal <- 15  #Signal to noise ratio during import.  Low intensities commonly lead to high recalibration errors(>1ppm). Lower this value if you have trouble finding recalibrants.
RecalRangeCheck <- TRUE #require that the recalibration covers a specified mass range (set below)
recalRange <- c(150,450) #Sample must be recalibrated within this entire range, else the sample is reported as failed recalibration.
AutoUpperLimit <- FALSE #should the upper limit of recalRange be set automatically to the m/z where a set percentage of the cumulative intensity is reached?
AutoUpperLimitTarget <- 0.7
Nitrogen_recal <- 1 #Number of nitrogen atoms allowed during homologous series search. Set to zero for CHO only.
max_loop <- 10 # maximum number of retries during recalibration

##ASSIGNMENT SETTINGS
# ELEMENT SETTINGS####
set_Nx <- 5 #10
set_Sx <-  2 #4
set_Px <- 1 #4
set_S34x <- 0
set_N15x <- 0
set_Dx <- 0
set_Ex <- 0 
set_Clx <- 0 #2 
set_Cl37x <- 0 
set_Fx <- 0
set_Brx <- 0 #2 
set_Br81x <- 0 
set_Ix <- 0 #2
set_Mx <- 0
set_NH4x <- 0
set_Zx <- 1 
set_Sval <- 2
set_Nval <- 3
set_S34val <- 2
set_N15val <- 3
set_Pval <- 5
set_Ox <- 30
set_O_Cmin <- 0.1
set_O_Cmax <- 2.3
set_H_Cmin <- 0.3
set_H_Cmax <- 3
set_DBEOmin <- -13
set_DBEOmax <- 13
set_Omin <- 1

# DETAIL SETTINGS####
set_POEx <- 0
set_NOEx <- 0
set_max_def <- 0.9
set_min_def <- 0.5
set_HetCut <- "off"
set_NMScut <- "on"
set_DeNovo <- 300
set_nLoop <- 10
set_SulfCheck <- "on"
set_Ambig <- "off"
set_MSMS <- "off"
set_S34_abund <- 30
set_C13_abund <- 60
set_N3corr <- "on"



if (parallel) {
  plan(multisession) #use plan(multicore) under linux
} else {
  plan(sequential)
}

handlers(handler_progress(format="[:bar] :percent :eta :message")) #formatting of the progress bar

# IMPORT
{
  time1 <- Sys.time()
  with_progress({
    p <- progressor(along = 1:length(filePaths))
    import <- future_lapply(1:length(filePaths), FUN = function(x) {
      p("importing spectra...")
      
      # load an averaged and peak picked mzml file
      mzint <- MALDIquantForeign::importMzMl(filePaths[x],
                                             centroided = TRUE)[[1]] %>%
        MALDIquant::as.matrix() %>%
        as.data.frame()
      noise_sample_specific <- MFAssignR::KMDNoise(mzint)$Noise
      #detect (mono)isotope peaks
      MonoIso <- try(MFAssignR::IsoFiltR(mzint,
                                         SN = noise_sample_specific  * SN_import),
                     silent = TRUE)
      
      # if no isotope peaks are detected, we set them to
      # "none" and continue without them   
      peakList <- list()
      if (class(MonoIso) == "try-error"){
        MonoIso <- list()
        MonoIso$Mono <- mzint
        MonoIso$Iso <- "none" 
      }
      
      tmp <- list(MonoIso, noise_sample_specific, path_ext_remove(basename(filePaths[x])))
      return(tmp)
    }#,
    #future.chunk.size = 5  #during parallel processing, this controls how many tasks(samples) are sent to each worker at a time
    )
  })
  time2 <- Sys.time()
  difftime(time2,time1)
}

#SAVE YOUR OUTPUT

save(import, file = paste0(outputPath,"import_","SN_",SN_import,"_", ionizationMode, ".Rdata"))

#IF CONTINUING PREVIOUS ANALYSES, LOAD PREVIOUS OUTPUT
#load("YOUR INPUT FILE PATH HERE")

#   RECALIBRATION
{
  time1 <- Sys.time()
  with_progress({
    p <- progressor(along = import)
    recal_result <- future_lapply(import, FUN = function(x) {
      p("recalibrating...")
      # Assign CHO sum formulas. This is needed as an
      # input for the identification of the recalibration
      # series.
      
      assignCHO <- try(MFAssign_RMD(peaks = x[[1]]$Mono, 
                                    isopeaks = x[[1]]$Iso, ionMode = ionizationMode, 
                                    lowMW = LowMassLimit, highMW = HighMassLimit, 
                                    ppm_err = CHO_assign_error,
                                    SN = SN_recal * x[[2]],
                                    Nx  = Nitrogen_recal
      )$Unambig,
      silent = FALSE)  
      df_recal <- try(RecalList(assignCHO),
                      silent = FALSE)
      
      #check if any homologous series were found
      if (class(df_recal) != "try-error" &&
          nrow(df_recal) > 0) {
        
        #extract their mass ranges (used later for ordering the series)
        
        df_recal2 <- df_recal %>%
          mutate(min_mass = as.numeric(str_extract(df_recal$`Mass Range`, pattern = "^[[:digit:]]+\\.[[:digit:]]+")),
                 max_mass = as.numeric(str_extract(df_recal$`Mass Range`, pattern = "[[:digit:]]+\\.[[:digit:]]+$"))) %>%
          arrange(min_mass)
        
        recal_series_initial <- c() #initialize
        start_mass_series_from <- c(0,25,50,100,150,200,250,300)
        start_mass_series_to <- c(25,50,100,150,200,250,300,400)
        for (i in 1:10){    #select series section-wise over the mass range, ca 50Da per section
          recal_series_initial[i] <- ifelse(any(between(df_recal2$min_mass, LowMassLimit + start_mass_series_from[i], LowMassLimit + start_mass_series_to[i])),
                                            df_recal2 %>%
                                              filter(between(df_recal2$`Tall Peak`, LowMassLimit + start_mass_series_from[i], LowMassLimit + start_mass_series_to[i])) %>%
                                              slice_min(order_by = `Series Index`, n = 1) %>%
                                              select("Series") %>%
                                              unlist(),
                                            NA)
        }
        recal_series <- recal_series_initial[!is.na(recal_series_initial)] 
        while (length(recal_series) < 10 &
               length(recal_series) < nrow({df_recal2 %>% filter(min_mass < 400)})){
          recal_series[length(recal_series) + 1] <- df_recal2 %>%
            filter(!Series %in% recal_series) %>%
            slice_min(order_by = `Series Index`, n = 1) %>%
            pull(Series)
        }
        
        
        #initialize an empty error object for the recalibration
        recal_output <- "placeholder for object initialization"
        class(recal_output) <- "try-error"
        
        #HERE THE ACTUAL RECALIBRATION STARTS
        #we go through all combinations of the top10 homologous series from df_recal
        #we start with all ten, then reduce to 9, 8, 7...
        #on each level we go through all possible combinations starting with low-mass combinations first
        
        
        #find first stable recalibration solution
        i <- length(recal_series)
        while (i >= 1){
          recal_series_combinations <- combn(recal_series, m = i)
          k <- 1
          while (k <= ncol(recal_series_combinations)) {
            recal_output <- try(
              Recal(
                df = assignCHO,
                peaks = x[[1]]$Mono,
                isopeaks = x[[1]]$Iso,
                mode = ionizationMode,
                series1 = ifelse(length(recal_series_combinations[,1])>=1,recal_series_combinations[1,k],NA),
                series2 = ifelse(length(recal_series_combinations[,1])>=2,recal_series_combinations[2,k],NA),
                series3 = ifelse(length(recal_series_combinations[,1])>=3,recal_series_combinations[3,k],NA),
                series4 = ifelse(length(recal_series_combinations[,1])>=4,recal_series_combinations[4,k],NA),
                series5 = ifelse(length(recal_series_combinations[,1])>=5,recal_series_combinations[5,k],NA),
                series6 = ifelse(length(recal_series_combinations[,1])>=6,recal_series_combinations[6,k],NA),
                series7 = ifelse(length(recal_series_combinations[,1])>=7,recal_series_combinations[7,k],NA),
                series8 = ifelse(length(recal_series_combinations[,1])>=8,recal_series_combinations[8,k],NA),
                series9 = ifelse(length(recal_series_combinations[,1])>=9,recal_series_combinations[9,k],NA),
                series10 = ifelse(length(recal_series_combinations[,1])>=10,recal_series_combinations[10,k],NA),
                SN = SN_recal * x[[2]]
              )[2:4],
              silent = TRUE)
            if (class(recal_output) == "list"){ if(str_detect(x[[3]],"[[B,b]]lank")) break
              if (RecalRangeCheck & (min(recal_output$RecalList$exp_mass) <= {recalRange[1]}) & (max(recal_output$RecalList$exp_mass) >= {
                ifelse(AutoUpperLimit,
                       x[[1]]$Mono %>%
                       mutate(cumsum = cumsum(abundance),
                              sum = sum(abundance)) %>%
                       filter(cumsum <= AutoUpperLimitTarget * sum(abundance)) %>%
                       pull(exp_mass) %>%
                       max(),
                       recalRange[2])
              })) break else {
                k <- k+1
                class(recal_output) <- "try-error"
              }
            } else {
              k <- k+1
            }
            if (k >= max_loop) break #maximum number of retries
          }
          if (class(recal_output) == "list") break
          if (i == 1) break
          i <- i-1
        }
        
        if (class(recal_output) == "list") {
          howdiditgo <- paste("Successfully recalibrated.")
        }
        
        if (class(recal_output) == "try-error") {
          howdiditgo <- paste("No stable recalibration solution in range.")
        }
      }
      
      if (class(df_recal) == "try-error") {
        howdiditgo <- paste("The recalibration step was skipped because no homologous series were found.")
        recal_output <- list(Mono = NA, Iso = NA, RecalList = NA)
      }
      
      if (class(recal_output) == "try-error"){
        recal_output <- "try-error" }
      
      tmp <- list(recal_output,
                  x[[2]],
                  howdiditgo,
                  x[[3]])
      return(tmp)
    }#,
    #future.chunk.size = 5  #during parallel processing, this controls how many tasks(samples) are sent to each worker at a time
    )
  })
  time2 <- Sys.time()
  difftime(time2,time1)
}

if (parallel) {
  plan(sequential) #if parallel processing was used, close worker instances.
}

save(recal_result, file = paste0(outputPath,"recal_result_withpotentialerrors_Jun25", ionizationMode, ".Rdata"))

#IF CONTINUING FROM PREVIOUS ATTEMPT, LOAD YOUR RECAL RESULT HERE
load("C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/mzml/SESO_HAI/output/recal_result_withpotentialerrors_Jun25neg.Rdata")

#CHECK IF ALL SAMPLES COULD BE RECALIBRATED
recal_results_report <- c(1:length(recal_result))
for (i in 1:length(recal_results_report)){
  recal_results_report[i] = paste0("Sample ",i, " ", recal_result[[i]][4],  ": ", recal_result[[i]][3]) 
}
recal_results_report <- as.data.frame(recal_results_report)
write.csv(recal_results_report,  file = "C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/mzml/SESO_HAI/recal_result_report.csv")
##save recal results
#CHECK THE NUMBER OF MASSES AT YOUR SELECTED NOISELEVEL

for (i in 1:length(import)){
  print(paste("Sample",i, import[[i]][[3]],"contains", import[[i]][[1]]$Mono %>% nrow(), "monoisotopic masses above SN =", SN_import,"and noiselevel =", round(import[[i]][[2]])))
  print(paste("Sample",i,":",recal_result[[i]][[3]]))
}

#WERE ANY UNSUCCESSFUL?
ifelse(any(sapply(recal_result, function(x){x[[3]]}) != "Successfully recalibrated."),
       "Some samples contain recalibration errors, please check and remove before continuing.",
       "All samples okay to conitnue.")

#exclude any samples that threw errors during recalibration.
recal_result <- recal_result[which(sapply(recal_result, function(x){x[[3]]}) == "Successfully recalibrated.")]

#check width of recalibrated mass range
minmz <- data.frame(sampleName = character(),
                    sampleNr = numeric(),
                    min_recal_mz = numeric(),
                    max_recal_mz = numeric())

for (i in 1:length(recal_result)){
  print(paste0("SampleNr", i,": recalibrated from ",recal_result[[i]][[1]]$RecalList$exp_mass %>% min() %>% round()," to ",
               recal_result[[i]][[1]]$RecalList$exp_mass %>% max() %>% round()," m/z"))
  minmz[i,1] <- recal_result[[i]][[4]]
  minmz[i,2] <- i
  minmz[i,3] <- recal_result[[i]][[1]]$RecalList$exp_mass %>% min()
  minmz[i,4] <- recal_result[[i]][[1]]$RecalList$exp_mass %>% max()}

#OPTIONAL MANUAL CHECKS
#DETAILED QUALITY ASSESSMENT WITH RECALIBRATION MASS ERROR PLOTS
#resulting mass_error should be less than 1ppm
sample_nr <-102 #select sample of interest (change this manually)
ggplot(data.frame(mass = recal_result[[sample_nr]][[1]]$RecalList$exp_mass,
                  intensity = recal_result[[sample_nr]][[1]]$RecalList$abundance), aes(x=mass, y = intensity))+
  geom_col(color="black")+
  ggtitle(recal_result[[sample_nr]][[4]])+
  scale_x_continuous(limits = c(100,1000),
                     breaks = seq(100,1000, by = 50))+
  theme(axis.title.x = element_blank(),
        panel.grid.minor = element_blank())+
  ggplot(data.frame(mass_error = recal_result[[sample_nr]][[1]]$RecalList$recal_err,
                    mass = recal_result[[sample_nr]][[1]]$RecalList$exp_mass), aes(x=mass, y = mass_error))+
  geom_point(shape = 1)+
  geom_smooth(color="red")+
  scale_x_continuous(limits = c(100,1000),
                     breaks = seq(100,1000, by = 50))+
  theme(panel.grid.minor = element_blank())+
  patchwork::plot_layout(ncol = 1)

#SAVE YOUR OUTPUT AS .RDATA FILE
save(recal_result,file = paste0(outputPath,"recal_result_clean_Jun25", ionizationMode, ".RData"))
load("C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/mzml/SESO_HAI/outputrecal_result_cleanneg.RData")

#JUNCTION ALGORITHM ADAPTED FROM https://github.com/JulianMerder/ICBM-OCEAN--Supplementary-Code-Modules/blob/master/fastjoin.R
junction_input <-  bind_rows(
  lapply(recal_result, function(x){
    colnames(x[[1]]$Mono)[1:2] <- c("mz", "I")
    colnames(x[[1]]$Iso)[1:2] <- c("mz", "I")
    rbind(
      x[[1]]$Mono[1:2],
      x[[1]]$Iso[1:2]) %>%
      as.data.frame()%>%
      arrange(`mz`) %>%
      mutate(`MDL_2.5` = 10,
             `ResPow` = 480000,
             index = fs::path_ext_remove(basename(x[[4]])))}))

tolp<-0.5 #your tolerance (usually no changes are needed here. should be less than assignment_errorrange)
mdlname<-names(junction_input)[3]
names(junction_input)[3]<-"MDL"

junction_output <- dplyr::arrange(junction_input,mz) %>%
  # Start with the smallest mass and check if the next mass is the tolerance range defined by 'tolp'. If Yes group togetherinto a mass cluster, if not start a new group (mass cluster). 
  dplyr::mutate(mz1=mz,group = cumsum( abs(((mz - dplyr::lag(mz, default = mz[1], order_by = mz) )/dplyr::lag(mz, default = mz[1], order_by = mz))*10^6) > tolp)) %>%
  # Group by found mass cluster and check for duplicate samples within a mass cluster
  dplyr::group_by(group) %>%
  dplyr::mutate(max=abs(mz[which.max(I)]-mz)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group, index) %>% dplyr::arrange(max)%>% dplyr::mutate(gjk=as.character(ifelse(dplyr::row_number()==1,paste0(group,"x",dplyr::row_number()),paste0(group,"x","x",dplyr::row_number()))))%>%            #dplyr::summarize(mz=(mz[which.min(max)]),I=(I[which.min(max)]),MDL=(MDL[which.min(max)]),ResPow=(ResPow[which.min(max)]),mz1=mz1[which.min(max)],refe=refe[which.min(max)],m1=m1[which.min(max)])%>%
  dplyr::ungroup() %>% 
  dplyr::select(c(-max,-group))%>%
  dplyr::group_by(gjk) %>%
  dplyr::mutate(
    # calculate weighted average mass within each mass group or use arithmethic mean with mean()
    mz = weighted.mean(mz, sqrt(I), na.rm=T),MDL=mean(MDL), ResPow=mean(ResPow),SE=sd(mz1)/mz*1e6,mz1=mean(mz1)) %>%  
  # Take out the  group assignment and columns not needed
  dplyr::ungroup() %>%
  dplyr::select(c(-mz1,-gjk)) %>%
  # Add the "B" leader to the index
  dplyr::mutate(index = paste0("Sample", index)) %>%
  #intensity filter
  #dplyr::filter(I > (SN_import * noise_abs_intensity)) %>%
  #remove singlets
  dplyr::group_by(mz) %>%
  dplyr::mutate(mzCount = n()) %>%
  dplyr::filter(mzCount > 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-mzCount) %>%
  # Put intensity values of each sample into separate columns
  tidyr::pivot_wider(names_from = index, values_from = I, values_fn = max) %>%
  #remove duplicate masses (weird error)
  dplyr::filter(!duplicated(mz)) %>%
  dplyr::arrange(mz)

save(junction_output,file = paste0(outputPath,"junction_output_", ionizationMode, "Jun25.RData"))
load("C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/mzml/SESO_HAI/output/junction_output_negJun25.RData")

#memory cleanup
gc()

#SUM FORMULA ASSIGNMENT 
## HAVE TO BE ADJUSTED TO HOW THE AVERAGED DATA LOOKS LIKE

#identify isotope peaks with max measured intensities
#for large data sets, this step can cause memory issues. 
#if this happens, restart your pc and continue with only the minimal number of objects and packages loaded in the environment
  
assignMonoIso <- MFAssignR::IsoFiltR(peaks = data.frame(mz=junction_output$mz,
                                                        abundance = junction_output %>%
                                                          select(starts_with(c("mz","Sample"))) %>%
                                                          tidyr::pivot_longer(cols = starts_with("Sample"), values_drop_na = TRUE) %>%
                                                          group_by(mz) %>% summarize(abundance_max = max(value)) %>% pull(abundance_max)))

#assignMonoIso <- MFAssignR::IsoFiltR(peaks = data.frame(mz=averaged$mz,
#                                                        abundance = averaged$Intensity))
#save(assignMonoIso,file = paste0(outputPath,"assignmonoiso_", ionizationMode, ".RData"))

assigned <- MFAssign(peaks = assignMonoIso$Mono, isopeaks = assignMonoIso$Iso, 
                     ionMode = ionizationMode, lowMW = LowMassLimit, 
                     highMW = HighMassLimit, 
                     ppm_err = assignment_errorrange,
                     POEx = set_POEx,
                     NOEx = set_NOEx,
                     Nx = set_Nx,
                     Sx = set_Sx,
                     Px = set_Px,
                     S34x = set_S34x,
                     N15x = set_N15x,
                     Dx = set_Dx,
                     Ex = set_Ex,
                     Clx = set_Clx,
                     Cl37x = set_Cl37x,
                     Fx = set_Fx,
                     Brx = set_Brx,
                     Br81x = set_Br81x,
                     Ix = set_Ix,
                     Mx = set_Mx,
                     NH4x = set_NH4x,
                     Zx = set_Zx,
                     Sval = set_Sval,
                     Nval = set_Nval,
                     S34val = set_S34val,
                     N15val = set_N15val,
                     Pval = set_Pval,
                     Ox = set_Ox,
                     O_Cmin = set_O_Cmin,
                     O_Cmax = set_O_Cmax,
                     H_Cmin = set_H_Cmin,
                     H_Cmax = set_H_Cmax,
                     DBEOmin = set_DBEOmin,
                     DBEOmax = set_DBEOmax,
                     Omin = set_Omin,
                     max_def = set_max_def,
                     min_def = set_min_def,
                     HetCut = set_HetCut,
                     NMScut = set_NMScut,
                     DeNovo = set_DeNovo,
                     nLoop = set_nLoop,
                     SulfCheck = set_SulfCheck,
                     Ambig = set_Ambig,
                     MSMS = set_MSMS,
                     S34_abund = set_S34_abund,
                     C13_abund = set_C13_abund,
                     N3corr = set_N3corr)

save(assigned,file = paste0(outputPath,"/assigned_Jun25", ionizationMode, ".RData"))
load("C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/mzml/SESO_HAI/output/assigned_Jun25neg.RData")

#how many sum formulas were assigned?
nrow(assigned$Unambig)
nrow(assigned$Ambig)
nrow(assigned$None)


#join the sum formulas with the intensities
join <- inner_join(x=assigned$Unambig, y=junction_output, by = c("exp_mass" = "mz"))

join_assigned <- join%>%
  select(starts_with("Sample")) %>%
  summarise(across(everything(), ~sum(!is.na(.)))) %>%
  tidyr::pivot_longer(everything(), 
                      names_to = "Sample", 
                      values_to = "assigned_mz") %>%
  filter(!str_detect(Sample, "[B,b]lank"),
         !str_detect(Sample, "reference"))%>%
  mutate(wells = str_remove_all(str_extract(Sample, 
                                            pattern = "H[12345][:punct:]?[1234][:punct:]{0,1}|S[1,2]_?[a,b]|[OBGP][:digit:]-[:digit:]{2}"),
                                pattern = "[[:punct:]]|[[a,b]]$"))%>%
  group_by(wells)%>%
  summarise(
    mean_assigned = mean(assigned_mz)
  )

QC_sampleNames <- colnames(join)[str_detect(colnames(join),"^Sample")] #select a sample name
#how much of the original intensity remains after sum formula assignment?
MF_assignment_report <- as.data.frame(matrix(nrow = 0, ncol = 0))

for (i in 1: length(QC_sampleNames)){
  MF_assignment_report[i, 1] = paste0(QC_sampleNames[i], " post-assignment intensity ratio: ")
  MF_assignment_report[i, 2] = as.numeric(
    round(sum(join %>% pull(QC_sampleNames[i]), na.rm = T) / sum(junction_output %>% pull(QC_sampleNames[i]), na.rm = T) * 100))
}


  sum(str_count(colnames(join), pattern = fixed("Sample"))) #391 Samples in join
  sum(str_count(colnames(junction_output), pattern = fixed("Sample"))) #392 samples in junction output... but no SampleH41_A_PNK137 present
y <- cbind(sort(colnames(join)[str_detect(colnames(join),"^Sample")]), 
           sort(colnames(junction_output)[str_detect(colnames(junction_output),"^Sample")]))

MF_assignment_report%>%dplyr::filter(grepl("SESO",V1))%>%pull(V2) %>%mean()

##average assigned intensity ratio is:
filtered_MF_assignment_report <- MF_assignment_report[!grepl("Blank", MF_assignment_report$V1), ]
AAIR <- mean(filtered_MF_assignment_report$V2)
sdAAIR <- sd(filtered_MF_assignment_report$V2)

##remove samples with assigned intensity ratio below 60%
low_AIR_samples <- filtered_MF_assignment_report[filtered_MF_assignment_report$V2<50, ]%>%
  mutate(samples = str_split(string = V1, simplify = TRUE, 
                             pattern = " ")[, 1])

##new average assigned intensity ratio after removing low assigned samples
filtered_MF_assignment_report <- MF_assignment_report[!grepl("Blank", MF_assignment_report$V1), ]
 join_new <- join_new[,!names(join_new) %in% low_AIR_samples$samples]


#in rare cases different masses can match to the same sum formula. let's remove them
if(any(duplicated(join$formula))){
  join <-
    join %>%
    distinct(formula, .keep_all = TRUE)
}

#save final output
save(join,file = paste0(outputPath,"assigned_join_Jun25_", ionizationMode, ".RData"))

#load joined data frame
load("C:/Users/aivanova/Documents/Orbitrap Data/aquadiva data 2024/mzml/SESO_HAI/output/assigned_join_Jun25_neg.RData")

### blank removal
blankformulas <-join %>%
  dplyr::select(formula,contains("blank")) %>%
  tidyr::pivot_longer(-formula) %>%
  filter(!is.na(value)) %>%
  arrange(desc(value)) %>%
  mutate(cumInt = cumsum(value),
         relcumSum = cumInt / sum(cumInt)) %>%
  filter(relcumSum < 0.95) %>%
  pull(formula) %>%
  unique()

join_noblank <- join[!join$formula %in% blankformulas, ]
join_noblank <- dplyr::select(join_noblank, !contains("Blank"))

join_for_PCA <- join_noblank%>%
  pivot_longer(cols = starts_with("Sample"),
                              values_to = "Intensity",
                              names_to = "Sample",
                              values_drop_na = TRUE) %>%
                 filter(Intensity > 0,
                        !str_detect(Sample, "[B,b]lank"))%>%
                 rename(mz = theor_mass) %>%
                 mutate(wells = str_remove_all(str_extract(Sample,
                                               pattern = "H[12345][:punct:]?[1234][:punct:]{0,1}|S[1,2]_?[a,b]|[OBGP][:digit:]-[:digit:]{2}"),
                                               pattern = "[[:punct:]]|[[a,b]]$"),
                        PNK = as.numeric(str_extract(str_extract(Sample, pattern = "PNK[[:digit:]]{2,3}|SESO[[:digit:]]{1,2}"), "\\d+")),
                        site = case_when(grepl("^H\\d{2}$", wells) ~ "Hainich",
                                         grepl("^S\\d{1}$", wells) ~ "SESO",
                                         TRUE ~ NA_character_))

MF_with_ref_hel <- join_noblank %>% 
  dplyr::select(!starts_with("NA")) %>% 
  dplyr::select(starts_with("Sample")) %>%
  replace(list = is.na(.), values = 0) %>%
  filter(rowSums(.)!=0) %>%
  decostand(method = "total",
            MARGIN = 2) %>%
  t()%>%
  hellinger()

PCA <- prcomp(x = MF_with_ref_hel,
              center = TRUE,
              scale. = TRUE)
summary(PCA)

PCA_variance <- PCA$sdev^2 / sum(PCA$sdev^2)
PCA_variance_df <- as.data.frame(PCA_variance)

PCA_scores <- PCA$x %>%
  as.data.frame() %>%
  mutate(sample_name = rownames(.),
         samplecode = sub(x = sample_name, 
                          pattern = "Sample", , replacement = ""),
         well = ifelse(str_detect(samplecode, "reference"), "reference", str_extract(samplecode, "^[^_]+")),
         PNK = ifelse(str_detect(samplecode, "reference"), NA, str_extract(samplecode, "(?<=PNK)\\d+|(?<=SESO)\\d+")),
         PNK = as.numeric(PNK))%>%
  left_join(dataref_new, by = "PNK")%>%
  mutate(site = case_when(grepl("^H\\d{2}$", well) ~ "Hainich",
                          grepl("^S\\d{1}$", well) ~ "SESO",
  ))

PCA_scores <- PCA$x %>%
  as.data.frame() %>%
  mutate(sample_name = rownames(.),
         samplecode = sub(x = sample_name, 
                          pattern = "Sample", , replacement = ""),
         samplecode = str_remove(samplecode, "-qb$"),
         well_raw = case_when(
           str_detect(samplecode, "PNK146_H14B") ~ "reference_A",
           str_detect(samplecode, "reference") ~ "reference_B",
           TRUE ~ str_extract(samplecode, "[SH]\\d+[A-Za-z]?")),
         well_type = case_when(
           well_raw %in% c("reference_A", "reference_B") ~ "reference",
           TRUE ~ "groundwater sample"))%>%
  dplyr::filter(!well_raw%in%c("H13", "H42", "H31") & !is.na(well_raw))

colors_type <- c( "H14" = "#871818", "H32" = "#cd5c5c", "H41" = "#ff6833", "H43" = "#1eb879",
                  "H51" = "#f5e23b", "H52" = "#2395d6", "H53" = "#233ed6", "S1" = "#d678dc", "S2" = "#da2f94", 
                  "reference_A" = "black", "reference_B" = "darkgray")

PC1_2 <- ggplot(PCA_scores, aes(x = PC1, y = PC2))+
  geom_point(aes(color = well_raw, shape = well_type, alpha = well_type), size = 3)+
  theme_minimal()+
  scale_color_manual(values = colors_type)+
  scale_alpha_manual(values = c(0.6, 1))

PC2_3 <- ggplot(PCA_scores, aes(x = PC2, y = PC3))+
  geom_point(aes(color = well_raw, shape = well_type, alpha = well_type), size = 3)+
  theme_minimal()+
  scale_color_manual(values = colors_type)+
  scale_alpha_manual(values = c(0.6, 1))

##supplementary figure 1:
(PC1_2 | PC2_3) + plot_layout(guides = "collect")

averaged_MF <- join_for_PCA %>%
  filter(!wells %in% c("H13","H31","H42"))%>%
  group_by(wells,formula, PNK) %>%
  add_count(name = "n_repl") %>%
  filter(!(wells %in% c("H14","H32","H41","H43","H51","H52","H53","S1","S2") & n_repl != 2)) %>%
  mutate(mean_int = mean(Intensity)) %>%
  ungroup() %>%
  select(-Intensity)%>%
  rename(Intensity = mean_int) %>%
  distinct(formula, wells, 
           PNK, .keep_all = TRUE) %>%
  group_by(wells,
           PNK) %>%
  mutate(sum_int = sum(Intensity)) %>%
  ungroup() %>%
  mutate(relInt = Intensity / sum_int,
         Sample_name = paste0(PNK, "_", wells)) %>%
  select(-c(Intensity,sum_int)) %>%
  rename(Intensity = relInt)

save(averaged_MF, file = "C:/Users/aivanova/Nextcloud/PhD/WP2/averaged_MF_Jun25.RData")
load("C:/Users/aivanova/Nextcloud/PhD/WP2/averaged_MF_Jun25.RData")
dataref_new <- read.csv("C:/Users/aivanova/Nextcloud/PhD/WP2/dateref_new.csv")

filtered_MF <- averaged_MF %>%
  filter(!is.na(Intensity) & Intensity > 0) %>%               
  group_by(formula) %>%
  summarise(sample_count = n_distinct(Sample_name)) %>%
  filter(sample_count >= 2) %>% #remove singlets
  inner_join(averaged_MF, by = "formula") %>%  
  group_by(Sample_name) %>%
  mutate(sum_int = sum(Intensity)) %>%
  ungroup() %>%
  mutate(relInt = Intensity / sum_int) %>%
  select(-c(Intensity,sum_int)) %>%
  rename(Intensity = relInt)

save(filtered_MF, file = "C:/Users/aivanova/Nextcloud/PhD/WP2/averaged_MF_no_singlets.RData")
