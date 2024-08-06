## Clear environment
rm(list = ls())
gc()
getwd()


EXP<- "40x magnified data"
FIRST_TIMEID <- 1
TIME_STEP <- 1.67 # Time interval between measurements in hours

library(tidyverse)
library(readr)
library(pacman)
p_load(multidplyr)
p_load(pheatmap)
p_load(ggplot2)
p_load(grid)
p_load("viridis")
p_load(DBI)
p_load(tidyverse)
p_load(ggnewscale)
p_load(readr)
install.packages("ggpubr")
library(ggpubr)
install.packages("zoo")
library(zoo)
library(dplyr)
install.packages("processx")
library(processx)
install.packages("gtable")
library("gtable")

# FUNCTIONS

make_full_tracks <- function(data, ftid){
  ### Fill timepoint t1
  tracks_t1 <- data %>% filter(timeID == ftid) %>% pull(alt_uid)
  df <- tibble(TrackName = tracks_t1,
               !!paste0("t",ftid) := tracks_t1)
  tracks_prev <- tracks_t1
  prev_col_name <- paste0("t",ftid)
  
  ### Fill following timepoints
  # Get the number of timepoints
  tps <- unique(data %>% pull(timeID))
  for (i in tps[2:length(tps)]) {
    tracks_next <- data %>% filter(timeID == i) %>% pull(alt_uid)
    next_col_name <- paste0("t",i)
    existing <- tracks_next[which(tracks_next %in% tracks_prev)]
    branch <- tracks_next[which(!tracks_next %in% tracks_prev & grepl("\\.",tracks_next))]
    new <- tracks_next[which(!tracks_next %in% tracks_prev & !grepl("\\.",tracks_next))]
    
    # Check whether all tracks are used and unique
    if (!all(sort(c(new,existing,branch)) == sort(tracks_next))) {warning("Not all new track points are taken along or some are double.")}
    
    # Make tibbles from existing, branch and new tracks
    existing <- tibble(!!prev_col_name := existing,
                       !!next_col_name := existing)
    new <- tibble(TrackName = new,
                  !!next_col_name := new)
    branch <- tibble(!!prev_col_name := sapply(branch,function(x){ 
      paste0(str_split(x,"\\.")[[1]][1:length(str_split(x,"\\.")[[1]])-1],collapse = ".")}),
      !!next_col_name := branch)
    
    # Add the existing track points to t2
    if (!dim(branch)[1] == 0) {
      df <- left_join(df, bind_rows(branch,existing), by = prev_col_name)
    } else {
      df <- left_join(df, existing, by = prev_col_name)
    }
    df <- bind_rows(df, new)
    
    # Rename previous colname and previous track
    tracks_prev <- tracks_next
    prev_col_name <- next_col_name
  }
  
  # Change double track names
  df <- df %>% mutate(TrackNameUnique = seq(1:nrow(df)))
  
  return(df)
}

# Calculate moving average
ma <- function(x, n = 5, side = 2){
  stats::filter(x, rep(1 / n, n), sides = side)
}

getPCN <- function(x) {
  # Get the unique values and how often they are repeated
  y <- rle(x)
  # Bind the values in blocks of three with a dash in between
  sequences <- paste0(y$values,"-",c(y$values[-1],NA),"-",c(y$values[-1:-2],NA,NA))
  # Repeat every sequence of 3 according to the length of the chunk
  out <- rep(c(NA,sequences[-length(sequences)]),y$length)
  out
}

getPrev <- function(x) {
  # Get the unique values and how often they are repeated
  y <- rle(x)
  # Bind the values in blocks of three with a dash in between
  prev <- y$values
  # Repeat every sequence of 3 according to the length of the chunk
  out <- rep(c(NA,prev[-length(prev)]),y$length)
  out
}
# Make track df
get_tracks <- function(well_id, data, ftid){
  # Get subset of data for one well location
  data_tmp <- data %>% filter(well_name_p == well_id)
  
  # Initiate temporary output df
  make_full_tracks(data_tmp, ftid = ftid) %>% mutate(well_name_p = well_id) %>% relocate(c("TrackName","TrackNameUnique","well_name_p"))
  #return(data_tmp)
}

# Make a cluster
if(!exists("cluster")) {
  cluster <- new_cluster(35)
  cluster_library(cluster,"tidyverse")
  cluster_send(cluster,
               ma <- function(x, n = 5, side = 2){
                 stats::filter(x, rep(1 / n, n), sides = side)
               })
  cluster_send(cluster,
               getPCN <- function(x) {
                 # Get the unique values and how often they are repeated
                 y <- rle(x)
                 # Bind the values in blocks of three with a dash in between
                 sequences <- paste0(y$values,"-",c(y$values[-1],NA),"-",c(y$values[-1:-2],NA,NA))
                 # Repeat every sequence of 3 according to the length of the chunk
                 out <- rep(c(NA,sequences[-length(sequences)]),y$length)
                 out
               })
  cluster_send(cluster,
               getPrev <- function(x) {
                 # Get the unique values and how often they are repeated
                 y <- rle(x)
                 # Bind the values in blocks of three with a dash in between
                 prev <- y$values
                 # Repeat every sequence of 3 according to the length of the chunk
                 out <- rep(c(NA,prev[-length(prev)]),y$length)
                 out
               })}


#read data

#Laden Data

DMSO_tracks_data <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/Database/fixed_tracks_DMSO.csv")
View(fixed_tracks)
Metadata_DMSO <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/Metadata_DMSO.csv")
View(Metadata_DMSO)
DMSO_layout <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/layout_DMSO.csv")

DMEM_tracks_data <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/Database/fixed_tracks_DMEM.csv")
View(DMEM_tracks_data)
Metadata_DMEM <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/Metadata_DMEM.csv")
View(Metadata_DMEM)
DMEM_layout <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/layout_DMEM.csv")

Cisplatin_tracks_data <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/Database/fixed_tracks_Cisplatin.csv")
View(Cisplatin_tracks_data)
Metadata_Cisplatin <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/Metadata_Cisplatin.csv")
View(Metadata_Cisplatin)
Cisplatin_layout <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/layout_Cisplatin.csv")

Etoposide_tracks_data <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/Database/fixed_tracks_Etoposide.csv")
View(Etoposide_tracks_data)
Metadata_Etoposide <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/Metadata_Etoposide.csv")
View(Metadata_Etoposide)
Etoposide_layout <- read_csv("D:/40x images/TD54_HepG2_FUCCI_tracking_384 well plate 1/layout_Etoposide.csv")

#Data corrigeren DMSO
Data_DMSO = DMSO_tracks_data %>%
  select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number","cid","uid","alt_uid","Nuclei_Intensity_IntegratedIntensity_Cy3","Nuclei_Intensity_IntegratedIntensity_dTomato","Nuclei_Intensity_IntegratedIntensity_GFP","Nuclei_Intensity_IntegratedIntensity_Hoechst","Nuclei_Intensity_MeanIntensity_Cy3","Nuclei_Intensity_MeanIntensity_dTomato","Nuclei_Intensity_MeanIntensity_GFP","Nuclei_Intensity_MeanIntensity_Hoechst","Nuclei_AreaShape_Area"))
view(Data_DMSO)

colnames(Data_DMSO)
Data_DMSO <- Data_DMSO %>% dplyr::rename(Cy3 = Nuclei_Intensity_MeanIntensity_Cy3,
                                                             GFP = Nuclei_Intensity_MeanIntensity_GFP,
                                                             dTomato=Nuclei_Intensity_MeanIntensity_dTomato,
                                                             Hoechst=Nuclei_Intensity_MeanIntensity_Hoechst, 
                                                             IntegratedHoechst = Nuclei_Intensity_IntegratedIntensity_Hoechst,
                                                             IntegratedGFP = Nuclei_Intensity_IntegratedIntensity_GFP,
                                                             IntegratedCy3 = Nuclei_Intensity_IntegratedIntensity_Cy3,
                                                             IntegrateddTomato = Nuclei_Intensity_IntegratedIntensity_dTomato,
                                                             nucleiSize = Nuclei_AreaShape_Area)
colnames(Metadata_DMSO)
Metadata_DMSO <- Metadata_DMSO %>% dplyr::rename(Image_Group_Number=loc)
Data_DMSO<- left_join(Data_DMSO,Metadata_DMSO,by="Image_Group_Number") %>% 
  relocate(c("row","col","pos","well_name","well_name_p",
             "cid","uid","alt_uid"))

DMSO_layout <- DMSO_layout %>% select(-c("...1"))
Data_DMSO <- left_join(Data_DMSO,DMSO_layout,by="well_name")%>% 
  relocate(c("treatment","condition","dose_uM","cell_line",
             "cid","uid","alt_uid"))

#Data corrigeren DMEM
Data_DMEM = DMEM_tracks_data %>%
  select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number","cid","uid","alt_uid","Nuclei_Intensity_IntegratedIntensity_Cy3","Nuclei_Intensity_IntegratedIntensity_dTomato","Nuclei_Intensity_IntegratedIntensity_GFP","Nuclei_Intensity_IntegratedIntensity_Hoechst","Nuclei_Intensity_MeanIntensity_Cy3","Nuclei_Intensity_MeanIntensity_dTomato","Nuclei_Intensity_MeanIntensity_GFP","Nuclei_Intensity_MeanIntensity_Hoechst","Nuclei_AreaShape_Area"))
view(Data_DMEM)

colnames(Data_DMEM)
Data_DMEM <- Data_DMEM %>% dplyr::rename(Cy3 = Nuclei_Intensity_MeanIntensity_Cy3,
                                         GFP = Nuclei_Intensity_MeanIntensity_GFP,
                                         dTomato=Nuclei_Intensity_MeanIntensity_dTomato,
                                         Hoechst=Nuclei_Intensity_MeanIntensity_Hoechst, 
                                         IntegratedHoechst = Nuclei_Intensity_IntegratedIntensity_Hoechst,
                                         IntegratedGFP = Nuclei_Intensity_IntegratedIntensity_GFP,
                                         IntegratedCy3 = Nuclei_Intensity_IntegratedIntensity_Cy3,
                                         IntegrateddTomato = Nuclei_Intensity_IntegratedIntensity_dTomato,
                                         nucleiSize = Nuclei_AreaShape_Area)
colnames(Metadata_DMEM)
Metadata_DMEM <- Metadata_DMEM %>% dplyr::rename(Image_Group_Number=loc)
Data_DMEM<- left_join(Data_DMEM,Metadata_DMEM,by="Image_Group_Number") %>% 
  relocate(c("row","col","pos","well_name","well_name_p",
             "cid","uid","alt_uid"))

DMEM_layout <- DMEM_layout %>% select(-c("...1"))
Data_DMEM <- left_join(Data_DMEM,DMEM_layout,by="well_name")%>% 
  relocate(c("treatment","condition","dose_uM","cell_line",
             "cid","uid","alt_uid"))

#Data corrigeren Cisplatin
Data_Cisplatin = Cisplatin_tracks_data %>%
  select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number","cid","uid","alt_uid","Nuclei_Intensity_IntegratedIntensity_Cy3","Nuclei_Intensity_IntegratedIntensity_dTomato","Nuclei_Intensity_IntegratedIntensity_GFP","Nuclei_Intensity_IntegratedIntensity_Hoechst","Nuclei_Intensity_MeanIntensity_Cy3","Nuclei_Intensity_MeanIntensity_dTomato","Nuclei_Intensity_MeanIntensity_GFP","Nuclei_Intensity_MeanIntensity_Hoechst","Nuclei_AreaShape_Area"))
view(Data_Cisplatin)

colnames(Data_Cisplatin)
Data_Cisplatin <- Data_Cisplatin %>% dplyr::rename(Cy3 = Nuclei_Intensity_MeanIntensity_Cy3,
                                         GFP = Nuclei_Intensity_MeanIntensity_GFP,
                                         dTomato=Nuclei_Intensity_MeanIntensity_dTomato,
                                         Hoechst=Nuclei_Intensity_MeanIntensity_Hoechst, 
                                         IntegratedHoechst = Nuclei_Intensity_IntegratedIntensity_Hoechst,
                                         IntegratedGFP = Nuclei_Intensity_IntegratedIntensity_GFP,
                                         IntegratedCy3 = Nuclei_Intensity_IntegratedIntensity_Cy3,
                                         IntegrateddTomato = Nuclei_Intensity_IntegratedIntensity_dTomato,
                                         nucleiSize = Nuclei_AreaShape_Area)
colnames(Metadata_Cisplatin)
Metadata_Cisplatin <- Metadata_Cisplatin %>% dplyr::rename(Image_Group_Number=loc)
Data_Cisplatin<- left_join(Data_Cisplatin,Metadata_Cisplatin,by="Image_Group_Number") %>% 
  relocate(c("row","col","pos","well_name","well_name_p",
             "cid","uid","alt_uid"))

Cisplatin_layout <- Cisplatin_layout %>% select(-c("...1"))
Data_Cisplatin <- left_join(Data_Cisplatin,Cisplatin_layout,by="well_name")%>% 
  relocate(c("treatment","condition","dose_uM","cell_line",
             "cid","uid","alt_uid"))

#Data corrigeren Etoposide
Data_Etoposide = Etoposide_tracks_data %>%
  select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number","cid","uid","alt_uid","Nuclei_Intensity_IntegratedIntensity_Cy3","Nuclei_Intensity_IntegratedIntensity_dTomato","Nuclei_Intensity_IntegratedIntensity_GFP","Nuclei_Intensity_IntegratedIntensity_Hoechst","Nuclei_Intensity_MeanIntensity_Cy3","Nuclei_Intensity_MeanIntensity_dTomato","Nuclei_Intensity_MeanIntensity_GFP","Nuclei_Intensity_MeanIntensity_Hoechst","Nuclei_AreaShape_Area"))
view(Data_Etoposide)

colnames(Data_Etoposide)
Data_Etoposide <- Data_Etoposide %>% dplyr::rename(Cy3 = Nuclei_Intensity_MeanIntensity_Cy3,
                                                   GFP = Nuclei_Intensity_MeanIntensity_GFP,
                                                   dTomato=Nuclei_Intensity_MeanIntensity_dTomato,
                                                   Hoechst=Nuclei_Intensity_MeanIntensity_Hoechst, 
                                                   IntegratedHoechst = Nuclei_Intensity_IntegratedIntensity_Hoechst,
                                                   IntegratedGFP = Nuclei_Intensity_IntegratedIntensity_GFP,
                                                   IntegratedCy3 = Nuclei_Intensity_IntegratedIntensity_Cy3,
                                                   IntegrateddTomato = Nuclei_Intensity_IntegratedIntensity_dTomato,
                                                   nucleiSize = Nuclei_AreaShape_Area)
colnames(Metadata_Etoposide)
Metadata_Etoposide <- Metadata_Etoposide %>% dplyr::rename(Image_Group_Number=loc)
Data_Etoposide<- left_join(Data_Etoposide,Metadata_Etoposide,by="Image_Group_Number") %>% 
  relocate(c("row","col","pos","well_name","well_name_p",
             "cid","uid","alt_uid"))

Etoposide_layout <- Etoposide_layout %>% select(-c("...1"))
Data_Etoposide <- left_join(Data_Etoposide,Etoposide_layout,by="well_name")%>% 
  relocate(c("treatment","condition","dose_uM","cell_line",
             "cid","uid","alt_uid"))

# Do min-max normalisation DMSO

Data_DMSO <- Data_DMSO %>% mutate(IntegratedHoechst = ((IntegratedHoechst - min(IntegratedHoechst, na.rm = T))/
                                               (max(IntegratedHoechst, na.rm = T) - min(IntegratedHoechst, na.rm = T))),
                        IntegratedGFP = ((IntegratedGFP - min(IntegratedGFP, na.rm = T))/
                                            (max(IntegratedGFP, na.rm = T) - min(IntegratedGFP, na.rm = T))),
                        IntegratedCy3 = ((IntegratedCy3 - min(IntegratedCy3, na.rm = T))/
                                           (max(IntegratedCy3, na.rm = T) - min(IntegratedCy3, na.rm = T))),
                        IntegrateddTomato = ((IntegrateddTomato - min(IntegrateddTomato, na.rm = T))/
                                             (max(IntegrateddTomato, na.rm = T) - min(IntegrateddTomato, na.rm = T))),
                        Hoechst = ((Hoechst- min(Hoechst, na.rm = T))/
                                     (max(Hoechst, na.rm = T) - min(Hoechst, na.rm = T))),
                        GFP = ((GFP- min(GFP, na.rm = T))/
                                  (max(GFP, na.rm = T) - min(GFP, na.rm = T))),
                        Cy3 = ((Cy3- min(Cy3, na.rm = T))/
                                 (max(Cy3, na.rm = T) - min(Cy3, na.rm = T))),
                        dTomato = ((dTomato- min(dTomato, na.rm = T))/
                                   (max(dTomato, na.rm = T) - min(dTomato, na.rm = T))),
                        nucleiSize = ((nucleiSize- min(nucleiSize, na.rm = T))/
                                        (max(nucleiSize, na.rm = T) - min(nucleiSize, na.rm = T))))

# Do min-max normalisation DMEM
Data_DMEM <- Data_DMEM %>% mutate(IntegratedHoechst = ((IntegratedHoechst - min(IntegratedHoechst, na.rm = T))/
                                                         (max(IntegratedHoechst, na.rm = T) - min(IntegratedHoechst, na.rm = T))),
                                  IntegratedGFP = ((IntegratedGFP - min(IntegratedGFP, na.rm = T))/
                                                     (max(IntegratedGFP, na.rm = T) - min(IntegratedGFP, na.rm = T))),
                                  IntegratedCy3 = ((IntegratedCy3 - min(IntegratedCy3, na.rm = T))/
                                                     (max(IntegratedCy3, na.rm = T) - min(IntegratedCy3, na.rm = T))),
                                  IntegrateddTomato = ((IntegrateddTomato - min(IntegrateddTomato, na.rm = T))/
                                                         (max(IntegrateddTomato, na.rm = T) - min(IntegrateddTomato, na.rm = T))),
                                  Hoechst = ((Hoechst- min(Hoechst, na.rm = T))/
                                               (max(Hoechst, na.rm = T) - min(Hoechst, na.rm = T))),
                                  GFP = ((GFP- min(GFP, na.rm = T))/
                                           (max(GFP, na.rm = T) - min(GFP, na.rm = T))),
                                  Cy3 = ((Cy3- min(Cy3, na.rm = T))/
                                           (max(Cy3, na.rm = T) - min(Cy3, na.rm = T))),
                                  dTomato = ((dTomato- min(dTomato, na.rm = T))/
                                               (max(dTomato, na.rm = T) - min(dTomato, na.rm = T))),
                                  nucleiSize = ((nucleiSize- min(nucleiSize, na.rm = T))/
                                                  (max(nucleiSize, na.rm = T) - min(nucleiSize, na.rm = T))))
# Do min-max normalisation Cisplatin
Data_Cisplatin <- Data_Cisplatin %>% mutate(IntegratedHoechst = ((IntegratedHoechst - min(IntegratedHoechst, na.rm = T))/
                                                         (max(IntegratedHoechst, na.rm = T) - min(IntegratedHoechst, na.rm = T))),
                                  IntegratedGFP = ((IntegratedGFP - min(IntegratedGFP, na.rm = T))/
                                                     (max(IntegratedGFP, na.rm = T) - min(IntegratedGFP, na.rm = T))),
                                  IntegratedCy3 = ((IntegratedCy3 - min(IntegratedCy3, na.rm = T))/
                                                     (max(IntegratedCy3, na.rm = T) - min(IntegratedCy3, na.rm = T))),
                                  IntegrateddTomato = ((IntegrateddTomato - min(IntegrateddTomato, na.rm = T))/
                                                         (max(IntegrateddTomato, na.rm = T) - min(IntegrateddTomato, na.rm = T))),
                                  Hoechst = ((Hoechst- min(Hoechst, na.rm = T))/
                                               (max(Hoechst, na.rm = T) - min(Hoechst, na.rm = T))),
                                  GFP = ((GFP- min(GFP, na.rm = T))/
                                           (max(GFP, na.rm = T) - min(GFP, na.rm = T))),
                                  Cy3 = ((Cy3- min(Cy3, na.rm = T))/
                                           (max(Cy3, na.rm = T) - min(Cy3, na.rm = T))),
                                  dTomato = ((dTomato- min(dTomato, na.rm = T))/
                                               (max(dTomato, na.rm = T) - min(dTomato, na.rm = T))),
                                  nucleiSize = ((nucleiSize- min(nucleiSize, na.rm = T))/
                                                  (max(nucleiSize, na.rm = T) - min(nucleiSize, na.rm = T))))
# Do min-max normalisation Etoposide
Data_Etoposide <- Data_Etoposide %>% mutate(IntegratedHoechst = ((IntegratedHoechst - min(IntegratedHoechst, na.rm = T))/
                                                                   (max(IntegratedHoechst, na.rm = T) - min(IntegratedHoechst, na.rm = T))),
                                            IntegratedGFP = ((IntegratedGFP - min(IntegratedGFP, na.rm = T))/
                                                               (max(IntegratedGFP, na.rm = T) - min(IntegratedGFP, na.rm = T))),
                                            IntegratedCy3 = ((IntegratedCy3 - min(IntegratedCy3, na.rm = T))/
                                                               (max(IntegratedCy3, na.rm = T) - min(IntegratedCy3, na.rm = T))),
                                            IntegrateddTomato = ((IntegrateddTomato - min(IntegrateddTomato, na.rm = T))/
                                                                   (max(IntegrateddTomato, na.rm = T) - min(IntegrateddTomato, na.rm = T))),
                                            Hoechst = ((Hoechst- min(Hoechst, na.rm = T))/
                                                         (max(Hoechst, na.rm = T) - min(Hoechst, na.rm = T))),
                                            GFP = ((GFP- min(GFP, na.rm = T))/
                                                     (max(GFP, na.rm = T) - min(GFP, na.rm = T))),
                                            Cy3 = ((Cy3- min(Cy3, na.rm = T))/
                                                     (max(Cy3, na.rm = T) - min(Cy3, na.rm = T))),
                                            dTomato = ((dTomato- min(dTomato, na.rm = T))/
                                                         (max(dTomato, na.rm = T) - min(dTomato, na.rm = T))),
                                            nucleiSize = ((nucleiSize- min(nucleiSize, na.rm = T))/
                                                            (max(nucleiSize, na.rm = T) - min(nucleiSize, na.rm = T))))

# Reset ImageNumber to timeID

Data_DMSO <- Data_DMSO %>% group_by(well_name_p,treatment,condition,dose_uM,cell_line) %>% 
  mutate(timeID = ImageNumber - min(ImageNumber)+FIRST_TIMEID) %>% 
  ungroup %>% relocate(c("well_name","well_name_p","timeID"))

Data_DMEM <- Data_DMEM %>% group_by(well_name_p,treatment,condition,dose_uM,cell_line) %>% 
  mutate(timeID = ImageNumber - min(ImageNumber)+FIRST_TIMEID) %>% 
  ungroup %>% relocate(c("well_name","well_name_p","timeID"))

Data_Cisplatin <- Data_Cisplatin %>% group_by(well_name_p,treatment,condition,dose_uM,cell_line) %>% 
  mutate(timeID = ImageNumber - min(ImageNumber)+FIRST_TIMEID) %>% 
  ungroup %>% relocate(c("well_name","well_name_p","timeID"))

Data_Etoposide <- Data_Etoposide %>% group_by(well_name_p,treatment,condition,dose_uM,cell_line) %>% 
  mutate(timeID = ImageNumber - min(ImageNumber)+FIRST_TIMEID) %>% 
  ungroup %>% relocate(c("well_name","well_name_p","timeID"))

# Remove quickly appearing and disappearing cells, based on the uid:
# If the uid only exists for 5 frames or less, this is removed
# Only remove track dead ends

#DMSO
Data_DMSO <- Data_DMSO %>% group_by(well_name_p, uid) %>% 
  mutate(fragment_length = length(uid)) %>% ungroup() %>%
  group_by(well_name_p) %>% 
  mutate(has_daughter = sapply(alt_uid, 
                               function(x){ifelse(paste0(x,".1") %in% alt_uid |
                                                    paste0(x,".2") %in% alt_uid,T,F)})) %>% 
  filter(!(fragment_length < 5 & !has_daughter)) %>%
  collect()

#DMEM
Data_DMEM <- Data_DMEM %>% group_by(well_name_p, uid) %>% 
  mutate(fragment_length = length(uid)) %>% ungroup() %>%
  group_by(well_name_p) %>% 
  mutate(has_daughter = sapply(alt_uid, 
                               function(x){ifelse(paste0(x,".1") %in% alt_uid |
                                                    paste0(x,".2") %in% alt_uid,T,F)})) %>% 
  filter(!(fragment_length < 5 & !has_daughter)) %>%
  collect()

#Cisplatin

Data_Cisplatin <- Data_Cisplatin %>% group_by(well_name_p, uid) %>% 
  mutate(fragment_length = length(uid)) %>% ungroup() %>%
  group_by(well_name_p) %>% 
  mutate(has_daughter = sapply(alt_uid, 
                               function(x){ifelse(paste0(x,".1") %in% alt_uid |
                                                    paste0(x,".2") %in% alt_uid,T,F)})) %>% 
  filter(!(fragment_length < 5 & !has_daughter)) %>%
  collect()

#Etoposide

Data_Etoposide <- Data_Etoposide %>% group_by(well_name_p, uid) %>% 
  mutate(fragment_length = length(uid)) %>% ungroup() %>%
  group_by(well_name_p) %>% 
  mutate(has_daughter = sapply(alt_uid, 
                               function(x){ifelse(paste0(x,".1") %in% alt_uid |
                                                    paste0(x,".2") %in% alt_uid,T,F)})) %>% 
  filter(!(fragment_length < 5 & !has_daughter)) %>%
  collect()

#Make full tracks per well_id

#wellIDs_DMSO <- unique(Data_DMSO %>% pull(well_name_p))
names(wellIDs_DMSO) <- wellIDs_DMSO
tracks_df_list_DMSO <- lapply(wellIDs_DMSO, get_tracks, Data_DMSO, ftid = FIRST_TIMEID)
names(tracks_df_list_DMSO)

# Bind the track data frames
#tracks_df_DMSO <- tracks_df_list_DMSO[[1]]
for (i in seq(2,length(tracks_df_list_DMSO))) {
  tracks_df_DMSO <- bind_rows(tracks_df_DMSO,tracks_df_list_DMSO[[i]])
}

# Collect the GFP expression info per track
#df_DMSO <- pivot_longer(tracks_df_DMSO,cols = colnames(tracks_df_DMSO)[grepl("t[1-9]",colnames(tracks_df_DMSO))], 
                   names_to = "time",values_to = "alt_uid") %>% mutate(timeID = as.numeric(substring(time, 2)))
track_info_DMSO <- left_join(df_DMSO, Data_DMSO, by = c("timeID","alt_uid","well_name_p"),relationship = "many-to-many")
track_info_DMSO <- track_info_DMSO %>% mutate(time_h = (timeID * TIME_STEP) - TIME_STEP)

Data_DMSO <- Data_DMSO %>%
  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)

Data_DMEM <- Data_DMEM %>%
  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)

Data_Cisplatin <- Data_Cisplatin %>%
  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)

Data_Etoposide <- Data_Etoposide %>%
  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)

#rolling mean

Data_DMSO_ma= Data_DMSO %>%
  group_by(well_name_p,alt_uid) %>%
  mutate(ma_Hoechst = rollmean(Hoechst, k=2, fill=NA, align='center'),
         ma_GFP = rollmean(GFP, k=2, fill=NA, align='center'),
         ma_Cy3 = rollmean(Cy3, k=2, fill=NA, align='center'),
         ma_dTomato = rollmean(dTomato, k=2, fill=NA, align='center'),
         ma_IntegratedHoechst = rollmean(IntegratedHoechst, k=2, fill=NA, align='center'),
         ma_IntegratedGFP = rollmean(IntegratedGFP, k=2, fill=NA, align='center'),
         ma_IntegratedCy3 = rollmean(IntegratedCy3, k=2, fill=NA, align='center'),
         ma_IntegrateddTomato = rollmean(IntegrateddTomato, k=2, fill=NA, align='center'),
         ma_nucleiSize= rollmean(nucleiSize, k=2, fill=NA, align='center'),
         track_break = ifelse(c(FALSE, alt_uid[-1L]!= alt_uid[-length(alt_uid)]),time_h,NA)) %>% 
  collect()

Data_DMEM_ma= Data_DMEM %>%
  group_by(well_name_p,alt_uid) %>%
  mutate(ma_Hoechst = rollmean(Hoechst, k=2, fill=NA, align='center'),
         ma_GFP = rollmean(GFP, k=2, fill=NA, align='center'),
         ma_Cy3 = rollmean(Cy3, k=2, fill=NA, align='center'),
         ma_dTomato = rollmean(dTomato, k=2, fill=NA, align='center'),
         ma_IntegratedHoechst = rollmean(IntegratedHoechst, k=2, fill=NA, align='center'),
         ma_IntegratedGFP = rollmean(IntegratedGFP, k=2, fill=NA, align='center'),
         ma_IntegratedCy3 = rollmean(IntegratedCy3, k=2, fill=NA, align='center'),
         ma_IntegrateddTomato = rollmean(IntegrateddTomato, k=2, fill=NA, align='center'),
         ma_nucleiSize= rollmean(nucleiSize, k=2, fill=NA, align='center'),
         track_break = ifelse(c(FALSE, alt_uid[-1L]!= alt_uid[-length(alt_uid)]),time_h,NA)) %>% 
  collect()

Data_Cisplatin_ma= Data_Cisplatin %>%
  group_by(well_name_p,alt_uid) %>%
  mutate(ma_Hoechst = rollmean(Hoechst, k=2, fill=NA, align='center'),
         ma_GFP = rollmean(GFP, k=2, fill=NA, align='center'),
         ma_Cy3 = rollmean(Cy3, k=2, fill=NA, align='center'),
         ma_dTomato = rollmean(dTomato, k=2, fill=NA, align='center'),
         ma_IntegratedHoechst = rollmean(IntegratedHoechst, k=2, fill=NA, align='center'),
         ma_IntegratedGFP = rollmean(IntegratedGFP, k=2, fill=NA, align='center'),
         ma_IntegratedCy3 = rollmean(IntegratedCy3, k=2, fill=NA, align='center'),
         ma_IntegrateddTomato = rollmean(IntegrateddTomato, k=2, fill=NA, align='center'),
         ma_nucleiSize= rollmean(nucleiSize, k=2, fill=NA, align='center'),
         track_break = ifelse(c(FALSE, alt_uid[-1L]!= alt_uid[-length(alt_uid)]),time_h,NA)) %>% 
  collect()

Data_Etoposide_ma= Data_Etoposide %>%
  group_by(well_name_p,alt_uid) %>%
  mutate(ma_Hoechst = rollmean(Hoechst, k=2, fill=NA, align='center'),
         ma_GFP = rollmean(GFP, k=2, fill=NA, align='center'),
         ma_Cy3 = rollmean(Cy3, k=2, fill=NA, align='center'),
         ma_dTomato = rollmean(dTomato, k=2, fill=NA, align='center'),
         ma_IntegratedHoechst = rollmean(IntegratedHoechst, k=2, fill=NA, align='center'),
         ma_IntegratedGFP = rollmean(IntegratedGFP, k=2, fill=NA, align='center'),
         ma_IntegratedCy3 = rollmean(IntegratedCy3, k=2, fill=NA, align='center'),
         ma_IntegrateddTomato = rollmean(IntegrateddTomato, k=2, fill=NA, align='center'),
         ma_nucleiSize= rollmean(nucleiSize, k=2, fill=NA, align='center'),
         track_break = ifelse(c(FALSE, alt_uid[-1L]!= alt_uid[-length(alt_uid)]),time_h,NA)) %>% 
  collect()


Data_DMSO_ma_filter <- Data_DMSO_ma   %>% group_by(alt_uid, well_name_p) %>% 
  mutate(track_length = sum(!is.na(Hoechst))) %>% 
  ungroup %>% 
  filter(track_length > 5,
         !is.na(cid))

Data_DMEM_ma_filter <- Data_DMEM_ma   %>% group_by(alt_uid, well_name_p) %>% 
  mutate(track_length = sum(!is.na(Hoechst))) %>% 
  ungroup %>% 
  filter(track_length > 5,
         !is.na(cid))

Data_Cisplatin_ma_filter <- Data_Cisplatin_ma   %>% group_by(alt_uid, well_name_p) %>% 
  mutate(track_length = sum(!is.na(Hoechst))) %>% 
  ungroup %>% 
  filter(track_length > 5,
         !is.na(cid))

Data_Etoposide_ma_filter <- Data_Etoposide_ma   %>% group_by(alt_uid, well_name_p) %>% 
  mutate(track_length = sum(!is.na(Hoechst))) %>% 
  ungroup %>% 
  filter(track_length > 5,
         !is.na(cid))

# kijken naar threshold for fase
print(max(Data_DMSO_ma_filter$ma_Cy3, na.rm=TRUE))*0.05
print(max(Data_DMSO_ma_filter$ma_GFP,na.rm=TRUE))*0.05

Data_DMSO_ma_filter1=Data_DMSO_ma_filter

Data_DMSO_ma_filter1= Data_DMSO_ma_filter1 %>%mutate(phase= ifelse(ma_GFP>0.04571677& ma_Cy3<0.04780901,"G2",
                                                                   ifelse(ma_Cy3>0.04780901& ma_GFP<0.04571677,"G1",
                                                                          ifelse(ma_Cy3>0.04780901& ma_GFP>0.04571677,"G1",
                                                                                 ifelse(ma_Cy3<0.04780901& ma_GFP<0.04571677,"Early G1",NA)))))
Data_DMSO_ma_filter1$phase= as.factor(Data_DMSO_ma_filter1$phase)
Data_DMSO_ma_filter1=Data_DMSO_ma_filter1%>%filter(!is.na(phase))

ggplot_DMSO_threshold_determination= ggplot()+
  geom_point(data=Data_DMSO_ma_filter1, aes(x=ma_GFP,y=ma_Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.04571677 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.04780901  , color = "black", linetype = "dotted")+
  ylab("Notmalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  ggtitle("DMSO")

print(max(Data_DMEM_ma_filter$ma_Cy3, na.rm=TRUE))*0.05
print(max(Data_DMEM_ma_filter$ma_GFP,na.rm=TRUE))*0.05

Data_DMEM_ma_filter1=Data_DMEM_ma_filter

Data_DMEM_ma_filter1= Data_DMEM_ma_filter1 %>%mutate(phase= ifelse(ma_GFP>0.0462688& ma_Cy3<0.04635716,"G2",
                                                                   ifelse(ma_Cy3>0.04635716& ma_GFP<0.0462688,"G1",
                                                                          ifelse(ma_Cy3>0.04635716& ma_GFP>0.0462688,"G1",
                                                                                 ifelse(ma_Cy3<0.0462688& ma_GFP<0.04635716,"Early G1",NA)))))
Data_DMEM_ma_filter1$phase= as.factor(Data_DMEM_ma_filter1$phase)
Data_DMEM_ma_filter1=Data_DMEM_ma_filter1%>%filter(!is.na(phase))

ggplot_DMEM_threshold_determination= ggplot()+
  geom_point(data=Data_DMEM_ma_filter1, aes(x=ma_GFP,y=ma_Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.0462688 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.04635716   , color = "black", linetype = "dotted")+
  ylab("Notmalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  ggtitle("DMEM")

print(max(Data_Cisplatin_ma_filter$ma_Cy3, na.rm=TRUE))*0.05
print(max(Data_Cisplatin_ma_filter$ma_GFP,na.rm=TRUE))*0.05

Data_Cisplatin_ma_filter1=Data_Cisplatin_ma_filter

Data_Cisplatin_ma_filter1= Data_Cisplatin_ma_filter1 %>%mutate(phase= ifelse(ma_GFP>0.04999831& ma_Cy3<0.04988634,"G2",
                                                                   ifelse(ma_Cy3>0.04988634& ma_GFP<0.04999831,"G1",
                                                                          ifelse(ma_Cy3>0.04988634& ma_GFP>0.04999831,"G1",
                                                                                 ifelse(ma_Cy3<0.04988634& ma_GFP<0.04999831,"Early G1",NA)))))
Data_Cisplatin_ma_filter1$phase= as.factor(Data_Cisplatin_ma_filter1$phase)
Data_Cisplatin_ma_filter1=Data_Cisplatin_ma_filter1%>%filter(!is.na(phase))

ggplot_Cisplatin_threshold_determination= ggplot()+
  geom_point(data=Data_Cisplatin_ma_filter1, aes(x=ma_GFP,y=ma_Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.04999831 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.04988634   , color = "black", linetype = "dotted")+
  ylab("Notmalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  ggtitle("Cisplatin")

print(max(Data_Etoposide_ma_filter$ma_Cy3, na.rm=TRUE))*0.05
print(max(Data_Etoposide_ma_filter$ma_GFP,na.rm=TRUE))*0.05

Data_Etoposide_ma_filter1=Data_Etoposide_ma_filter

Data_Etoposide_ma_filter1= Data_Etoposide_ma_filter1 %>%mutate(phase= ifelse(ma_GFP>0.04998375& ma_Cy3<0.04998375,"G2",
                                                                   ifelse(ma_Cy3>0.04998375& ma_GFP<0.04998375,"G1",
                                                                          ifelse(ma_Cy3>0.04998375& ma_GFP>0.04998375,"G1",
                                                                                 ifelse(ma_Cy3<0.04998375& ma_GFP<0.04998375,"Early G1",NA)))))
Data_Etoposide_ma_filter1$phase= as.factor(Data_Etoposide_ma_filter1$phase)
Data_Etoposide_ma_filter1=Data_Etoposide_ma_filter1%>%filter(!is.na(phase))

ggplot_Etoposide_threshold_determination= ggplot()+
  geom_point(data=Data_Etoposide_ma_filter1, aes(x=ma_GFP,y=ma_Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.04999959 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.04998375   , color = "black", linetype = "dotted")+
  ylab("Notmalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  ggtitle("Etoposide")


#phase detemination 
#DMSO
Data_DMSO_ma_filter<- Data_DMSO_ma_filter %>% group_by(alt_uid, well_name_p) %>% mutate(diff_hoechst = ma(abs(Hoechst - ma_Hoechst), n = 2),
                                                                                           diff_size = c(NA,diff(nucleiSize)),
                                                                                           diff_GFP = c(NA,diff(ma(GFP, n = 2))),
                                                                                           diff_Cy3 = c(NA,ma(diff(ma_Cy3)))) %>%
  mutate(gem_over_cdt1 = ifelse(ma_GFP == 0, -100,
                                ifelse(ma_Cy3 < 0.04780901 & ma_GFP < 0.04571677, -100,
                                       ifelse(ma_Cy3 == 0 & ma_GFP > 0.04571677, 100,
                                              log(ma_GFP/(ma_Cy3))))),
         peaks_Hoechst_binary = diff_hoechst > 0.03,
         peaks_Hoechst_time = ifelse(peaks_Hoechst_binary,time_h,NA),
         peak_size_binary = diff_size > 0.05,
         peak_size_binshift = c(NA,peak_size_binary[-length(peak_size_binary)]),
         dip_size_binary = diff_size < -0.05,
         peaks_size_time = ifelse(dip_size_binary & (!peak_size_binshift),time_h,NA),
         dip_GFP_binary = diff_GFP < -0.5*sd(GFP, na.rm = T),
         peaks_GFP_time = ifelse(dip_GFP_binary,time_h,NA),
         Cy3_decrease = diff_Cy3 < -0.15*sd(Cy3, na.rm = T)
  ) %>%
  collect()

Data_DMSO_ma_filter_break <- Data_DMSO_ma_filter %>% group_by(alt_uid, well_name_p, time_h) %>%
  mutate(Consistency = sum(c(dip_GFP_binary,(dip_size_binary & !peak_size_binshift),peaks_Hoechst_binary), na.rm = T),
         Break = dip_GFP_binary & Consistency > 0,
         Division = ifelse(dip_GFP_binary & Consistency == 3,"High",
                           ifelse(dip_GFP_binary & Consistency == 2,"Medium",
                                  ifelse(dip_GFP_binary & Consistency == 1,"Low",NA)))) %>%
  collect()

Data_DMSO_ma_filter_break <- Data_DMSO_ma_filter_break %>% group_by(alt_uid, well_name_p, time_h)%>%
  mutate(Break_filtered = ifelse(c(diff(Break),NA) == -1,T,F),
         Division_time = ifelse(Break_filtered, time_h,NA),
         Cycle_length = replace(Division_time,which(!is.na(Division_time)),c(NA,diff(na.omit(Division_time)))),
         Phase = ifelse(gem_over_cdt1 == -100, "G1",
                        ifelse(gem_over_cdt1 == 100, "G2",
                               ifelse(abs(gem_over_cdt1) < 2,"G1/S",
                                      ifelse(ma_GFP > ma_Cy3, "G2", "G1")))))%>%  
  select(-c(peaks_Hoechst_binary,peak_size_binary,peak_size_binshift,diff_hoechst,diff_size,diff_GFP,diff_Cy3)) %>%
  collect()

Data_DMSO_ma_filter_break <- Data_DMSO_ma_filter_break %>%filter(sum(Break_filtered, na.rm = T) < 5) %>% ungroup()
Data_DMSO_ma_filter_break <- Data_DMSO_ma_filter_break %>%group_by(alt_uid, well_name_p) %>%
  mutate(PrevCurNext = getPCN(Phase),
         Prev = getPrev(Phase),
         PhaseLength = rep(rle(Phase)$length,rle(Phase)$length)) %>%
  mutate(Phase_corrected = ifelse(PrevCurNext == "G2-G1/S-G1", "G1", 
                                  ifelse(PrevCurNext %in% c("G2-G1-G2",
                                                            "G2-G1/S-G2",
                                                            "G1-G2-G1",
                                                            "G1-G1/S-G1",
                                                            "G1/S-G1-G1/S",
                                                            "G1/S-G2-G1/S") & PhaseLength < 5,
                                         Prev, Phase))) %>%
  collect()

phase_stats <- Data_DMSO_ma_filter_break %>% group_by(alt_uid, well_name_p) %>%
  mutate(phaseLength = unlist(sapply(rle(Phase_corrected)$length,function(x){seq(1,x)}, simplify = T)),
         phaseChange = c(!Phase_corrected[-1] == Phase_corrected[-length(Phase_corrected)],NA),
         timeUntilPhaseChange = ifelse(phaseChange,time_h,NA),
         timeUntilPhaseChange = replace(timeUntilPhaseChange,which(!is.na(timeUntilPhaseChange)),
                                        c(na.omit(timeUntilPhaseChange)[1],diff(na.omit(timeUntilPhaseChange))))) %>%
  filter(phaseChange, timeUntilPhaseChange > 2) %>% 
  mutate(phaseID = row_number()) %>% 
  ungroup() 
phase_stats <- phase_stats %>% group_by(cid,alt_uid,well_name_p, well_name, Phase_corrected, phaseID,
                                        treatment,condition,dose_uM,cell_line) %>%
  summarise(phaseLength = mean(timeUntilPhaseChange, na.rm = T)) %>%
  ungroup()

#DMEM

Data_DMEM_ma_filter<- Data_DMEM_ma_filter %>% group_by(alt_uid, well_name_p) %>% mutate(diff_hoechst = ma(abs(Hoechst - ma_Hoechst), n = 2),
                                                                                        diff_size = c(NA,diff(nucleiSize)),
                                                                                        diff_GFP = c(NA,diff(ma(GFP, n = 2))),
                                                                                        diff_Cy3 = c(NA,ma(diff(ma_Cy3)))) %>%
  mutate(gem_over_cdt1 = ifelse(ma_GFP == 0, -100,
                                ifelse(ma_Cy3 < 0.04635716  & ma_GFP < 0.0462688, -100,
                                       ifelse(ma_Cy3 == 0 & ma_GFP > 0.0462688, 100,
                                              log(ma_GFP/(ma_Cy3))))),
         peaks_Hoechst_binary = diff_hoechst > 0.03,
         peaks_Hoechst_time = ifelse(peaks_Hoechst_binary,time_h,NA),
         peak_size_binary = diff_size > 0.05,
         peak_size_binshift = c(NA,peak_size_binary[-length(peak_size_binary)]),
         dip_size_binary = diff_size < -0.05,
         peaks_size_time = ifelse(dip_size_binary & (!peak_size_binshift),time_h,NA),
         dip_GFP_binary = diff_GFP < -0.5*sd(GFP, na.rm = T),
         peaks_GFP_time = ifelse(dip_GFP_binary,time_h,NA),
         Cy3_decrease = diff_Cy3 < -0.15*sd(Cy3, na.rm = T)
  ) %>%
  collect()

Data_DMEM_ma_filter_break <- Data_DMEM_ma_filter %>% group_by(alt_uid, well_name_p, time_h) %>%
  mutate(Consistency = sum(c(dip_GFP_binary,(dip_size_binary & !peak_size_binshift),peaks_Hoechst_binary), na.rm = T),
         Break = dip_GFP_binary & Consistency > 0,
         Division = ifelse(dip_GFP_binary & Consistency == 3,"High",
                           ifelse(dip_GFP_binary & Consistency == 2,"Medium",
                                  ifelse(dip_GFP_binary & Consistency == 1,"Low",NA)))) %>%
  collect()

Data_DMEM_ma_filter_break <- Data_DMEM_ma_filter_break %>% group_by(alt_uid, well_name_p, time_h)%>%
  mutate(Break_filtered = ifelse(c(diff(Break),NA) == -1,T,F),
         Division_time = ifelse(Break_filtered, time_h,NA),
         Cycle_length = replace(Division_time,which(!is.na(Division_time)),c(NA,diff(na.omit(Division_time)))),
         Phase = ifelse(gem_over_cdt1 == -100, "G1",
                        ifelse(gem_over_cdt1 == 100, "G2",
                               ifelse(abs(gem_over_cdt1) < 2,"G1/S",
                                      ifelse(ma_GFP > ma_Cy3, "G2", "G1")))))%>%  
  select(-c(peaks_Hoechst_binary,peak_size_binary,peak_size_binshift,diff_hoechst,diff_size,diff_GFP,diff_Cy3)) %>%
  collect()

Data_DMEM_ma_filter_break <- Data_DMEM_ma_filter_break %>%filter(sum(Break_filtered, na.rm = T) < 5) %>% ungroup()
Data_DMEM_ma_filter_break <- Data_DMEM_ma_filter_break %>%group_by(alt_uid, well_name_p) %>%
  mutate(PrevCurNext = getPCN(Phase),
         Prev = getPrev(Phase),
         PhaseLength = rep(rle(Phase)$length,rle(Phase)$length)) %>%
  mutate(Phase_corrected = ifelse(PrevCurNext == "G2-G1/S-G1", "G1", 
                                  ifelse(PrevCurNext %in% c("G2-G1-G2",
                                                            "G2-G1/S-G2",
                                                            "G1-G2-G1",
                                                            "G1-G1/S-G1",
                                                            "G1/S-G1-G1/S",
                                                            "G1/S-G2-G1/S") & PhaseLength < 5,
                                         Prev, Phase))) %>%
  collect()

phase_stats_DMEM <- Data_DMEM_ma_filter_break %>% group_by(alt_uid, well_name_p) %>%
  mutate(phaseLength = unlist(sapply(rle(Phase_corrected)$length,function(x){seq(1,x)}, simplify = T)),
         phaseChange = c(!Phase_corrected[-1] == Phase_corrected[-length(Phase_corrected)],NA),
         timeUntilPhaseChange = ifelse(phaseChange,time_h,NA),
         timeUntilPhaseChange = replace(timeUntilPhaseChange,which(!is.na(timeUntilPhaseChange)),
                                        c(na.omit(timeUntilPhaseChange)[1],diff(na.omit(timeUntilPhaseChange))))) %>%
  filter(phaseChange, timeUntilPhaseChange > 2) %>% 
  mutate(phaseID = row_number()) %>% 
  ungroup() 
phase_stats_DMEM <- phase_stats_DMEM %>% group_by(cid,alt_uid,well_name_p, well_name, Phase_corrected, phaseID,
                                        treatment,condition,dose_uM,cell_line) %>%
  summarise(phaseLength = mean(timeUntilPhaseChange, na.rm = T)) %>%
  ungroup()

#Cisplatin

Data_Cisplatin_ma_filter<- Data_Cisplatin_ma_filter %>% group_by(alt_uid, well_name_p) %>% mutate(diff_hoechst = ma(abs(Hoechst - ma_Hoechst), n = 2),
                                                                                        diff_size = c(NA,diff(nucleiSize)),
                                                                                        diff_GFP = c(NA,diff(ma(GFP, n = 2))),
                                                                                        diff_Cy3 = c(NA,ma(diff(ma_Cy3)))) %>%
  mutate(gem_over_cdt1 = ifelse(ma_GFP == 0, -100,
                                ifelse(ma_Cy3 < 0.04988634  & ma_GFP < 0.04999831, -100,
                                       ifelse(ma_Cy3 == 0 & ma_GFP > 0.04999831, 100,
                                              log(ma_GFP/(ma_Cy3))))),
         peaks_Hoechst_binary = diff_hoechst > 0.03,
         peaks_Hoechst_time = ifelse(peaks_Hoechst_binary,time_h,NA),
         peak_size_binary = diff_size > 0.05,
         peak_size_binshift = c(NA,peak_size_binary[-length(peak_size_binary)]),
         dip_size_binary = diff_size < -0.05,
         peaks_size_time = ifelse(dip_size_binary & (!peak_size_binshift),time_h,NA),
         dip_GFP_binary = diff_GFP < -0.5*sd(GFP, na.rm = T),
         peaks_GFP_time = ifelse(dip_GFP_binary,time_h,NA),
         Cy3_decrease = diff_Cy3 < -0.15*sd(Cy3, na.rm = T)
  ) %>%
  collect()

Data_Cisplatin_ma_filter_break <- Data_Cisplatin_ma_filter %>% group_by(alt_uid, well_name_p, time_h) %>%
  mutate(Consistency = sum(c(dip_GFP_binary,(dip_size_binary & !peak_size_binshift),peaks_Hoechst_binary), na.rm = T),
         Break = dip_GFP_binary & Consistency > 0,
         Division = ifelse(dip_GFP_binary & Consistency == 3,"High",
                           ifelse(dip_GFP_binary & Consistency == 2,"Medium",
                                  ifelse(dip_GFP_binary & Consistency == 1,"Low",NA)))) %>%
  collect()

Data_Cisplatin_ma_filter_break <- Data_Cisplatin_ma_filter_break %>% group_by(alt_uid, well_name_p, time_h)%>%
  mutate(Break_filtered = ifelse(c(diff(Break),NA) == -1,T,F),
         Division_time = ifelse(Break_filtered, time_h,NA),
         Cycle_length = replace(Division_time,which(!is.na(Division_time)),c(NA,diff(na.omit(Division_time)))),
         Phase = ifelse(gem_over_cdt1 == -100, "G1",
                        ifelse(gem_over_cdt1 == 100, "G2",
                               ifelse(abs(gem_over_cdt1) < 2,"G1/S",
                                      ifelse(ma_GFP > ma_Cy3, "G2", "G1")))))%>%  
  select(-c(peaks_Hoechst_binary,peak_size_binary,peak_size_binshift,diff_hoechst,diff_size,diff_GFP,diff_Cy3)) %>%
  collect()

Data_Cisplatin_ma_filter_break <- Data_Cisplatin_ma_filter_break %>%filter(sum(Break_filtered, na.rm = T) < 5) %>% ungroup()
Data_Cisplatin_ma_filter_break <- Data_Cisplatin_ma_filter_break %>%group_by(alt_uid, well_name_p) %>%
  mutate(PrevCurNext = getPCN(Phase),
         Prev = getPrev(Phase),
         PhaseLength = rep(rle(Phase)$length,rle(Phase)$length)) %>%
  mutate(Phase_corrected = ifelse(PrevCurNext == "G2-G1/S-G1", "G1", 
                                  ifelse(PrevCurNext %in% c("G2-G1-G2",
                                                            "G2-G1/S-G2",
                                                            "G1-G2-G1",
                                                            "G1-G1/S-G1",
                                                            "G1/S-G1-G1/S",
                                                            "G1/S-G2-G1/S") & PhaseLength < 5,
                                         Prev, Phase))) %>%
  collect()

phase_stats_Cisplatin <- Data_Cisplatin_ma_filter_break %>% group_by(alt_uid, well_name_p) %>%
  mutate(phaseLength = unlist(sapply(rle(Phase_corrected)$length,function(x){seq(1,x)}, simplify = T)),
         phaseChange = c(!Phase_corrected[-1] == Phase_corrected[-length(Phase_corrected)],NA),
         timeUntilPhaseChange = ifelse(phaseChange,time_h,NA),
         timeUntilPhaseChange = replace(timeUntilPhaseChange,which(!is.na(timeUntilPhaseChange)),
                                        c(na.omit(timeUntilPhaseChange)[1],diff(na.omit(timeUntilPhaseChange))))) %>%
  filter(phaseChange, timeUntilPhaseChange > 2) %>% 
  mutate(phaseID = row_number()) %>% 
  ungroup() 
phase_stats <- phase_stats %>% group_by(cid,alt_uid,well_name_p, well_name, Phase_corrected, phaseID,
                                        treatment,condition,dose_uM,cell_line) %>%
  summarise(phaseLength = mean(timeUntilPhaseChange, na.rm = T)) %>%
  ungroup()

#Etoposide

Data_Etoposide_ma_filter<- Data_Etoposide_ma_filter %>% group_by(alt_uid, well_name_p) %>% mutate(diff_hoechst = ma(abs(Hoechst - ma_Hoechst), n = 2),
                                                                                                  diff_size = c(NA,diff(nucleiSize)),
                                                                                                  diff_GFP = c(NA,diff(ma(GFP, n = 2))),
                                                                                                  diff_Cy3 = c(NA,ma(diff(ma_Cy3)))) %>%
  mutate(gem_over_cdt1 = ifelse(ma_GFP == 0, -100,
                                ifelse(ma_Cy3 < 0.04998375  & ma_GFP < 0.04999959, -100,
                                       ifelse(ma_Cy3 == 0 & ma_GFP > 0.04999959, 100,
                                              log(ma_GFP/(ma_Cy3))))),
         peaks_Hoechst_binary = diff_hoechst > 0.03,
         peaks_Hoechst_time = ifelse(peaks_Hoechst_binary,time_h,NA),
         peak_size_binary = diff_size > 0.05,
         peak_size_binshift = c(NA,peak_size_binary[-length(peak_size_binary)]),
         dip_size_binary = diff_size < -0.05,
         peaks_size_time = ifelse(dip_size_binary & (!peak_size_binshift),time_h,NA),
         dip_GFP_binary = diff_GFP < -0.5*sd(GFP, na.rm = T),
         peaks_GFP_time = ifelse(dip_GFP_binary,time_h,NA),
         Cy3_decrease = diff_Cy3 < -0.15*sd(Cy3, na.rm = T)
  ) %>%
  collect()

Data_Etoposide_ma_filter_break <- Data_Etoposide_ma_filter %>% group_by(alt_uid, well_name_p, time_h) %>%
  mutate(Consistency = sum(c(dip_GFP_binary,(dip_size_binary & !peak_size_binshift),peaks_Hoechst_binary), na.rm = T),
         Break = dip_GFP_binary & Consistency > 0,
         Division = ifelse(dip_GFP_binary & Consistency == 3,"High",
                           ifelse(dip_GFP_binary & Consistency == 2,"Medium",
                                  ifelse(dip_GFP_binary & Consistency == 1,"Low",NA)))) %>%
  collect()

Data_Etoposide_ma_filter_break <- Data_Etoposide_ma_filter_break %>% group_by(alt_uid, well_name_p, time_h)%>%
  mutate(Break_filtered = ifelse(c(diff(Break),NA) == -1,T,F),
         Division_time = ifelse(Break_filtered, time_h,NA),
         Cycle_length = replace(Division_time,which(!is.na(Division_time)),c(NA,diff(na.omit(Division_time)))),
         Phase = ifelse(gem_over_cdt1 == -100, "G1",
                        ifelse(gem_over_cdt1 == 100, "G2",
                               ifelse(abs(gem_over_cdt1) < 2,"G1/S",
                                      ifelse(ma_GFP > ma_Cy3, "G2", "G1")))))%>%  
  select(-c(peaks_Hoechst_binary,peak_size_binary,peak_size_binshift,diff_hoechst,diff_size,diff_GFP,diff_Cy3)) %>%
  collect()

Data_Etoposide_ma_filter_break <- Data_Etoposide_ma_filter_break %>%filter(sum(Break_filtered, na.rm = T) < 5) %>% ungroup()
Data_Etoposide_ma_filter_break <- Data_Etoposide_ma_filter_break %>%group_by(alt_uid, well_name_p) %>%
  mutate(PrevCurNext = getPCN(Phase),
         Prev = getPrev(Phase),
         PhaseLength = rep(rle(Phase)$length,rle(Phase)$length)) %>%
  mutate(Phase_corrected = ifelse(PrevCurNext == "G2-G1/S-G1", "G1", 
                                  ifelse(PrevCurNext %in% c("G2-G1-G2",
                                                            "G2-G1/S-G2",
                                                            "G1-G2-G1",
                                                            "G1-G1/S-G1",
                                                            "G1/S-G1-G1/S",
                                                            "G1/S-G2-G1/S") & PhaseLength < 5,
                                         Prev, Phase))) %>%
  collect()

phase_stats_Etoposide <- Data_Etoposide_ma_filter_break %>% group_by(alt_uid, well_name_p) %>%
  mutate(phaseLength = unlist(sapply(rle(Phase_corrected)$length,function(x){seq(1,x)}, simplify = T)),
         phaseChange = c(!Phase_corrected[-1] == Phase_corrected[-length(Phase_corrected)],NA),
         timeUntilPhaseChange = ifelse(phaseChange,time_h,NA),
         timeUntilPhaseChange = replace(timeUntilPhaseChange,which(!is.na(timeUntilPhaseChange)),
                                        c(na.omit(timeUntilPhaseChange)[1],diff(na.omit(timeUntilPhaseChange))))) %>%
  filter(phaseChange, timeUntilPhaseChange > 2) %>% 
  mutate(phaseID = row_number()) %>% 
  ungroup() 
phase_stats <- phase_stats %>% group_by(cid,alt_uid,well_name_p, well_name, Phase_corrected, phaseID,
                                        treatment,condition,dose_uM,cell_line) %>%
  summarise(phaseLength = mean(timeUntilPhaseChange, na.rm = T)) %>%
  ungroup()

#boxplot phase length
count = Data_DMSO_ma_filter_break %>% group_by(alt_uid,well_name_p,dose_uM,condition)%>%count(PhaseLength,Phase)%>%filter(!is.na(Phase))
count_corrected = count %>% group_by(alt_uid,well_name_p,Phase,dose_uM,condition) %>% summarize(total=sum(n, na.rm = TRUE))

count_filter1_DMSO = Data_DMSO_ma_filter1 %>% group_by(alt_uid,well_name_p,dose_uM,condition)%>%count(track_length,phase)%>%filter(!is.na(phase))
count_corrected_filter_DMSO = count_filter1_DMSO  %>% group_by(alt_uid,well_name_p,phase,dose_uM,condition) %>% summarize(total=sum(n, na.rm = TRUE))

plot_DMSO_metG1/S<- ggplot()+
  geom_boxplot(data = count_corrected, aes(x = condition, y = total*1.67, fill = Phase)) +
  ylab("Phase length (h)") + xlab("condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","G1/S" = "gold")) + 
  theme_classic()

plot_DMSO_metEarlyG1<- ggplot()+
  geom_boxplot(data = count_corrected_filter_DMSO, aes(x = condition, y = total*1.67, fill = phase)) +
  ylab("Phase length (h)") + xlab("condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","Early G1" = "grey")) +
  theme_classic()

count_DMEM = Data_DMEM_ma_filter_break %>% group_by(alt_uid,well_name_p,dose_uM,condition)%>%count(PhaseLength,Phase)%>%filter(!is.na(Phase))
count_corrected_DMEM = count_DMEM %>% group_by(alt_uid,well_name_p,Phase,dose_uM,condition) %>% summarize(total=sum(n, na.rm = TRUE))

count_filter1_DMEM = Data_DMEM_ma_filter1 %>% group_by(alt_uid,well_name_p,dose_uM,condition)%>%count(track_length,phase)%>%filter(!is.na(phase))
count_corrected_filter_DMEM = count_filter1_DMEM  %>% group_by(alt_uid,well_name_p,phase,dose_uM,condition) %>% summarize(total=sum(n, na.rm = TRUE))

plot_DMEM_metG1/S<- ggplot()+
  geom_boxplot(data = count_corrected_DMEM, aes(x = condition, y = total*1.67, fill = Phase)) +
  ylab("Phase length (h)") + xlab("condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","G1/S" = "gold")) + 
  theme_classic()

plot_DMEM_metEarlyG1<- ggplot()+
  geom_boxplot(data = count_corrected_filter_DMEM, aes(x = condition, y = total*1.67, fill = phase)) +
  ylab("Phase length (h)") + xlab("condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","Early G1" = "grey")) +
  theme_classic()

count_Cisplatin = Data_Cisplatin_ma_filter_break %>% group_by(alt_uid,well_name_p,dose_uM,condition)%>%count(PhaseLength,Phase)%>%filter(!is.na(Phase))
count_corrected_Cisplatin = count_Cisplatin %>% group_by(alt_uid,well_name_p,Phase,dose_uM,condition) %>% summarize(total=sum(n, na.rm = TRUE))

count_filter1_Cisplatin = Data_Cisplatin_ma_filter1 %>% group_by(alt_uid,well_name_p,dose_uM,condition)%>%count(track_length,phase)%>%filter(!is.na(phase))
count_corrected_filter_Cisplatin = count_filter1_Cisplatin  %>% group_by(alt_uid,well_name_p,phase,dose_uM,condition) %>% summarize(total=sum(n, na.rm = TRUE))

plot_Cisplatin_metG1/S<- ggplot()+
  geom_boxplot(data = count_corrected_Cisplatin, aes(x = condition, y = total*1.67, fill = Phase)) +
  ylab("Phase length (h)") + xlab("Cisplatin") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","G1/S" = "gold")) + 
  theme_classic()

belangrijk_figuur_Cisplatin = ggplot()+
  geom_boxplot(data = count_corrected, aes(x = condition, y = total*1.67, fill = Phase))+
  geom_boxplot(data = count_corrected_DMEM, aes(x = condition, y = total*1.67, fill = Phase)) +
  geom_boxplot(data = count_corrected_Cisplatin, aes(x = condition, y = total*1.67, fill = Phase)) +
  ylab("Phase length (h)") + xlab("Condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","G1/S" = "gold")) + 
  theme_classic()

plot_Cisplatin_metEarlyG1<- ggplot()+
  geom_boxplot(data = count_corrected_filter_Cisplatin, aes(x = condition, y = total*1.67, fill = phase)) +
  ylab("Phase length (h)") + xlab("condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","Early G1" = "grey")) +
  theme_classic()

plot_Cistplatin_DMSO_earlyg1 = ggplot()+
  geom_boxplot(data = count_corrected_filter_DMSO, aes(x = condition, y = total*1.67, fill = phase))+
  geom_boxplot(data = count_corrected_filter_DMEM, aes(x = condition, y = total*1.67, fill = phase)) +
  geom_boxplot(data = count_corrected_filter_Cisplatin, aes(x = condition, y = total*1.67, fill = phase)) +
  ylab("Phase length (h)") + xlab("condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","Early G1" = "grey")) +
  theme_classic()


count_Etoposide = Data_Etoposide_ma_filter_break %>% group_by(alt_uid,well_name_p,dose_uM,condition)%>%count(PhaseLength,Phase)%>%filter(!is.na(Phase))
count_corrected_Etoposide = count_Etoposide %>% group_by(alt_uid,well_name_p,Phase,dose_uM,condition) %>% summarize(total=sum(n, na.rm = TRUE))

count_filter1_Etoposide = Data_Etoposide_ma_filter1 %>% group_by(alt_uid,well_name_p,dose_uM,condition)%>%count(track_length,phase)%>%filter(!is.na(phase))
count_corrected_filter_Etoposide = count_filter1_Etoposide  %>% group_by(alt_uid,well_name_p,phase,dose_uM,condition) %>% summarize(total=sum(n, na.rm = TRUE))

plot_Etoposide_metG1/S<- ggplot()+
  geom_boxplot(data = count_corrected_Etoposide, aes(x = condition, y = total*1.67, fill = Phase)) +
  ylab("Phase length (h)") + xlab("Etoposide") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","G1/S" = "gold")) + 
  theme_classic()

belangrijk_figuur_Etoposide = ggplot()+
  geom_boxplot(data = count_corrected, aes(x = condition, y = total*1.67, fill = Phase))+
  geom_boxplot(data = count_corrected_DMEM, aes(x = condition, y = total*1.67, fill = Phase)) +
  geom_boxplot(data = count_corrected_Etoposide, aes(x = condition, y = total*1.67, fill = Phase)) +
  ylab("Phase length (h)") + xlab("Condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","G1/S" = "gold")) + 
  theme_classic()

plot_Cisplatin_metEarlyG1<- ggplot()+
  geom_boxplot(data = count_corrected_filter_Etoposide, aes(x = condition, y = total*1.67, fill = phase)) +
  ylab("Phase length (h)") + xlab("condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","Early G1" = "grey")) +
  theme_classic()

plot_Cistplatin_DMSO_earlyg1 = ggplot()+
  geom_boxplot(data = count_corrected_filter_DMSO, aes(x = condition, y = total*1.67, fill = phase))+
  geom_boxplot(data = count_corrected_filter_DMEM, aes(x = condition, y = total*1.67, fill = phase)) +
  geom_boxplot(data = count_corrected_filter_Etoposide, aes(x = condition, y = total*1.67, fill = phase)) +
  ylab("Phase length (h)") + xlab("condition") +
  scale_fill_manual(values = c("G1" = "red", "G2" = "green","Early G1" = "grey")) +
  theme_classic()
#plot
X = 1
filtered_track_ma_2 = Data_DMEM_ma  %>%
  filter(cid == X, Image_Group_Number==2)

X = 1
filtered_track = Data_DMEM %>%
  filter(cid == X, Image_Group_Number==1)

ma_mean_intensity_2<- ggplot()+
  geom_line(data = filtered_track, aes(x = timeID, y = Cy3, linetype = alt_uid), color = "red") + # must include argument label "data"
  geom_line(data = filtered_track, aes(x = timeID, y = GFP, linetype = alt_uid),color = "green")+
  geom_line(data = filtered_track, aes(x = timeID, y = dTomato, linetype = alt_uid),color = "darkred")+
  geom_line(data = filtered_track, aes(x = timeID, y = Hoechst, linetype = alt_uid),color = "blue")+
  geom_line(data = filtered_track, aes(x = timeID, y = nucleiSize, linetype = alt_uid),color = "black")+
  labs(x ="timeID", y = "Mean Intensity")

ma_DMSO <- ggplot()+
  geom_line(data = filtered_track_ma_2, aes(x = timeID, y = ma_Cy3, linetype = alt_uid), color = "red") + # must include argument label "data"
  geom_line(data = filtered_track_ma_2, aes(x = timeID, y = ma_GFP, linetype = alt_uid),color = "green")+
  geom_line(data = filtered_track_ma_2, aes(x = timeID, y = ma_dTomato, linetype = alt_uid),color = "darkred")+
  geom_line(data = filtered_track_ma_2, aes(x = timeID, y = ma_Hoechst, linetype = alt_uid),color = "blue")+
  geom_line(data = filtered_track_ma_2, aes(x = timeID, y = ma_nucleiSize, linetype = alt_uid),color = "black")+
  labs(x ="timeID", y = "Mean Intensity")

ggplot()+
  geom_line(data = filtered_track_ma_2, aes(x = timeID, y = ma_Cy3, linetype = alt_uid), color = "red")+
  geom_line(data = filtered_track_ma_2, aes(x = timeID, y = ma_GFP, linetype = alt_uid),color = "green")

#staistical tests

#Cisplatin
G1_Data_Cisplatin = count_corrected_Cisplatin %>% filter(Phase == "G1")
G2_Data_Cisplatin = count_corrected_Cisplatin %>% filter(Phase == "G2")
geel_Data_Cisplatin = count_corrected_Cisplatin %>% filter(Phase == "G1/S")


t.test_G1_Cisplatin = print(pairwise.t.test(G1_Data_Cisplatin %>% pull (total), G1_Data_Cisplatin %>% pull (condition),
                      p.adjust.method = "bonferroni"))

t.test__geel_Cisplatin = print(pairwise.t.test(geel_Data_Cisplatin %>% pull (total), geel_Data_Cisplatin %>% pull (condition),
                                            p.adjust.method = "bonferroni"))

t.test_G2_Cisplatin = print(pairwise.t.test(G2_Data_Cisplatin %>% pull (total), G2_Data_Cisplatin %>% pull (condition),
                                              p.adjust.method = "bonferroni"))

#Etoposide
G1_Data_Etoposide = count_corrected_Etoposide %>% filter(Phase == "G1")
G2_Data_Etoposide = count_corrected_Etoposide %>% filter(Phase == "G2")
geel_Data_Etoposide = count_corrected_Etoposide %>% filter(Phase == "G1/S")

t.test_G1_Etoposide = print(pairwise.t.test(G1_Data_Etoposide %>% pull (total), G1_Data_Etoposide %>% pull (condition),
                                            p.adjust.method = "bonferroni"))

t.test__geel_Etoposide = print(pairwise.t.test(geel_Data_Etoposide %>% pull (total), geel_Data_Etoposide %>% pull (condition),
                                               p.adjust.method = "bonferroni"))

t.test_G2_Etoposide = print(pairwise.t.test(G2_Data_Etoposide %>% pull (total), G2_Data_Etoposide %>% pull (condition),
                                            p.adjust.method = "bonferroni"))

#Vergelijken met DMSO en Cisplatin en DMEM phase length with G1-S

combinded_DMSO_Cisplatin_DMEM= rbind(count_corrected,count_corrected_Cisplatin,count_corrected_DMEM)

combinded_DMSO_Cisplatin_DMEM__G1 = combinded_DMSO_Cisplatin_DMEM %>% filter(Phase == "G1")
combinded_DMSO_Cisplatin_DMEM__geel = combinded_DMSO_Cisplatin_DMEM %>% filter(Phase == "G1/S")
combinded_DMSO_Cisplatin_DMEM__G2= combinded_DMSO_Cisplatin_DMEM %>% filter(Phase == "G2")

t.test_G1_vergelijken = print(pairwise.t.test(combinded_DMSO_Cisplatin_DMEM__G1 %>% pull (total), combinded_DMSO_Cisplatin_DMEM__G1%>% pull (condition),
                                            p.adjust.method = "bonferroni"))

t.test_geel_vergelijken = print(pairwise.t.test(combinded_DMSO_Cisplatin_DMEM__geel %>% pull (total), combinded_DMSO_Cisplatin_DMEM__geel%>% pull (condition),
                                              p.adjust.method = "bonferroni"))
t.test_G2_vergelijken = print(pairwise.t.test(combinded_DMSO_Cisplatin_DMEM__G2 %>% pull (total), combinded_DMSO_Cisplatin_DMEM__G2%>% pull (condition),
                                                p.adjust.method = "bonferroni"))



#Vergelijken met DMSO en Etoposide en DMEM phase length with G1-S

combinded_DMSO_Etoposide_DMEM= rbind(count_corrected,count_corrected_Etoposide,count_corrected_DMEM)

combinded_DMSO_Etoposide_DMEM__G1 = combinded_DMSO_Etoposide_DMEM %>% filter(Phase == "G1")
combinded_DMSO_Etoposide_DMEM__geel = combinded_DMSO_Etoposide_DMEM %>% filter(Phase == "G1/S")
combinded_DMSO_Etoposide_DMEM__G2= combinded_DMSO_Etoposide_DMEM %>% filter(Phase == "G2")

t.test_G1_vergelijken_Etoposide = print(pairwise.t.test(combinded_DMSO_Etoposide_DMEM__G1 %>% pull (total), combinded_DMSO_Etoposide_DMEM__G1%>% pull (condition),
                                              p.adjust.method = "bonferroni"))

t.test_geel_vergelijken_Etoposide = print(pairwise.t.test(combinded_DMSO_Etoposide_DMEM__geel %>% pull (total), combinded_DMSO_Etoposide_DMEM__geel%>% pull (condition),
                                                p.adjust.method = "bonferroni"))
t.test_G2_vergelijken_Etoposide = print(pairwise.t.test(combinded_DMSO_Etoposide_DMEM__G2 %>% pull (total), combinded_DMSO_Etoposide_DMEM__G2%>% pull (condition),
                                              p.adjust.method = "bonferroni"))

#Vergelijken met DMSO en Cisplatin en DMEM phase length with Early G1

combinded_DMSO_Cisplatin_DMEM_earlyG1= rbind(count_corrected_filter_DMSO,count_corrected_filter_Cisplatin,count_corrected_filter_DMEM)

combinded_DMSO_Cisplatin_DMEM_earlyG1__G1 = combinded_DMSO_Cisplatin_DMEM_earlyG1 %>% filter(phase == "G1")
combinded_DMSO_Cisplatin_DMEM_earlyG1__grey = combinded_DMSO_Cisplatin_DMEM_earlyG1 %>% filter(phase == "Early G1")
combinded_DMSO_Cisplatin_DMEM_earlyG1__G2 = combinded_DMSO_Cisplatin_DMEM_earlyG1 %>% filter(phase == "G2")

t.test_G1_vergelijken_Cisplatin_EarlyG1 = print(pairwise.t.test(combinded_DMSO_Cisplatin_DMEM_earlyG1__G1 %>% pull (total), combinded_DMSO_Cisplatin_DMEM_earlyG1__G1%>% pull (condition),
                                                        p.adjust.method = "bonferroni"))

t.test_geel_vergelijken_Cisplatin_EarlyG1 = print(pairwise.t.test(combinded_DMSO_Cisplatin_DMEM_earlyG1__grey %>% pull (total), combinded_DMSO_Cisplatin_DMEM_earlyG1__grey%>% pull (condition),
                                                          p.adjust.method = "bonferroni"))
t.test_G2_vergelijken_Cisplatin_EerlyG1 = print(pairwise.t.test(combinded_DMSO_Cisplatin_DMEM_earlyG1__G2 %>% pull (total), combinded_DMSO_Cisplatin_DMEM_earlyG1__G2%>% pull (condition),
                                                        p.adjust.method = "bonferroni"))

#Vergelijken met DMSO en Etoposide en DMEM phase length with Early G1

combinded_DMSO_Etoposide_DMEM_earlyG1= rbind(count_corrected_filter_DMSO,count_corrected_filter_Etoposide,count_corrected_filter_DMEM)

combinded_DMSO_Etoposide_DMEM_earlyG1__G1 = combinded_DMSO_Etoposide_DMEM_earlyG1 %>% filter(phase == "G1")
combinded_DMSO_Etoposide_DMEM_earlyG1__grey = combinded_DMSO_Etoposide_DMEM_earlyG1 %>% filter(phase == "Early G1")
combinded_DMSO_Etoposide_DMEM_earlyG1__G2 = combinded_DMSO_Etoposide_DMEM_earlyG1 %>% filter(phase == "G2")

t.test_G1_vergelijken_Etoposide_EarlyG1 = print(pairwise.t.test(combinded_DMSO_Etoposide_DMEM_earlyG1__G1 %>% pull (total), combinded_DMSO_Etoposide_DMEM_earlyG1__G1%>% pull (condition),
                                                                p.adjust.method = "bonferroni"))

t.test_geel_vergelijken_Etoposide_EarlyG1 = print(pairwise.t.test(combinded_DMSO_Etoposide_DMEM_earlyG1__grey %>% pull (total), combinded_DMSO_Etoposide_DMEM_earlyG1__grey%>% pull (condition),
                                                                  p.adjust.method = "bonferroni"))
t.test_G2_vergelijken_Etoposide_EerlyG1 = print(pairwise.t.test(combinded_DMSO_Etoposide_DMEM_earlyG1__G2 %>% pull (total), combinded_DMSO_Etoposide_DMEM_earlyG1__G2%>% pull (condition),
                                                                p.adjust.method = "bonferroni"))

#count number of cells per image

#DMSO
cellen_tellen_DMSO= DMSO_tracks_data %>% select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number"))
cellen_tellen_DMSO_1 = cellen_tellen_DMSO %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
cellen_tellen_DMSO_2 = cellen_tellen_DMSO_1 %>% group_by(Image_Group_Index)%>%summarize(total=sum(n, na.rm = TRUE)/24)
cellen_tellen_DMSO_2 = cellen_tellen_DMSO_2 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)
cellen_tellen_DMSO_1 = cellen_tellen_DMSO_1 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)
cellen_tellen_DMSO_3 = cellen_tellen_DMSO_1 %>% mutate(n_normalized = ((n - min(cellen_tellen_DMSO_1$n, na.rm = T))/
                                                                 (max(cellen_tellen_DMSO_1$n, na.rm = T) - min(cellen_tellen_DMSO_1$n, na.rm = T))))
cellen_tellen_DMSO_4 = cellen_tellen_DMSO_3 %>% group_by(Image_Group_Index)%>%summarize(total_n_normalized=sum(n_normalized, na.rm = TRUE)/24)
cellen_tellen_DMSO_4 = cellen_tellen_DMSO_4 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)
cellen_tellen_DMSO_5= cellen_tellen_DMSO_2 %>% mutate(total_normalized = ((total - min(total, na.rm = T))/
                                                    (max(total, na.rm = T) - min(total, na.rm = T))))


plot_cellcount_DMSO= ggplot()+
  geom_line(data = cellen_tellen_DMSO_2, aes(x = time_h, y = total,))+
  ylab("Number of Cells")+xlab("Time (h)")+
theme_classic()

plot_cellcount_DMSO= ggplot()+
  geom_line(data = cellen_tellen_DMSO_1, aes(x = time_h, y = n))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_DMSO= ggplot()+
  geom_line(data = cellen_tellen_DMSO_4, aes(x = time_h, y = total_n_normalized))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_DMSO= ggplot()+
  geom_line(data = cellen_tellen_DMSO_5, aes(x = time_h, y = total_normalized))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()

#DMEM
cellen_tellen_DMEM= DMEM_tracks_data %>% select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number"))
cellen_tellen_DMEM_1 = cellen_tellen_DMEM %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
cellen_tellen_DMEM_2 = cellen_tellen_DMEM_1 %>% group_by(Image_Group_Index)%>%summarize(total=sum(n, na.rm = TRUE)/24)
cellen_tellen_DMEM_2 = cellen_tellen_DMEM_2 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)
cellen_tellen_DMEM_3 = cellen_tellen_DMEM_1 %>% mutate(n_normalized = ((n - min(n, na.rm = T))/
                                                                         (max(n, na.rm = T) - min(n, na.rm = T))))
cellen_tellen_DMEM_4 = cellen_tellen_DMEM_3 %>% group_by(Image_Group_Index)%>%summarize(total_n_normalized=sum(n_normalized, na.rm = TRUE)/24)
cellen_tellen_DMEM_4 = cellen_tellen_DMEM_4 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)
cellen_tellen_DMEM_5= cellen_tellen_DMEM_2 %>% mutate(total_normalized = ((total - min(total, na.rm = T))/
                                                                            (max(total, na.rm = T) - min(total, na.rm = T))))
cellen_tellen_DMEM_1 = cellen_tellen_DMEM_1 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)
cellen_tellen_DMEM_2 <- cellen_tellen_DMEM_2[-nrow(cellen_tellen_DMEM_2), ]
cellen_tellen_DMEM_4 <- cellen_tellen_DMEM_4[-nrow(cellen_tellen_DMEM_4), ]

cellen_tellen_DMEM_6 = cellen_tellen_DMEM_1%>%mutate(n_normalized = ((n - min(cellen_tellen_DMEM_1$n, na.rm = T))/
                                                                    (max(cellen_tellen_DMEM_1$n, na.rm = T) - min(cellen_tellen_DMEM_1$n, na.rm = T))))
cellen_tellen_DMEM_7 = cellen_tellen_DMEM_6 %>% group_by(Image_Group_Index)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/24)
cellen_tellen_DMEM_7 = cellen_tellen_DMEM_7 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP) - TIME_STEP)
cellen_tellen_DMEM_7 <- cellen_tellen_DMEM_7[-nrow(cellen_tellen_DMEM_7), ]

plot_cellcount_DMEM= ggplot()+
  geom_line(data = cellen_tellen_DMEM_2, aes(x = time_h, y = total,))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_DMEM= ggplot()+
  geom_line(data = cellen_tellen_DMEM_1, aes(x = time_h, y = n))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_DMSO= ggplot()+
  geom_line(data = cellen_tellen_DMEM_4, aes(x = time_h, y = total_n_normalized))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_DMSO= ggplot()+
  geom_line(data = cellen_tellen_DMEM_5, aes(x = time_h, y = total_normalized))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_DMSO= ggplot()+
  geom_line(data = cellen_tellen_DMEM_7, aes(x = time_h, y = total_normalized))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()


#cisplatin
cellen_tellen_Cisplatin= Cisplatin_tracks_data %>% select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number"))
cellen_tellen_Cisplatin_1 = cellen_tellen_Cisplatin %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)%>%ungroup()
Metedata_Cisplatin_2 = Metadata_Cisplatin%>% select(c("well_name","Image_Group_Number","well_name_p"))
cellen_tellen_Cisplatin_2<- left_join(cellen_tellen_Cisplatin_1,Metedata_Cisplatin_2,by="Image_Group_Number")
Cisplatin_layout_2= Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
cellen_tellen_Cisplatin_3<- left_join(cellen_tellen_Cisplatin_2,Cisplatin_layout_2,by="well_name")
cellen_tellen_Cisplatin_4= cellen_tellen_Cisplatin_3 %>% group_by(Image_Group_Index,condition)%>%summarize(total=sum(n, na.rm = TRUE)/8)
cellen_tellen_Cisplatin_5 = cellen_tellen_Cisplatin_4 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP)- TIME_STEP)
cellen_tellen_Cisplatin_6 = cellen_tellen_Cisplatin_5[-((nrow(cellen_tellen_Cisplatin_5)-1):nrow(cellen_tellen_Cisplatin_5)), ]
cellen_tellen_Cisplatin_7 = cellen_tellen_Cisplatin_3%>%  mutate(time_h = (Image_Group_Index * TIME_STEP)- TIME_STEP)
cellen_tellen_Cisplatin_8 = cellen_tellen_Cisplatin_7%>% filter(condition== "2500nM")
cellen_tellen_Cisplatin_9 = cellen_tellen_Cisplatin_7%>% filter(condition== "5000nM")
cellen_tellen_Cisplatin_10 = cellen_tellen_Cisplatin_7%>% filter(condition== "1000nM")
cellen_tellen_Cisplatin_11 = cellen_tellen_Cisplatin_7%>% filter(condition== "100000nM")

cellen_tellen_Cisplatin_12 = cellen_tellen_Cisplatin_1 %>% mutate(n_normalized = ((n - min(cellen_tellen_Cisplatin_1$n, na.rm = T))/
                                                                         (max(cellen_tellen_Cisplatin_1$n, na.rm = T) - min(cellen_tellen_Cisplatin_1$n, na.rm = T))))
Metedata_Cisplatin_2 = Metadata_Cisplatin%>% select(c("well_name","Image_Group_Number","well_name_p"))
cellen_tellen_Cisplatin_13<- left_join(cellen_tellen_Cisplatin_12,Metedata_Cisplatin_2,by="Image_Group_Number")
Cisplatin_layout_2= Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
cellen_tellen_Cisplatin_14<- left_join(cellen_tellen_Cisplatin_13,Cisplatin_layout_2,by="well_name")
cellen_tellen_Cisplatin_15= cellen_tellen_Cisplatin_14 %>% group_by(Image_Group_Index,condition)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/8)
cellen_tellen_Cisplatin_16 = cellen_tellen_Cisplatin_15 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP)- TIME_STEP)
cellen_tellen_Cisplatin_16 = cellen_tellen_Cisplatin_16[-((nrow(cellen_tellen_Cisplatin_16)-1):nrow(cellen_tellen_Cisplatin_16)), ]
cellen_tellen_Cisplatin_17 = cellen_tellen_Cisplatin_16%>% group_by(condition) %>% mutate(ma_total_normalized = rollmean(total_normalized, k=5, fill=NA, align='center'))
cellen_tellen_Cisplatin_18 = cellen_tellen_Cisplatin_6%>% group_by(condition) %>% mutate(ma_total_normalized = rollmean(total, k=5, fill=NA, align='center'))

print(max(cellen_tellen_Cisplatin_1$n, na.rm=TRUE))
print(min(cellen_tellen_Cisplatin_1$n, na.rm=TRUE))

plot_cellcount_Cisplatin= ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_6, aes(x = time_h, y = total,color= condition))+
  ylab("Number of Cells")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "red", "5000nM" = "green", "2500nM" = "blue","1000nM" = "black" ),
                     name = "Dose Cisplatin")+ 
  theme_classic()

plot_cellcount_Cisplatin_normalized =ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_16, aes(x = time_h, y = total_normalized,color= condition))+
  ylab("Normalized cell conuts")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "red", "5000nM" = "green", "2500nM" = "blue","1000nM" = "black" ),
                     name = "Dose Cisplatin")+ 
  theme_classic()

plot_cellcount_Cisplatin_normalized =ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_17, aes(x = time_h, y = ma_total_normalized,color= condition))+
  ylab("Normalized cell conuts")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "red", "5000nM" = "green", "2500nM" = "blue","1000nM" = "black" ),
                     name = "Dose Cisplatin")+ 
  theme_classic()

plot_cellcount_Cisplatin_normalized =ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_18, aes(x = time_h, y = ma_total_normalized,color= condition))+
  ylab("Number of cells")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "red", "5000nM" = "green", "2500nM" = "blue","1000nM" = "black" ),
                     name = "Dose Cisplatin")+ 
  theme_classic()



plot_cellcount_Cisplatin_2500nM= ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_8, aes(x = time_h, y = n,color= Image_Group_Number))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_Cisplatin_5000nM= ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_9, aes(x = time_h, y = n,color= Image_Group_Number))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_Cisplatin_1000nM= ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_10, aes(x = time_h, y = n,color= Image_Group_Number))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_Cisplatin_100000nM= ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_11, aes(x = time_h, y = n,color= Image_Group_Number))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

#etoposide
cellen_tellen_Etoposide= Etoposide_tracks_data %>% select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number"))
cellen_tellen_Etoposide_1 = cellen_tellen_Etoposide %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Metedata_Etoposide_2 = Metadata_Etoposide%>% select(c("well_name","Image_Group_Number","well_name_p"))
cellen_tellen_Etoposide_2<- left_join(cellen_tellen_Etoposide_1,Metedata_Etoposide_2,by="Image_Group_Number")
Etoposide_layout_2= Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
cellen_tellen_Etoposide_3<- left_join(cellen_tellen_Etoposide_2,Etoposide_layout_2,by="well_name")
cellen_tellen_Etoposide_4= cellen_tellen_Etoposide_3 %>% group_by(Image_Group_Index,condition)%>%summarize(total=sum(n, na.rm = TRUE)/8)
cellen_tellen_Etoposide_5 = cellen_tellen_Etoposide_4 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP)- TIME_STEP)
cellen_tellen_Etoposide_6 = cellen_tellen_Etoposide_5[-((nrow(cellen_tellen_Etoposide_5)-1):nrow(cellen_tellen_Etoposide_5)), ]
cellen_tellen_Etoposide_7 = cellen_tellen_Etoposide_3%>%  mutate(time_h = (Image_Group_Index * TIME_STEP)- TIME_STEP)
cellen_tellen_Etoposide_8 = cellen_tellen_Etoposide_7%>% filter(condition== "2500nM")
cellen_tellen_Etoposide_9 = cellen_tellen_Etoposide_7%>% filter(condition== "5000nM",Image_Group_Number==11)
cellen_tellen_Etoposide_10 = cellen_tellen_Etoposide_7%>% filter(condition== "1000nM")
cellen_tellen_Etoposide_11 = cellen_tellen_Etoposide_7%>% filter(condition== "10000nM",Image_Group_Number==16)
cellen_tellen_Etoposide_12 = cellen_tellen_Etoposide_6%>% group_by(condition) %>% mutate(ma_total = rollmean(total, k=5, fill=NA, align='center'))

cellen_tellen_Etoposide_13 = cellen_tellen_Etoposide_1 %>% mutate(n_normalized = ((n - min(cellen_tellen_Etoposide_1$n, na.rm = T))/
                                                                                    (max(cellen_tellen_Etoposide_1$n, na.rm = T) - min(cellen_tellen_Etoposide_1$n, na.rm = T))))
Metedata_Etoposide_2 = Metadata_Etoposide%>% select(c("well_name","Image_Group_Number","well_name_p"))
cellen_tellen_Etoposide_14<- left_join(cellen_tellen_Etoposide_13,Metedata_Etoposide_2,by="Image_Group_Number")
Etoposide_layout_2= Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
cellen_tellen_Etoposide_15<- left_join(cellen_tellen_Etoposide_14,Etoposide_layout_2,by="well_name")
cellen_tellen_Etoposide_16= cellen_tellen_Cisplatin_15 %>% group_by(Image_Group_Index,condition)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/8)
cellen_tellen_Etoposide_17 = cellen_tellen_Cisplatin_16 %>%  mutate(time_h = (Image_Group_Index * TIME_STEP)- TIME_STEP)
cellen_tellen_Etoposide_18 = cellen_tellen_Etoposide_17[-((nrow(cellen_tellen_Etoposide_17)-1):nrow(cellen_tellen_Etoposide_17)), ]
cellen_tellen_Etoposde_18 = cellen_tellen_Etoposide_18%>% group_by(condition) %>% mutate(ma_total_normalized = rollmean(total_normalized, k=5, fill=NA, align='center'))

plot_cellcount_Eroposide_normal= ggplot()+
  geom_line(data = cellen_tellen_Etoposide_6, aes(x = time_h, y = total,color= condition))+
  ylab("Number of Cells")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "red", "5000nM" = "green", "2500nM" = "blue","1000nM" = "black" ),
                     name = "Dose Etoposide")+ 
  theme_classic()

plot_cellcount_Eroposide_normal_rolling_mean= ggplot()+
  geom_line(data = cellen_tellen_Etoposde_18, aes(x = time_h, y =total_normalized,color= condition))+
  ylab("Number of Cells")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "red", "5000nM" = "green", "2500nM" = "blue","1000nM" = "black" ),
                     name = "Dose Etoposide")+ 
  theme_classic()

plot_cellcount_Eroposide_normal_rolling_mean= ggplot()+
  geom_line(data = cellen_tellen_Etoposde_18, aes(x = time_h, y = ma_total_normalized,color= condition))+
  ylab("Normalized cell conut")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "red", "5000nM" = "green", "2500nM" = "blue","1000nM" = "black" ),
                     name = "Dose Etoposide")+ 
  theme_classic()

plot_cellcount_Etoposide_2500nM= ggplot()+
  geom_line(data = cellen_tellen_Etoposide_8, aes(x = time_h, y = n,color= Image_Group_Number))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_Cisplatin_5000nM= ggplot()+
  geom_line(data = cellen_tellen_Etoposide_9, aes(x = time_h, y = n,color= Image_Group_Number))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_Cisplatin_1000nM= ggplot()+
  geom_line(data = cellen_tellen_Etoposide_10, aes(x = time_h, y = n,color= Image_Group_Number))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()

plot_cellcount_Cisplatin_10000nM= ggplot()+
  geom_line(data = cellen_tellen_Cisplatin_11, aes(x = time_h, y = n,color= Image_Group_Number))+
  ylab("Number of Cells")+xlab("Time (h)")+
  theme_classic()


# DMSO en Cisplatin

view(cellen_tellen_DMSO_3)
view(cellen_tellen_Cisplatin_12)
view(cellen_tellen_DMEM_3)

stap1= cellen_tellen_DMSO_3 %>% select(-c(n_normalized,time_h))
stap1_2 =cellen_tellen_DMEM_3%>% select(-c(n_normalized,time_h))
stap2= cellen_tellen_Cisplatin_12 %>% select(-c(n_normalized))
stap3= left_join(stap2,Metedata_Cisplatin_2,by="Image_Group_Number")
stap4 = left_join(stap3,Cisplatin_layout_2,by="well_name")
stap5 = DMSO_layout %>% select(c("well_name","condition","dose_uM"))
stap6 = Metadata_DMSO%>% select(c("well_name","Image_Group_Number","well_name_p"))
stap7 = left_join(stap1,stap6,by="Image_Group_Number")
stap8 = left_join(stap7,stap5,by="well_name")
stap1_2_3= DMEM_layout %>% select(c("well_name","condition","dose_uM"))
stap1_2_3_4= Metadata_DMEM%>% select(c("well_name","Image_Group_Number","well_name_p"))
stap1_2_3_4_5 = left_join(stap1_2,stap1_2_3_4,by="Image_Group_Number")
stap1_2_3_4_5_6 = left_join(stap1_2_3_4_5 ,stap1_2_3,by="well_name")
stap9 = full_join(stap4,stap8)
stap9= full_join(stap9,stap1_2_3_4_5_6)
stap10= stap9 %>%  mutate(n_normalized = ((n - min(stap9$n, na.rm = T))/
                                     (max(stap9$n, na.rm = T) - min(stap9$n, na.rm = T))))

stap11= stap10 %>% filter(condition=="DMSO") %>% group_by(Image_Group_Index,condition)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/24)
stap11_1= stap10 %>% filter(condition=="DMEM") %>% group_by(Image_Group_Index,condition)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/24)
stap12= stap10 %>% filter(!condition=="DMSO") %>% filter(!condition=="DMEM")%>% group_by(Image_Group_Index,condition)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/8)
stap13= stap12[-((nrow(stap12)-1):nrow(stap12)), ]
stap11_1_2= stap11_1[-nrow(stap11_1), ]
stap14= full_join(stap11,stap13)
stap14=full_join(stap14,stap11_1_2)
stap15= stap14%>% group_by(condition) %>% mutate(ma_total_normalized = rollmean(total_normalized, k=5, fill=NA, align='center'))
stap16= stap15%>%mutate(time_h = (Image_Group_Index * TIME_STEP)- TIME_STEP)
stap17= ggplot()+
  geom_line(data = stap16, aes(x = time_h, y =ma_total_normalized,color= condition))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+ 
  theme_classic()

filter1000 = stap16%>% filter(condition== "1000nM")
filter2500 = stap16%>% filter(condition== "2500nM")
filter5000 = stap16%>% filter(condition== "5000nM")
filter100000 = stap16%>% filter(condition== "100000nM")

DMSO_stijging=print(max(stap11$total_normalized, na.rm=TRUE))-print(min(stap11$total_normalized, na.rm=TRUE))#0.3446078
Cisplatin1000_stijging=print(max(filter1000$total_normalized, na.rm=TRUE))-print(min(filter1000$total_normalized, na.rm=TRUE))#0.3985294
Cisplatin2500_stijging=print(max(filter2500$total_normalized, na.rm=TRUE))-print(min(filter2500$total_normalized, na.rm=TRUE))#0.2492647
Cisplatin5000_stijging=print(max(filter5000$total_normalized, na.rm=TRUE))-print(min(filter5000$total_normalized, na.rm=TRUE))#0.1757353
Cisplatin100000_stijging=print(max(filter100000$total_normalized, na.rm=TRUE))-print(min(filter100000$total_normalized, na.rm=TRUE))#0.05367647
DMEM_stijging= print(max(stap11_1_2$total_normalized, na.rm=TRUE))-print(min(stap11_1_2$total_normalized, na.rm=TRUE))#0.3215686

stap18=filter(stap16,!is.na(ma_total_normalized))
stap18$time_h <- as.factor(stap18$time_h)
print(pairwise.t.test(stap18 %>% pull (ma_total_normalized), stap18 %>% pull (condition,time_h),
                      p.adjust.method = "bonferroni"))

# DMSO,DMEM en Etoposide

view(cellen_tellen_Etoposide_3)

stap1_Etoposide= full_join(stap1_2_3_4_5_6,cellen_tellen_Etoposide_3)
stap2_Etoposide= full_join(stap1_Etoposide,stap8)
stap3_Etoposide= stap2_Etoposide%>%  mutate(n_normalized = ((n - min(stap2_Etoposide$n, na.rm = T))/
                                               (max(stap2_Etoposide$n, na.rm = T) - min(stap2_Etoposide$n, na.rm = T))))

stap4_Etoposide= stap3_Etoposide %>% filter(condition=="DMSO") %>% group_by(Image_Group_Index,condition)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/24)
stap5_Etoposide= stap3_Etoposide %>% filter(condition=="DMEM") %>% group_by(Image_Group_Index,condition)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/24)
stap6_Etoposide= stap3_Etoposide %>% filter(!condition=="DMSO") %>% filter(!condition=="DMEM")%>% group_by(Image_Group_Index,condition)%>%summarize(total_normalized=sum(n_normalized, na.rm = TRUE)/8)
stap7_Etoposide= stap6_Etoposide[-((nrow(stap6_Etoposide)-1):nrow(stap6_Etoposide)), ]
stap8_Etoposide= stap5_Etoposide[-nrow(stap5_Etoposide), ]
stap9_Etoposide= full_join(stap8_Etoposide,stap7_Etoposide)
stap10_Etoposide=full_join(stap9_Etoposide,stap4_Etoposide)
stap11_Etoposide=stap10_Etoposide %>%group_by(condition) %>% mutate(ma_total_normalized = rollmean(total_normalized, k=5, fill=NA, align='center'))
stap12_Etoposide= stap11_Etoposide%>%mutate(time_h = (Image_Group_Index * TIME_STEP)- TIME_STEP)

stap13_Etoposide= ggplot()+
  geom_line(data = stap12_Etoposide, aes(x = time_h, y =ma_total_normalized,color= condition))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+ 
  theme_classic()

filter1000_E =stap12_Etoposide %>% filter(condition== "1000nM")
filter2500_E = stap12_Etoposide%>% filter(condition== "2500nM")
filter5000_E = stap12_Etoposide%>% filter(condition== "5000nM")
filter10000_E = stap12_Etoposide%>% filter(condition== "10000nM")

DMSO_stijging_E=print(max(stap4_Etoposide$total_normalized, na.rm=TRUE))-print(min(stap4_Etoposide$total_normalized, na.rm=TRUE))#0.3446078
E1000_stijging=print(max(filter1000_E$total_normalized, na.rm=TRUE))-print(min(filter1000_E$total_normalized, na.rm=TRUE))#0.2529412
E2500_stijging=print(max(filter2500_E$total_normalized, na.rm=TRUE))-print(min(filter2500_E$total_normalized, na.rm=TRUE))#0.2544118
E5000_stijging=print(max(filter5000_E$total_normalized, na.rm=TRUE))-print(min(filter5000_E$total_normalized, na.rm=TRUE))#0.2683824
E100000_stijging=print(max(filter10000_E$total_normalized, na.rm=TRUE))-print(min(filter10000_E$total_normalized, na.rm=TRUE))#0.1963235
DMEM_stijging_E= print(max(stap11_1_2$total_normalized, na.rm=TRUE))-print(min(stap11_1_2$total_normalized, na.rm=TRUE))#0.3215686

stap13_Etoposide=filter(stap12_Etoposide,!is.na(ma_total_normalized))
stap13_Etoposide$time_h <- as.factor(stap13_Etoposide$time_h)
print(pairwise.t.test(stap13_Etoposide %>% pull (ma_total_normalized), stap13_Etoposide %>% pull (condition,time_h),
                      p.adjust.method = "bonferroni"))




#normalizeren cell count van 1 DMSO

view(cellen_tellen_DMSO)
cellen_tellen_DMSO_1 = cellen_tellen_DMSO %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
DMSOfilter1= cellen_tellen_DMSO_1 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
DMSO1norm= cellen_tellen_DMSO_1%>% left_join(DMSOfilter1 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(n_norm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_DMSO = DMSO_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_DMSO = Metadata_DMSO%>% select(c("well_name","Image_Group_Number","well_name_p"))
DMSO1norm= left_join(DMSO1norm,Wellnames_DMSO,by="Image_Group_Number")
DMSO1norm=left_join(DMSO1norm,conditions_DMSO,by="well_name")
DMSO1norm_total=DMSO1norm %>% group_by(condition,Image_Group_Index.x)%>%summarize(total_nnorm=sum(n_norm, na.rm = TRUE)/24)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)%>%mutate(ma_total_nnorm = rollmean(total_nnorm, k=2, fill=NA, align='center'))
DMSO1norm_total=DMSO1norm_total%>% relocate(c("Image_Group_Index.x","condition"))

DMSO1norm_total_ggplot= ggplot()+
  geom_line(data =DMSO1norm_total , aes(x = time_h, y =ma_total_nnorm))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()

#normalizeren cell count van 1 DMEM
view(cellen_tellen_DMEM)
cellen_tellen_DMEM_1 = cellen_tellen_DMEM %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
DMEMfilter1= cellen_tellen_DMEM_1 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
DMEM1norm= cellen_tellen_DMEM_1%>% left_join(DMEMfilter1 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(n_norm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)%>%filter(!Image_Group_Index.x==41)
conditions_DMEM = DMEM_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_DMEM = Metadata_DMEM%>% select(c("well_name","Image_Group_Number","well_name_p"))
DMEM1norm= left_join(DMEM1norm,Wellnames_DMEM,by="Image_Group_Number")
DMEM1norm=left_join(DMEM1norm,conditions_DMEM,by="well_name")
DMEM1norm_total=DMEM1norm %>% group_by(condition,Image_Group_Index.x)%>%summarize(total_nnorm=sum(n_norm, na.rm = TRUE)/24)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)%>%mutate(ma_total_nnorm = rollmean(total_nnorm, k=2, fill=NA, align='center'))
DMEM1norm_total=DMEM1norm_total%>% relocate(c("Image_Group_Index.x","condition"))


DMEM1norm_total_ggplot= ggplot()+
  geom_line(data =DMEM1norm_total , aes(x = time_h, y =ma_total_nnorm))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()

#normalizeren cell count van 1 cisplatin
view(cellen_tellen_Cisplatin)
cellen_tellen_Cisplatin_1 = cellen_tellen_Cisplatin %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Cisplatinfilter1= cellen_tellen_Cisplatin_1 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
Cisplatin1norm= cellen_tellen_Cisplatin_1%>% left_join(Cisplatinfilter1 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(n_norm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)%>%filter(!Image_Group_Index.x==41)
conditions_cisplatin = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_cisplatin = Metadata_Cisplatin%>% select(c("well_name","Image_Group_Number","well_name_p"))
Cisplatin1norm= left_join(Cisplatin1norm,Wellnames_cisplatin,by="Image_Group_Number")
Cisplatin1norm=left_join(Cisplatin1norm,conditions_cisplatin,by="well_name")
Cipslatin1norm_5000= Cisplatin1norm%>%filter(condition=="5000nM")%>%filter(!well_name_p=="M06_2"&!well_name_p=="M11_1")
Cisplatin1norm= Cisplatin1norm%>%filter(!condition=="5000nM")
Cisplatin1norm= Cisplatin1norm%>%group_by(Image_Group_Index.x,condition)%>%summarize(total_nnorm=sum(n_norm, na.rm = TRUE)/8)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Cisplatin1norm_5000= Cipslatin1norm_5000%>%group_by(Image_Group_Index.x,condition)%>%summarize(total_nnorm=sum(n_norm, na.rm = TRUE)/6)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Cisplatin1norm=  full_join(Cisplatin1norm,Cisplatin1norm_5000)
Cisplatin1nor_ma= Cisplatin1norm%>%group_by(condition)%>% mutate(ma_total_nnorm = rollmean(total_nnorm, k=2, fill=NA, align='center'))
Cisplatin_DMSO= full_join(Cisplatin1nor_ma,DMSO1norm_total)
Cisplatin_DMSO_DMEM= full_join(Cisplatin_DMSO,DMEM1norm_total)


Cisplatin1norm_total_ggplot= ggplot()+
  geom_line(data = Cisplatin_DMSO_DMEM , aes(x = time_h, y =ma_total_nnorm,color= condition),size = 1.1)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Cisplatin_DMSO_DMEM,aes(x = time_h, ymax = ma_total_nnorm+0.2, ymin = ma_total_nnorm-0.2, color = condition,fill=condition), alpha = 0.2)+
  scale_fill_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()

#etoposide

view(cellen_tellen_Etoposide)
cellen_tellen_Etoposide_1 = cellen_tellen_Etoposide %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Etoposidefilter1= cellen_tellen_Etoposide_1 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
Etoposide1norm= cellen_tellen_Etoposide_1%>% left_join(Etoposidefilter1 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(n_norm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)%>%filter(!Image_Group_Index.x==41)
conditions_Etoposide = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_Etoposide = Metadata_Etoposide%>% select(c("well_name","Image_Group_Number","well_name_p"))
Etoposide1norm= left_join(Etoposide1norm,Wellnames_Etoposide,by="Image_Group_Number")
Etoposide1norm=left_join(Etoposide1norm,conditions_Etoposide,by="well_name")
Etoposide1norm_10000= Etoposide1norm%>%filter(condition=="10000nM")%>%filter(!well_name_p=="C06_2")
Etoposide1norm= Etoposide1norm%>%filter(!condition=="10000nM")
Etoposide1norm= Etoposide1norm%>%group_by(Image_Group_Index.x,condition)%>%summarize(total_nnorm=sum(n_norm, na.rm = TRUE)/8)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Etoposide1norm_10000= Etoposide1norm_10000%>%group_by(Image_Group_Index.x,condition)%>%summarize(total_nnorm=sum(n_norm, na.rm = TRUE)/7)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Etoposide1norm=  full_join(Etoposide1norm,Etoposide1norm_10000)
Etoposide1nor_ma= Etoposide1norm%>%group_by(condition)%>% mutate(ma_total_nnorm = rollmean(total_nnorm, k=2, fill=NA, align='center'))
Etoposide_DMSO= full_join(Etoposide1nor_ma,DMSO1norm_total)
Etoposide_DMSO_DMEM= full_join(Etoposide_DMSO,DMEM1norm_total)

#Etoposide1norm= Etoposide1norm%>%group_by(Image_Group_Index.x,condition)%>%summarize(total_nnorm=sum(n_norm, na.rm = TRUE)/8)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
#Etoposide1norm_ma= Etoposide1norm%>%group_by(condition)%>% mutate(ma_total_nnorm = rollmean(total_nnorm, k=2, fill=NA, align='center'))

Etoposide1norm_total_ggplot= ggplot()+
  geom_line(data = Etoposide_DMSO_DMEM , aes(x = time_h, y =ma_total_nnorm,color= condition))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Etoposide_DMSO_DMEM,aes(x = time_h, ymax = ma_total_nnorm+0.1, ymin = ma_total_nnorm-0.1, color = condition,fill=condition), alpha = 0.2)+
  scale_fill_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()

#opnieuw analyse op cell population niveau
#Etoposide
cellen_tellen_Etoposide= Etoposide_tracks_data %>% select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number"))
cellen_tellen_Etoposide_1 = cellen_tellen_Etoposide %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Etoposidefilter1= cellen_tellen_Etoposide_1 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
Etoposide1norm= cellen_tellen_Etoposide_1%>% left_join(Etoposidefilter1 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(n_norm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)%>%filter(!Image_Group_Index.x==41)
conditions_Etoposide = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_Etoposide = Metadata_Etoposide%>% select(c("well_name","loc","well_name_p"))
Wellnames_Etoposide <- Wellnames_Etoposide  %>% dplyr::rename(Image_Group_Number=loc)
Etoposide1norm= left_join(Etoposide1norm,Wellnames_Etoposide,by="Image_Group_Number")
Etoposide1norm=left_join(Etoposide1norm,conditions_Etoposide,by="well_name")
view(Etoposide1norm)

Etoposide1norm_10000nM= Etoposide1norm%>% filter(condition=="10000nM")
Etoposide1norm_10000nM_filter= Etoposide1norm_10000nM %>% filter(!well_name_p=="B06_2",!well_name_p=="B11_2",!well_name_p=="C06_2")# minder dan 15 cellen bij tijdspunt 0
Etoposide1normmean_sd2= Etoposide1norm_10000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
Etoposide1normmean_sd3= Etoposide1norm_10000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
Etoposide1normmmean_sd4= left_join(Etoposide1normmean_sd2,Etoposide1normmean_sd3,by = c("well_name", "Image_Group_Index.x"))
Etoposide1normmmean_sd4=Etoposide1normmmean_sd4%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Etoposide1normmmean_sd4=left_join(Etoposide1normmmean_sd4,conditions_Etoposide, by= "well_name")
Etoposide1normmean_sd5= Etoposide1normmmean_sd4  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Etoposdie1normmean_sd6= Etoposide1normmmean_sd4  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Etoposide1normmean_sd7= left_join(Etoposide1normmean_sd5,Etoposdie1normmean_sd6,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

Etoposide1norm_5000nM= Etoposide1norm%>% filter(condition=="5000nM")
Etoposide1norm_5000nM_filter= Etoposide1norm_5000nM %>% filter(!well_name_p=="D06_2",!well_name_p=="E06_2")# minder dan 15 cellen bij tijdspunt 0
Etoposide1normmean_sd2_5000= Etoposide1norm_5000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
Etoposide1normmean_sd3_5000= Etoposide1norm_5000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
Etoposide1normmmean_sd4_5000= left_join(Etoposide1normmean_sd2_5000,Etoposide1normmean_sd3_5000,by = c("well_name", "Image_Group_Index.x"))
Etoposide1normmmean_sd4_5000=Etoposide1normmmean_sd4_5000%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Etoposide1normmmean_sd4_5000=left_join(Etoposide1normmmean_sd4_5000,conditions_Etoposide, by= "well_name")
Etoposide1normmean_sd5_5000= Etoposide1normmmean_sd4_5000  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Etoposdie1normmean_sd6_5000= Etoposide1normmmean_sd4_5000  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Etoposide1normmean_sd7_5000= left_join(Etoposide1normmean_sd5_5000,Etoposdie1normmean_sd6_5000,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

Etoposide1norm_2500nM= Etoposide1norm%>% filter(condition=="2500nM")
Etoposide1normmean_sd2_2500= Etoposide1norm_2500nM  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
Etoposide1normmean_sd3_2500= Etoposide1norm_2500nM  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
Etoposide1normmmean_sd4_2500= left_join(Etoposide1normmean_sd2_2500,Etoposide1normmean_sd3_2500,by = c("well_name", "Image_Group_Index.x"))
Etoposide1normmmean_sd4_2500=Etoposide1normmmean_sd4_2500%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Etoposide1normmmean_sd4_2500=left_join(Etoposide1normmmean_sd4_2500,conditions_Etoposide, by= "well_name")
Etoposide1normmean_sd5_2500= Etoposide1normmmean_sd4_2500  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Etoposdie1normmean_sd6_2500= Etoposide1normmmean_sd4_2500  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Etoposide1normmean_sd7_2500= left_join(Etoposide1normmean_sd5_2500,Etoposdie1normmean_sd6_2500,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

Etoposide1norm_1000nM= Etoposide1norm%>% filter(condition=="1000nM")
Etoposide1normmean_sd2_1000= Etoposide1norm_1000nM  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
Etoposide1normmean_sd3_1000= Etoposide1norm_1000nM  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
Etoposide1normmmean_sd4_1000= left_join(Etoposide1normmean_sd2_1000,Etoposide1normmean_sd3_1000,by = c("well_name", "Image_Group_Index.x"))
Etoposide1normmmean_sd4_1000=Etoposide1normmmean_sd4_1000%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Etoposide1normmmean_sd4_1000=left_join(Etoposide1normmmean_sd4_1000,conditions_Etoposide, by= "well_name")
Etoposide1normmean_sd5_1000= Etoposide1normmmean_sd4_1000  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Etoposdie1normmean_sd6_1000= Etoposide1normmmean_sd4_1000  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Etoposide1normmean_sd7_1000= left_join(Etoposide1normmean_sd5_1000,Etoposdie1normmean_sd6_1000,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

fifa=ggplot()+
  geom_line(data = Etoposide1normmmean_sd4_1000 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("H06" = "darkred", "H11" = "red", "I06" = "blue","I11" = "green"),
                    name = "Well_name")+
  theme_classic()+
  ggtitle("1 µM Etoposide")
  geom_ribbon(data = Etoposide1normmmean_sd4_1000,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("H06" = "darkred", "H11" = "red", "I06" = "blue","I11" = "green" ))+
  guides(fill = "none")+

 

fifa=ggplot()+
  geom_line(data = Etoposide1normmmean_sd4_5000 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("D06" = "darkred", "D11" = "red", "E06" = "blue","E11" = "green"),
                     name = "Well_name")+
  theme_classic()+
  ggtitle("5 µM Etoposide")
  geom_ribbon(data = Etoposide1normmmean_sd4_5000,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("D06" = "darkred", "D11" = "red", "E06" = "blue","E11" = "green"))+
  guides(fill = "none")


fifa=ggplot()+
  geom_line(data = Etoposide1normmmean_sd4_2500 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("F06" = "darkred", "F11" = "red", "G06" = "blue","G11" = "green"),
                     name = "Well_name")+
  theme_classic()+
  ggtitle("2.5 µM Etoposide")
  geom_ribbon(data = Etoposide1normmmean_sd4_2500,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("F06" = "darkred", "F11" = "red", "G06" = "blue","G11" = "green"))+
  guides(fill = "none")+


fifa=ggplot()+
  geom_line(data = Etoposide1normmmean_sd4 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("B06" = "darkred", "B11" = "red", "C06" = "blue","C11" = "green"),
                     name = "Well_name")+
  theme_classic()+
  ggtitle("10 µM Etoposide")
  geom_ribbon(data = Etoposide1normmmean_sd4,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("B06" = "darkred", "B11" = "red", "C06" = "blue","C11" = "green"))+
  guides(fill = "none")



#DMSO
cellen_tellen_DMSO= DMSO_tracks_data %>% select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number"))
cellen_tellen_DMSO_1 = cellen_tellen_DMSO %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
DMSOfilter1= cellen_tellen_DMSO_1 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
DMSO1norm= cellen_tellen_DMSO_1%>% left_join(DMSOfilter1 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(n_norm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_DMSO = DMSO_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_DMSO = Metadata_DMSO%>% select(c("well_name","loc","well_name_p"))
Wellnames_DMSO <- Wellnames_DMSO  %>% dplyr::rename(Image_Group_Number=loc)
DMSO1norm= left_join(DMSO1norm,Wellnames_DMSO,by="Image_Group_Number")
DMSO1norm=left_join(DMSO1norm,conditions_DMSO,by="well_name")
view(DMSO1norm)
DMSO1norm_filter= DMSO1norm %>% filter(!well_name_p=="H02_1",!well_name_p=="H02_2",!well_name_p=="I02_1",!well_name_p=="I02_2",!well_name_p=="J02_1",!well_name_p=="J02_2",!well_name_p=="K07_2",!well_name_p=="L02_1",!well_name_p=="L02_2",!well_name_p=="M02_1",!well_name_p=="M02_2",!well_name_p=="K02_1",!well_name_p=="K02_2",!well_name_p=="L07_1")# minder dan 15 cellen bij tijdspunt 0
#DMSO1normmean_sd =DMSO1norm  %>% group_by(condition,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))%>%summarize(sd_nnorm=sd(Mean_nnorm))
DMSO1normmean_sd2= DMSO1norm_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
DMSO1normmean_sd3= DMSO1norm_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
DMSO1normmean_sd4= left_join(DMSO1normmean_sd2,DMSO1normmean_sd3,by = c("well_name", "Image_Group_Index.x"))
DMSO1normmean_sd4=DMSO1normmean_sd4%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_DMSO = DMSO_layout %>% select(c("well_name","condition","dose_uM"))
DMSO1normmean_sd4=left_join(DMSO1normmean_sd4,conditions_DMSO, by= "well_name")
DMSO1normmean_sd5= DMSO1normmean_sd4  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
DMSO1normmean_sd6= DMSO1normmean_sd4  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
DMSO1normmean_sd7= left_join(DMSO1normmean_sd5,DMSO1normmean_sd6,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

fifa=ggplot()+
  geom_line(data = DMSO1normmean_sd4 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("H07" = "darkred", "I07" = "red", "J07" = "blue","K07" = "green","L07"= "orange","M07"="black"),
                     name = "Well_name")+
  theme_classic()+
  ggtitle("DMSO")+
  geom_ribbon(data = DMSO1normmean_sd4,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("H07" = "darkred", "I07" = "red", "J07" = "blue","K07" = "green","L07"= "orange","M07"="black"))+
  guides(fill = "none")+
  





fifa=ggplot()+
  geom_line(data = DMSO1norm , aes(x = time_h, y =n_norm,color= well_name_p))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic() 
fifa=ggplot()+
  geom_line(data = DMSO1norm_filter , aes(x = time_h, y =n_norm,color= well_name_p))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic() 
fifa=ggplot()+
  geom_line(data =DMSO1normmean_sd4 , aes(x = Image_Group_Index.x, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()

fifa=ggplot()+
  geom_line(data =DMSO1normmean_sd7 , aes(x = time_h, y =Mean_nnorm,color= condition))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = DMSO1normmean_sd7,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = condition,fill=condition),alpha=0.5)+
  scale_fill_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()

#DMEM
cellen_tellen_DMEM= DMEM_tracks_data %>% select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number"))
cellen_tellen_DMEM_1 = cellen_tellen_DMEM %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
DMEMfilter1= cellen_tellen_DMEM_1 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
DMEM1norm= cellen_tellen_DMEM_1%>% left_join(DMEMfilter1 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(n_norm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)%>%filter(!Image_Group_Index.x==41)
conditions_DMEM = DMEM_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_DMEM = Metadata_DMEM%>% select(c("well_name","loc","well_name_p"))
Wellnames_DMEM <- Wellnames_DMEM  %>% dplyr::rename(Image_Group_Number=loc)
DMEM1norm= left_join(DMEM1norm,Wellnames_DMEM,by="Image_Group_Number")
DMEM1norm=left_join(DMEM1norm,conditions_DMEM,by="well_name")
view(DMEM1norm)
DMEM1norm_filter= DMEM1norm %>% filter(!well_name_p=="B02_1",!well_name_p=="C02_1",!well_name_p=="D02_1",!well_name_p=="D02_2",!well_name_p=="E02_1",!well_name_p=="F02_1",!well_name_p=="F07_2",!well_name_p=="G02_1",!well_name_p=="G02_2",!well_name_p=="F02_2")# minder dan 15 cellen bij tijdspunt 0
DMEM1normmean_sd2= DMEM1norm_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
DMEM1normmean_sd3= DMEM1norm_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
DMEM1normmean_sd4= left_join(DMEM1normmean_sd2,DMEM1normmean_sd3,by = c("well_name", "Image_Group_Index.x"))
DMEM1normmean_sd4=DMEM1normmean_sd4%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_DMEM = DMEM_layout %>% select(c("well_name","condition","dose_uM"))
DMEM1normmean_sd4=left_join(DMEM1normmean_sd4,conditions_DMEM, by= "well_name")
DMEM1normmean_sd5= DMEM1normmean_sd4  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
DMEM1normmean_sd6= DMEM1normmean_sd4  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
DMEM1normmean_sd7= left_join(DMEM1normmean_sd5,DMEM1normmean_sd6,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

fifa=ggplot()+
  geom_line(data = DMEM1normmean_sd4 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("B02" = "darkred", "B07" = "red", "C02" = "blue","C07" = "green","D07"= "orange","E02"="black","E07"="purple","F07"="brown","G07"="maroon"),
                     name = "Well_name")+
  ggtitle("DMEM")+
  theme_classic()+
  geom_ribbon(data = DMEM1normmean_sd4,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("B02" = "darkred", "B07" = "red", "C02" = "blue","C07" = "green","D07"= "orange","E02"="black","E07"="purple","F07"="brown","G07"="maroon"))+
  guides(fill = "none")
  



fifacis=ggplot()+
  geom_line(data = DMEM1norm_filter , aes(x = time_h, y =n_norm,color= well_name_p))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic()

#cisplatin
cellen_tellen_Cisplatin= Cisplatin_tracks_data %>% select(c("ImageNumber","Image_Group_Number","Image_Group_Index","Nuclei_Number_Object_Number"))
cellen_tellen_Cisplatin_1 = cellen_tellen_Cisplatin %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Cisplatinfilter1= cellen_tellen_Cisplatin_1 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
Cisplatin1norm= cellen_tellen_Cisplatin_1%>% left_join(Cisplatinfilter1 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(n_norm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)%>%filter(!Image_Group_Index.x==41)
conditions_cisplatin = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_cisplatin = Metadata_Cisplatin%>% select(c("well_name","loc","well_name_p"))
Wellnames_cisplatin <- Wellnames_cisplatin  %>% dplyr::rename(Image_Group_Number=loc)
Cisplatin1norm= left_join(Cisplatin1norm,Wellnames_cisplatin,by="Image_Group_Number")
Cisplatin1norm=left_join(Cisplatin1norm,conditions_cisplatin,by="well_name")
view(Cisplatin1norm)

Cisplatin1norm_100000nM= Cisplatin1norm%>% filter(condition=="100000nM")
Cisplatin1norm_100000nM_filter= Cisplatin1norm_100000nM %>% filter(!well_name_p=="J06_2",!well_name_p=="K06_1")# minder dan 15 cellen bij tijdspunt 0
Cisplatin1normmean_sd2= Cisplatin1norm_100000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
Cisplatin1normmean_sd3= Cisplatin1norm_100000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
Cisplatin1normmmean_sd4= left_join(Cisplatin1normmean_sd2,Cisplatin1normmean_sd3,by = c("well_name", "Image_Group_Index.x"))
Cisplatin1normmmean_sd4=Cisplatin1normmmean_sd4%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_cisplatin = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Cisplatin1normmmean_sd4=left_join(Cisplatin1normmmean_sd4,conditions_cisplatin, by= "well_name")
Cisplatin1normmean_sd5= Cisplatin1normmmean_sd4  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Cisplatin1normmean_sd6= Cisplatin1normmmean_sd4  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Cisplatin1normmean_sd7= left_join(Cisplatin1normmean_sd5,Cisplatin1normmean_sd6,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

Cisplatin1norm_5000nM= Cisplatin1norm%>% filter(condition=="5000nM")
Cisplatin1norm_5000nM_filter= Cisplatin1norm_5000nM %>% filter(!well_name_p=="M06_2",!well_name_p=="M11_1")# minder dan 15 cellen bij tijdspunt 0
Cisplatin1normmean_sd2_5000= Cisplatin1norm_5000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
Cisplatin1normmean_sd3_5000= Cisplatin1norm_5000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
Cisplatin1normmmean_sd4_5000= left_join(Cisplatin1normmean_sd2_5000,Cisplatin1normmean_sd3_5000,by = c("well_name", "Image_Group_Index.x"))
Cisplatin1normmmean_sd4_5000=Cisplatin1normmmean_sd4_5000%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_cisplatin = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Cisplatin1normmmean_sd4_5000=left_join(Cisplatin1normmmean_sd4_5000,conditions_cisplatin, by= "well_name")
Cisplatin1normmean_sd5_5000= Cisplatin1normmmean_sd4_5000  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Cisplatin1normmean_sd6_5000= Cisplatin1normmmean_sd4_5000  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Cisplatin1normmean_sd7_5000= left_join(Cisplatin1normmean_sd5_5000,Cisplatin1normmean_sd6_5000,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

Cisplatin1norm_2500nM= Cisplatin1norm%>% filter(condition=="2500nM")
Cisplatin1norm_2500nM_filter= Cisplatin1norm_2500nM # minder dan 15 cellen bij tijdspunt 0
Cisplatin1normmean_sd2_2500= Cisplatin1norm_2500nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
Cisplatin1normmean_sd3_2500= Cisplatin1norm_2500nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
Cisplatin1normmmean_sd4_2500= left_join(Cisplatin1normmean_sd2_2500,Cisplatin1normmean_sd3_2500,by = c("well_name", "Image_Group_Index.x"))
Cisplatin1normmmean_sd4_2500=Cisplatin1normmmean_sd4_2500%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_cisplatin = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Cisplatin1normmmean_sd4_2500=left_join(Cisplatin1normmmean_sd4_2500,conditions_cisplatin, by= "well_name")
Cisplatin1normmean_sd5_2500= Cisplatin1normmmean_sd4_2500  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Cisplatin1normmean_sd6_2500= Cisplatin1normmmean_sd4_2500  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Cisplatin1normmean_sd7_2500= left_join(Cisplatin1normmean_sd5_2500,Cisplatin1normmean_sd6_2500,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

Cisplatin1norm_1000nM= Cisplatin1norm%>% filter(condition=="1000nM")
Cisplatin1norm_1000nM_filter= Cisplatin1norm_1000nM%>% filter(!well_name_p=="D13_1") # minder dan 15 cellen bij tijdspunt 0
Cisplatin1normmean_sd2_1000= Cisplatin1norm_1000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(Mean_nnorm=mean(n_norm))
Cisplatin1normmean_sd3_1000= Cisplatin1norm_1000nM_filter  %>% group_by(well_name,Image_Group_Index.x)%>%summarize(sd_nnorm=sd(n_norm))
Cisplatin1normmmean_sd4_1000= left_join(Cisplatin1normmean_sd2_1000,Cisplatin1normmean_sd3_1000,by = c("well_name", "Image_Group_Index.x"))
Cisplatin1normmmean_sd4_1000=Cisplatin1normmmean_sd4_1000%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_cisplatin = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Cisplatin1normmmean_sd4_1000=left_join(Cisplatin1normmmean_sd4_1000,conditions_cisplatin, by= "well_name")
Cisplatin1normmean_sd5_1000= Cisplatin1normmmean_sd4_1000  %>% group_by(Image_Group_Index.x,condition)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Cisplatin1normmean_sd6_1000= Cisplatin1normmmean_sd4_1000  %>% group_by(Image_Group_Index.x,condition)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Cisplatin1normmean_sd7_1000= left_join(Cisplatin1normmean_sd5_1000,Cisplatin1normmean_sd6_1000,by = c("Image_Group_Index.x","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)

fifacis=ggplot()+
  geom_line(data = Cisplatin1normmmean_sd4 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("J06" = "darkred", "J11" = "red", "K06" = "blue","K11" = "green"),
                     name = "Well_name")+
  theme_classic()+
  ggtitle("100 µM Cisplatin")+
  geom_ribbon(data = Cisplatin1normmmean_sd4,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("J06" = "darkred", "J11" = "red", "K06" = "blue","K11" = "green"))+
  guides(fill = "none")
  

fifacis2=ggplot()+
  geom_line(data = Cisplatin1normmmean_sd4_1000 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("D12" = "darkred", "D13" = "red", "E12" = "blue","E13" = "green"),
                     name = "Well_name")+
  theme_classic()+
  ggtitle("1 µM Cisplatin")+
  geom_ribbon(data = Cisplatin1normmmean_sd4_1000,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("D12" = "darkred", "D13" = "red", "E12" = "blue","E13" = "green"))+
  guides(fill = "none")

fifacis3=ggplot()+
  geom_line(data = Cisplatin1normmmean_sd4_2500 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("B12" = "darkred", "B13" = "red", "C12" = "blue","C13" = "green"),
                     name = "Well_name")+
  theme_classic()+
  ggtitle("2.5 µM Cisplatin")+
  geom_ribbon(data = Cisplatin1normmmean_sd4_2500,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("B12" = "darkred", "B13" = "red", "C12" = "blue","C13" = "green"))+
  guides(fill = "none")

fifacis4=ggplot()+
  geom_line(data = Cisplatin1normmmean_sd4_5000 , aes(x = time_h, y =Mean_nnorm,color= well_name))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("L06" = "darkred", "L11" = "red", "M06" = "blue","M11" = "green"),
                     name = "Well_name")+
  theme_classic()+
  ggtitle("5 µM Cisplatin")+
  geom_ribbon(data = Cisplatin1normmmean_sd4_5000,aes(x = time_h, ymax = Mean_nnorm+sd_nnorm, ymin = Mean_nnorm-sd_nnorm,color = well_name,fill=well_name),alpha=0.05)+
  scale_fill_manual(values = c("L06" = "darkred", "L11" = "red", "M06" = "blue","M11" = "green"))+
  guides(fill = "none")

fifa=ggplot()+
  geom_line(data = DMSO1norm , aes(x = time_h, y =n_norm,color= well_name_p))+
  ylab("Normalized cell count")+xlab("Time (h)")+
  theme_classic() 

#plot Cisplatin
Cisplatin_full_join= full_join(Cisplatin1normmean_sd7,Cisplatin1normmean_sd7_1000)
Cisplatin_full_join=full_join(Cisplatin_full_join,Cisplatin1normmean_sd7_2500)
Cisplatin_full_join=full_join(Cisplatin_full_join,Cisplatin1normmean_sd7_5000)
Cisplatin_DMSO= full_join(Cisplatin_full_join,DMSO1normmean_sd7)
Cisplatin_DMSO_DMEM= full_join(Cisplatin_DMSO,DMEM1normmean_sd7)
Cisplatin_DMSO_DMEM_ma= Cisplatin_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")

Cisplatin1000_DMSO= full_join(DMSO1normmean_sd7,Cisplatin1normmean_sd7_1000)
Cisplatin1000_DMSO_DMEM=full_join(Cisplatin1000_DMSO,DMEM1normmean_sd7)
Cisplatin1000_DMSO_DMEM_ma= Cisplatin1000_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")
Cisplatin2500_DMSO= full_join(DMSO1normmean_sd7,Cisplatin1normmean_sd7_2500)
Cisplatin2500_DMSO_DMEM=full_join(Cisplatin2500_DMSO,DMEM1normmean_sd7)
Cisplatin2500_DMSO_DMEM_ma= Cisplatin2500_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")
Cisplatin5000_DMSO= full_join(DMSO1normmean_sd7,Cisplatin1normmean_sd7_5000)
Cisplatin5000_DMSO_DMEM=full_join(Cisplatin5000_DMSO,DMEM1normmean_sd7)
Cisplatin5000_DMSO_DMEM_ma= Cisplatin5000_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")
Cisplatin100000_DMSO= full_join(DMSO1normmean_sd7,Cisplatin1normmean_sd7)
Cisplatin100000_DMSO_DMEM=full_join(Cisplatin100000_DMSO,DMEM1normmean_sd7)
Cisplatin100000_DMSO_DMEM_ma= Cisplatin100000_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")

Plot_Cisplatin_all=ggplot()+
  geom_line(data =Cisplatin_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Cisplatin_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.05)+
  scale_fill_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()

Plot_Cisplatin_100000=ggplot()+
  geom_line(data =Cisplatin100000_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Cisplatin100000_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.05)+
  scale_fill_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("100 µM")

Plot_Cisplatin_5000=ggplot()+
  geom_line(data =Cisplatin5000_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Cisplatin5000_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.05)+
  scale_fill_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("5 µM")

Plot_Cisplatin_2500=ggplot()+
  geom_line(data =Cisplatin2500_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Cisplatin2500_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.05)+
  scale_fill_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("2.5 µM")

Plot_Cisplatin_1000=ggplot()+
  geom_line(data =Cisplatin1000_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Cisplatin1000_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.05)+
  scale_fill_manual(values = c("100000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("1 µM")


#plot alle Etoposide
Etoposide_full_join= full_join(Etoposide1normmean_sd7,Etoposide1normmean_sd7_1000)
Etoposide_full_join=full_join(Etoposide_full_join,Etoposide1normmean_sd7_2500)
Etoposide_full_join=full_join(Etoposide_full_join,Etoposide1normmean_sd7_5000)
Etoposide_DMSO= full_join(Etoposide_full_join,DMSO1normmean_sd7)
Etoposide_DMSO_DMEM= full_join(Etoposide_DMSO,DMEM1normmean_sd7)
Etoposide_DMSO_DMEM_ma= Etoposide_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")

#per condition Etoposide
Etoposide1000_DMSO= full_join(DMSO1normmean_sd7,Etoposide1normmean_sd7_1000)
Etoposide1000_DMSO_DMEM=full_join(Etoposide1000_DMSO,DMEM1normmean_sd7)
Etoposide1000_DMSO_DMEM_ma= Etoposide1000_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")
Etoposide2500_DMSO=full_join(DMSO1normmean_sd7,Etoposide1normmean_sd7_2500)
Etoposide2500_DMSO_DMEM=full_join(Etoposide2500_DMSO,DMEM1normmean_sd7)
Etoposide2500_DMSO_DMEM_ma= Etoposide2500_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")
Etoposide5000_DMSO= full_join(DMSO1normmean_sd7,Etoposide1normmean_sd7_5000)
Etoposide5000_DMSO_DMEM= full_join(Etoposide5000_DMSO,DMEM1normmean_sd7)
Etoposide5000_DMSO_DMEM_ma=Etoposide5000_DMSO_DMEM   %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")
Etoposide10000_DMSO= full_join(Etoposide1normmean_sd7,DMSO1normmean_sd7)
Etoposide10000_DMSO_DMEM= full_join(Etoposide10000_DMSO,DMEM1normmean_sd7)
Etoposide10000_DMSO_DMEM_ma= Etoposide10000_DMSO_DMEM %>%group_by(condition) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="40")



Plot_Etoposide_all=ggplot()+
  geom_line(data =Etoposide_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Etoposide_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.05)+
  scale_fill_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()

Plot_Etoposide1000_DMSO=ggplot()+
  geom_line(data =Etoposide1000_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Etoposide1000_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.1)+
  scale_fill_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("1 µM")

Plot_Etoposide2500_DMSO_DMEM= ggplot()+
  geom_line(data =Etoposide2500_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Etoposide2500_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.1)+
  scale_fill_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("2.5 µM")

Plot_Etoposide5000_DMSO=ggplot()+
  geom_line(data =Etoposide5000_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Etoposide5000_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.1)+
  scale_fill_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("5 µM")

Plot_Etoposide10000_DMSO_DMEM=ggplot()+
  geom_line(data =Etoposide10000_DMSO_DMEM_ma , aes(x = time_h, y =ma_total_Mean_normalized,color= condition),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ),
                     name = "Condition")+
  geom_ribbon(data = Etoposide10000_DMSO_DMEM_ma,aes(x = time_h, ymax = ma_total_Mean_normalized+sd_nnorm, ymin = ma_total_Mean_normalized-sd_nnorm,color = condition,fill=condition),alpha=0.1)+
  scale_fill_manual(values = c("10000nM" = "darkred", "5000nM" = "red", "2500nM" = "blue","1000nM" = "green","DMSO"= "black","DMEM"="orange" ))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("10 µM")

#cell populatie nivea phase.
#DMSO
Data_DMSO1=Data_DMSO
print(max(Data_DMSO1$Cy3, na.rm=TRUE))*0.05
print(max(Data_DMSO1$GFP,na.rm=TRUE))*0.05

Data_DMSO1= Data_DMSO1 %>%mutate(phase= ifelse(GFP>0.05& Cy3<0.05,"G2",
                                                                   ifelse(Cy3>0.05& GFP<0.05,"G1",
                                                                          ifelse(Cy3>0.05& GFP>0.05,"G1",
                                                                                 ifelse(Cy3<0.05& GFP<0.05,"Early G1",NA)))))

Data_DMSO1$phase= as.factor(Data_DMSO1$phase)
Data_DMSO1=Data_DMSO1%>%filter(!is.na(phase))
Data_DMSO1_0hours= Data_DMSO1%>%filter(time_h==0)
Data_DMSO1_24hours= Data_DMSO1%>%filter(timeID==15)
Data_DMSO1_66hours= Data_DMSO1%>%filter(timeID==40)


ggplot_DMSO_threshold_determination1= ggplot()+
  geom_point(data=Data_DMSO1_0hours, aes(x=GFP,y=Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.05 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.05  , color = "black", linetype = "dotted")+
  ylab("Normalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  scale_y_continuous(limits = c(0, max(1)))+
  scale_x_continuous(limits = c(0, max(1)))+
  theme_classic()+
  ggtitle("0 hours")

ggplot_DMSO_threshold_determination1= ggplot()+
  geom_point(data=Data_DMSO1_24hours, aes(x=GFP,y=Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.05 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.05  , color = "black", linetype = "dotted")+
  ylab("Normalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  scale_y_continuous(limits = c(0, max(1)))+
  scale_x_continuous(limits = c(0, max(1)))+
  ggtitle("24 hours")

ggplot_DMSO_threshold_determination1= ggplot()+
  geom_point(data=Data_DMSO1_66hours, aes(x=GFP,y=Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.05 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.05  , color = "black", linetype = "dotted")+
  ylab("Normalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  scale_y_continuous(limits = c(0, max(1)))+
  scale_x_continuous(limits = c(0, max(1)))+
  ggtitle("66 hours")

Data_DMSO2= Data_DMSO1 %>% group_by(Image_Group_Number,Image_Group_Index)%>% count(phase)
Data_DMSO3 = Data_DMSO1 %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Data_DMSO3= Data_DMSO3 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
Data_DMSO4=Data_DMSO2 %>% left_join(Data_DMSO3 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(phasenorm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_DMSO = DMSO_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_DMSO = Metadata_DMSO%>% select(c("well_name","Image_Group_Number","well_name_p"))
#Wellnames_DMEM <- Wellnames_DMEM  %>% dplyr::rename(Image_Group_Number=loc)
Data_DMSO5= left_join(Data_DMSO4,Wellnames_DMSO,by="Image_Group_Number")
Data_DMSO5=left_join(Data_DMSO5,conditions_DMSO,by="well_name")
Data_DMSO6= Data_DMSO5 %>% filter(!well_name_p=="H02_1",!well_name_p=="H02_2",!well_name_p=="I02_1",!well_name_p=="I02_2",!well_name_p=="J02_1",!well_name_p=="J02_2",!well_name_p=="K07_2",!well_name_p=="L02_1",!well_name_p=="L02_2",!well_name_p=="M02_1",!well_name_p=="M02_2",!well_name_p=="K02_1",!well_name_p=="K02_2",!well_name_p=="L07_1")# minder dan 15 cellen bij tijdspunt 0
Data_DMSO7= Data_DMSO6  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_DMSO8= Data_DMSO6  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_DMSO9= left_join(Data_DMSO7,Data_DMSO8,by = c("well_name", "Image_Group_Index.x","phase"))
Data_DMSO9=Data_DMSO9%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_DMSO1 = DMSO_layout %>% select(c("well_name","condition","dose_uM"))
Data_DMSO9=left_join(Data_DMSO9,conditions_DMSO1, by= "well_name")
Data_DMSO10= Data_DMSO9  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_DMSO11= Data_DMSO9 %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_DMSO12= left_join(Data_DMSO10,Data_DMSO11,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_DMSO13= Data_DMSO12 %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="41")

plot_DMSO =ggplot()+
  geom_line(data = Data_DMSO13  , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_DMSO13,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
theme_classic()+
  ggtitle("0 µM")

#DMEM
Data_DMEM1=Data_DMEM
print(max(Data_DMEM1$Cy3, na.rm=TRUE))*0.05
print(max(Data_DMEM1$GFP,na.rm=TRUE))*0.05

Data_DMEM1= Data_DMEM1 %>%mutate(phase= ifelse(GFP>0.05& Cy3<0.05,"G2",
                                               ifelse(Cy3>0.05& GFP<0.05,"G1",
                                                      ifelse(Cy3>0.05& GFP>0.05,"G1",
                                                             ifelse(Cy3<0.05& GFP<0.05,"Early G1",NA)))))

Data_DMEM1$phase= as.factor(Data_DMEM1$phase)
Data_DMEM1=Data_DMEM1%>%filter(!is.na(phase))


ggplot_DMEM_threshold_determination1= ggplot()+
  geom_point(data=Data_DMEM1, aes(x=GFP,y=Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.05 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.05  , color = "black", linetype = "dotted")+
  ylab("Notmalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  ggtitle("DMEM")


Data_DMEM2= Data_DMEM1 %>% group_by(Image_Group_Number,Image_Group_Index)%>% count(phase)
Data_DMEM3 = Data_DMEM1 %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Data_DMEM3= Data_DMEM3 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
Data_DMEM4=Data_DMEM2 %>% left_join(Data_DMEM3 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(phasenorm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_DMEM = DMEM_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_DMEM = Metadata_DMEM%>% select(c("well_name","Image_Group_Number","well_name_p"))
#Wellnames_DMEM <- Wellnames_DMEM  %>% dplyr::rename(Image_Group_Number=loc)
Data_DMEM5= left_join(Data_DMEM4,Wellnames_DMEM,by="Image_Group_Number")
Data_DMEM5=left_join(Data_DMEM5,conditions_DMEM,by="well_name")
Data_DMEM6= Data_DMEM5 %>% filter(!well_name_p=="B02_1",!well_name_p=="C02_1",!well_name_p=="D02_1",!well_name_p=="D02_2",!well_name_p=="E02_1",!well_name_p=="F02_1",!well_name_p=="F07_2",!well_name_p=="G02_1",!well_name_p=="G02_2",!well_name_p=="F02_2")# minder dan 15 cellen bij tijdspunt 0
Data_DMEM7= Data_DMEM6  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_DMEM8= Data_DMEM6  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_DMEM9= left_join(Data_DMEM7,Data_DMEM8,by = c("well_name", "Image_Group_Index.x","phase"))
Data_DMEM9=Data_DMEM9%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_DMEM1 = DMEM_layout %>% select(c("well_name","condition","dose_uM"))
Data_DMEM9=left_join(Data_DMEM9,conditions_DMEM1, by= "well_name")
Data_DMEM10= Data_DMEM9  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_DMEM11= Data_DMEM9 %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_DMEM12= left_join(Data_DMEM10,Data_DMEM11,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_DMEM13= Data_DMEM12 %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>%filter(!Image_Group_Index.x=="41")

plot_DMSO=ggplot()+
  geom_line(data = Data_DMEM13  , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_DMEM13,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("0 µM DMEM")

#Cisplatin
Data_Cisplatin1=Data_Cisplatin
print(max(Data_Cisplatin1$Cy3, na.rm=TRUE))*0.05
print(max(Data_Cisplatin1$GFP,na.rm=TRUE))*0.05

Data_Cisplatin1= Data_Cisplatin1 %>%mutate(phase= ifelse(GFP>0.05& Cy3<0.05,"G2",
                                               ifelse(Cy3>0.05& GFP<0.05,"G1",
                                                      ifelse(Cy3>0.05& GFP>0.05,"G1",
                                                             ifelse(Cy3<0.05& GFP<0.05,"Early G1",NA)))))

Data_Cisplatin1$phase= as.factor(Data_Cisplatin1$phase)
Data_Cisplatin1=Data_Cisplatin1%>%filter(!is.na(phase))
Data_Cisplatin_0hours= Data_Cisplatin1%>%filter(timeID==1&dose_uM==5)
Data_Cisplatin_24hours= Data_Cisplatin1%>%filter(timeID==15&dose_uM==5)
Data_Cisplatin_66hours= Data_Cisplatin1%>%filter(timeID==40&dose_uM==5)

ggplot_Cisplatin_threshold_determination1= ggplot()+
  geom_point(data=Data_Cisplatin_0hours, aes(x=GFP,y=Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.05 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.05  , color = "black", linetype = "dotted")+
  ylab("Normalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  scale_y_continuous(limits = c(0, max(1)))+
  scale_x_continuous(limits = c(0, max(1)))+
  ggtitle("0 hours")

ggplot_Cisplatin_threshold_determination1= ggplot()+
  geom_point(data=Data_Cisplatin_24hours, aes(x=GFP,y=Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.05 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.05  , color = "black", linetype = "dotted")+
  ylab("Normalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  scale_y_continuous(limits = c(0, max(1)))+
  scale_x_continuous(limits = c(0, max(1)))+
  ggtitle("24 hours")

ggplot_Cisplatin_threshold_determination1= ggplot()+
  geom_point(data=Data_Cisplatin_66hours, aes(x=GFP,y=Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.05 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.05  , color = "black", linetype = "dotted")+
  ylab("Normalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  scale_y_continuous(limits = c(0, max(1)))+
  scale_x_continuous(limits = c(0, max(1)))+
  ggtitle("66 hours")


Data_Cisplatin2= Data_Cisplatin1 %>% group_by(Image_Group_Number,Image_Group_Index)%>% count(phase)
Data_Cisplatin3 = Data_Cisplatin1 %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Data_Cisplatin3= Data_Cisplatin3 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
Data_Cisplatin4=Data_Cisplatin2 %>% left_join(Data_Cisplatin3 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(phasenorm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Cisplatin = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_Cisplatin = Metadata_Cisplatin%>% select(c("well_name","Image_Group_Number","well_name_p"))
#Wellnames_DMEM <- Wellnames_DMEM  %>% dplyr::rename(Image_Group_Number=loc)
Data_Cisplatin5= left_join(Data_Cisplatin4,Wellnames_Cisplatin,by="Image_Group_Number")
Data_Cisplatin5=left_join(Data_Cisplatin5,conditions_cisplatin,by="well_name")

Data_Cisplatin_100000nM= Data_Cisplatin5%>% filter(condition=="100000nM")
Data_Cisplatin1_100000nM= Data_Cisplatin_100000nM %>% filter(!well_name_p=="J06_2",!well_name_p=="K06_1")# minder dan 15 cellen bij tijdspunt 0
Data_Cisplatin2_100000nM= Data_Cisplatin1_100000nM %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_Cisplatin3_100000nM= Data_Cisplatin1_100000nM  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_Cisplatin4_100000nM= left_join(Data_Cisplatin2_100000nM,Data_Cisplatin3_100000nM,by = c("well_name", "Image_Group_Index.x","phase"))
Data_Cisplatin4_100000nM=Data_Cisplatin4_100000nM%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Cisplatin1 = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Data_Cisplatin4_100000nM=left_join(Data_Cisplatin4_100000nM,conditions_Cisplatin1, by= "well_name")
Data_Cisplatin5_100000nM= Data_Cisplatin4_100000nM  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_Cisplatin6_100000nM= Data_Cisplatin4_100000nM %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_Cisplatin7_100000nM= left_join(Data_Cisplatin5_100000nM,Data_Cisplatin6_100000nM,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_Cisplatin8_100000nM= Data_Cisplatin7_100000nM %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>% filter(!Image_Group_Index.x=="40")

Plot_Cisplatin_100000nM=ggplot()+
  geom_line(data = Data_Cisplatin8_100000nM  , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_Cisplatin8_100000nM,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("100 µM Cisplatin")

Data_Cisplatin_5000nM= Data_Cisplatin5%>% filter(condition=="5000nM")
Data_Cisplatin1_5000nM= Data_Cisplatin_5000nM %>% filter(!well_name_p=="M06_2",!well_name_p=="M11_1")# minder dan 15 cellen bij tijdspunt 0
Data_Cisplatin2_5000nM= Data_Cisplatin1_5000nM %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_Cisplatin3_5000nM= Data_Cisplatin1_5000nM  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_Cisplatin4_5000nM= left_join(Data_Cisplatin2_5000nM,Data_Cisplatin3_5000nM,by = c("well_name", "Image_Group_Index.x","phase"))
Data_Cisplatin4_5000nM=Data_Cisplatin4_5000nM%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Cisplatin1 = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Data_Cisplatin4_5000nM=left_join(Data_Cisplatin4_5000nM,conditions_Cisplatin1, by= "well_name")
Data_Cisplatin5_5000nM= Data_Cisplatin4_5000nM  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_Cisplatin6_5000nM= Data_Cisplatin4_5000nM %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_Cisplatin7_5000nM= left_join(Data_Cisplatin5_5000nM,Data_Cisplatin6_5000nM,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_Cisplatin8_5000nM= Data_Cisplatin7_5000nM %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>% filter(!Image_Group_Index.x=="40")

Plot_Cisplatin_5000nM=ggplot()+
  geom_line(data = Data_Cisplatin8_5000nM  , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_Cisplatin8_5000nM,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("5 µM Cisplatin")

Data_Cisplatin_2500nM= Data_Cisplatin5%>% filter(condition=="2500nM")
Data_Cisplatin1_2500nM= Data_Cisplatin_2500nM # minder dan 15 cellen bij tijdspunt 0
Data_Cisplatin2_2500nM= Data_Cisplatin1_2500nM %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_Cisplatin3_2500nM= Data_Cisplatin1_2500nM  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_Cisplatin4_2500nM= left_join(Data_Cisplatin2_2500nM,Data_Cisplatin3_2500nM,by = c("well_name", "Image_Group_Index.x","phase"))
Data_Cisplatin4_2500nM=Data_Cisplatin4_2500nM%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Cisplatin1 = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Data_Cisplatin4_2500nM=left_join(Data_Cisplatin4_2500nM,conditions_Cisplatin1, by= "well_name")
Data_Cisplatin5_2500nM= Data_Cisplatin4_2500nM  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_Cisplatin6_2500nM= Data_Cisplatin4_2500nM %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_Cisplatin7_2500nM= left_join(Data_Cisplatin5_2500nM,Data_Cisplatin6_2500nM,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_Cisplatin8_2500nM= Data_Cisplatin7_2500nM %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>% filter(!Image_Group_Index.x=="41")

Plot_Cisplatin_2500nM=ggplot()+
  geom_line(data = Data_Cisplatin8_2500nM  , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_Cisplatin8_2500nM,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("2.5 µM Cisplatin")

Data_Cisplatin_1000nM= Data_Cisplatin5%>% filter(condition=="1000nM")
Data_Cisplatin1_1000nM= Data_Cisplatin_1000nM%>%filter(!well_name_p=="D13_1") # minder dan 15 cellen bij tijdspunt 0
Data_Cisplatin2_1000nM= Data_Cisplatin1_1000nM %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_Cisplatin3_1000nM= Data_Cisplatin1_1000nM  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_Cisplatin4_1000nM= left_join(Data_Cisplatin2_1000nM,Data_Cisplatin3_1000nM,by = c("well_name", "Image_Group_Index.x","phase"))
Data_Cisplatin4_1000nM=Data_Cisplatin4_1000nM%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Cisplatin1 = Cisplatin_layout %>% select(c("well_name","condition","dose_uM"))
Data_Cisplatin4_1000nM=left_join(Data_Cisplatin4_1000nM,conditions_Cisplatin1, by= "well_name")
Data_Cisplatin5_1000nM= Data_Cisplatin4_1000nM  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_Cisplatin6_1000nM= Data_Cisplatin4_1000nM %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_Cisplatin7_1000nM= left_join(Data_Cisplatin5_1000nM,Data_Cisplatin6_1000nM,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_Cisplatin8_1000nM= Data_Cisplatin7_1000nM %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>% filter(!Image_Group_Index.x=="40")

Plot_Cisplatin_1000nM=ggplot()+
  geom_line(data = Data_Cisplatin8_1000nM  , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_Cisplatin8_1000nM,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("1 µM Cisplatin")

#Etoposide
Data_Etoposide1=Data_Etoposide
print(max(Data_Etoposide1$Cy3, na.rm=TRUE))*0.05
print(max(Data_Etoposide1$GFP,na.rm=TRUE))*0.05

Data_Etoposide1= Data_Etoposide1 %>%mutate(phase= ifelse(GFP>0.05& Cy3<0.05,"G2",
                                                         ifelse(Cy3>0.05& GFP<0.05,"G1",
                                                                ifelse(Cy3>0.05& GFP>0.05,"G1",
                                                                       ifelse(Cy3<0.05& GFP<0.05,"Early G1",NA)))))

Data_Etoposide1$phase= as.factor(Data_Etoposide1$phase)
Data_Etoposide1=Data_Etoposide1%>%filter(!is.na(phase))


ggplot_Etoposide_threshold_determination1= ggplot()+
  geom_point(data=Data_Etoposide1, aes(x=GFP,y=Cy3, color=phase))+
  scale_colour_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  geom_vline(xintercept = 0.05 , color = "black", linetype = "dotted")+
  geom_hline(yintercept = 0.05  , color = "black", linetype = "dotted")+
  ylab("Notmalized Cy3 intensity")+xlab("Normalized GFP intensity")+
  theme_classic()+
  ggtitle("Etoposide")

Data_Etoposide2= Data_Etoposide1 %>% group_by(Image_Group_Number,Image_Group_Index)%>% count(phase)
Data_Etoposide3 = Data_Etoposide1 %>% group_by(Image_Group_Number)%>% count(Image_Group_Index)
Data_Etoposide3= Data_Etoposide3 %>%filter(Image_Group_Index==1)%>% rename(n0=n)
Data_Etoposide4=Data_Etoposide2 %>% left_join(Data_Etoposide3 ,by="Image_Group_Number")%>%select(-c(Image_Group_Index.y))%>%mutate(phasenorm=n/n0)%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Wellnames_Etoposide = Metadata_Etoposide%>% select(c("well_name","Image_Group_Number","well_name_p"))
#Wellnames_DMEM <- Wellnames_DMEM  %>% dplyr::rename(Image_Group_Number=loc)
Data_Etoposide5= left_join(Data_Etoposide4,Wellnames_Etoposide,by="Image_Group_Number")
Data_Etoposide5=left_join(Data_Etoposide5,conditions_Etoposide,by="well_name")

Data_Etoposide_10000nM= Data_Etoposide5%>% filter(condition=="10000nM")
Data_Etoposide1_10000nM= Data_Etoposide_10000nM %>% filter(!well_name_p=="B06_2",!well_name_p=="B11_2",!well_name_p=="C06_2")# minder dan 15 cellen bij tijdspunt 0
Data_Etoposide2_10000nM= Data_Etoposide1_10000nM %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_Etoposide3_10000nM= Data_Etoposide1_10000nM  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_Etoposide4_10000nM= left_join(Data_Etoposide2_10000nM,Data_Etoposide3_10000nM,by = c("well_name", "Image_Group_Index.x","phase"))
Data_Etoposide4_10000nM=Data_Etoposide4_10000nM%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide1 = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Data_Etoposide4_10000nM=left_join(Data_Etoposide4_10000nM,conditions_Etoposide1, by= "well_name")
Data_Etoposide5_10000nM= Data_Etoposide4_10000nM  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_Etoposide6_10000nM= Data_Etoposide4_10000nM %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_Etoposide7_10000nM= left_join(Data_Etoposide5_10000nM,Data_Etoposide6_10000nM,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_Etoposide8_10000nM= Data_Etoposide7_10000nM %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>% filter(!Image_Group_Index.x=="41")

Plot_Etoposide_10000nM=ggplot()+
  geom_line(data = Data_Etoposide8_10000nM , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_Etoposide8_10000nM,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("10 µM Etoposide")

Data_Etoposide_5000nM= Data_Etoposide5%>% filter(condition=="5000nM")
Data_Etoposide1_5000nM= Data_Etoposide_5000nM %>% filter(!well_name_p=="D06_2",!well_name_p=="E06_2")# minder dan 15 cellen bij tijdspunt 0
Data_Etoposide2_5000nM= Data_Etoposide1_5000nM %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_Etoposide3_5000nM= Data_Etoposide1_5000nM  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_Etoposide4_5000nM= left_join(Data_Etoposide2_5000nM,Data_Etoposide3_5000nM,by = c("well_name", "Image_Group_Index.x","phase"))
Data_Etoposide4_5000nM=Data_Etoposide4_5000nM%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide1 = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Data_Etoposide4_5000nM=left_join(Data_Etoposide4_5000nM,conditions_Etoposide1, by= "well_name")
Data_Etoposide5_5000nM= Data_Etoposide4_5000nM  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_Etoposide6_5000nM= Data_Etoposide4_5000nM %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_Etoposide7_5000nM= left_join(Data_Etoposide5_5000nM,Data_Etoposide6_5000nM,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_Etoposide8_5000nM= Data_Etoposide7_5000nM %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>% filter(!Image_Group_Index.x=="41")

Plot_Etoposide_5000nM=ggplot()+
  geom_line(data = Data_Etoposide8_5000nM , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_Etoposide8_5000nM,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("5 µM Etoposide")

Data_Etoposide_2500nM= Data_Etoposide5%>% filter(condition=="2500nM")
Data_Etoposide1_2500nM= Data_Etoposide_2500nM # minder dan 15 cellen bij tijdspunt 0
Data_Etoposide2_2500nM= Data_Etoposide1_2500nM %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_Etoposide3_2500nM= Data_Etoposide1_2500nM  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_Etoposide4_2500nM= left_join(Data_Etoposide2_2500nM,Data_Etoposide3_2500nM,by = c("well_name", "Image_Group_Index.x","phase"))
Data_Etoposide4_2500nM=Data_Etoposide4_2500nM%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide1 = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Data_Etoposide4_2500nM=left_join(Data_Etoposide4_2500nM,conditions_Etoposide1, by= "well_name")
Data_Etoposide5_2500nM= Data_Etoposide4_2500nM  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_Etoposide6_2500nM= Data_Etoposide4_2500nM %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_Etoposide7_2500nM= left_join(Data_Etoposide5_2500nM,Data_Etoposide6_2500nM,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_Etoposide8_2500nM= Data_Etoposide7_2500nM %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>% filter(!Image_Group_Index.x=="41")

Plot_Etoposide_2500nM=ggplot()+
  geom_line(data = Data_Etoposide8_2500nM , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_Etoposide8_2500nM,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("2.5 µM Etoposide")

Data_Etoposide_1000nM= Data_Etoposide5%>% filter(condition=="1000nM")
Data_Etoposide1_1000nM= Data_Etoposide_1000nM # minder dan 15 cellen bij tijdspunt 0
Data_Etoposide2_1000nM= Data_Etoposide1_1000nM %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(Mean_nnorm=mean(phasenorm))
Data_Etoposide3_1000nM= Data_Etoposide1_1000nM  %>% group_by(well_name,Image_Group_Index.x,phase)%>%summarize(sd_nnorm=sd(phasenorm))
Data_Etoposide4_1000nM= left_join(Data_Etoposide2_1000nM,Data_Etoposide3_1000nM,by = c("well_name", "Image_Group_Index.x","phase"))
Data_Etoposide4_1000nM=Data_Etoposide4_1000nM%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
conditions_Etoposide1 = Etoposide_layout %>% select(c("well_name","condition","dose_uM"))
Data_Etoposide4_1000nM=left_join(Data_Etoposide4_1000nM,conditions_Etoposide1, by= "well_name")
Data_Etoposide5_1000nM= Data_Etoposide4_1000nM  %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(Mean_nnorm=mean(Mean_nnorm))
Data_Etoposide6_1000nM= Data_Etoposide4_1000nM %>% group_by(Image_Group_Index.x,condition,phase)%>%summarize(sd_nnorm=sd(Mean_nnorm))
Data_Etoposide7_1000nM= left_join(Data_Etoposide5_1000nM,Data_Etoposide6_1000nM,by = c("Image_Group_Index.x","phase","condition"))%>%mutate(time_h = (Image_Group_Index.x * TIME_STEP)- TIME_STEP)
Data_Etoposide8_1000nM= Data_Etoposide7_1000nM %>%group_by(phase) %>% mutate(ma_total_Mean_normalized = rollmean(Mean_nnorm, k=2, fill=NA, align='center'))%>% mutate(ma_total_sd_normalized = rollmean(sd_nnorm, k=2, fill=NA, align='center'))%>% filter(!Image_Group_Index.x=="41")

Plot_Etoposide_1000nM=ggplot()+
  geom_line(data = Data_Etoposide8_1000nM , aes(x = time_h, y =ma_total_Mean_normalized,color= phase),size = 2)+
  ylab("Normalized cell count")+xlab("Time (h)")+
  scale_color_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"),
                     name = "Phases")+
  geom_ribbon(data = Data_Etoposide8_1000nM,aes(x = time_h, ymax = ma_total_Mean_normalized+ma_total_sd_normalized ,ymin = ma_total_Mean_normalized-ma_total_sd_normalized ,color = phase,fill=phase),alpha=0.05)+
  scale_fill_manual(values = c("G2" = "green", "G1" = "red", "Early G1" = "grey"))+
  guides(fill = "none")+
  theme_classic()+
  ggtitle("1 µM Etoposide")