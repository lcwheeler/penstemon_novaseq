library(tidyverse)

#define function to reorganize pixy data
#because tidyverse is *obviously* superior to base R
#and I'm definitely not just doing this because copy/paste
#(code from https://pixy.readthedocs.io/en/latest/plotting.html)
pixy_to_long <- function(pixy_files){
  
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
      
      
    } else{
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
      
    }
    
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  
}


#read in files
pixy_folder <- "~/project storage/project_dasanthera_novaseq/analysis/analysis_pixy_50kb/analysis_pixy_50kb"
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)


#interested only in newberryi and rupicola, so subset the df
#also replace single-species metrics of interest so we can plot those in addition
subset_pixy_df <- pixy_df %>%
  filter(is.na(pop2) | pop2 == "P_rupicola") %>%
  filter(pop1 == c("P_newberryi", "P_rupicola")) %>%
  mutate(statistic = replace(statistic, statistic =="avg_pi" & pop1 == "P_newberryi", "newberryi_avg_pi")) %>%
  mutate(statistic = replace(statistic, statistic =="avg_pi" & pop1 == "P_rupicola", "rupicola_avg_pi"))

  

# custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(newberryi_avg_pi = "pi[new]",
                               rupicola_avg_pi = "pi[rup]",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)



# plotting summary statistics along a single chromosome
#set scaffold list and taxon list
scaflist <- unique(subset_pixy_df$chromosome)


setwd('~/Desktop/hummingbird_dxy_fst/')
#set up pdf plotting
for (for_scaf in scaflist){
  pdf(file = paste("pixyplot_",for_scaf,".pdf", sep = ''))
  t <- subset_pixy_df %>%
    filter(chromosome == for_scaf) %>%
    filter(statistic %in% c("newberryi_avg_pi", "rupicola_avg_pi", "avg_dxy", "avg_wc_fst")) %>%
    mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000)
    print(
      ggplot(t, aes(x = chr_position, y = value, color = statistic))+
        geom_line(size = 0.25)+
        facet_grid(statistic ~ .,
                   scales = "free_y", switch = "x", space = "free_x",
                   labeller = labeller(statistic = pixy_labeller,
                                       value = label_value))+
        xlab("Position on Chromosome (Mb)")+
        ylab("Statistic Value")+
        ggtitle(paste("P. newberryi vs. P. rupicola", for_scaf, sep = ' '))+
        theme_bw()+
        theme(panel.spacing = unit(0.1, "cm"),
              strip.background = element_blank(),
              strip.placement = "outside",
              legend.position = "none")+
        scale_x_continuous(expand = c(0, 0))+
        scale_y_continuous(expand = c(0, 0))+
        scale_color_brewer(palette = "Set1")
    )
    dev.off()
}


