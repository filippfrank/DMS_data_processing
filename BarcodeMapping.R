library(ggplot2)
library(hrbrthemes)
library("dplyr")
library(RColorBrewer)
library(tidyr)
library(plyr)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(seqinr)
require(ggseqlogo)
library(tidyverse)
library(grid)

################################################################################
###Rename 
################################################################################
antibody <- "C706"

escape_counts <- read_tsv("C706_N05_found.txt")
input_counts <- read_tsv("220311_N05_found.txt")

escape_counts <- escape_counts[,2:3]
input_counts <- input_counts[,2:3]

#####################################################################
#Combine duplicated mutations (same mutation with multiple barcodes)
#####################################################################
escape_combined <- escape_counts %>% 
  group_by(mutation) %>% 
  summarise_all(funs(sum))

reference_combined <- input_counts %>% 
  group_by(mutation) %>% 
  summarise_all(funs(sum))

colnames(reference_combined) <- c("name", "depth")
colnames(escape_combined) <- c("name", "depth")

################################################################################
###Convert mutation "names" into: wildtype, site, and mutations
################################################################################
####################################################################
##Renaming for escape
####################################################################

fun_first <- function(x) {
  substring(x, 1,1)
}

fun_last <- function(x){
  substr(x, nchar(x), nchar(x))
}

fun_site <- function(x){
  substr(x, 2, nchar(x)-1)
}

wildtype <- data.frame(sapply(escape_combined[1], fun_first))
#colnames(wildtype) <- c("wildtype")
escape_combined$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(escape_combined[1], fun_last))
colnames(mutation) <- c("mutation")
escape_combined$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(escape_combined[1], fun_site))
colnames(site) <- c("site")
escape_combined$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
escape_combined <- escape_combined[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
escape_combined <- escape_combined[
  with(escape_combined, order(site, mutation)),
]

####################################################################
##Renaming for reference
####################################################################
wildtype <- data.frame(sapply(reference_combined[1], fun_first))
#colnames(wildtype) <- c("wildtype")
reference_combined$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(reference_combined[1], fun_last))
colnames(mutation) <- c("mutation")
reference_combined$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(reference_combined[1], fun_site))
colnames(site) <- c("site")
reference_combined$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
reference_combined <- reference_combined[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
reference_combined <- reference_combined[
  with(reference_combined, order(site, mutation)),
]

#Clean-up
remove(mutation)
remove(site)
remove(wildtype)

escape_trimmed <- escape_combined[which(escape_combined$mutation != "*"),]
reference_trimmed <- reference_combined[which(reference_combined$mutation != "*"),]

###############################################################################
###Add missing residues
###############################################################################
Wuhan <- read.fasta("N_Wuhan.fasta")  #Nucleocapsid sequence fasta file 

###First reference data
### If there are residues without any data, these have to be added manually
### so that plotting hte heatmap works
#1. Define the range of amino acids present
seq_range <- min(reference_trimmed$site):max(reference_trimmed$site)
missing_aa <- seq_range[!seq_range %in% unique(reference_trimmed$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Wuhan[["N_Wuhan"]][missing_aa]), missing_aa, toupper(Wuhan[["N_Wuhan"]][missing_aa])),
                           mutation = toupper(Wuhan[["N_Wuhan"]][missing_aa]), 
                           wildtype = toupper(Wuhan[["N_Wuhan"]][missing_aa]))
complete_reference <- rbind.fill(reference_trimmed, missing_data)

#Expand the dataset to include NA values for synonymous amino acids (wild type)
all <- complete_reference %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
reference_trimmed = complete_reference %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)

###Next escape data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(escape_trimmed$site):max(escape_trimmed$site)
missing_aa <- seq_range[!seq_range %in% unique(escape_trimmed$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Wuhan[["N_Wuhan"]][missing_aa]), missing_aa, toupper(Wuhan[["N_Wuhan"]][missing_aa])),
                           mutation = toupper(Wuhan[["N_Wuhan"]][missing_aa]), 
                           wildtype = toupper(Wuhan[["N_Wuhan"]][missing_aa]))
complete_escape <- rbind.fill(escape_trimmed, missing_data)

#Expand the dataset to include NA values for synonymous amino acids (wild types)
all <- complete_escape %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
escape_trimmed = complete_escape %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)


###Order reference_trimmed and escape_trimmed by site and mutation
escape_trimmed <- escape_trimmed[
  with(escape_trimmed, order(site, mutation)),
]
reference_trimmed <- reference_trimmed[
  with(reference_trimmed, order(site, mutation)),
]

###################################################################
### Calculate abundance: n/N
###################################################################
escape_trimmed$abundance <- escape_trimmed$depth/sum(escape_trimmed$depth, na.rm=TRUE)
reference_trimmed$abundance <- reference_trimmed$depth/sum(reference_trimmed$depth, na.rm=TRUE)

###Then adjust the lowest reference abundance values
#Use 95 percentile
#This helps remove exaggerated escape fractions caused by dividing by a very small number
cutoff  <- quantile(reference_trimmed$abundance, probs = .05, na.rm = T)
reference_trimmed$abundance <- ifelse(reference_trimmed$abundance<cutoff, cutoff, reference_trimmed$abundance)

##Fishers exact test
esc_fisher <- escape_trimmed[,1:5]
esc_fisher$ref_depth <- reference_trimmed$depth
esc_sum <- sum(na.omit(esc_fisher$depth))
ref_sum <- sum(na.omit(esc_fisher$ref_depth))

p_values <- data.frame(matrix(ncol = 1, nrow = nrow(esc_fisher)))
colnames(p_values) <- "p"
#Calculate p values for each comparison (using Fisher's exact test)
for (i in 1:nrow(esc_fisher)){
  if (!anyNA(esc_fisher[i,])){
    f_ <- matrix(unlist(c(esc_fisher[i,5], 
                          esc_fisher[i,6],
                          esc_sum-esc_fisher[i,5],
                          ref_sum-esc_fisher[i,6])),2)
    p_ <- fisher.test(f_)
    p_values$p[i] <-p_$p.value
  }
}

p_values$p_adj <- p.adjust(p_values$p, "bonferroni")
#Some p values are so small they are interpreted as 0; set these to 1e-300
p_values$p_adj <- ifelse(p_values$p_adj == 0, 1e-300, p_values$p_adj)


###################################################################
### Calculate escape fractions
### And add adjusted p values
###################################################################
escape_trimmed$esc_fraction <- (escape_trimmed$abundance/reference_trimmed$abundance)
escape_trimmed$p_adj <- p_values$p_adj
max(escape_trimmed$esc_fraction, na.rm=T)
min(escape_trimmed$esc_fraction, na.rm=T)
#Write out raw escape score before normalizing
escape_raw_out <- na.omit(escape_trimmed[,c(1, ncol(escape_trimmed)-2)])
colnames(escape_raw_out) <- c("Mutation", "EscapeScore")
write.csv(escape_raw_out, file = paste(antibody,"_raw_escape.csv", sep = ""), row.names = FALSE)


###Normalize the data between the 99 and 1 percentiles
max(escape_trimmed$esc_fraction, na.rm=T)
upper_limit  <- quantile(escape_trimmed$esc_fraction, probs = .99, na.rm = T)
lower_limit <- min(escape_trimmed$esc_fraction, na.rm=T)
escape_trimmed$esc_fraction <- ifelse(escape_trimmed$esc_fraction>upper_limit, upper_limit, escape_trimmed$esc_fraction)
escape_trimmed$esc_fraction <- (escape_trimmed$esc_fraction-lower_limit)/(upper_limit-lower_limit)

ggsave(filename = paste(antibody,"_EscapeHistogram_preArcSine.png", sep=""), 
       ggplot(escape_trimmed, aes(x=esc_fraction), title = "Escape Fractions")+
         #geom_histogram(aes(y=..density..), color = "black", fill = "grey", binwidth = 0.025) +
         geom_histogram(color = "black", fill = "grey") +
         xlim(0,quantile(escape_trimmed$esc_fraction, probs = .99, na.rm = T)) +
         xlab("Escape Fraction") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

###############################################################################
### Transform data using arcsine of the squareroot
### This transforms the data into something closer to a normal distribution
###############################################################################
escape_trimmed$esc_asin <- asin(sqrt(escape_trimmed$esc_fraction))
max(escape_trimmed$esc_asin, na.rm=T)
ggsave(filename = paste(antibody, "_EscapeHistogram_postArcSine.png", sep=""), 
       ggplot(escape_trimmed, aes(x=esc_asin), title = "Escape Fractions")+
         #geom_histogram(aes(y=..density..), color = "black", fill = "grey", binwidth = 0.025) +
         geom_histogram(color = "black", fill = "grey") +
         xlim(0,max(escape_trimmed$esc_asin, na.rm=T)) +
         xlab("Escape Fraction") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

################################################################################
###Z Normalization
################################################################################
m <- mean(escape_trimmed$esc_asin, na.rm=T)
s <- sd(escape_trimmed$esc_asin, na.rm=T)
escape_trimmed$esc_z <- (escape_trimmed$esc_asin-m)/s

ggsave(filename = paste(antibody, "_EscapeHistogram_postZnorm.png"), 
       ggplot(escape_trimmed, aes(x=esc_z), title = "Escape Fractions")+
         #geom_histogram(aes(y=..density..), color = "black", fill = "grey", binwidth = 0.025) +
         geom_histogram(color = "black", fill = "grey") +
         xlim(min(escape_trimmed$esc_z, na.rm=T),max(escape_trimmed$esc_z, na.rm=T)) +
         xlab("Escape Fraction") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')


################################################################################
### Generate a matrix for plotting sequence logos
### This matrix is also used to calculate total escape scores
################################################################################
logo_matrix <- matrix(ncol = nrow(escape_trimmed)/20,nrow=20)
row.names(logo_matrix) <- escape_trimmed$mutation[1:20]
colnames(logo_matrix) <- seq(2,ncol(logo_matrix)+1, 1)

for(i in 0:ncol(logo_matrix)-1){
  logo_matrix[,i+1] <- escape_trimmed$esc_z[seq(from = 20*i+1, to=20*i+20, by = 1)]
}

################################################################################
##Plot heat maps
################################################################################

#Generate a dataframe containing the wild type amino acids 
#This will be used to mark wild type in the tiled heat map
tmp <- data.frame(sapply(escape_trimmed[1], fun_first), 
                  sapply(escape_trimmed[1], fun_site), 
                  sapply(escape_trimmed[1], fun_first))
frames <- na.omit(distinct(tmp))
colnames(frames) <- c("mutation", "site", "wildtype")

#Change data type to integer for "site"
#Required for proper heatmap plotting
frames$site <- as.integer(frames$site)
escape_trimmed$site <- as.integer(escape_trimmed$site)

#Set the order for amino acids in the heatmap (by aa property)
polar <- c("H", "C", "S", "T", "N", "Q")
nonpolar <- c("G", "A", "V", "L", "I", "M", "P")
aromatic <- c("F", "Y", "W")
positive <- c("K", "R")
negative <- c("D", "E")

aa_order <- c(negative,positive, polar, nonpolar, aromatic)

#Residues 1-209
start = 0
end = 209
mut_range <- subset(escape_trimmed, site>start & site<end)

#Reverse the order so plotting will happen top to bottom:
#mut_range$mutation <- factor(mut_range$mutation,])
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)

ggsave(filename = paste(antibody,"_EscapeFraction_heatmap01.png", sep=""), 
       ggplot(mut_range, aes(site, mutation, size = -log(p_adj))) + 
         geom_tile(color = "white",
                   fill = '#FCF0F0', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 1) +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_tile(data=frames_range, 
                   size=0,
                   height = 1,
                   fill='white', 
                   colour="white") +
         geom_point(aes(colour = esc_z), alpha=1) +
         scale_colour_distiller(palette = "RdPu", direction = +1,
                                na.value = '#FCF0F0',
                                limits=c(0,max(escape_trimmed$esc_z, na.rm=T))) +
         scale_size(range = c(0, 2.5)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames_range,
                    aes(site, mutation),
                    shape=16,
                    size = 1.25,
                    colour="black") +
         xlab("Nucleocapsid Site") + ylab("Mutation") +
         scale_x_discrete(breaks = levels(mut_range$site)[c(rep(F,3),T,F)]) +
         theme(# Hide panel borders and remove grid lines
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           # Change axis line
           #axis.line = element_line(colour = "black"),
           legend.key.size = unit(0.15, 'inch')
         ),
       width = 16, height = 3, dpi = 300, units = "in", device='png')


#Residues 210-end
start = 210
end = 420
mut_range <- subset(escape_trimmed, site>start & site<end)
mut_range$site <- as.factor(mut_range$site)
frames_range <- subset(frames, site>start & site<end)
frames_range$site <- as.factor(frames_range$site)

ggsave(filename = paste(antibody,"_EscapeFraction_heatmap02.png", sep=""), 
       ggplot(mut_range, aes(site, mutation, size = -log(p_adj))) + 
         geom_tile(color = "white",
                   fill = '#FCF0F0', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 1) +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_tile(data=frames_range, 
                   size=0,
                   height = 1,
                   fill='grey95', 
                   colour="black") +
         geom_point(aes(colour = esc_z), alpha=1) +
         scale_colour_distiller(palette = "RdPu", direction = +1,
                                na.value = '#FCF0F0',
                                limits=c(0,max(escape_trimmed$esc_z, na.rm=T))) +
         scale_size(range = c(0, 2.5)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames_range,
                    aes(site, mutation),
                    shape=16,
                    size = 1.25,
                    colour="black") +
         xlab("Nucleocapsid Site") + ylab("Mutation") +
         scale_x_discrete(breaks = levels(mut_range$site)[c(rep(F,4),T)]) +
         theme(# Hide panel borders and remove grid lines
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           # Change axis line
           #axis.line = element_line(colour = "black"),
           legend.key.size = unit(0.15, 'inch')
         ),
       width = 16, height = 3, dpi = 300, units = "in", device='png')

################################################################################
###Plot a mini-heatmap 
################################################################################
ggsave(filename = paste(antibody,"_EscapeFraction_heatmap_mini.png", sep=""), 
       ggplot(escape_trimmed, aes(site, mutation, size = -log(p_adj))) + 
         geom_tile(color = "white",
                   fill = '#FCF0F0', #azure1 for blueish grey background
                   lwd = 0.1,
                   linetype = 1,
                   alpha = 0.5) +
         geom_tile(data=frames, 
                   size=0,
                   height = 1,
                   fill='white', 
                   colour="black") +
         #scale_y_discrete(limits=rev) +
         ylim(rev(aa_order)) +
         geom_point(aes(colour = esc_z), alpha=1) +
         scale_colour_distiller(palette = "RdPu", direction = +1,
                                na.value = '#FCF0F0',
                                limits=c(0,max(escape_trimmed$esc_z, na.rm=T))) +
         scale_size(range = c(0, 1.5)) +
         #The next section adds lightgrey tiles and black points for 
         #the wild type residues
         geom_point(inherit.aes = FALSE, 
                    data = frames,
                    aes(site, mutation),
                    shape=16,
                    size = 0.75,
                    colour="black") +
         scale_x_discrete(breaks = levels(escape_trimmed$site)[c(rep(F,8),T,F)]) +
         theme_void() +
         theme(legend.position = "none") +
         scale_fill_distiller(palette = "RdPu", direction = +1, 
                              na.value = '#FCF0F0',
                              limits=c(0,max(escape_trimmed$esc_z, na.rm=T))),
       width = 24, height = 1.6, dpi = 300, units = "in", device='png')

################################################################################
### Plot Total Escape histogram
################################################################################
total_escape <- as.data.frame(colSums(logo_matrix, na.rm=T))
total_escape$site <- seq(2,419,1)
colnames(total_escape) <- c("escape", "site")

ggsave(filename = paste(antibody,"_Total_escape.png", sep=""), 
       ggplot(total_escape, aes(site, escape)) + 
         geom_bar(stat="identity") + 
         xlim(2,419) +
         theme_minimal() + 
         theme_void(),
       width = 12, height = .4, dpi = 300, units = "in", device='png')


################################################################################
###Plot histogram of escape scores
################################################################################
ggsave(filename = paste(antibody, "_TotalEscapeHistogram.png", sep=""), 
       ggplot(total_escape, aes(x=escape), title = "Escape Fractions")+
         #geom_histogram(aes(y=..density..), color = "black", fill = "grey", binwidth = 0.25) +
         geom_histogram(color = "black", fill = "grey") +
         #xlim(0,1) +
         geom_vline(aes(xintercept=m_total + 4*s_total),
                    linetype="dashed", color = "red") + 
         #geom_density(alpha=.1, fill="#FF6666") +
         #geom_text(aes(label = mu), x = 3, y = 3, vjust = "inward", hjust = "inward") +
         #ggtitle("Escape Fractions") +
         xlab("Total Escape") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')
