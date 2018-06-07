#!/usr/bin/env Rscript

#  Copyright (c) 2014-2018, Centre for Genomic Regulation (CRG).
#  Copyright (c) 2014-2018, Jose Espinosa-Carrasco and the respective authors.
#
#  This file is part of Pergola.
#
#  Pergola is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Pergola is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################
### Jose Espinosa-Carrasco NPMMD/CB-CRG Group. May 2018                  ###
#############################################################################
### Mean values for each group of mice                                    ###
### Using bed files with bouts intersected with experimental              ###
### phases or other phases                                                ###
#############################################################################

## Loading libraries
library ("ggplot2")

## Loading parameters for plots
source("https://gist.githubusercontent.com/JoseEspinosa/307430167b05e408ac071b8724701abf/raw/06b26f764953ceea53d334190d7b736308e5141d/param_plot_publication.R")

# Whole experiment
path2files <- "/Users/jespinosa/2018_chrom_hmm/test_feeding/input/inputdir/"
setwd(path2files)
files <- list.files(pattern=paste("tr_.*.bed$", sep=""))

# Habituation
# path2files <- "/Users/jespinosa/git/mouse_chrom_hmm/results_habituation_binning/dir_bed_hab/"
# setwd(path2files)
# files <- list.files(pattern=paste("habituation.tr_.*.bed$", sep=""))

pwd <- getwd()

data.frame_bed <- NULL

for (bed_file in files) {
  
  info = file.info (bed_file)
  if (info$size == 0) { next }
  df <- read.csv(bed_file, header=F, sep="\t")
  phenotype <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[1])
  genotype <- unlist(strsplit(phenotype, split="_",fixed=T))[1]
  mouse <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[2])
  data_type <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[3])
  data_type <- gsub ("food_fat", "HF", data_type)
  data_type <- gsub ("food_sc", "SC", data_type)
  phase <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[4])
  exp_phase <- gsub ("tr_", "", unlist(strsplit(bed_file, split=".",fixed=T))[5])
  df$phenotype <- phenotype
  df$phenotype <- factor(df$phenotype, levels=c("wt_saline", "wt_nicotine", "KO_cb1_saline", "KO_cb1_nicotine"),
                         labels=c("wt_saline", "wt_nicotine", "KO_cb1_saline", "KO_cb1_nicotine"))                                                  
  df$mouse <- mouse
  df$genotype <- genotype
  df$genotype <- factor(df$genotype, levels=c("wt", "KO"), labels=c("wt", "KO"))
  df$mouse <- mouse
  df$data_type <- data_type
  df$phase <- phase
  df$exp_phase <- gsub("_", " ", exp_phase)
  df$group2plot <- paste (phase, data_type)
  data.frame_bed <- rbind(data.frame_bed, df)
}

setwd(pwd)

data.frame_bed$length <- data.frame_bed$V3 - data.frame_bed$V2
x <- data.frame_bed$length 
# data.frame_bed_lt_1000 <- data.frame_bed [which(data.frame_bed$length > 1000), ]

h<-hist (data.frame_bed$length)
# xfit<-seq(min(x),max(x),length=40) 
# yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
# yfit <- yfit*diff(h$mids[1:2])*length(x) 
# lines(xfit, yfit, col="blue", lwd=2)

x_lt_10min <- x[which (x<600)]
hist(x_lt_10min, breaks=seq(0,600,by=60), main="Breaks is vector of breakpoints")
hist(x_lt_10min, breaks=seq(0,600,by=10), main="Breaks is vector of breakpoints")
max(x)

hist(x, breaks=seq(0, 10000,by=100))
hist_bed <- hist(x, breaks=seq(0, 10000,by=10), plot=FALSE)
# hist_bed$counts

plot(hist_bed$counts, main="Breaks is vector of breakpoints",type='h')
# plot(hist_bed$counts, main="Breaks is vector of breakpoints",type='h', lwd=10, lend=2)

x_lt_10min <- data.frame(x=x_lt_10min)

ggplot(x_lt_10min, aes(x=x))+ geom_histogram(breaks = c(0, 30, 120, 600)) +
  stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5) 

## Histogram of values
# head(data.frame_bed)

value <- data.frame_bed$V5
v_max <- 8

value_lt_4 <- value[which (value < v_max)]
value <- data.frame(value=value_lt_4)
# head(value)

## https://stackoverflow.com/questions/24646594/how-to-improve-the-aspect-of-ggplot-histograms-with-log-scales-and-discrete-valu?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
ggplot(value, aes(x=value))+ geom_histogram(breaks = c(0, 0.08, 1.2, 8)) +
  stat_bin(aes(y=..count.., label=..count..), binwidth=1, geom="text", vjust=-.5) 

### plot distribution habituation
## no logarithmic
ggplot(data.frame_bed, aes(x=length)) +
  stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.3) +
  # scale_x_continuous(breaks=c(0, 30, 120,600), expand=c(0,0)) +
  scale_y_continuous(breaks=c(0,125,250,375,500,625,750), expand=c(0,0)) +
  theme_bw()

# log of length
data.frame_bed$length_log <- log(data.frame_bed$length)

ggplot(data.frame_bed, aes(x=length)) +
  stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.1) +
  scale_x_continuous(breaks=c(0, 30, 120,600), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0, 25000, 5000), expand=c(0,0)) +
  theme_bw()

hist_br_100 <- hist(data.frame_bed$length_log, breaks=100 )
# hist(data.frame_bed$length_log, breaks=100 )$mids
# hist(data.frame_bed$length_log, breaks=100 )$counts
# hist_br_100$density
value_min = min(hist_br_100$counts[ which (hist_br_100$mid > 3 & hist_br_100$mid < 4 ) ])
min_values <- hist_br_100$mids[ which (hist_br_100$mid > 3 & 
                                       hist_br_100$mid < 4 &
                                       hist_br_100$counts == value_min) ]
exp(min_values)

ggplot(data.frame_bed, aes(x=length_log)) +
  stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.1) +
  scale_x_continuous(breaks=c(0:9), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0, 25000, 5000), expand=c(0,0)) +
  theme_bw() + geom_vline(xintercept=min_values, linetype="dashed")
  
# exp(3.375) # aprox 30
v_line_intercept <- exp(min_values)

ggplot(data.frame_bed, aes(x=length)) +
  stat_density(aes(y=..count..), color="black", fill="blue", alpha=0.3) +
  scale_x_continuous(breaks=c(0, 30, 120,600), trans="log1p", expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0, 25000, 5000), expand=c(0,0)) +
  theme_bw() + geom_vline(xintercept=v_line_intercept, linetype="dotted")

# hist(data.frame_bed$length, log="x", breaks=200, xlim=c(100,10000))
# quantile(data.frame_bed$length_log)
# 0%      25%      50%      75%     100% 
# 2.079442 4.890349 5.736567 6.359574 8.585226

# exp(2.079442)
# exp(4.890349) # 133
# exp(5.736567) # 309
# exp(6.359574) # 578

## Terciles habituation data
# quantile(data.frame_bed$length_log, probs = seq(0, 1, 0.33))
# exp(2.079442) 
# exp(5.204007) # 182.0001
# exp(6.167516) # 476.9998
# exp(7.557229) # 1914

## Ahora hacer lo mismo con los datos del development
## Terciles all data
# exp(2.079442)
# exp(3.663562) # 182.0001
# exp(5.049856) # 476.9998
# exp(7.557229) # 1914

# ggsave (file=name_file, width=plot_width, height=plot_height, dpi=300)
