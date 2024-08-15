## Libraries
library(tidyverse)
library(gghalves)
library(gridExtra)

## Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Read in data
data <- read.csv('Data_Output/simulationresults.csv')


## Plot themes
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

raincloud_theme = theme(
  text = element_text(size = 10, colour="black"),
  axis.text = element_text(size = 10, colour="black"),
  legend.position = "none",
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line(colour = 'black', size=0.5, linetype='solid'),
  plot.title = element_text(vjust = - 6, hjust = 0.05, face='bold'))

plotheight <- 4
plotwidth <- 4
plotresolution <- 600


## Raincloud plots for initial optimized solutions (Figure 2)
initial_data <- data %>%
  filter(Joint == "Wrist" & Joint_Pert == 0) %>%
  select(-Joint, -Joint_Pert, -X)

fig2_data <- initial_data %>%
  select(Subject, Sex, Reach_Height, RW_Moment, RE_Moment, RS_Moment, B_Moment)

fig2_data$Reach_Height <- as.factor(fig2_data$Reach_Height)
fig2_data <- gather(fig2_data, Joint, Moment, RW_Moment, RE_Moment, RS_Moment, B_Moment)

fig2_data_Mwrist <- fig2_data %>% filter(Joint == "RW_Moment" & Sex == "M")
fig2_data_Melbow <- fig2_data %>% filter(Joint == "RE_Moment" & Sex == "M")
fig2_data_Mshoulder <- fig2_data %>% filter(Joint == "RS_Moment" & Sex == "M")
fig2_data_Mback <- fig2_data %>% filter(Joint == "B_Moment" & Sex == "M")
fig2_data_Fwrist <- fig2_data %>% filter(Joint == "RW_Moment" & Sex == "F")
fig2_data_Felbow <- fig2_data %>% filter(Joint == "RE_Moment" & Sex == "F")
fig2_data_Fshoulder <- fig2_data %>% filter(Joint == "RS_Moment" & Sex == "F")
fig2_data_Fback <- fig2_data %>% filter(Joint == "B_Moment" & Sex == "F")

fig2a <- ggplot(data = fig2_data_Mwrist, aes(x = Reach_Height, y = Moment, fill=Reach_Height)) +
  geom_flat_violin(aes(colour = Reach_Height), position = position_nudge(x = 0.1, y = 0), alpha = .8, adjust=1) +
  geom_point(aes(x = as.numeric(Reach_Height)-0.15, y = Moment, color = Reach_Height), position = position_jitter(width = .05), size = 0.2, alpha = 0.1) +
  geom_boxplot(aes(x = as.numeric(Reach_Height)), width = .1, outlier.shape = NA, alpha = 0.3, colour = "black") +
  expand_limits(x = 5) +
  scale_y_continuous(limits = c(-40, 160), breaks=c(-40, 0, 40, 80, 120, 160)) +
  xlab("Reach Height (cm)") + 
  ylab("Moment (Nm)") +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  theme_bw(base_line_size = 0.25) +
  raincloud_theme +
  coord_flip()

fig2a <- fig2a + ggtitle('Male - Wrist')
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2-Mwrist.tiff", height=plotheight, width=plotwidth, dpi=plotresolution, compression="lzw")

fig2b <- fig2a %+% fig2_data_Melbow + ggtitle('Male - Elbow')
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2-Melbow.tiff", fig2b, height=plotheight, width=plotwidth, dpi=plotresolution, compression="lzw")

fig2c <- fig2a %+% fig2_data_Mshoulder + ggtitle('Male - Shoulder')
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2-Mshoulder.tiff", fig2c, height=plotheight, width=plotwidth, dpi=plotresolution, compression="lzw")

fig2d <- fig2a %+% fig2_data_Mback + ggtitle('Male - Back')
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2-Mback.tiff", fig2d, height=plotheight, width=plotwidth, dpi=plotresolution, compression="lzw")

fig2e <- fig2a %+% fig2_data_Fwrist + ggtitle('Female - Wrist')
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2-Fwrist.tiff", fig2e, height=plotheight, width=plotwidth, dpi=plotresolution, compression="lzw")

fig2f <- fig2a %+% fig2_data_Felbow + ggtitle('Female - Elbow')
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2-Felbow.tiff", fig2f, height=plotheight, width=plotwidth, dpi=plotresolution, compression="lzw")

fig2g <- fig2a %+% fig2_data_Fshoulder + ggtitle('Female - Shoulder')
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2-Fshoulder.tiff", fig2g, height=plotheight, width=plotwidth, dpi=plotresolution, compression="lzw")

fig2h <- fig2a %+% fig2_data_Fback + ggtitle('Female - Back')
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2-Fback.tiff", fig2h, height=plotheight, width=plotwidth, dpi=plotresolution, compression="lzw")

fig2 <- arrangeGrob(fig2a, fig2b, fig2c, fig2d, fig2e, fig2f, fig2g, fig2h, nrow = 2)
ggsave("Figures/Fig2-InitialOptimalSolutions/Fig2.tiff", fig2, height=plotheight*1.5, width=plotwidth*3, dpi=plotresolution*2, compression="lzw")




### Plots for showing changes with joint perturbation

figure3function <- function(data, moment_column, sex, reachheight, joint, color1, color2, color3, 
                            ymin, ymax, yinterval, plottitle, yaxislabel) {
  fig3data <- data %>% 
    filter(Sex == sex & Reach_Height == reachheight & Joint == joint)
  
  #fig3data$Joint_Pert_jitter <- fig3data$Joint_Pert
  #fig3data$Joint_Pert_jitter[fig3data$Joint_Pert_jitter == -5] <- -5 + runif(sum(fig3data$Joint_Pert_jitter == -5, na.rm = TRUE), -0.3, 0.3)
  #fig3data$Joint_Pert_jitter[fig3data$Joint_Pert_jitter == 0] <- 0 + runif(sum(fig3data$Joint_Pert_jitter == 0, na.rm = TRUE), -0.3, 0.3)
  #fig3data$Joint_Pert_jitter[fig3data$Joint_Pert_jitter == 5] <- 5 + runif(sum(fig3data$Joint_Pert_jitter == 5, na.rm = TRUE), -0.3, 0.3)
  fig3data$Joint_Pert_jitter <- jitter(fig3data$Joint_Pert, amount = 0.3)
  
  fig <- ggplot(data = fig3data, aes_string(y = moment_column)) +
    geom_line(aes(x = Joint_Pert_jitter, group = Subject), color = 'lightgray', alpha = 0.1) +
    geom_point(data = fig3data %>% filter(Joint_Pert==-5), aes(x=Joint_Pert_jitter), color = color1, size = 2, alpha = .1) +
    geom_point(data = fig3data %>% filter(Joint_Pert==0), aes(x=Joint_Pert_jitter), color = color2, size = 2, alpha = .1) +
    geom_point(data = fig3data %>% filter(Joint_Pert==5), aes(x=Joint_Pert_jitter), color = color3, size = 2,  alpha = .1) +
    #geom_half_violin(data = fig3data %>% filter(Joint_Pert==-5), aes(x = Joint_Pert), position = position_nudge(x = 11), side = "r", fill = color1, alpha = .5, color = color1, trim = TRUE) +
    #geom_half_violin(data = fig3data %>% filter(Joint_Pert==0), aes(x = Joint_Pert), position = position_nudge(x = 6), side = "r", fill = color2, alpha = .5, color = color2, trim = TRUE) +
    #geom_half_violin(data = fig3data %>% filter(Joint_Pert==5), aes(x = Joint_Pert), position = position_nudge(x = 1), side = "r", fill = color3, alpha = .5, color = color3, trim = TRUE) +
    scale_x_continuous(breaks=c(-5, 0, 5), labels=c("-5", "0", "5"), limits=c(-7, 7)) +
    scale_y_continuous(breaks=seq(ymin, ymax, by = yinterval), limits=c(ymin, ymax)) +
    xlab("Joint Alteration (deg)") + ylab(yaxislabel) +
    theme_bw(base_line_size = 0.25) +
    ggtitle(plottitle) +
    raincloud_theme +
    theme(plot.title = element_text(vjust = 0, hjust = 0.5))
  
  fig
}

data_anglesfiltered <- data %>%
  filter(Wrist_Angle >= -25.1 & Wrist_Angle <= 40.1) %>%
  filter(Elbow_Angle >= 29.9 & Elbow_Angle <= 180.1) %>%
  filter(Shoulder_Angle >= -180.1 & Shoulder_Angle <= 30.1) %>%
  filter(Trunk_Angle >= -30.1 & Trunk_Angle <=180.1)

fig3a <- figure3function(data, 'RE_Moment', 'M', 80, 'Wrist', 'green3', 'black', 'darkgreen', 
                         -6, 1, 1, '(A) Male - Wrist Alteration', 'Elbow Moment (Nm)')
fig3b <- figure3function(data, 'RE_Moment', 'M', 80, 'Elbow', 'green3', 'black',' darkgreen', 
                         -6, 1, 1, '(B) Male - Elbow Alteration', 'Elbow Moment (Nm)')
fig3c <- figure3function(data, 'RE_Moment', 'M', 80, 'Shoulder', 'green3', 'black',' darkgreen', 
                         -6, 1, 1, '(C) Male - Shoulder Alteration', 'Elbow Moment (Nm)')
fig3d <- figure3function(data, 'RS_Moment', 'M', 80, 'Wrist', 'firebrick1', 'black', 'firebrick4', 
                         -24, -12, 4, '', 'Shoulder Moment (Nm)')
fig3e <- figure3function(data, 'RS_Moment', 'M', 80, 'Elbow', 'firebrick1', 'black',' firebrick4', 
                         -24, -12, 4, '', 'Shoulder Moment (Nm)')
fig3f <- figure3function(data, 'RS_Moment', 'M', 80, 'Shoulder', 'firebrick1', 'black',' firebrick4', 
                         -24, -12, 4, '', 'Shoulder Moment (Nm)')
fig3g <- figure3function(data, 'B_Moment', 'M', 80, 'Wrist', 'dodgerblue', 'black', 'dodgerblue4', 
                         -20, 220, 40, '', 'Back Moment (Nm)')
fig3h <- figure3function(data, 'B_Moment', 'M', 80, 'Elbow', 'dodgerblue', 'black',' dodgerblue4', 
                         -20, 220, 40, '', 'Back Moment (Nm)')
fig3i <- figure3function(data, 'B_Moment', 'M', 80, 'Shoulder', 'dodgerblue', 'black',' dodgerblue4', 
                         -20, 220, 40, '', 'Back Moment (Nm)')
fig3male <- arrangeGrob(fig3a, fig3b, fig3c, fig3d, fig3e, fig3f, fig3g, fig3h, fig3i, nrow = 3)
#ggsave("Figures/Fig3_6-PerturbedSolutions/fig3-male-80.tiff", fig3male, height=plotheight*2, width=plotwidth*3, dpi=plotresolution, compression="lzw")


fig3j <- figure3function(data, 'RE_Moment', 'F', 80, 'Wrist', 'green3', 'black', 'darkgreen', 
                         -6, 1, 1, '(D) Female - Wrist Alteration', 'Elbow Moment (Nm)')
fig3k <- figure3function(data, 'RE_Moment', 'F', 80, 'Elbow', 'green3', 'black',' darkgreen', 
                         -6, 1, 1, '(E) Female - Elbow Alteration', 'Elbow Moment (Nm)')
fig3l <- figure3function(data, 'RE_Moment', 'F', 80, 'Shoulder', 'green3', 'black',' darkgreen', 
                         -6, 1, 1, '(F) Female - Shoulder Alteration', 'Elbow Moment (Nm)')
fig3m <- figure3function(data, 'RS_Moment', 'F', 80, 'Wrist', 'firebrick1', 'black', 'firebrick4', 
                         -24, -12, 4, '', 'Shoulder Moment (Nm)')
fig3n <- figure3function(data, 'RS_Moment', 'F', 80, 'Elbow', 'firebrick1', 'black',' firebrick4', 
                         -24, -12, 4, '', 'Shoulder Moment (Nm)')
fig3o <- figure3function(data, 'RS_Moment', 'F', 80, 'Shoulder', 'firebrick1', 'black',' firebrick4', 
                         -24, -12, 4, '', 'Shoulder Moment (Nm)')
fig3p <- figure3function(data, 'B_Moment', 'F', 80, 'Wrist', 'dodgerblue', 'black', 'dodgerblue4', 
                         -20, 220, 40, '', 'Back Moment (Nm)')
fig3q <- figure3function(data, 'B_Moment', 'F', 80, 'Elbow', 'dodgerblue', 'black',' dodgerblue4', 
                         -20, 220, 40, '', 'Back Moment (Nm)')
fig3r <- figure3function(data, 'B_Moment', 'F', 80, 'Shoulder', 'dodgerblue', 'black',' dodgerblue4', 
                         -20, 220, 40, '', 'Back Moment (Nm)')
fig3female <- arrangeGrob(fig3j, fig3k, fig3l, fig3m, fig3n, fig3o, fig3p, fig3q, fig3r, nrow = 3)
#ggsave("Figures/Fig3_6-PerturbedSolutions/fig3-female-80.tiff", fig3female, height=plotheight*2, width=plotwidth*3, dpi=plotresolution, compression="lzw")


fig3 <- grid.arrange(fig3male, fig3female, nrow = 1)
ggsave("Figures/Fig3_6-PerturbedSolutions/fig3.tiff", fig3, height=plotheight*2, width=plotwidth*6, dpi=plotresolution, compression="lzw")



fig4a <- figure3function(data, 'RE_Moment', 'M', 100, 'Wrist', 'green3', 'black', 'darkgreen', 
                         -6, 3, 3, '(A) Male - Wrist Alteration', 'Elbow Moment (Nm)')
fig4b <- figure3function(data, 'RE_Moment', 'M', 100, 'Elbow', 'green3', 'black',' darkgreen', 
                         -6, 3, 3, '(B) Male - Elbow Alteration', 'Elbow Moment (Nm)')
fig4c <- figure3function(data, 'RE_Moment', 'M', 100, 'Shoulder', 'green3', 'black',' darkgreen', 
                         -6, 3, 3, '(C) Male - Shoulder Alteration', 'Elbow Moment (Nm)')
fig4d <- figure3function(data, 'RS_Moment', 'M', 100, 'Wrist', 'firebrick1', 'black', 'firebrick4', 
                         -24, 4, 4, '', 'Shoulder Moment (Nm)')
fig4e <- figure3function(data, 'RS_Moment', 'M', 100, 'Elbow', 'firebrick1', 'black',' firebrick4', 
                         -24, 4, 4, '', 'Shoulder Moment (Nm)')
fig4f <- figure3function(data, 'RS_Moment', 'M', 100, 'Shoulder', 'firebrick1', 'black',' firebrick4', 
                         -24, 4, 4, '', 'Shoulder Moment (Nm)')
fig4g <- figure3function(data, 'B_Moment', 'M', 100, 'Wrist', 'dodgerblue', 'black', 'dodgerblue4', 
                         -40, 120, 40, '', 'Back Moment (Nm)')
fig4h <- figure3function(data, 'B_Moment', 'M', 100, 'Elbow', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 120, 40, '', 'Back Moment (Nm)')
fig4i <- figure3function(data, 'B_Moment', 'M', 100, 'Shoulder', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 120, 40, '', 'Back Moment (Nm)')
fig4male <- arrangeGrob(fig4a, fig4b, fig4c, fig4d, fig4e, fig4f, fig4g, fig4h, fig4i, nrow = 3)
#ggsave("Figures/Fig3_6-PerturbedSolutions/fig4-male-100.tiff", fig4male, height=plotheight*2, width=plotwidth*3, dpi=plotresolution, compression="lzw")


fig4j <- figure3function(data, 'RE_Moment', 'F', 100, 'Wrist', 'green3', 'black', 'darkgreen', 
                         -6, 3, 3, '(D) Female - Wrist Alteration', 'Elbow Moment (Nm)')
fig4k <- figure3function(data, 'RE_Moment', 'F', 100, 'Elbow', 'green3', 'black',' darkgreen', 
                         -6, 3, 3, '(E) Female - Elbow Alteration', 'Elbow Moment (Nm)')
fig4l <- figure3function(data, 'RE_Moment', 'F', 100, 'Shoulder', 'green3', 'black',' darkgreen', 
                         -6, 3, 3, '(F) Female - Shoulder Alteration', 'Elbow Moment (Nm)')
fig4m <- figure3function(data, 'RS_Moment', 'F', 100, 'Wrist', 'firebrick1', 'black', 'firebrick4', 
                         -24, 4, 4, '', 'Shoulder Moment (Nm)')
fig4n <- figure3function(data, 'RS_Moment', 'F', 100, 'Elbow', 'firebrick1', 'black',' firebrick4', 
                         -24, 4, 4, '', 'Shoulder Moment (Nm)')
fig4o <- figure3function(data, 'RS_Moment', 'F', 100, 'Shoulder', 'firebrick1', 'black',' firebrick4', 
                         -24, 4, 4, '', 'Shoulder Moment (Nm)')
fig4p <- figure3function(data, 'B_Moment', 'F', 100, 'Wrist', 'dodgerblue', 'black', 'dodgerblue4', 
                         -40, 120, 40, '', 'Back Moment (Nm)')
fig4q <- figure3function(data, 'B_Moment', 'F', 100, 'Elbow', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 120, 40, '', 'Back Moment (Nm)')
fig4r <- figure3function(data, 'B_Moment', 'F', 100, 'Shoulder', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 120, 40, '', 'Back Moment (Nm)')
fig4female <- arrangeGrob(fig4j, fig4k, fig4l, fig4m, fig4n, fig4o, fig4p, fig4q, fig4r, nrow = 3)
#ggsave("Figures/Fig3_6-PerturbedSolutions/fig4-female-100.tiff", fig4female, height=plotheight*2, width=plotwidth*3, dpi=plotresolution, compression="lzw")

fig4 <- grid.arrange(fig4male, fig4female, nrow = 1)
ggsave("Figures/Fig3_6-PerturbedSolutions/fig4.tiff", fig4, height=plotheight*2, width=plotwidth*6, dpi=plotresolution, compression="lzw")



fig5a <- figure3function(data, 'RE_Moment', 'M', 120, 'Wrist', 'green3', 'black', 'darkgreen', 
                         -4, 16, 4, '(A) Male - Wrist Alteration', 'Elbow Moment (Nm)')
fig5b <- figure3function(data, 'RE_Moment', 'M', 120, 'Elbow', 'green3', 'black',' darkgreen', 
                         -4, 16, 4, '(B) Male - Elbow Alteration', 'Elbow Moment (Nm)')
fig5c <- figure3function(data, 'RE_Moment', 'M', 120, 'Shoulder', 'green3', 'black',' darkgreen', 
                         -4, 16, 4, '(C) Male - Shoulder Alteration', 'Elbow Moment (Nm)')
fig5d <- figure3function(data, 'RS_Moment', 'M', 120, 'Wrist', 'firebrick1', 'black', 'firebrick4', 
                         -12, 12, 4, '', 'Shoulder Moment (Nm)')
fig5e <- figure3function(data, 'RS_Moment', 'M', 120, 'Elbow', 'firebrick1', 'black',' firebrick4', 
                         -12, 12, 4, '', 'Shoulder Moment (Nm)')
fig5f <- figure3function(data, 'RS_Moment', 'M', 120, 'Shoulder', 'firebrick1', 'black',' firebrick4', 
                         -12, 12, 4, '', 'Shoulder Moment (Nm)')
fig5g <- figure3function(data, 'B_Moment', 'M', 120, 'Wrist', 'dodgerblue', 'black', 'dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig5h <- figure3function(data, 'B_Moment', 'M', 120, 'Elbow', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig5i <- figure3function(data, 'B_Moment', 'M', 120, 'Shoulder', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig5male <- arrangeGrob(fig5a, fig5b, fig5c, fig5d, fig5e, fig5f, fig5g, fig5h, fig5i, nrow = 3)
#ggsave("Figures/Fig3_6-PerturbedSolutions/fig5-male-120.tiff", fig5male, height=plotheight*2, width=plotwidth*3, dpi=plotresolution, compression="lzw")


fig5j <- figure3function(data, 'RE_Moment', 'F', 120, 'Wrist', 'green3', 'black', 'darkgreen', 
                         -4, 16, 4, '(D) Female - Wrist Alteration', 'Elbow Moment (Nm)')
fig5k <- figure3function(data, 'RE_Moment', 'F', 120, 'Elbow', 'green3', 'black',' darkgreen', 
                         -4, 16, 4, '(E) Female - Elbow Alteration', 'Elbow Moment (Nm)')
fig5l <- figure3function(data, 'RE_Moment', 'F', 120, 'Shoulder', 'green3', 'black',' darkgreen', 
                         -4, 16, 4, '(F) Female - Shoulder Alteration', 'Elbow Moment (Nm)')
fig5m <- figure3function(data, 'RS_Moment', 'F', 120, 'Wrist', 'firebrick1', 'black', 'firebrick4', 
                         -12, 12, 4, '', 'Shoulder Moment (Nm)')
fig5n <- figure3function(data, 'RS_Moment', 'F', 120, 'Elbow', 'firebrick1', 'black',' firebrick4', 
                         -12, 12, 4, '', 'Shoulder Moment (Nm)')
fig5o <- figure3function(data, 'RS_Moment', 'F', 120, 'Shoulder', 'firebrick1', 'black',' firebrick4', 
                         -12, 12, 4, '', 'Shoulder Moment (Nm)')
fig5p <- figure3function(data, 'B_Moment', 'F', 120, 'Wrist', 'dodgerblue', 'black', 'dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig5q <- figure3function(data, 'B_Moment', 'F', 120, 'Elbow', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig5r <- figure3function(data, 'B_Moment', 'F', 120, 'Shoulder', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig5female <- arrangeGrob(fig5j, fig5k, fig5l, fig5m, fig5n, fig5o, fig5p, fig5q, fig5r, nrow = 3)
#ggsave("Figures/Fig3_6-PerturbedSolutions/fig5-female-120.tiff", fig5female, height=plotheight*2, width=plotwidth*3, dpi=plotresolution, compression="lzw")


fig5 <- grid.arrange(fig5male, fig5female, nrow = 1)
ggsave("Figures/Fig3_6-PerturbedSolutions/fig5.tiff", fig5, height=plotheight*2, width=plotwidth*6, dpi=plotresolution, compression="lzw")




fig6a <- figure3function(data, 'RE_Moment', 'M', 140, 'Wrist', 'green3', 'black', 'darkgreen', 
                         -2, 14, 2, '(A) Male - Wrist Alteration', 'Elbow Moment (Nm)')
fig6b <- figure3function(data, 'RE_Moment', 'M', 140, 'Elbow', 'green3', 'black',' darkgreen', 
                         -2, 14, 2, '(B) Male - Elbow Alteration', 'Elbow Moment (Nm)')
fig6c <- figure3function(data, 'RE_Moment', 'M', 140, 'Shoulder', 'green3', 'black',' darkgreen', 
                         -2, 14, 2, '(C) Male - Shoulder Alteration', 'Elbow Moment (Nm)')
fig6d <- figure3function(data, 'RS_Moment', 'M', 140, 'Wrist', 'firebrick1', 'black', 'firebrick4', 
                         -10, 25, 5, '', 'Shoulder Moment (Nm)')
fig6e <- figure3function(data, 'RS_Moment', 'M', 140, 'Elbow', 'firebrick1', 'black',' firebrick4', 
                         -10, 25, 5, '', 'Shoulder Moment (Nm)')
fig6f <- figure3function(data, 'RS_Moment', 'M', 140, 'Shoulder', 'firebrick1', 'black',' firebrick4', 
                         -10, 25, 5, '', 'Shoulder Moment (Nm)')
fig6g <- figure3function(data, 'B_Moment', 'M', 140, 'Wrist', 'dodgerblue', 'black', 'dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig6h <- figure3function(data, 'B_Moment', 'M', 140, 'Elbow', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig6i <- figure3function(data, 'B_Moment', 'M', 140, 'Shoulder', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig6male <- arrangeGrob(fig6a, fig6b, fig6c, fig6d, fig6e, fig6f, fig6g, fig6h, fig6i, nrow = 3)
#ggsave("Figures/Fig3_6-PerturbedSolutions/fig6-male-140.tiff", fig6male, height=plotheight*2, width=plotwidth*3, dpi=plotresolution, compression="lzw")


fig6j <- figure3function(data, 'RE_Moment', 'F', 140, 'Wrist', 'green3', 'black', 'darkgreen', 
                          -2, 14, 2, '(D) Female - Wrist Alteration', 'Elbow Moment (Nm)')
fig6k <- figure3function(data, 'RE_Moment', 'F', 140, 'Elbow', 'green3', 'black',' darkgreen', 
                          -2, 14, 2, '(E) Female - Elbow Alteration', 'Elbow Moment (Nm)')
fig6l <- figure3function(data, 'RE_Moment', 'F', 140, 'Shoulder', 'green3', 'black',' darkgreen', 
                          -2, 14, 2, '(F) Female - Shoulder Alteration', 'Elbow Moment (Nm)')
fig6m <- figure3function(data, 'RS_Moment', 'F', 140, 'Wrist', 'firebrick1', 'black', 'firebrick4', 
                         -10, 25, 5, '', 'Shoulder Moment (Nm)')
fig6n <- figure3function(data, 'RS_Moment', 'F', 140, 'Elbow', 'firebrick1', 'black',' firebrick4', 
                         -10, 25, 5, '', 'Shoulder Moment (Nm)')
fig6o <- figure3function(data, 'RS_Moment', 'F', 140, 'Shoulder', 'firebrick1', 'black',' firebrick4', 
                         -10, 25, 5, '', 'Shoulder Moment (Nm)')
fig6p <- figure3function(data, 'B_Moment', 'F', 140, 'Wrist', 'dodgerblue', 'black', 'dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig6q <- figure3function(data, 'B_Moment', 'F', 140, 'Elbow', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig6r <- figure3function(data, 'B_Moment', 'F', 140, 'Shoulder', 'dodgerblue', 'black',' dodgerblue4', 
                         -40, 40, 10, '', 'Back Moment (Nm)')
fig6female <- arrangeGrob(fig6j, fig6k, fig6l, fig6m, fig6n, fig6o, fig6p, fig6q, fig6r, nrow = 3)
#ggsave("Figures/Fig3_6-PerturbedSolutions/fig6-female-140.tiff", fig6female, height=plotheight*2, width=plotwidth*3, dpi=plotresolution, compression="lzw")


fig6 <- grid.arrange(fig6male, fig6female, nrow = 1)
ggsave("Figures/Fig3_6-PerturbedSolutions/fig6.tiff", fig6, height=plotheight*2, width=plotwidth*6, dpi=plotresolution, compression="lzw")


## EXPLORATORY CHECK (140 cm reach height) -- bimodal response
testwrist <- data %>% filter(Reach_Height == 140 & Joint == 'Wrist' & Joint_Pert == -5)
testwrist_mass <- ggplot(data = testwrist, aes(x = Subject_Mass, y = Norm_B_Moment, color = Sex, group = Sex)) + geom_point()
testwrist_height <- ggplot(data = testwrist, aes(x = Subject_Height, y = Norm_B_Moment, color = Sex, group = Sex)) + geom_point()

testelbow <- data %>% filter(Reach_Height == 140 & Joint == 'Elbow' & Joint_Pert == -5)
testelbow_mass <- ggplot(data = testelbow, aes(x = Subject_Mass, y = Norm_B_Moment, color = Sex, group = Sex)) + geom_point()
testelbow_height <- ggplot(data = testelbow, aes(x = Subject_Height, y = Norm_B_Moment, color = Sex, group = Sex)) + geom_point()

testshoulder <- data %>% filter(Reach_Height == 140 & Joint == 'Shoulder' & Joint_Pert == -5)
testshoulder_mass <- ggplot(data = testshoulder, aes(x = Subject_Mass, y = Norm_B_Moment, color = Sex, group = Sex)) + geom_point()
testshoulder_height <- ggplot(data = testshoulder, aes(x = Subject_Height, y = Norm_B_Moment, color = Sex, group = Sex)) + geom_point()
