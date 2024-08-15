## Libraries
library(car)
library(gdata)
library(gghalves)
library(gridExtra)
library(tidyverse)
library(performance)

## Set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Read in data
data <- read.csv('Data_Output/simulationresults.csv')

## Descriptive Statistics for initial optimized solutions
initial_data <- data %>%
  filter(Joint == "Wrist" & Joint_Pert == 0) %>%
  select(-Joint, -Joint_Pert, -X)

initial_means <- initial_data %>%
  group_by(Sex, Reach_Height) %>%
  summarise(across(everything(), mean))

initial_sds <- initial_data %>%
  group_by(Sex, Reach_Height) %>%
  summarise(across(everything(), sd))

initial_means
initial_sds

write.csv(initial_means, 'StatisticalAnalysis/initial_means.csv')
write.csv(initial_sds, 'StatisticalAnalysis/initial_sds.csv')

rm(initial_means, initial_sds)



## Equality of Variance tests
"
Analysis to answer: 
  - Does there exist an inequity of variance in optimal moments as reach height changes

Test: Levene's test
  - Asses equality of variances for a variable (optimal joint moments) calculated for two or more groups (reach height)
  - Brown-forythe was chosen due to the utilization of the median instead of the mean (like levenes test) due 
  to the robustness of the utilization of the median against non-normal data (as can be seen in the elbow for instance)
  - Interpretation: using alpha of 0.05 (bonferroni corrected), anything less, the variances should not be assumed to be equal
  - Performed overall across all reach heights. If significant, proceeded with pairwise comparisons between pairs of reach heights.

"

# Placeholder storage
variancetest.sex <- character()
variancetest.moment <- character()
variancetest.num <- vector()
variancetest.den <- vector()
variancetest.F <- vector()
variancetest.p <- vector()
variancetest.sig <- vector()

variancetest_pw.sex <- character()
variancetest_pw.moment <- character()
variancetest_pw.ht1 <- vector()
variancetest_pw.ht2 <- vector()
variancetest_pw.num <- vector()
variancetest_pw.den <- vector()
variancetest_pw.F <- vector()
variancetest_pw.p <- vector()

# comparisons and conditions
initial_data$Reach_Height <- as.factor(initial_data$Reach_Height)  # Setting reach height as a factor
sexes <- c('M', 'F')
momentcolumns <- which(colnames(initial_data) %in% c("Norm_RW_Moment", "Norm_RE_Moment", "Norm_RS_Moment", "Norm_B_Moment"))
nvariancetests <- length(sexes)*length(momentcolumns)
reachheights <- unique(initial_data$Reach_Height)
pairwisecomps <- combn(reachheights, 2)
npairwisecomps <- length(pairwisecomps)/2

# Levene's test
for (sex in sexes) {
  initial_data_bysex <- initial_data %>% filter(Sex == sex)
  for (moment in momentcolumns) {
    
    # Overall levene's test across all reach heights (homogeneity of variance)
    tempstatsoutput <- leveneTest(initial_data_bysex[[moment]], initial_data_bysex$Reach_Height, center=median)
    variancetest.sex <- append(variancetest.sex, sex)
    variancetest.moment <- append(variancetest.moment, colnames(initial_data_bysex)[moment])
    variancetest.num <- append(variancetest.num, tempstatsoutput$Df[1])
    variancetest.den <- append(variancetest.den, tempstatsoutput$Df[2])
    variancetest.F <- append(variancetest.F, tempstatsoutput$`F value`[1])
    variancetest.p <- append(variancetest.p, tempstatsoutput$`Pr(>F)`[1])
    
    # Perform pairwise Levene's test between pairs of reach heights if overall Levene's test is significant (w/ corrected alpha)
    if (tempstatsoutput$`Pr(>F)`[1] < 0.05 / nvariancetests) {
      
      variancetest.sig <- append(variancetest.sig, 1)

      for (comparison in 1:npairwisecomps) {
        initial_data_bysexandheight <- initial_data %>% 
          filter(Sex == sex) %>% 
          filter(Reach_Height %in% c(pairwisecomps[1, comparison], pairwisecomps[2, comparison]))
        tempstatsoutput_pw <- leveneTest(initial_data_bysexandheight[[moment]], initial_data_bysexandheight$Reach_Height, center=median)
        variancetest_pw.sex <- append(variancetest_pw.sex, sex)
        variancetest_pw.moment <- append(variancetest_pw.moment, colnames(initial_data_bysex)[moment])
        variancetest_pw.ht1 <- append(variancetest_pw.ht1, pairwisecomps[1, comparison])
        variancetest_pw.ht2 <- append(variancetest_pw.ht2, pairwisecomps[2, comparison])
        variancetest_pw.num <- append(variancetest_pw.num, tempstatsoutput_pw$Df[1])
        variancetest_pw.den <- append(variancetest_pw.den, tempstatsoutput_pw$Df[2])
        variancetest_pw.F <- append(variancetest_pw.F, tempstatsoutput_pw$`F value`[1])
        variancetest_pw.p <- append(variancetest_pw.p, tempstatsoutput_pw$`Pr(>F)`[1])
      }
    }
    else {
      variancetest.sig <- append(variancetest.sig, 0)
    }
  }
}

variancetest <- data.frame(sex = variancetest.sex,
                           moment = variancetest.moment,
                           num = variancetest.num,
                           den = variancetest.den,
                           F = variancetest.F,
                           p = variancetest.p,
                           significant = variancetest.sig)
variancetest

# Bonferroni correction of pairwise tests
npairwisetests <- length(variancetest_pw.p)
alpha_levene_pw <- 0.05 / npairwisetests
variancetest_pw.sig <- numeric(npairwisetests)
variancetest_pw.sig[ which(variancetest_pw.p < alpha_levene_pw) ] <- 1

variancetest_pw <- data.frame(sex = variancetest_pw.sex,
                              moment = variancetest_pw.moment,
                              height1 = variancetest_pw.ht1,
                              height2 = variancetest_pw.ht2,
                              num = variancetest_pw.num,
                              den = variancetest_pw.den,
                              F = variancetest_pw.F,
                              p = variancetest_pw.p,
                              significant = variancetest_pw.sig)
variancetest_pw

write.csv(variancetest, 'StatisticalAnalysis/levene_overall.csv')
write.csv(variancetest_pw, 'StatisticalAnalysis/levene_pairwise.csv')

# Individual analysis check:
# testdata <- initial_data %>% filter(Sex == 'F') %>% filter(Reach_Height %in% c(100, 140))
# leveneTest(testdata$Norm_RS_Moment, testdata$Reach_Height, center=median)

gdata::keep(data, sure = TRUE)



## Descriptive statistics for changes in normalized joint moments at -5/5 degree perturbation
changeinjointmoments <- data %>% 
  group_by(Subject, Reach_Height, Joint) %>%
  mutate(diff_RE_Moment = RE_Moment - RE_Moment[Joint_Pert == 0],
         diff_RS_Moment = RS_Moment - RS_Moment[Joint_Pert == 0],
         diff_B_Moment = B_Moment - B_Moment[Joint_Pert == 0],
         diff_Norm_RE_Moment = Norm_RE_Moment - Norm_RE_Moment[Joint_Pert == 0],
         diff_Norm_RS_Moment = Norm_RS_Moment - Norm_RS_Moment[Joint_Pert == 0],
         diff_Norm_B_Moment = Norm_B_Moment - Norm_B_Moment[Joint_Pert == 0])

changeinjointmoments <- changeinjointmoments %>%
  filter(Joint_Pert %in% c(-5, 5)) %>%
  group_by(Sex, Joint_Pert, Reach_Height, Joint) %>%
  summarize(meandiff_RE_Moment = mean(diff_RE_Moment),
            meandiff_RS_Moment = mean(diff_RS_Moment),
            meandiff_B_Moment = mean(diff_B_Moment),
            stddiff_RE_Moment = sd(diff_RE_Moment),
            stddiff_RS_Moment = sd(diff_RS_Moment),
            stddiff_B_Moment = sd(diff_B_Moment),
            meandiff_Norm_RE_Moment = mean(diff_Norm_RE_Moment),
            meandiff_Norm_RS_Moment = mean(diff_Norm_RS_Moment),
            meandiff_Norm_B_Moment = mean(diff_Norm_B_Moment),
            stddiff_Norm_RE_Moment = sd(diff_Norm_RE_Moment),
            stddiff_Norm_RS_Moment = sd(diff_Norm_RS_Moment),
            stddiff_Norm_B_Moment = sd(diff_Norm_B_Moment))

write.csv(changeinjointmoments, 'StatisticalAnalysis/changeinjointmoments_descriptives.csv')

rm(changeinjointmoments)



## Linear models for evaluating significant changes in normalized joint moments from optimal

linearmodel.sex <- character()
linearmodel.joint <- character()
linearmodel.reachheight <- vector()
linearmodel.deviation <- vector()
linearmodel.moment <- character()
linearmodel.F <- vector()
linearmodel.numdf <- vector()
linearmodel.dendf <- vector()
linearmodel.adjrsquared <- vector()
linearmodel.residualstderror <- vector()
linearmodel.slope_coef <- vector()
linearmodel.slope_stde <- vector()
linearmodel.slope_tvalue <- vector()
linearmodel.slope_pvalue <- vector()

sexes <- unique(data$Sex)
reachheights <- unique(data$Reach_Height)
joints <- unique(data$Joint)
momentcolumns <- which(colnames(data) %in% c("Norm_RE_Moment", "Norm_RS_Moment", "Norm_B_Moment"))
deviations <- c(-1, 1)

for (sex in sexes) {
  for (reachheight in reachheights) {
    for (joint in joints) {
      for (deviation in deviations) {
        
        modeldata <- data %>%
          filter(Sex == sex) %>% filter(Reach_Height == reachheight) %>% filter(Joint == joint)
        
        
        if (deviation < 0) {modeldata <- modeldata %>% filter(Joint_Pert <= 0)}
        else {modeldata <- modeldata %>% filter(Joint_Pert >= 0)}
        
        for (moment in momentcolumns) {
          
          model <- lm(as.formula(paste(colnames(modeldata)[moment], "~", paste(colnames(modeldata)[8]), sep = "" )), data = modeldata)
          summarymodel <- summary(model)
          
          linearmodel.sex <- append(linearmodel.sex, sex)
          linearmodel.joint <- append(linearmodel.joint, joint)
          linearmodel.reachheight <- append(linearmodel.reachheight, reachheight)
          linearmodel.deviation <- append(linearmodel.deviation, deviation)
          linearmodel.moment <- append(linearmodel.moment, colnames(modeldata)[moment])
          linearmodel.F <- append(linearmodel.F, summarymodel$fstatistic[1])
          linearmodel.numdf <- append(linearmodel.numdf, summarymodel$fstatistic[2])
          linearmodel.dendf <- append(linearmodel.dendf, summarymodel$fstatistic[3])
          linearmodel.adjrsquared <- append(linearmodel.adjrsquared, summarymodel$adj.r.squared)
          linearmodel.residualstderror <- append(linearmodel.residualstderror, summarymodel$sigma)
          linearmodel.slope_coef <- append(linearmodel.slope_coef, summarymodel$coefficients[2,1])
          linearmodel.slope_stde <- append(linearmodel.slope_stde, summarymodel$coefficients[2,2])
          linearmodel.slope_tvalue <- append(linearmodel.slope_tvalue, summarymodel$coefficients[2,3])
          linearmodel.slope_pvalue <- append(linearmodel.slope_pvalue, summarymodel$coefficients[2,4])
          
          #check_model(model)

        }
      }
    }
  }
}

linearmodelresults <- data.frame(sex = linearmodel.sex,
                                 joint = linearmodel.joint,
                                 reachheight = linearmodel.reachheight,
                                 deviation = linearmodel.deviation,
                                 moment = linearmodel.moment,
                                 Fstat = linearmodel.F,
                                 numdf = linearmodel.numdf,
                                 dendf = linearmodel.dendf,
                                 adjrsquared = linearmodel.adjrsquared,
                                 residualstderror = linearmodel.residualstderror,
                                 slope_coefficient = linearmodel.slope_coef,
                                 slope_standarderror = linearmodel.slope_stde,
                                 slope_tvalue = linearmodel.slope_tvalue,
                                 slope_pvalue = linearmodel.slope_pvalue)

write.csv(linearmodelresults, 'StatisticalAnalysis/linearmodeloutput.csv')


# Individual analysis check:
datacheck <- data %>% filter(Sex == 'F' & Reach_Height == 80 & Joint == 'Wrist' & Joint_Pert >= 0)
modelcheck <- lm(Norm_RE_Moment ~ Joint_Pert, data = datacheck)
summary(modelcheck)
