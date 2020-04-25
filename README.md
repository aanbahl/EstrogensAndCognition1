EstrogensAndCognition1R

This project analyses data from the estrogens and cognition study (E&C). The aim of this project is to analyse data to discern if there is a gene x environment interaction that affects memory in women who have undergone a bilateral salpingo oophorectomy. Specifically, we are looking to see how the body's estrogen levels and Val/Met polymorphisms of the COMT (Val158Met) and BDNF (Val66Met) genes impact memory. 

Step 1: Set the working directory 
```
getwd()
#setwd("C:/Users/RA/Desktop/Aanya Data")
setwd("/Volumes/EC_RAs/Aanya - COMTBDNF/Aanya Data/")
dir()
genotype_data <- read.csv("Genotype Data.csv")
neuro_data <- read.csv("SOLVED Neuro Data January 22 2020.csv")
dataset <- merge(genotype_data, neuro_data, by="ID")
```

Step 2: Install packages
I find that installing packages in the beginning of your code makes it cleaner 
```
# Install packages -----
install.packages("carData")
install.packages ("car")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("FSA")
install.packages("ggpubr")
install.packages("HardyWeinberg")
install.packages("afex")
install.packages("knitr")
install.packages("arsenal")
install.packages("magrittr")
```
Don't forget to load your packages for each new session!
```
# Load packages ---- 
library(carData)
library(car)
library(dplyr)
library(tidyverse)
library(FSA)
library(emmeans)
library(ggpubr)
library(HardyWeinberg)
library(afex)
library(knitr)
library(arsenal) 
library(magrittr) 
```
Step 3: Assign your group names! 

The BSO category should include women who have undergone a BSO and have never taken estradiol therapy
```
dataset$group_membership[dataset$Group == 1 & dataset$HRTever == 0] <- "BSO"
```
The BSO-E2 category should include women who have undergone a BSO, have ever taken hormone therapy, and are taking it now
```
dataset$group_membership[dataset$Group == 1 & dataset$HRTever == 1 & dataset$E2now == 1] <- "BSO-E2"
```
The aged matched controls (AMC) group should iclude women who are age matched controls and have never taken hormone therapy
```
dataset$group_membership[dataset$Group == 3 & dataset$HRTever == 0] <- "AMC"
```
However, people who have taken hormone therapy in the past should be excluded from teh AMC and BSO categories. Create a new group for them.
```
dataset$group_membership[dataset$Group == 1 & dataset$HRTever == 1 & dataset$HRTnow == 0] <- "BSO_HRT_Past"
```
Use the subset() function to select and exclude variables. The ! operator means 'is not equal to'. The following line of code is used to look for all the data in the data set (by use of subset()) that is not equal to women who have undergone hormone therapy in the past (BSO_HRT_Past). You are then creating a new dataset from these values, called dataset_new
```
dataset_new <- subset(dataset, !group_membership=="BSO_HRT_Past")
```
Remove all the N/A cells from the data with the !is.na() function.
```
dataset_new <- dataset_new[!is.na(dataset_new$group_membership),]
dataset_new <- dataset_new[!is.na(dataset_new$ID),]
dataset_new <- dataset_new[!is.na(dataset_new$SPWM_WME_T1),]
View(dataset_new)
```
If necessary, you can also set your white spaces into N/As.
```
dataset_new[dataset_new == ""] <- NA
View(dataset_new)
```

# Adding BDNF Met Carriers to the dataset ----
dataset_new$BDNF_combined[dataset_new$BDNF == "Val"] <- "Val"
dataset_new$BDNF_combined[dataset_new$BDNF == "Het" | dataset_new$BDNF == "Met"] <- "Met_carrier"

#use clean_data as a dataset

clean_data <- read.csv("clean_data.csv")

clean_data <-  dataset_new 

# Group Counts for demographics table ---- 

# Use the %>% (pipe operator) to simplify code which can otherwise have a lot of parentheses!

clean_data %>% group_by(group_membership) %>% tally()
clean_data %>% group_by(group_membership, COMT) %>% tally()
# get a tally for BDNF numbers 
clean_data %>% group_by(group_membership, BDNF_combined) %>% tally()
# get a tally for the APOE #s
clean_data %>% group_by(group_membership, e4) %>% tally()
#get a tally for the PdG numbers 
clean_data %>% group_by(group_membership, PdG) %>% tally()
#get a tally for the cancner numbers
clean_data %>% group_by(group_membership, Cancer_Treatments) %>% tally()


# Standard error = SD / sqrt(n-1)
# Introducing the sem calculation into the summarize funtion makes your life a lot easier - eliminates many lines of code!

View(clean_data)
clean_data %>% group_by(group_membership) %>% summarize(mean=mean(E1G, na.rm = T),
                                                        sd=sd(E1G, na.rm = T),
                                                        n=length(Group),
                                                        sem=sd(E1G, na.rm = T)/sqrt(n-1)) 


clean_data %>% group_by(group_membership) %>% summarize(mean=mean(MenopauseAge, na.rm = T), 
                                                        sd=sd(MenopauseAge, na.rm = T),
                                                        n=length(Group),
                                                        sem=sd(MenopauseAge, na.rm = T)/sqrt(n-1))

clean_data %>% group_by(group_membership) %>% summarize(mean = mean(EduYears, na.rm = T),
                                                        sd =sd(EduYears, na.rm = T), 
                                                        n=length(Group),
                                                        sem=sd(EduYears, na.rm = T)/sqrt(n-1)) 

clean_data %>% group_by(group_membership) %>% summarize(mean=mean(Age, na.rm = T),
                                                        sd=sd(Age, na.rm = T),
                                                        n=length(Group),
                                                        sem=sd(Age, na.rm = T)/sqrt(n-1)) 

clean_data %>% group_by(group_membership) %>% summarize(mean=mean(BMI, na.rm = T),
                                                        sd=sd(BMI, na.rm = T),
                                                        n=length(Group),
                                                        sem=sd(BMI, na.rm = T)/sqrt(n-1)) 


clean_data %>% group_by(group_membership) %>% summarize(mean(Cancer_Treatments, na.rm = T), sd(Cancer_Treatments, na.rm = T)) 
clean_data %>% group_by(group_membership) %>% summarize(mean=mean(PdG, na.rm = T),
                                                        sd=sd(PdG, na.rm = T),
                                                        n=length(Group),
                                                        sem=sd(PdG, na.rm = T)/sqrt(n-1)) 

# Group comparison tests for demo ---- 
# One factorial ANOVA

#Age
attach(clean_data)
Age.aov<-aov(Age~group_membership)
summary(Age.aov)
leveneTest(Age~group_membership) #assess the equality of the variances

# Extract the residuals
Age_aovresiduals <- residuals(object = Age.aov )
shapiro.test(x = Age_aovresiduals )
#Does pass Levene's Test and Shapiro

#No. With History of Cancer Treatment
Cancer_Treatment.aov <- aov(Cancer_Treatments~group_membership, data = clean_data)
summary(Cancer_Treatment.aov)
leveneTest(Cancer_Treatments~group_membership) #Assess the quality of the variances

#Extract the residuals
Cancer_Treatments_aovresiduals <- residuals(object = Cancer_Treatment.aov)
shapiro.test(x = Cancer_Treatments_aovresiduals)

#Is significant so run a Kurskal-Wallis Rank Sum Test
kruskal.test(Cancer_Treatments~group_membership)

# Post-hoc for KW is the Dunn Test 
clean_data$group_membership = factor(clean_data$group_membership,
                                     levels=c("AMC", "BSO", "BSO-E2"))
### Dunn test
PT = dunnTest(Cancer_Treatments ~ group_membership,
              data=clean_data,
              method="bonferroni")    
PT


#MenopauseAge test
Menopause.aov<-aov(MenopauseAge~group_membership)
summary(Menopause.aov)
leveneTest(MenopauseAge~group_membership) #assess the equality of the variances

# Extract the residuals
Menopause_aovresiduals <- residuals(object = Menopause.aov )
shapiro.test(x = Menopause_aovresiduals )

#Years of Education
EduYears.aov<-aov(EduYears~group_membership)
summary(EduYears.aov)
leveneTest(EduYears~group_membership) #assess the equality of the variances

# Extract the residuals
EduYears_aovresiduals <- residuals(object = EduYears.aov )
shapiro.test(x = EduYears_aovresiduals )

#BMI
BMI.aov<-aov(BMI~group_membership)
summary(BMI.aov)
leveneTest(BMI~group_membership) #assess the equality of the variances

# Extract the residuals
BMI_aovresiduals <- residuals(object = BMI.aov )
shapiro.test(x = BMI_aovresiduals )

#P value less than 0.05 so perform Kruskal-Wallis Rank Sum Test
kruskal.test(BMI~group_membership)

#PdG
PdG.aov <- aov(PdG~group_membership)
summary(PdG.aov)

#E1G
E1G.aov<-aov(E1G~group_membership)
summary(E1G.aov)
leveneTest(E1G~group_membership) #assess the equality of the variances

# Extract the residuals
E1G_aovresiduals <- residuals(object = E1G.aov )
shapiro.test(x = E1G_aovresiduals )

#P value smaller than 0.05 so run a Kruskal-Wallis Rank Sum Test
kruskal.test(E1G~group_membership)

# Post-hoc for KW is the Dunn Test 
clean_data$group_membership = factor(clean_data$group_membership,
                                     levels=c("AMC", "BSO", "BSO-E2"))
### Dunn test
PT = dunnTest(E1G ~ group_membership,
              data=clean_data,
              method="bonferroni")    
PT
#Make graphs

#Compute ANOVA between two independent variables, genotype and group membership
#Since the samples sizes in the groups are uneven, you ned to do an unbalanced two way ANOVA

#When you run a linear model, there are some missing values in the COMT column. Remove those missing values
clean_data$COMT2 <- factor(clean_data$COMT, exclude = "", levels = c("Met","Het","Val"))

#Remove missing values for in the BDNF column too
clean_data$BDNF2 <- factor(clean_data$BDNF, exclude = "", c("Met_carrier", "Val"))

#First, check the contrasts. What are the contrasts? It shows how one group is compared to another group. 
#Each type of contrast will give you a different sum of squares. The default setting is contr.treatment(for unordered data) and contr.poly(for ordered data)
#when you open R and you check the contrasts, it'll show you treatment and poly. 
#You need to set it to contr.sum and contr.poly. You can't use treatment because the tests wont give you sensible results. This is in the help file(According to the person who made it in car)

#The data is already sorted intro reference groups. The reference group is what the other groups are compared to. Use the following code to change the reference groups
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("AMC","BSO-E2","BSO")) #to make AMC the reference group
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("BSO","BSO-E2","AMC")) #to make BSO the reference group


knitr::kable(nice(ancova_RAVLTlearn))%>%  # this code allows you to save out your table directly in a word doc 
  write2word(file="RAVLTlearn.doc", quiet = TRUE) # see this website, https://cran.r-project.org/web/packages/arsenal/vignettes/write2.html#introduction


# Two factorial ANOVAs -----

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(SPWM_WME_T1 ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data), type = 3)

#Since results are significant, run post hoc test
#The TukeyHSD function is only for aov, it cannot be used on the Anova function
#Run a linear regression to see where the differences in the groups is
#Don't forget to change the contrasts back to treatment and poly! (the default)

options(contrasts = c("contr.treatment", "contr.poly"))
COMT_SPWM <- lm(SPWM_WME_T1 ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data)
summary(COMT_SPWM)

#Results: 
#COMT2Val: The difference between Val and Het is significant in AMCs
#group_membershipBSO-E2: The difference between BSO-E2 and AMCs is significant for Het 
#COMT2Val:group_membershipBSO-E2: The interaction term, the difference of the differences

#My interpretation: difference between Vals in BSO and BSO-E2 (in this case the constant is the val)

#Het and AMCs are the reference groups here. To check reference groups, use contrasts function
contrasts(clean_data$group_membership)
clean_data$COMT2 <- factor(clean_data$COMT2, levels = c("Het","Met","Val"))
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("AMC","BSO-E2","BSO"))
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("BSO","BSO-E2","AMC"))

#Since the value for AMC is 0 in both BSO and BSO-E2, it is the reference group
contrasts(clean_data$COMT2)
#Since the value for Het is 0 in both Met and Val, it is the reference group 
# Makes more sense if the reference group is Met (we want to compare the Met to Val)

#I want to know if DigitsForward scores have an effect on COMT genotype and group membership

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(DigitsForward ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data), type = 3)

#Not signficant

#I want to know if DigitsBackward scores have an effect on COMT genotype and group membership

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(DigitsBackward ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data), type = 3)

#Not significant

#I want to know if DigitsTotal scores have an effect on COMT genotype and group membership

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(DigitsTotal ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data), type = 3)

#Not significant, but there's a trend in the COMT2 scores 
#Run a linear model maybe

options(contrasts = c("contr.treatment", "contr.poly"))
COMT_DigitsTotal <- lm(DigitsTotal ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data)
summary(COMT_DigitsTotal)

#I want to know if RAVLT_learn scores have an effect on BDNF genotype and group membership

#To change reference groups
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("AMC","BSO-E2","BSO"))
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("BSO","BSO-E2","AMC"))

#To check reference groups
contrasts(clean_data$BDNF_combined)
contrasts(clean_data$group_membership)

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(RAVLT_Learn ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data), type = 3)
#contr.sum(clean_data$BDNF_combined)

options(contrasts = c("contr.treatment", "contr.poly"))
BDNF_RAVLT<- lm(RAVLT_Learn ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data)
summary(BDNF_RAVLT)

#Logical Memory Immediate
options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(LMA_Imm_Verbatim ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data), type = 3)

options(contrasts = c("contr.treatment", "contr.poly"))
BDNF_LMA_Imm_RAVLT<- lm(LMA_Imm_Verbatim ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data)
summary(BDNF_LMA_Imm_RAVLT)

#Logical Memory Delayed
options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(LMA_Del_Verbatim ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data), type = 3)

options(contrasts = c("contr.treatment", "contr.poly"))
BDNF_LMA_Del_RAVLT<- lm(LMA_Del_Verbatim ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data)
summary(BDNF_LMA_Del_RAVLT)

#To change the reference group
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("AMC","BSO-E2","BSO"))
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("BSO","BSO-E2","AMC"))


#making  graphs ----


#Make a graph for the E1G levels
#make box plot for E1G levels
#First make a graph theme that can be applied to all your graphs
E1GTheme <- theme(plot.title= element_text(family = "Palatino", face = "bold")) #determining font

E1GPlot <- ggplot(clean_data, aes(X=group_membership, y = E1G, fill=group_membership)) + geom_boxplot(outlier.size = .8, lwd=1)
print(E1GPlot + E1GTheme + labs(title="17-β-Estradiol Levels in Different Women Groups", 
                                y = "Estrone-3-Glucoronide (ng/ml)", x="Group Membership",
                                fill = "Group Membership"))

#Box plot with jittered points
#Change outline colours by groups
#Use a custom colour palette
E1G_Boxplot <- ggboxplot(clean_data, x = "group_membership", y = "E1G",
                         color = "group_membership", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                         add = "jitter", shape = "group_membership", legend = "none")
E1GTheme <- theme(plot.title= element_text(family = "Helvetica", size = (13), hjust = 0.5))

my_comparisons <- list( c("BSO", "BSO-E2"), c("BSO-E2", "AMC"), c("BSO", "AMC") ) #if you want to include Kruskal-Wallis values 

print(E1G_Boxplot
      + E1GTheme 
      + stat_compare_means(comparisons = my_comparisons) #if you want to include Kruskal-Wallis values
      + stat_compare_means(label.y = 165) #if you want to include Kruskal Wallis values
      + labs(title="17-β-Estradiol Levels in Different Women Groups", 
             y = "Estrone-3-Glucoronide (ng/ml)", x="Group Membership",
             fill = "Group Membership"))

#Violin plot
E1GViolinplot <- ggviolin(clean_data, x = "group_membership", y = "E1G", fill = "group_membership",
                          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                          add.params = list(fill = "white"),
                          add = "boxplot", legend = "none")

E1GTheme <- theme(plot.title= element_text(family = "Helvetica", size = (13), hjust = 0.5))

my_comparisons <- list( c("BSO", "BSO-E2"), c("BSO-E2", "AMC"), c("BSO", "AMC") ) #if you want to include Kruskal-Wallis values 

print(E1GViolinplot 
      + E1GTheme
      + stat_compare_means(comparisons = my_comparisons) #if you want to include Kruskal-Wallis values
      + stat_compare_means(label.y = 165) #if you want to include Kruskal Wallis values
      + labs(title="17-β-Estradiol Levels in Different Women Groups", 
             y = "Estrone-3-Glucoronide (ng/ml)", x="Group Membership",
             fill = "Group Membership"))

#Make a bar graph for COMT and SPWM

BarTheme <- theme(plot.title= element_text(family = "Helvetica", size = (13),hjust = 0.5),
                  panel.background = element_rect(fill = "white", size = 4),
                  axis.line = element_line(size = 0.5, colour = "black"),
                  axis.ticks = element_line(colour = "black", size = 0.5))


SPWMBar <- ggplot(clean_data, aes(x=group_membership, y=SPWM_WME_T1, fill=COMT2)) +
  stat_summary(fun=mean,  geom="bar", position=position_dodge(width=.9)) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge(width=.9), width=.3) + # need to set width=.9, or else will squish all the errorbars together
  labs(title = "Spatial Working Memory Errors in COMT Genotypes", x = "Group Membership", y = "Working Memory Errors") +
  scale_fill_discrete(name = "COMT Genotype", labels = c("Val/Met", "Met/Met", "Val/Val"))

print(SPWMBar + BarTheme)

#?ggpar

#Make a bar graph for COMT and DigitsForward
DigitsForwardBar <- ggplot(clean_data, aes(x=group_membership, y=DigitsForward, fill=COMT2)) +
  stat_summary(fun=mean,  geom="bar", position=position_dodge(width=.9)) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge(width=.9), width=.3) + # need to set width=.9, or else will squish all the errorbars together
  labs( x = "Group Membership", y = "Mean Score") +
  scale_fill_discrete(name = "COMT Genotype", labels = c("Val/Met", "Met/Met", "Val/Val"))

print(DigitsForwardBar + BarTheme)

#Make a bar graph for COMT and DigitsBackwards
DigitsBackwardBar <- ggplot(clean_data, aes(x=group_membership, y=DigitsBackward, fill=COMT2)) +
  stat_summary(fun=mean,  geom="bar", position=position_dodge(width=.9)) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge(width=.9), width=.3) + # need to set width=.9, or else will squish all the errorbars together
  labs(x = "Group Membership", y = "Mean Score") +
  scale_fill_discrete(name = "COMT Genotype", labels = c("Val/Met", "Met/Met", "Val/Val"))

print(DigitsBackwardBar + BarTheme)



#Make a bar graph for BDNF and RAVLT_Learn
RAVLT_LearnBar <-  ggplot(clean_data, aes(x=group_membership, y=RAVLT_Learn, fill=BDNF_combined)) +
  stat_summary(fun=mean,  geom="bar", position=position_dodge(width=.9)) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge(width=.9), width=.3) + # need to set width=.9, or else will squish all the errorbars together
  labs(x = "Group Membership", y = "RAVLT Percent Learn (A5/A6 x 100%)") +
  scale_fill_discrete(name = "BDNF Genotype", labels = c("Met Carrier", "Val/Val"))

print(RAVLT_LearnBar + BarTheme)

#Make a bar graph for BDNF and LM_Imm_Verbatim
LMA_Imm_VerbatimBar <-  ggplot(clean_data, aes(x=group_membership, y=LMA_Imm_Verbatim, fill=BDNF_combined)) +
  stat_summary(fun=mean,  geom="bar", position=position_dodge(width=.9)) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge(width=.9), width=.3) + # need to set width=.9, or else will squish all the errorbars together
  labs( x = "Group Membership", y = "Immediate Verbatim Mean Score") +
  scale_fill_discrete(name = "BDNF Genotype", labels = c("Met Carrier", "Val/Val"))

print(LMA_Imm_VerbatimBar + BarTheme)

#Make a bar graph for BDNF and LM_Imm_Verbatim
LMA_Del_VerbatimBar <-  ggplot(clean_data, aes(x=group_membership, y=LMA_Del_Verbatim, fill=BDNF_combined)) +
  stat_summary(fun=mean,  geom="bar", position=position_dodge(width=.9)) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge(width=.9), width=.3) + # need to set width=.9, or else will squish all the errorbars together
  labs(x = "Group Membership", y = "Delayed Verbatim Mean Score") +
  scale_fill_discrete(name = "BDNF Genotype", labels = c("Met Carrier", "Val/Val"))

print(LMA_Del_VerbatimBar + BarTheme)

#Hardy Weinberg Test for COMT AMCs
x <- c(MM=4,VM=11,VV=7) 
HW.test <- HWChisq(x, verbose=TRUE)
#P value 0.79 therefore population in HWE, X2 = 0.0689

#Hardy Weinberg Test for COMT BSOs
x <- c(MM=6,VM=6,VV=5) 
HW.test <- HWChisq(x, verbose=TRUE)
#P value 0.39 therefore population in HWE, X2 = 0.7155

#Hardy Weinberg Test for COMT BSO-E2s
x <- c(MM=4,VM=7,VV=4) 
HW.test <- HWChisq(x, verbose=TRUE)
#P value 0.85 therefore population in HWE, X2 = 0.03

#Hardy Weinberg Test for BDNF AMCs
x <- c(MM=1,VM=3,VV=9) 
HW.test <- HWChisq(x, verbose=TRUE)
#P value 0.78 therefore population in HWE, X2 = 0.07

#Hardy Weinberg Test for BDNF BSOs
x <- c(MM=2,VM=5,VV=15) 
HW.test <- HWChisq(x, verbose=TRUE)
#P value 0.3 therefore population in HWE, X2 = 0.77

#Hardy Weinberg Test for BDNF BSOE2s
x <- c(MM=0,VM=2,VV=11) 
HW.test <- HWChisq(x, verbose=TRUE)
#P value is 0.12 therefore population in HWE, X2 = 2.4

help(HardyWeinberg)

citation("car")
citation("HardyWeinberg")

