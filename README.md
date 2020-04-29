# Estrogens_And_Cognition_1_R

This project analyses data to discern if there is a gene x environment interaction that affects memory in women who have undergone a bilateral salpingo oophorectomy. Specifically, we are looking to see how the body's estrogen levels and Val/Met polymorphisms of the COMT (Val158Met) and BDNF (Val66Met) genes impact memory. 

Step 1: Set the working directory 
```
getwd()
#setwd("C:/Users/RA/Desktop/Aanya Data")
setwd("/Volumes/EC_RAs/Aanya - COMTBDNF/Aanya Data/")
dir()
```
Step 2: read your files and merge them to create the dataset you want.
```
genotype_data <- read.csv("Genotype Data.csv")
neuro_data <- read.csv("SOLVED Neuro Data January 22 2020.csv")
dataset <- merge(genotype_data, neuro_data, by="ID")
```

Step 3: Install packages

(I find that installing packages in the beginning of your code makes it cleaner) 
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
Step 4: Assign your group names! 

The BSO category should include women who have undergone a BSO and have never taken estradiol therapy.
```
dataset$group_membership[dataset$Group == 1 & dataset$HRTever == 0] <- "BSO"
```
The BSO-E2 category should include women who have undergone a BSO, have ever taken hormone therapy, and are taking it now.
```
dataset$group_membership[dataset$Group == 1 & dataset$HRTever == 1 & dataset$E2now == 1] <- "BSO-E2"
```
The aged matched controls (AMC) group should iclude women who are age matched controls and have never taken hormone therapy.
```
dataset$group_membership[dataset$Group == 3 & dataset$HRTever == 0] <- "AMC"
```
However, people who have taken hormone therapy in the past should be excluded from the AMC and BSO categories. Create a new group for them.
```
dataset$group_membership[dataset$Group == 1 & dataset$HRTever == 1 & dataset$HRTnow == 0] <- "BSO_HRT_Past"
```
Use the ```subset()``` function to select and exclude variables. The ```!``` operator means 'is not equal to'. 

The following line of code looks for all the data in the dataset (by use of subset()) that is not equal to women who have undergone hormone therapy in the past (```BSO_HRT_Past```). You are then creating a new dataset from these values, called ```dataset_new```.
```
dataset_new <- subset(dataset, !group_membership=="BSO_HRT_Past")
```
Remove all the N/A cells from the data with the ```!is.na()``` function.
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
Step 5: The BDNF data size is not that large, so merge the Heterozygotes and the Met Homozygotes into a single Met carrier category. Use the ```|``` operator to code for 'or'
```
# Adding BDNF Met Carriers to the dataset ----
dataset_new$BDNF_combined[dataset_new$BDNF == "Val"] <- "Val"
dataset_new$BDNF_combined[dataset_new$BDNF == "Het" | dataset_new$BDNF == "Met"] <- "Met_carrier"
```
Rename ```dataset_new``` as ```clean_data```
```
clean_data <- dataset_new
clean_data <- read.csv("clean_data.csv")
```

Step 6: Get counts for demographic data table

Operators and Functions:

Use the ```%>%``` (pipe operator) to simplify code which can otherwise have a lot of parentheses!

Use the ```group_by``` function to group the data by a single category/variable

Use the ```tally()``` function to count the number of data points within the particular group you have created
```
clean_data %>% group_by(group_membership) %>% tally()
#get a tally for the COMT numbers
clean_data %>% group_by(group_membership, COMT) %>% tally()
#get a tally for BDNF numbers 
clean_data %>% group_by(group_membership, BDNF_combined) %>% tally()
#get a tally for the APOE #s
clean_data %>% group_by(group_membership, e4) %>% tally()
#get a tally for the PdG numbers 
clean_data %>% group_by(group_membership, PdG) %>% tally()
#get a tally for the cancner numbers
clean_data %>% group_by(group_membership, Cancer_Treatments) %>% tally()
```
Step 7: Calculate the mean and standard error for the demographic groups.

Standard error = SD / sqrt(n-1)

Introducing the ```sem``` calculation into the summarize funtion makes your life a lot easier - eliminates many lines of code!

```
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
```
Step 8: Group comparison tests for the demographic data. Create a one factorial ANOVA for each group to see if there are any significant differences between groups. 
```
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
```
If the differences are significant, run a Kruskal-Wallis Rank Sum Test to clarify that the differences are significant.
```
#Is significant so run a Kurskal-Wallis Rank Sum Test
kruskal.test(Cancer_Treatments~group_membership)
```
The post-hoc test for the Kruskal-Wallis test is the Dunn test.
```
# Post-hoc for KW is the Dunn Test 
clean_data$group_membership = factor(clean_data$group_membership,
                                     levels=c("AMC", "BSO", "BSO-E2"))
### Dunn test
PT = dunnTest(Cancer_Treatments ~ group_membership,
              data=clean_data,
              method="bonferroni")    
PT
```
Continue the group comparison tests.
```
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
```
Step 9: Compute a two factorial ANOVA between two independent variables, genotype and group membership.

Since the samples sizes in the groups are uneven, you need to do an unbalanced two way ANOVA on a linear model. 

But when you run a linear model ```lm()```, there are some missing values in the COMT column. Remove those missing values.
```
clean_data$COMT2 <- factor(clean_data$COMT, exclude = "", levels = c("Met","Het","Val"))
```
Remove missing values for in the BDNF column too.
```
clean_data$BDNF2 <- factor(clean_data$BDNF, exclude = "", c("Met_carrier", "Val"))
```
Run a type 3 two factorial ANOVA using the ```Anova()``` function from the ```car``` package. (There are 3 types of two way ANOVAs, type 3 is the easiest). However, you need to change the default settings before starting.

First, check the contrasts. What are the contrasts? Contrasts show how one group is compared to another group. 

Each type of contrast will give you a different sum of squares. The default setting is ```contr.treatment```(for unordered data) and ```contr.poly```(for ordered data). When you open R and you check the contrasts, it'll show you treatment (```contr.treatment```) and poly (```contr.poly```). 
You need to set it to ```contr.sum``` and ```contr.poly```. (This has been done below for each test). According to the creator of the ```car``` package, you can't use the treatment setting because the tests wont give you sensible results. This information can be found in the help file of the ```car``` package.

Secondly, the data is already sorted into reference groups. The reference group is what the other groups are compared to. Use the following code to change the reference groups:
```
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("AMC","BSO-E2","BSO")) #to make AMC the reference group
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("BSO","BSO-E2","AMC")) #to make BSO the reference group
```
If you'd like to save the results from your ANOVA tables onto a work document, use the following code:
```
knitr::kable(nice(COMT_SPWM))%>%   
  write2word(file="RAVLTlearn.doc", quiet = TRUE)
```  
For more information, see this website, https://cran.r-project.org/web/packages/arsenal/vignettes/write2.html#introduction

In the code below, an ANOVA is being run on a linear model. Since there was a significant difference in the proportion of study participants with a history of cancer treatments, this variable has been added to the ANOVA as a predictor variable.
```
# Two factorial ANOVAs -----

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(SPWM_WME_T1 ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data), type = 3)
```

If results are significant, run post hoc test. Unfortunately, the ```TukeyHSD``` function can only be applied after the ```aov()``` function, it cannot be used after the ```Anova()``` function. Therefore, run a linear regression to see where the differences in the groups is.

Don't forget to change the contrasts back to treatment and poly! (the default)

```
options(contrasts = c("contr.treatment", "contr.poly"))
COMT_SPWM <- lm(SPWM_WME_T1 ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data)
summary(COMT_SPWM)
```

Het and AMCs are the reference groups here. To check reference groups, use contrasts function
```
contrasts(clean_data$group_membership)
clean_data$COMT2 <- factor(clean_data$COMT2, levels = c("Het","Met","Val"))
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("AMC","BSO-E2","BSO"))
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("BSO","BSO-E2","AMC"))
```
Continue conducting comparisons for all the other test scores and genotypes.
```
#I want to know if DigitsForward scores are affected by COMT genotype and group membership

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(DigitsForward ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data), type = 3)

#I want to know if DigitsBackward scores have an effect on COMT genotype and group membership

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(DigitsBackward ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data), type = 3)

#I want to know if DigitsTotal scores are affected by COMT genotype and group membership

options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(DigitsTotal ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data), type = 3)

options(contrasts = c("contr.treatment", "contr.poly"))
COMT_DigitsTotal <- lm(DigitsTotal ~ COMT2 * group_membership + Cancer_Treatments, data = clean_data)
summary(COMT_DigitsTotal)

#I want to know if RAVLT_learn scores are affected by BDNF genotype and group membership

#Remember, to change reference groups:
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

#I want to know if Logical Memory Immediate scores are affected by BDNF genotype and group membership
options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(LMA_Imm_Verbatim ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data), type = 3)

options(contrasts = c("contr.treatment", "contr.poly"))
BDNF_LMA_Imm_RAVLT<- lm(LMA_Imm_Verbatim ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data)
summary(BDNF_LMA_Imm_RAVLT)

#I want to know if Logical Memory Delayed scores are affected by BDNF genotype and group membership
options()$contrasts #to check the contrasts
options(contrasts = c("contr.sum","contr.poly"))
Anova(lm(LMA_Del_Verbatim ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data), type = 3)

options(contrasts = c("contr.treatment", "contr.poly"))
BDNF_LMA_Del_RAVLT<- lm(LMA_Del_Verbatim ~ BDNF_combined * group_membership + Cancer_Treatments, data = clean_data)
summary(BDNF_LMA_Del_RAVLT)

#To change the reference group
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("AMC","BSO-E2","BSO"))
clean_data$group_membership <- factor(clean_data$group_membership, levels = c("BSO","BSO-E2","AMC"))
```
Step 10: Make visualisations for the E1G levels among the groups. You can make a variety of graphs!

First, make a graph for the E1G levels (a measure of how much estrogen is in a woman's body)

Make a graph theme that can be applied to all your graphs. Use the ```theme()``` function to determine the theme, inluding the ```element_text()``` function to decide the font for the graph text.
```
E1GTheme <- theme(plot.title= element_text(family = "Palatino", face = "bold")) #determining font
```
Next, make a boxplot using the ```ggplot()``` function. Use the```aes()``` function to make 'aesthetic mappings', which configure the visual settings (aesthetics) of your graph. Use the ```geom_boxplot()``` function to make the boxplot, and finally print both the boxplot and the theme together. 
```
E1GPlot <- ggplot(clean_data, aes(X=group_membership, y = E1G, fill=group_membership)) + geom_boxplot(outlier.size = .8, lwd=1)
print(E1GPlot + E1GTheme + labs(title="17-β-Estradiol Levels in Different Women Groups", 
                                y = "Estrone-3-Glucoronide (ng/ml)", x="Group Membership",
                                fill = "Group Membership"))
```
You can also make a boxplot with jittered points. You can even change the colours if you want.
```
E1G_Boxplot <- ggboxplot(clean_data, x = "group_membership", y = "E1G",
                         color = "group_membership", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                         add = "jitter", shape = "group_membership", legend = "none")
E1GTheme <- theme(plot.title= element_text(family = "Helvetica", size = (13), hjust = 0.5))
```
If you want to compare Kruskal-Wallis p values, you can create a variable called ```my_comparisons``` and use the ``` list()``` function to denote which groups you are comparing.
```
my_comparisons <- list( c("BSO", "BSO-E2"), c("BSO-E2", "AMC"), c("BSO", "AMC") ) 
```
Finally, use the ```print()``` function to add all the graph elements together.
```
print(E1G_Boxplot
      + E1GTheme 
      + stat_compare_means(comparisons = my_comparisons) #if you want to include Kruskal-Wallis values
      + stat_compare_means(label.y = 165) #if you want to include Kruskal Wallis values
      + labs(title="17-β-Estradiol Levels in Different Women Groups", 
             y = "Estrone-3-Glucoronide (ng/ml)", x="Group Membership",
             fill = "Group Membership"))
```
You can also create a violin plot if you'd like with the ```ggplot()``` function.
```
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
```
Step 11: Make visualisations for your neuropsych test results. 

```
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
```
Step 12: Test whether the genotype frequencies are in accordance with Hardy-Weinberg equilibrium with the ```HWChisq()``` function. MM means Met/Met, VM means Val/Met, VV means Val/Val. The number next to each genotype denotes how many samples of that genotype we have.
```
#Hardy Weinberg Test for COMT AMCs
x <- c(MM=4,VM=11,VV=7) 
HW.test <- HWChisq(x, verbose=TRUE)

#Hardy Weinberg Test for COMT BSOs
x <- c(MM=6,VM=6,VV=5) 
HW.test <- HWChisq(x, verbose=TRUE)

#Hardy Weinberg Test for COMT BSO-E2s
x <- c(MM=4,VM=7,VV=4) 
HW.test <- HWChisq(x, verbose=TRUE)

#Hardy Weinberg Test for BDNF AMCs
x <- c(MM=1,VM=3,VV=9) 
HW.test <- HWChisq(x, verbose=TRUE)

#Hardy Weinberg Test for BDNF BSOs
x <- c(MM=2,VM=5,VV=15) 
HW.test <- HWChisq(x, verbose=TRUE)

#Hardy Weinberg Test for BDNF BSOE2s
x <- c(MM=0,VM=2,VV=11) 
HW.test <- HWChisq(x, verbose=TRUE)
```
For more information, use the ```help()``` function to locate the information page. 
```
help(HardyWeinberg)
```
Step 13: To cite the packages you used, use the ```citation()``` function.
```
citation("car")
citation("HardyWeinberg")
```
