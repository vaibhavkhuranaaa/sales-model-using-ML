rm(list=ls())

# Loading the packages
library(Hmisc)
library(dplyr)
library(haven)
library(tidyr)
library(ggplot2)
library(fastDummies)


# Load dataset
data.raw = read.csv("C:/Users/Vaibhav Khurana/Downloads/S10(1).csv")  
View(data.raw)
str(data.raw)
nrow(data.raw)
ncol(data.raw)
#arranging the variables and forming anew data set
data <- data.raw %>% dplyr::select(SeqNum,MOSTYP, MOSHOO, MGEMOM, MAANTH, MGEMLE, everything())
View(data)
str(data)

#mutating data before prepping for logistic regression analysis
#FOR L0 and L2

#L0 <- data %>% mutate(master= if_else(.$MOSTYP==2,1,0),
# phd   = if_else(.$MOSTYP==3,1,0)
#)

L02 <- dummy_cols(data, select_columns = c('MOSTYP', 'MOSHOO'),remove_selected_columns = TRUE)
View(L02)
str(L02)

#for L1

L1 <- L02 %>% mutate(across(MGEMLE, ~ case_when(
  (.x == 0 ~ 25),
  (.x == 1 ~ 35),
  (.x == 2 ~ 45),
  (.x == 3 ~ 55),
  (.x == 5 ~ 65),
  (.x == 6 ~ 75),
  (TRUE ~ -99)
)))

View(L1)
str(L1)

#mutating data on the basis of L3
#Alternate code
L3 <- L1 %>% mutate(across(MOPLHO:MAUT0, ~ case_when(
  (.x == 0 ~ 0),
  (.x == 1 ~ 5.5),
  (.x == 2 ~ 17),
  (.x == 3 ~ 30),
  (.x == 4 ~ 43),
  (.x == 5 ~ 56),
  (.x == 6 ~ 69),
  (.x == 7 ~ 82),
  (.x == 8 ~ 94),
  (.x == 9 ~ 100),
  (TRUE ~ -99)
)))

#ALT for L3
'L3 <- L2 %>% mutate(across(PPERSA:PWAPAR,~if_else(.x==0,   0,
                                                   if_else(.x==1,  5.5,
                                                           if_else(.x==2,   17,
                                                                   if_else(.x==3,  30,
                                                                           if_else(.x==4,  43,
                                                                                   if_else(.x==5, 56,
                                                                                           if_else(.x==6, 69,
                                                                                                   if_else(.x==7, 82,
                                                                                                           if_else(.x==8, 94,
                                                                                                                   if_else(.x==9, 100, -99))))))))))))'

View(L3)
str(L3)

#For L4, naming it as before split

beforesplit <- L3 %>% mutate(across(PPERSA:AWAPAR, ~ case_when(
  (.x == 0 ~ 0),
  (.x == 1 ~ 25),
  (.x == 2 ~ 75),
  (.x == 3 ~ 150),
  (.x == 5 ~ 350),
  (.x == 5 ~ 750),
  (.x == 6 ~ 3000),
  (.x == 7 ~ 7500),
  (.x == 8 ~ 15000),
  (.x == 9 ~ 30000),
  (TRUE ~ -99)
)))

View(beforesplit)
str(beforesplit)


#Splitting into training and testing data
set.seed(1)
train <- beforesplit %>% dplyr::sample_frac(0.70)
test  <- dplyr::anti_join(beforesplit, train, by = "SeqNum" )

#remove seqnum and train logistic model

train<-subset(train,select=-c(SeqNum))
test <-subset(test,select=-c(SeqNum))

fullModel = glm(Resp ~ ., family = 'binomial', data = train) # model with all variables
nullModel = glm(Resp ~ 1, family = 'binomial', data = train) # model with intercept only

#we need this library for stepAIC
library(MASS)
interim<-summary(stepAIC(nullModel, # start with a model containing no variables
                         direction = 'forward', # run forward selection
                         scope = list(upper = fullModel, # the maximum to consider is a model with all variables
                                      lower = nullModel), # the minimum to consider is a model with no variables
                         trace = 0)) # do not show the step-by-step process of model selection

coef<-data.frame(interim[['coefficients']])
final<-coef[coef$Pr...z..<0.05,]
print(final)

#Retyping list of variables into final model build
varnames<-rownames(final)
varnames<-varnames[2:length(varnames)]
finalmodel<-glm(Resp ~ PPERSA+MAANTH+MGODPR+MGEMOM+MOSTYP_34+MOSTYP_21+MAUT2, family = 'binomial', data = train)


#Evaluating performance on test data
test$pred<-predict(finalmodel,newdata=test,type="response")
test<-test[order(-test$pred),]
test$one<-1
test$cumprospects<-cumsum(test$one)
test$cumresp    <-cumsum(test$Resp)

Perf<-subset(test,select=c(pred,cumprospects,cumresp))
Perf$PctProspect<-Perf$cumprospects/nrow(Perf)
Perf$PctResp    <-Perf$cumresp/max(Perf$cumresp)

cutpoint<-subset(Perf,PctProspect>0.745 & PctProspect<0.755)
cutpoint


ForMemo <-lm(Resp ~ PPERSA+MAANTH+MGODPR+MGEMOM+MOSTYP_34+MOSTYP_21+MAUT2, data = train)
summary(ForMemo)

write.csv(Perf, "C:/Users/Vaibhav Khurana/Downloads/perf.csv")

Perf$Lift<-Perf$PctResp-Perf$PctProspect
MaxLift<-Perf[Perf$Lift==max(Perf$Lift),]
MaxLift
