#Libraries
library(ggplot2)
library(dplyr)
library(AICcmodavg)
library(MASS)
library(lme4)
library(gamlss)
library(reshape2)

#Datasets
florida <- read.table("Florida.csv", sep=",", header = T)
states <- read.table("FloridaWisconsinVirginia.csv", sep=",", header = T)

#Analysis
  #Calculation of blood Hg from feather Hg - formula from Eagles-Smith et al. (2008)
  states$Blood.Hg <- exp(0.673*log(states$Hg)-1.673)

  #Histograns
  ggplot(data=states, aes(x=Hg))+
    geom_histogram(fill="#222222", size=1, bins=20)+
    theme_bw()+
    xlab("Concentration of THg (µg/g)")+
    ylab("Number of Samples")

  ggplot(data=florida, aes(x=Cort))+
    geom_histogram(fill="#222222", size=1, bins=20)+
    theme_bw()+
    xlab("Concentration of Corticosterone (pg/mg)")+
    ylab("Number of Samples")

  #Hg
    #Model selection
    Hg.all<-lm(Hg~State+Sex, data=states) #Created for stepAIC
    Hg.step<-stepAIC(Hg.all, trace = F) #Model selected by stepAIC was Null
    Hg.state <- lm(Hg~State, data=states)
    Hg.sex <- lm(Hg~Sex, data=states)

    cand <- list("State"=Hg.state,
                 "Sex"=Hg.sex,
                 "Null"=Hg.step)
    
    aictab(cand.set = cand) #Null model selected
    
    #Plot of non-selected model using breeding location ("State" variable) as predictor. 
    ##Plot intended for visualization of non-correlation only. 
    ggplot(data=states, aes(x=State, y=Hg))+
      geom_boxplot(size = 1)+
      theme_bw()+
      xlab("Breeding Location")+
      ylab("Concentration of THg (µg/g)")+
      scale_y_continuous(limits = c(1, 9), breaks = seq(1, 9, by = 2))
    
  #Cort
    #Model selection
    cort.all<-lm(Cort~State+Hg+Sex, data=states) #Created for stepAIC
    cort.step<-stepAIC(cort.all, trace = F) #Model selected by stepAIC was State + Hg
    cort.state <- lm(Cort~State, data=states)
    cort.hg <- lm(Cort~Hg, data=states)
    cort.null <- lm(Cort~1, data=states)
    cort.sex <- lm(Cort~Sex, data=states)
    
    cand <- list("State"=cort.state,
                 "Hg"=cort.hg,
                 "Sex"=cort.sex,
                 "State+Hg"=cort.step,
                 "Null"=cort.null)
    
    aictab(cand.set = cand) #Model using breeding location as sole predictor selected
    
    #Summary of selected model
    summary(cort.state)
    states %>% group_by(State) %>% summarise(mean=mean(Cort))
    
    #Plot of selected model. Log() function used to facilitate visualization
    ggplot(data=states, aes(y=log(Cort), x=State))+
      geom_boxplot(fill="white", size=1)+
      theme_bw()+
      xlab("Breeding Location")+
      ylab("Concentration of Corticosterone (ln pg/mg)")
    
  #Mass
    #Model selection
    mass.null<-lm(Mass~1, data=florida)
    mass.all<-lm(Mass~log(Hg)+Age+Sex, data=florida) #Created for stepAIC
    mass.step<-stepAIC(mass.all, trace = F) #Model selected by stepAIC was log(Hg)+Sex
    mass.hg<-lm(Mass~log(Hg), data=florida) #log(Hg) used to facilitate visualization
    mass.age<-lm(Mass~Age, data=florida)
    mass.sex<-lm(Mass~Sex, data=florida)
    
    cand <- list("Null"=mass.null,
                 "ln(Hg)"=mass.hg,
                 "Age"=mass.age,
                 "Sex"=mass.sex,
                 "ln(Hg)+Sex"=mass.step)
    
    aictab(cand.set = cand) #Model using log(Hg) as sole predictor selected
    
    #Plot of selected model. Log() function used to facilitate visualization
    ggplot(data=florida, aes(y=Mass, x=log(Hg)))+
      geom_point(size=1)+
      geom_line(aes(y=predict(mass.hg, type="response")), size=1, color="black")+
      theme_bw()+
      xlab("Concentration of THg (ln µg/g)")+
      ylab("Mass (g)")

  
  #Fat Score
    #Model selection
    florida$Fat.Score<-factor(florida$Fat.Score, levels=c(4,3,2,1)) #polr requires response to be a factor variable
    fat.all<-polr(Fat.Score~Hg+Sex+Age, data=florida, Hess = TRUE) #Created for stepAIC
    fat.step<-stepAIC(fat.all, trace=F) #Model selected by stepAIC was Hg + Sex + Age. Same model as fat.all
    fat.null<-polr(Fat.Score~1, data=florida, Hess = TRUE)
    fat.hg<-polr(Fat.Score~Hg, data=florida, Hess = TRUE)
    fat.sex<-polr(Fat.Score~Sex, data=florida, Hess = TRUE)
    fat.age<-polr(Fat.Score~Age, data=florida, Hess = TRUE)
    fat.hgage<-polr(Fat.Score~Age+Hg, data=florida, Hess = TRUE) #Removed from article to decrease table size. Not selected by aictab
    fat.hgsex<-polr(Fat.Score~Sex+Hg, data=florida, Hess = TRUE) #Removed from article to decrease table size. Not selected by aictab
    fat.sexage<-polr(Fat.Score~Sex+Age, data=florida, Hess = TRUE) #Removed from article to decrease table size. Not selected by aictab
        
    cand <- list("Null"=fat.null,
                 "Hg"=fat.hg,
                 "Age"=fat.age,
                 "Sex"=fat.sex,
                 "Hg+Age"=fat.hgage,
                 "Hg+Sex"=fat.hgsex,
                 "Sex+Age"=fat.sexage,
                 "Hg+Sex+Age"=fat.all)
    
    aictab(cand.set = cand)
    
    #Summary of selected model
    summary(fat.hg_sex_age)
    
    #Creation of data.frame from selected model for visualization. 
    ##Code adapted from https://stats.oarc.ucla.edu/r/dae/ordinal-logistic-regression/
    newdat <- data.frame(
      Age = rep(1:7, len=630),
      Hg = rep(1:9, len=630),
      Sex = rep(c("F","M"), each=315))
    newdat <- cbind(newdat, predict(fat.hg_sex_age, newdat, type="probs"))
    
    lnewdat <- melt(newdat, id.vars = c("Age","Hg", "Sex"),
                    variable.name = "Fat.Score", value.name="Probability")
    
    lnewdat <- subset(lnewdat, Age %in% 2:7)
    
    #Plot of selected model.
    ggplot(data=lnewdat, aes(x=Hg, y=Probability, fill=Fat.Score))+
      geom_bar(position=position_fill(), stat="identity")+
      facet_grid(Age~Sex, labeller = label_both)+
      xlab("THg (ug/g)")+
      ylab("Probability")+
      labs(fill="Fat Score")+
      scale_fill_manual(values=c("#55aaff","#6fca6f","#ffb86c","#ff5555"))+
      theme_bw()
    