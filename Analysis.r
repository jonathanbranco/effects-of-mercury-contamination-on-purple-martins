#Libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(AICcmodavg)
library(MASS)
library(lme4)
library(gamlss)
library(reshape2)

#Datasets ----
florida <- read.table("Florida.csv", sep=",", header = T)
states <- read.table("FloridaWisconsinVirginia.csv", sep=",", header = T)

#Analysis ----
  #Calculation of blood Hg from feather Hg - formula from Eagles-Smith et al. (2008)
  states$Blood.Hg <- exp(0.673*log(states$Hg)-1.673)

  #Histograns
  ggplot(data=states, aes(x=Hg))+
    geom_histogram(fill="#222222", size=1, bins=20)+
    theme_bw()+
    xlab("Concentration of THg (ug/g)")+
    ylab("Number of Samples") -> hg_hist

  ggplot(data=states, aes(x=Cort))+
    geom_histogram(fill="#222222", size=1, bins=20)+
    theme_bw()+
    xlab("Concentration of Corticosterone (pg/mg)")+
    ylab("Number of Samples") -> cort_hist
  
    #Joins histogram plot
    hist_plot <- grid.arrange(hg_hist, cort_hist, ncol=1)
    
  #Statistical Models ----
    #Hg ----
      #Model selection
      Hg_state_sex <- lm(Hg~State+Sex, data=states)
      Hg_null <- lm(Hg~1, data=states)
      Hg_state <- lm(Hg~State, data=states)
      Hg_sex <- lm(Hg~Sex, data=states)
  
      cand_Hg <- list("Breeding Location"=Hg_state,
                      "Sex"= Hg_sex,
                      "Breeding Location + Sex" = Hg_state_sex,
                      "Null"=Hg_null)
      
      aictab_hg <- aictab(cand.set = cand_Hg) #Null model selected
      
      #Plot of non-selected model using breeding location ("State" variable) as predictor. 
      #Plot intended for visualization of non-correlation only. 
      ggplot(data=states, aes(x=State, y=Hg))+
        geom_boxplot(size = 1)+
        theme_bw()+
        xlab("Breeding Location")+
        ylab("Concentration of THg (ug/g)")+
        scale_y_continuous(limits = c(1, 9), breaks = seq(1, 9, by = 2)) -> hg_plot
      
    #Cort ----
      #Model selection
      cort_all <- lm(Cort~State+Hg+Sex, data=states)
      cort_state_sex <- lm(Cort~State+Sex, data=states) 
      cort_state <- lm(Cort~State, data=states)
      cort_hg <- lm(Cort~Hg, data=states)
      cort_null <- lm(Cort~1, data=states)
      cort_sex <- lm(Cort~Sex, data=states)
      
      cand_cort <- list("Breeding Location" = cort_state,
                        "THg" = cort_hg,
                        "Sex" = cort_sex,
                        "Breeding Location + THg" = cort_state_sex,
                        "Breeding Location + THg + Sex" = cort_all,
                        "Null" = cort_null)
      
      aictab_cort <- aictab(cand.set = cand_cort) #Model using breeding location as sole predictor selected
      
      #Plot of selected model. Log10() function used to facilitate visualization
      ggplot(data=states, aes(y=log10(Cort), x=State))+
        geom_boxplot(fill="white", size=1)+
        theme_bw()+
        xlab("Breeding Location")+
        ylab("Concentration of Corticosterone (log10 pg/mg)") -> cort_plot
      
        #Joins Cort and Hg boxplot
        boxplots <- grid.arrange(hg_plot, cort_plot, ncol=1)
      
    #Mass ----
      #Model selection
      mass_null <- lm(Mass~1, data=florida)
      mass_all <- lm(Mass~log10(Hg)+Age+Sex, data=florida) 
      mass_hg_age <- lm(Mass~log10(Hg)+Age, data=florida) 
      mass_hg_sex <- lm(Mass~log10(Hg)+Sex, data=florida)
      mass_hg <- lm(Mass~log10(Hg), data=florida) #Selected model. Log10(THg) used to facilitate visualization
      mass_age_sex <- lm(Mass~Age+Sex, data=florida)
      mass_age <- lm(Mass~Age, data=florida)
      mass_sex <- lm(Mass~Sex, data=florida)
      
      cand_mass <- list("Null"=mass_null,
                        "log10(THg) + Age + Sex"=mass_all,
                        "log10(THg) + Age"=mass_hg_age,
                        "log10(THg) + Sex"=mass_hg_sex,
                        "log10(THg)"=mass_hg,
                        "Age + Sex"=mass_age_sex,
                        "Age"=mass_age,
                        "Sex"=mass_sex)
      
      aictab_mass <- aictab(cand.set = cand_mass) #Model using log10(THg) as sole predictor selected
      
      #Plot of selected model. Log10() function used to facilitate visualization
      ggplot(data=florida, aes(y=Mass, x=log10(Hg)))+
        geom_point(size=1)+
        geom_line(aes(y=predict(mass_hg, type="response")), size=1, color="black")+
        theme_bw()+
        xlab("Concentration of THg (log10 ug/g)")+
        ylab("Mass (g)") -> mass_plot
  
    #Fat Score ----
      #Model selection
      florida$Fat.Score<-factor(florida$Fat.Score, levels=c(4,3,2,1)) #polr requires response to be a factor variable
      
      fat_all<-polr(Fat.Score~Hg+Sex+Age, data=florida, Hess = TRUE)
      fat_null<-polr(Fat.Score~1, data=florida, Hess = TRUE)
      fat_hg<-polr(Fat.Score~Hg, data=florida, Hess = TRUE)
      fat_sex<-polr(Fat.Score~Sex, data=florida, Hess = TRUE)
      fat_age<-polr(Fat.Score~Age, data=florida, Hess = TRUE)
      fat_hg_age<-polr(Fat.Score~Age+Hg, data=florida, Hess = TRUE) #Removed from article to decrease table size. Not selected by aictab
      fat_hg_sex<-polr(Fat.Score~Sex+Hg, data=florida, Hess = TRUE) #Removed from article to decrease table size. Not selected by aictab
      fat_sex_age<-polr(Fat.Score~Sex+Age, data=florida, Hess = TRUE) #Removed from article to decrease table size. Not selected by aictab
          
      cand_fat_score <- list("Null"=fat_null,
                             "THg"=fat_hg,
                             "Age"=fat_age,
                             "Sex"=fat_sex,
                             "THg + Age"=fat_hg_age,
                             "THg + Sex"=fat_hg_sex,
                             "Sex + Age"=fat_sex_age,
                             "THg + Sex + Age"=fat_all)
      
      aictab_fat_score <- aictab(cand.set = cand_fat_score)

      #Creation of data.frame from selected model for visualization. 
      ##Code adapted from https://stats.oarc.ucla.edu/r/dae/ordinal-logistic-regression/
      newdat <- data.frame(Age = rep(1:7, len=630),
                           Hg = rep(1:9, len=630),
                           Sex = rep(c("F","M"), each=315))
      newdat <- cbind(newdat, predict(fat_all, newdat, type="probs"))
      
      lnewdat <- melt(newdat, id.vars = c("Age","Hg", "Sex"),
                      variable.name = "Fat.Score", value.name="Probability")
      
      lnewdat <- subset(lnewdat, Age %in% 2:7)
      
      #Plot of selected model.
      ggplot(data=lnewdat, aes(x=Hg, y=Probability, fill=Fat.Score))+
        geom_bar(position=position_fill(), stat="identity")+
        facet_grid(Age~Sex, labeller = label_both)+
        xlab("Concentration of THg (ug/g)")+
        ylab("Probability")+
        labs(fill="Fat Score")+
        scale_fill_manual(values=c("#55aaff","#6fca6f","#ffb86c","#ff5555"))+
        theme_bw() -> fat_score_plot
      
#Exporting results ----
  #Prepares and exports aictabs
  for(tab in c("aictab_hg","aictab_cort","aictab_mass","aictab_fat_score")){
    assign(tab,
           as_tibble(get(tab))) #Turns aictabs into tibbles to facilitate changes
    assign(tab,
           cbind(get(tab)[1:2], #Gets columns 1 and 2 as normal
                 round(as_tibble(get(tab))[3:8],2))) #Rounds columns 3 to 8 up to 2 decimal points
    assign(tab,
           rename(get(tab), 'Model predictor' = Modnames)) #Renames Modnames column
    write.csv(get(tab), paste("aictabs/",tab,".csv",sep=""), row.names=F) #Exports aictabs
  }
    
  #Exports plots
  ggsave("graphs/hist_plot.png", plot=hist_plot, width = 6, height=8, dpi="print")
  ggsave("graphs/boxplots.png", plot=boxplots, width = 6, height=8, dpi="print")
  ggsave("graphs/mass_plot.png", plot=mass_plot, width = 6, height=5, dpi="print")
  ggsave("graphs/fat_score_plot.png", plot=fat_score_plot, width = 6, height=9, dpi="print")
    
