# Modelise dose-response curve and extract EC50
# Extract plate data, modelize and calculate EC50

#you need to export data from Omega in table format with absorbance and stimulant concentration
#in the table excel file, replace absorbance value you do not use (in grey) by NA 

setwd("~/MOrDOr_project/Phenotype/")
require(readxl)
require(drc)

### Import dataset ----

abs_data <- read_excel("test.xlsx")

colnames(abs_data) <- abs_data[12,] #set column names
abs_data <- abs_data[c(13:108),] #remove the unecessary rows
abs_data$`Standard Concentrations [M]`[abs_data$Content == "Blank B"] <- 1e-13 #set blank concentration to 1e-13
abs_data$Group[abs_data$Content == "Blank B"] <- "Blank" #set blank group
abs_data$absorbance <- as.numeric(abs_data$`Difference based on Raw Data and Raw Data (calculated)`) #set absorbance as numeric variables
abs_data$concentration <- as.numeric(abs_data$`Standard Concentrations [M]`) #set concentration as numeric variables
select <- complete.cases(abs_data) #select complete rows
abs_data <- abs_data[select,] #remove empty well


### Subset by molecules ----

# set all tested moelcules in that plate /!\ /!\ /!\
molecules=c("2PEITC","racGR24")

abs_data$molecule <- abs_data$Group
if (length(molecules) == 2) {
  abs_data$molecule[abs_data$Group == "A"] <- molecules[1]
  abs_data$molecule[abs_data$Group == "B"] <- molecules[2]
} else {
  abs_data$molecule[abs_data$Group == "A"] <- molecules[1]
  abs_data$molecule[abs_data$Group == "B"] <- molecules[2]
  abs_data$molecule[abs_data$Group == "C"] <- molecules[3]
  abs_data$molecule[abs_data$Group == "D"] <- molecules[4]
}

#choose one molecule to  modelize /!\ /!\ /!\
molec = "2PEITC"

mod_data <- subset(abs_data, abs_data$molecule %in% c(molec,"Blank"))

### Calculate relative absorbance ----

Tmax <- mean(mod_data$absorbance[mod_data$concentration>= 1e-7],
             na.rm = T)
Tmin <- mean(mod_data$absorbance[mod_data$concentration <= 1e-13],
             na.rm = T)

for (i in 1:nrow(mod_data)) {
  mod_data[i,"rel_abs"] <- ((mod_data[i,"absorbance"]-Tmin)/(Tmax-Tmin))
}

### Modelize dose-response curve and extract EC50 ----

#modelize (if Convergence failed, see section below)
mod1<- drm(mod_data$rel_abs~mod_data$concentration,
           #plot_data$sample,
           data=mod_data,
           fct = LL.4())
#plot model
plot(mod1,  xlab=" GS moles", ylab="ABS percent of germination",
     ylim=c(-0.1,2), # ylim depend of abs max
     xlim=c(0, 0.2),type= "average", # type can be "bars" or "none"
     main ="GS stimulation",
     cex=1.2, cex.axis=1.2, lwd=2, legendPos = c(1e-12,1.5))# EQUALS NUMBER OF SAMPLE 
#EXTRACT EC50
ED(mod1, 50)

### If convergence failed ----

#add positive control 
mod_data <- rbind(mod_data,c("NA","NA","NA","NA","NA",
                             max(mod_data$absorbance, na.rm = T),
                             1e-5,"X",
                             max(mod_data$rel_abs, na.rm = T)))
#add negative control
mod_data <- rbind(mod_data,c("NA","NA","NA","NA","NA",
                             0,0,"X",0))
