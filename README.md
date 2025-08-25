# Insecticide_Resistance_Kwale_County
The code below describes the association btn phenotypic and genotypic insecticide resistance monitoring in Anopheles mosquitoes in Kwale County
###setting work directory
insecticde<-read.csv("C:/DATA/Msc Research Docs/Thesis/Thesis documents/Reviewed Thesis/Marta/Insecticide_data_Kwale.csv")
getwd("C:/DATA/Msc Research Docs/Thesis/Thesis documents/Reviewed Thesis/Marta")
names(insecticde)
View(insecticde)
#### load library for data cleaning
library(psych)
library(DescTools)
library(Hmisc)
library(compareGroups)
library(janitor)
library(sjPlot)
describe(insecticde)
tabyl(insecticde, Insecticide, Synergist.Assay)

##Remove empty rows/spaces
insecticde<-insecticde[grepl("without synergist",insecticde$Synergist.Assay),]
insecticde<-insecticde[grepl("Deltamethrin|Permethrin",insecticde$Insecticide),]

###Categorizing variables
insecticde$CYP6Pb.detection_Second.PCR.Results<-factor(insecticde$CYP6Pb.detection_Second.PCR.Results)
insecticde$status<-factor(insecticde$status)
insecticde$Lab_complex<-factor(insecticde$Lab_complex)
insecticde$sibling.Species_PCR<-factor(insecticde$sibling.Species_PCR)
insecticde$L119F.GSTe2<-factor(insecticde$L119F.GSTe2)
insecticde$KDR.Detection<-factor(insecticde$KDR.Detection)
insecticde$CYP6Pa_Second.PCR.Results<-factor(insecticde$CYP6Pa_Second.PCR.Results)
insecticde$CYP6Pb.detection_Second.PCR.Results<-factor(insecticde$CYP6Pb.detection_Second.PCR.Results)
insecticde$X6.5kb.SV<-factor(insecticde$X6.5kb.SV)

##Checking for distribution of Mosquitoes
table(insecticde$Lab_complex)

### Filtering funestus by complexes
af<-insecticde[grepl("AN. FUNESTUS",insecticde$Lab_complex),]
af
#### Filtering by resistance status CYP
b<-af[grepl("Resistant|Susceptible",af$CYP6Pb.detection_Second.PCR.Results),]
a<-af[grepl("Resistant|Susceptible|Heterozygote",af$CYP6Pa_Second.PCR.Results),]
sv<-af[grepl("Resistant|Susceptible|Heterozygote",af$X6.5kb.SV),]
gst<-af[grepl("Resistant|Susceptible|Heterozygote",af$L119F.GSTe2),]

#### Factoring variables 
b$CYP6Pb.detection_Second.PCR.Results<-factor(b$CYP6Pb.detection_Second.PCR.Results)
b$Lab_complex<-factor(b$Lab_complex)
b$sibling.Species_PCR<-factor(b$sibling.Species_PCR)
sv$sibling.Species_PCR<-factor(sv$sibling.Species_PCR)
sv$X6.5kb.SV<-factor(sv$X6.5kb.SV)
gst$sibling.Species_PCR<-factor(gst$sibling.Species_PCR)
gst$L119F.GSTe2<-factor(gst$L119F.GSTe2)

tabyl(factor(gst$L119F.GSTe2))
table(gst$sibling.Species_PCR, gst$L119F.GSTe2)

####Perfom chi-square for all genes associated with An. funestus

cypb<-compareGroups(formula = status~CYP6Pb.detection_Second.PCR.Results,
                    data = b,
                    max.xlev = 20,
                    max.ylev = 30)

cypb_table<-createTable(x = cypb,show.n = TRUE,
                        show.p.ratio = TRUE, show.ratio = TRUE)
summary(cypb_table)
export2csv(x = cypb_table,file = "CYP6P9b.csv")

cypa<-compareGroups(formula = status~CYP6Pa_Second.PCR.Results,
                    data = a,
                    max.xlev = 20,
                    max.ylev = 30)

cypa_table<-createTable(x = cypa,show.n = TRUE,
                        show.p.ratio = TRUE, show.ratio = TRUE)
summary(cypa_table)
export2csv(x = cypa_table,file = "CYP6P9a.csv")


svar<-compareGroups(formula = status~X6.5kb.SV,
                    data = sv,
                    max.xlev = 20,
                    max.ylev = 30)

sv_table<-createTable(x = svar,show.n = TRUE,
                      show.p.ratio = TRUE, show.ratio = TRUE)
summary(sv_table)
export2csv(x = sv_table,file = "structural variation.csv")

gste<-compareGroups(formula = status~L119F.GSTe2,
                    data = gst,
                    max.xlev = 20,
                    max.ylev = 30)

gste_table<-createTable(x = gste,show.n = TRUE,
                        show.p.ratio = TRUE, show.ratio = TRUE)
summary(gste_table)
export2csv(x = gste_table,file = "GSTE.csv")


##### Chi-square for KDR genes on Anopeheles gambiae 
ag<-insecticde[grepl("AN. GAMBIAE",insecticde$Lab_complex),]

kdr<-ag[grepl("Heterozygote-WT/KDr East|Heterozygote-WT/KDr West|Resistant-KDr East|Resistant KDr West|Susceptible-WT",ag$KDR.Detection),]

kdrCG<-compareGroups(formula = status~KDR.Detection,
                     data = kdr,
                     max.xlev = 20,
                     max.ylev = 30)

kdrCG_table<-createTable(x = kdrCG,show.n = TRUE,
                         show.p.ratio = TRUE, show.ratio = TRUE)
summary(kdrCG_table)
export2csv(x = kdrCG_table,file = "kdr.csv")
