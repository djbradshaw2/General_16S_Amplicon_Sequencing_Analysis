#Load libraries#
library(phyloseq)
library(indicspecies)
library(vegan)
library(ggplot2)
library(plyr)
library(tidyverse)
library(FSA)
library(reshape)
library(ggplot2)
library(Hmisc)
library(gdata)
library(tidyverse)

###LOAD UP SUMMARYSE FUNCTION###
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

####UPLOADING AND GENERAL EDITS####
#Import input files from QIIME2 if snakefile worked to phyloseq.biom#
BIOM <- import_biom(file.choose()) #phyloseq.biom
TREE =  read_tree(file.choose()) #file = tree -> tree.nwk
META <- import_qiime_sample_data(file.choose()) #map.txt

#Merge data#
LWSSdata <- merge_phyloseq (BIOM,TREE,META)

#Alternative Uploading if snakefile did not work to phyloseq.biom#
asv_table = read.delim(file.choose(), row.names=1, header = T) #file = "feature-table.tsv    #open file in Excel/TextEdit and remove 1st row in file before import
taxa_table = read.delim(file.choose(), row.names = 1) #file= "taxonomy -> taxonomy2.tsv"   #need to delete Confidence column and separate Taxonomy column using "Text to Columns" (Data tab) in Excel by "semicolon" and fill column headers with taxonomic ranks
taxa_table = as.matrix(taxa_table)
DottedMETA = read.delim(file.choose(), row.names=1) #mapping file
dotted_meta = sample_data(DottedMETA) 
ASV = otu_table(asv_table, taxa_are_rows = TRUE)
TAX = tax_table(taxa_table)
TREE =  read_tree(file.choose()) #file = tree -> tree.nwk

LWSSdata <- merge_phyloseq(ASV, TAX, dotted_meta, TREE)

#Change taxa names if needed#
colnames(tax_table(LWSSdata))
colnames(tax_table(LWSSdata))=c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
colnames(tax_table(LWSSdata))

#Check number of taxa, remove taxa that occur less than 11 times, recheck taxa#
ntaxa(LWSSdata) #110575#
FLWSSDataA = filter_taxa(LWSSDataA, function(x) sum(x) >10, TRUE)
ntaxa(FLWSSDataA) #16027#

####GENERAL SEQUENCE STATISTICS####
#Make table of trimmed sequences per ASV#
LWSSDataA_seqs_per_ASV <- as.data.frame(taxa_sums(LWSSDataA))
colnames(LWSSDataA_seqs_per_ASV) <- c("Sequences")
sum(LWSSDataA_seqs_per_ASV$Sequences) #1857744#
LWSSDataA_seqs_per_ASV = cbind(as(LWSSDataA_seqs_per_ASV, "data.frame"), as(tax_table(LWSSDataA)[rownames(LWSSDataA_seqs_per_ASV), ], "matrix"))
write.table(data.frame("OTUID" =rownames(LWSSDataA_seqs_per_ASV), LWSSDataA_seqs_per_ASV) , "Full LWSS Sequences Per Taxa.txt", sep="\t", row.names=FALSE)

#Make table of trimmed sequences per ASV#
FLWSSDataA_seqs_per_ASV <- as.data.frame(taxa_sums(FLWSSDataA))
colnames(FLWSSDataA_seqs_per_ASV) <- c("Sequences")
sum(FLWSSDataA_seqs_per_ASV$Sequences) #1598653#
FLWSSDataA_seqs_per_ASV = cbind(as(FLWSSDataA_seqs_per_ASV, "data.frame"), as(tax_table(FLWSSDataA)[rownames(FLWSSDataA_seqs_per_ASV), ], "matrix"))
write.table(data.frame("OTUID" =rownames(FLWSSDataA_seqs_per_ASV), FLWSSDataA_seqs_per_ASV) , "Trimmed LWSS Sequences Per Taxa.txt", sep="\t", row.names=FALSE)

#Make table of full and trimmed sequences per sample#
LWSSDataA_seqs_per_sample <- as.data.frame(sample_sums(LWSSDataA))
colnames(LWSSDataA_seqs_per_sample) <- c("Full_Sequences")
sum(LWSSDataA_seqs_per_sample$Full_Sequences) #1857744#
FLWSSDataA_seqs_per_sample <- as.data.frame(sample_sums(FLWSSDataA))
colnames(FLWSSDataA_seqs_per_sample) <- c("Trimmed_Sequences")
sum(FLWSSDataA_seqs_per_sample$Trimmed_Sequences) #1598653#
ComLWSSDataA_seqs_per_sample = cbind(as(LWSSDataA_seqs_per_sample, "data.frame"), as(FLWSSDataA_seqs_per_sample, "data.frame"))
write.table(data.frame("OTUID" =rownames(ComLWSSDataA_seqs_per_sample), ComLWSSDataA_seqs_per_sample) , "LWSS Sequences Per Sample.txt", sep="\t", row.names=FALSE)

####EXPORT DATA FOR USE IN PRIMER7####
# Extract abundance matrix from the phyloseq object
LWSS_Bio = as(otu_table(FLWSSDataA), "matrix")
# Coerce to data.frame
LWSS_Bio = as.data.frame(LWSS_Bio)
#Add in taxanomy#
LWSS_Bio_Tax = cbind(as(LWSS_Bio, "data.frame"), as(tax_table(FLWSSDataA)[rownames(LWSS_Bio), ], "matrix"))
#Make a csv file
write.csv(LWSS_Bio_Tax, file = 'LWSS_Bio_Tax2.csv')

####CREATE FIGURES FROM NWS WEATHER DATA####
###MAKE FIGURE S1A###
#Upload table containing the year or years in the first column, the monthly means from Jan to Dec in columns 2-13, and the last column containing the sensor location#
NWS_Temp <- read.csv(file.choose())

#Summarize by melting the months into one column based upon year and location#
NWS_Temp_melted <- melt(NWS_Temp, id.vars=c("Year","Location"))

#Write a csv to manipulate the data by adding a column indicating what sampling season the value was in or if it was it was not part of desired range, label Not Tested, then reupload#
write.csv(NWS_Temp_melted, "NWS Temp Melted.csv")
Adjusted_NWS_Temp <- read.csv(file.choose())

#Summarize data by Sampling Season#
NWS_Temp3 <- summarySE(data=Adjusted_NWS_Temp, measurevar="value", groupvars=c("Sampling.Season"), na.rm = FALSE, conf.interval = .95, .drop = TRUE)

#Summarize by Month#
NWS_Temp4 <- summarySE(data=Adjusted_NWS_Temp, measurevar="value", groupvars=c("variable"), na.rm = FALSE, conf.interval = .95, .drop = TRUE)

#Copy and paste two tables above to Excel and manipulate data to summarize by sampling season and the appropriate year, then reupload, and change order for graphing#
Temp_SS <- read.csv(file.choose())
Temp_SS$Sampling.Period <- factor(Temp_SS$Sampling.Period, levels = c("Aug-Sep", "Mar-Apr", "Oct-Nov", "Apr"))
Temp_SS <- Temp_SS[order(Temp_SS$Sampling.Period),]

#Make the graph#
ggplot(Temp_SS, aes(x=Sampling.Period, y=value, fill=factor(Year, levels =c("2016", "2017a", "2017b", "2018", "1990-2018"))))+
  geom_bar(position=position_dodge(), stat="identity", colour="black", size=0.3)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.3, width=0.2, position=position_dodge(0.9))+
  xlab("Sampling Period")+
  ylab("Temperature (°C)")+
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ggtitle("IRL Wide NWS Average Temperature (°C) by Sampling Period") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim=c(18,30))+
  scale_fill_manual(name="Year(s)", values=(c("chocolate", "seagreen", "chocolate1", "seagreen1", "grey20")), breaks=c("2016", "2017a", "2017b", "2018", "1990-2018"))

###MAKE FIGURE S1B###
#Upload table containing the year or years in the first column, the monthly means from Jan to Dec in columns 2-13, and the last column containing the sensor location#
NWS_Rain <- read.csv(file.choose())

#Summarize by melting the months into one column based upon year and location#
NWS_Rain_melted <- melt(NWS_Rain, id.vars=c("Year","Location"))

#Write a csv to manipulate the data by adding a column indicating what sampling season the value was in or if it was it was not part of desired range, label Not Tested, then reupload#
write.csv(NWS_Rain_melted, "NWS Rain Melted.csv")
Adjusted_NWS_Rain <- read.csv(file.choose())

#Summarize data by Sampling Season#
NWS_Rain3 <- summarySE(data=Adjusted_NWS_Rain, measurevar="value", groupvars=c("Sampling.Season"), na.rm = FALSE, conf.interval = .95, .drop = TRUE)

#Summarize data by Month#
NWS_Rain4 <- summarySE(data=Adjusted_NWS_Rain, measurevar="value", groupvars=c("variable"), na.rm = FALSE, conf.interval = .95, .drop = TRUE)

#Copy and paste two tables above to Excel and manipulate data to summarize by sampling season and the appropriate year, then reupload, and change order for graphing#
Rain_SS <- read.csv(file.choose())
Rain_SS$Sampling.Period <- factor(Rain_SS$Sampling.Period, levels = c("Aug-Sep", "Mar-Apr", "Oct-Nov", "Apr"))
Rain_SS <- Rain_SS[order(Rain_SS$Sampling.Period),]

#Make the graph#
ggplot(Rain_SS, aes(x=Sampling.Period, y=value, fill=factor(Year, levels =c("2016", "2017a", "2017b", "2018", "1990-2018"))))+
  geom_bar(position=position_dodge(), stat="identity", colour="black", size=0.3)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.3, width=0.2, position=position_dodge(0.9))+
  xlab("Sampling Period")+
  ylab("Rain Sum (mm)")+
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  ggtitle("IRL Wide NWS Average Rain Sum (mm) by Sampling Period") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(name="Year(s)", values=(c("chocolate", "seagreen", "chocolate1", "seagreen1", "grey20")), breaks=c("2016", "2017a", "2017b", "2018", "1990-2018"))

####SUMMARIZE STREAMFLOW INFO FOR TABLE S2####
#Upload streamflow information with the first column listing the Stream/Canal, the second listing one of the periods being tested or Not Tested for a monthly average that will not be used, and column three with the monlty average from USGS or SFWMD#
Streamflowsa <- read.csv(file.choose())
#Summarize by Period and Stream/Canal#
summstreamflowsa <- summarySE(data=Streamflowsa, measurevar="Monthly.Mean", groupvars=c("Stream.Canal","Period"), na.rm = FALSE, conf.interval = .95, .drop = TRUE)
#Summarize by Stream/Canal across all Periods#
summstreamflows2a <- summarySE(data=Streamflowsa, measurevar="Monthly.Mean", groupvars=c("Stream.Canal"), na.rm = FALSE, conf.interval = .95, .drop = TRUE)
#Summarize by Period across all Canals/Streams
summstreamflows3a <- summarySE(data=Streamflowsa, measurevar="Monthly.Mean", groupvars=c("Period"), na.rm = FALSE, conf.interval = .95, .drop = TRUE)
tgc$variable  <- factor(tgc$variable , levels = c("Water.Content","X..fines", "Loss.On.Ignition.aka.Total.Organic.Matter"))
write.csv(as.data.frame(summstreamflows), file="StreamFlowsSumm.csv")
write.csv(as.data.frame(summstreamflows2), file="StreamFlowsSummbyStream.csv")
write.csv(as.data.frame(summstreamflows3), file="StreamFlowsSummbyPeriod.csv")
#Manipulate data into format seen in paper in Excel#

####ENVIRONMENTAL DATA####
#Upload environmental data#
LWSS_Enviro <- read.csv(file.choose())

###TABLE S3 ENVIRONMENTAL STATISICAL DATA###
##Test for environmental SLE data for differences between sampling seasons##
#Subset SLE data#
SLE_LWSS_Enviro <- subset(LWSS_Enviro, Estuary == "SLE")
#Test porewater salinity for total significance#
kruskal.test(Porewater.salinity ~ Estuary.Sampling, data = SLE_LWSS_Enviro)
#Test porewater salinity for pairwise significance#
dunnTest(Porewater.salinity ~ Estuary.Sampling, data = SLE_LWSS_Enviro, method="bh")
#Test sediment temperature for total significance#
kruskal.test(Sediment.Temperature ~ Estuary.Sampling, data = SLE_LWSS_Enviro)
#Test porewater salinity for pairwise significance#
dunnTest(Sediment.Temperature ~ Estuary.Sampling, data = SLE_LWSS_Enviro, method="bh")

##Test for environmental IRL data for differences between sampling seasons##
#Subset IRL data#
IRL_LWSS_Enviro <- subset(LWSS_Enviro, Estuary == "IRL")
#Test porewater salinity for total significance#
kruskal.test(Porewater.salinity ~ Estuary.Sampling, data = IRL_LWSS_Enviro)
#Test porewater salinity for pairwise significance#
dunnTest(Porewater.salinity ~ Estuary.Sampling, data = IRL_LWSS_Enviro, method="bh") 
#Test sediment temperature for total significance#
kruskal.test(Sediment.Temperature ~ Estuary.Sampling, data = IRL_LWSS_Enviro)
#Test sediment temperature for pairwise significance#
dunnTest(Sediment.Temperature ~ Estuary.Sampling, data = IRL_LWSS_Enviro, method="bh")

##Test for environmental data for differences between locations##
#Test porewater salinity for total significance#
kruskal.test(Porewater.salinity ~ Location, data = LWSS_Enviro)
#Test porewater salinity for pairwise significance#
dunnTest(Porewater.salinity ~ Location, data = LWSS_Enviro, method="bh") 
#Test sediment temperature for total significance#
kruskal.test(Sediment.Temperature ~ Location, data = LWSS_Enviro)
#Test sediment temperature for pairwise significance#
dunnTest(Sediment.Temperature ~ Location, data = LWSS_Enviro, method="bh")

###MAKE FIGURE 2A###
#Subset and summarize porewater salinity and sediment temperature by Estuary by Sampling Period (melt)#
meltPWST_EbyS <- subset(LWSS_Enviro, select=c("Estuary.Sampling","Porewater.salinity", "Sediment.Temperature"))
meltPWST_EbyS <- melt(meltPWST_EbyS, id=c("Estuary.Sampling"))
#Create labels for use in graph#
PWST_EbySlabels <- c(Porewater.salinity = "Porewater Salinity", Sediment.Temperature = "Sediment Temperature")

#Plot graph#
ggplot(meltPWST_EbyS, aes(x=Estuary.Sampling, y=value, color=Estuary.Sampling)) +
  geom_boxplot() +
  facet_wrap(~variable, labeller=labeller(variable=PWST_EbySlabels)) +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_discrete(name="Estuary by Sampling Period", breaks=c("IRL-W16", "IRL-D17", "IRL-W17", "IRL-D18", "SLE-W16", "SLE-D17", "SLE-W17", "SLE-D18"), labels=c("IRL Aug/Sept 2016", "IRL Mar/Apr 2017", "IRL Oct/Nov 2017", "IRL Apr 2018", "SLE Aug/Sept 2016", "SLE Mar/Apr 2017", "SLE Oct/Nov 2017", "SLE Apr 2018")) +
  xlab("Estuary by Sampling Period") +
  scale_color_manual(values=c("seagreen4", "seagreen3", "seagreen2", "seagreen1", "chocolate4", "chocolate3", "chocolate2", "chocolate1"), limits=c("IRL-W16", "IRL-D17", "IRL-W17", "IRL-D18", "SLE-W16", "SLE-D17", "SLE-W17", "SLE-D18"), labels=c("IRL Aug/Sept 2016", "IRL Mar/Apr 2017", "IRL Oct/Nov 2017", "IRL Apr 2018", "SLE Aug/Sept 2016", "SLE Mar/Apr 2017", "SLE Oct/Nov 2017", "SLE Apr 2018")) +
  ggtitle("Porewater Salinty and Sediment Temperature by Estuary by Sampling Period") +
  ylab("Porewater Salinity (ppt) / Sediment Temperature(°C)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete( limits=c("IRL-W16", "IRL-D17", "IRL-W17", "IRL-D18", "SLE-W16", "SLE-D17", "SLE-W17", "SLE-D18"), labels=c("IRL Aug/Sept 2016", "IRL Mar/Apr 2017", "IRL Oct/Nov 2017", "IRL Apr 2018", "SLE Aug/Sept 2016", "SLE Mar/Apr 2017", "SLE Oct/Nov 2017", "SLE Apr 2018"))

#Get summary statistics for porewater salinity by Estuary by Sampling Period#
ddply(LWSS_Enviro, .(Estuary.Sampling), summarise, mean=mean(Porewater.salinity), sd=sd(Porewater.salinity), median=median(Porewater.salinity), IQR=IQR(Porewater.salinity))
#Get summary statistics for porewater salintiy by Site#
ddply(LWSS_Enviro, .(Site), summarise, mean=mean(Porewater.salinity), sd=sd(Porewater.salinity), median=median(Porewater.salinity), IQR=IQR(Porewater.salinity))

#Get summary statistics for sediment temperature by Estuary by Sampling Period#
ddply(LWSS_Enviro, .(Estuary.Sampling), summarise, mean=mean(Sediment.Temperature), sd=sd(Sediment.Temperature), median=median(Sediment.Temperature), IQR=IQR(Sediment.Temperature))

###MAKE FIGURE 2B###
#Subset and summarize porewater salinity and sediment temperature by Estuary by Sampling Period (melt)#
meltPWST_Loc <- subset(LWSS_Enviro, select=c("Location","Porewater.salinity", "Sediment.Temperature"))
meltPWST_Loc <- melt(meltPWST_Loc, id=c("Location"))
#Make the graph#
ggplot(meltPWST_Loc, aes(x=Location, y=value, color=Location)) +
  geom_boxplot() +
  facet_wrap(~variable, labeller=labeller(variable=PWSTlabels)) +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  scale_color_discrete(name="Location", breaks=c("North", "North Central", "South Central", "South", "SLE")) +
  xlab("Location") +
  scale_color_manual(values=c("seagreen4", "seagreen3", "seagreen2", "seagreen1", "chocolate4"), limits=c("North", "North Central", "South Central", "South", "SLE")) +
  ggtitle("Porewater Salinty and Sediment Temperature by Location") +
  ylab("Porewater Salinity (ppt) / Sediment Temperature(°C)") +
  scale_x_discrete( limits=c("North", "North Central", "South Central", "South", "SLE")) +
  theme(plot.title = element_text(hjust = 0.5))

#Get summary statistics for porewater salinity by Location#
ddply(LWSS_Enviro, .(Location), summarise, mean=mean(Porewater.salinity), sd=sd(Porewater.salinity), median=median(Porewater.salinity), IQR=IQR(Porewater.salinity))
#Get summary statistics for porewater salinity by Location#
ddply(LWSS_Enviro, .(Location), summarise, mean=mean(Sediment.Temperature), sd=sd(Sediment.Temperature), median=median(Sediment.Temperature), IQR=IQR(Sediment.Temperature))

###MAKE FIGURE 3###
#Subset and summarize muck characteristics by Site (melt)#
meltMuck <- subset(LWSS_Enviro, select=c("Site","X..fines", "Water.Content", "Loss.On.Ignition.aka.Total.Organic.Matter"))
meltMuck <- melt(meltMuck, id=c("Site")) 

#Create labels for graph#
Mucklabels<- c(X..fines = "Percent Fines", Water.Content = "Water Content", Loss.On.Ignition.aka.Total.Organic.Matter = "Total Organic Matter")

#Make melted table a dataframe#
meltMuck <- as.data.frame(meltMuck)

#Get summary statistics for each muck characteristic by Site#
tgc <- summarySE(data=meltMuck, measurevar="value", groupvars=c("variable","Site"), na.rm = FALSE, conf.interval = .95, .drop = TRUE)
tgc <-as.data.frame(tgc)
tgc$variable  <- factor(tgc$variable , levels = c("Water.Content","X..fines", "Loss.On.Ignition.aka.Total.Organic.Matter"))

#Make the graph#
ggplot(tgc, aes(x=Site, y=value, fill=variable))+
  geom_bar(position=position_dodge(), stat="identity", colour="black", size=0.3)+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), size=0.3, width=0.2, position=position_dodge(0.9))+
  xlab("Site")+
  ylab("Percentage")+
  theme(axis.text.x  = element_text(angle=45, hjust=1))+
  scale_fill_discrete(name="Muck Characteristic", breaks=c("Water.Content", "X..fines", "Loss.On.Ignition.aka.Total.Organic.Matter"), labels=c("Water Content", "Percent Fines", "Total Organic Matter")) + 
  geom_hline(yintercept=60, linetype="dashed", color = "darkgreen", size = 1) +
  geom_hline(yintercept=75, linetype="dashed", color = "red", size = 1) +
  geom_hline(yintercept=10, linetype="dashed", color = "dodgerblue", size = 1) +
  ggtitle("Summary of Muck Characteristics by Site") +
  theme(plot.title = element_text(hjust = 0.5))

###MAKE FIGURE 4###
ggplot(data = LWSS_Enviro, aes( x = Loss.On.Ignition.aka.Total.Organic.Matter, y = Cu, color=Site, shape=Muck)) +
  geom_point(size=2) +
  geom_hline(yintercept=65, linetype="dashed", color = "red") +
  geom_vline(xintercept=10, linetype="dashed", color = "blue") +
  xlab("Total Organic Matter (%)") +
  ylab("Copper Concentration (µg/g)") +
  ggtitle("Copper Concentration vs. Total Organic Matter") +
  scale_colour_manual(values=c("midnightblue", "blue", "red", "purple", "yellow", "violet", "maroon", "darkgreen", "dodgerblue", "gold",  "moccasin",  "darkgray", "powderblue", "green", "lightgray", "black", "slateblue", "greenyellow", "lightseagreen")) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_shape_discrete(name="Number of Muck Characteristics", limits=c("Muck", "Mucky", "Muckish", "Not"), labels=c("3", "2", "1", "0"))

###ALPHA DIVERSITY###
"Estimate alpha diversity and export it"
FLWSS_alpha_diversity = estimate_richness(FLWSSDataA, split = TRUE)
write.csv(FLWSS_alpha_diversity, file='LWSS_alpha_div.csv')

#Add sample data to alpha diversity table#
RFLWSS_alpha_diversity = merge(FLWSS_alpha_diversity, dotted_sample_data, by = 0) 
RFLWSS_alpha_diversity = as.data.frame(RFLWSS_alpha_diversity)

#Test for normality in all alpha diversity statistics#
shapiro.test(RFLWSS_alpha_diversity$Shannon)
shapiro.test(RFLWSS_alpha_diversity$Observed)
shapiro.test(RFLWSS_alpha_diversity$Fisher)
shapiro.test(RFLWSS_alpha_diversity$Simpson)
shapiro.test(RFLWSS_alpha_diversity$Chao1)

#Get summary statistics for all alpha diversity metrics#
ddply(RFLWSS_alpha_diversity, .(Medium), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))
ddply(RFLWSS_alpha_diversity, .(Medium), summarise, mean=mean(Observed), sd=sd(Observed), median=median(Observed), IQR=IQR(Observed))
ddply(RFLWSS_alpha_diversity, .(Medium), summarise, mean=mean(Fisher), sd=sd(Fisher), median=median(Fisher), IQR=IQR(Fisher))
ddply(RFLWSS_alpha_diversity, .(Medium), summarise, mean=mean(Simpson), sd=sd(Simpson), median=median(Simpson), IQR=IQR(Simpson))
ddply(RFLWSS_alpha_diversity, .(Medium), summarise, mean=mean(Chao1), sd=sd(Chao1), median=median(Chao1), IQR=IQR(Chao1))

#Test for correlations between Shannon diversity and other metrics#
cor.test(RFLWSS_alpha_diversity$Shannon, RFLWSS_alpha_diversity$Observed, method = "spearman", exact = FALSE)
cor.test(RFLWSS_alpha_diversity$Shannon, RFLWSS_alpha_diversity$Fisher, method = "spearman", exact = FALSE)
cor.test(RFLWSS_alpha_diversity$Shannon, RFLWSS_alpha_diversity$Simpson, method = "spearman", exact = FALSE)
cor.test(RFLWSS_alpha_diversity$Shannon, RFLWSS_alpha_diversity$Chao1, method = "spearman", exact = FALSE)

#Multiple test adjustment for correlation tests#
alpha_div_cor_p.value=c(2.2e-16,2.2e-16,2.2e-16,2.2e-16)
alpha_div_cor_p.value=data.frame(alpha_div_cor_p.value)
alpha_div_cor_p.value$padj <- p.adjust(alpha_div_cor_p.value$alpha_div_cor_p.value, method = "BH")
alpha_div_cor_p.value

###GET STATISTICS FOR TABLE S3 ALPHA DIVERSITY SECTION###
#Test Shannon Diverstiy at the Muck level#
#Test for overall significance#
kruskal.test(Shannon ~ Muck, data = RFLWSS_alpha_diversity)
#Test for pairwise significance#
dunnTest(Shannon ~ Muck, data = RFLWSS_alpha_diversity, method="bh")
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(Muck), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

#Test Shannon Diversity at Estuary level#
#Test for pairwise significance#
wilcox.test(Shannon ~ Estuary, data = RFLWSS_alpha_diversity)
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(Estuary), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

#Test Shannon Diversity at LOI.Cu level#
#Test for overall significance#
kruskal.test(Shannon ~ LOI.Cu, data = RFLWSS_alpha_diversity)
#Test for pairwise significance#
dunnTest(Shannon ~ LOI.Cu, data = RFLWSS_alpha_diversity, method="bh")
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(LOI.Cu), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

#Test Shannon Diversity at Site level#
#Test for overal significance#
kruskal.test(Shannon ~ Site, data = RFLWSS_alpha_diversity)
#Test for pairwise significance#
dunnTest(Shannon ~ Site, data = RFLWSS_alpha_diversity, method="bh")
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(Site), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

#Test Shannon Diversity at Location level#
#Test for overal significance#
kruskal.test(Shannon ~ Location, data = RFLWSS_alpha_diversity)
#Test for pairwaise significance#
dunnTest(Shannon ~ Location, data = RFLWSS_alpha_diversity, method="bh")
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(Location), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

#Test Shannon Diversity at Season level
#Test for total significance#
wilcox.test(Shannon ~ Season, data = RFLWSS_alpha_diversity)
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(Season), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))


#Test Shannon Diversity at Estuary.Season level#
#Test for total significance#
kruskal.test(Shannon ~ Estuary.Season, data = RFLWSS_alpha_diversity)
#Test for pairwise significance#
dunnTest(Shannon ~ Estuary.Season, data = RFLWSS_alpha_diversity, method="bh") 
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(Estuary.Season), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

#Test Shannon Diversity at Sampling level#
#Test for total significance#
kruskal.test(Shannon ~ AbbSeasonAbbYear, data = RFLWSS_alpha_diversity)
#Test for pairwise significance#
dunnTest(Shannon ~ AbbSeasonAbbYear, data = RFLWSS_alpha_diversity, method="bh") 
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(AbbSeasonAbbYear), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

#Create IRL and SLE focused tables#
SLE_RFLWSS_alpha_diversity <- subset(RFLWSS_alpha_diversity, Estuary == "SLE")
IRL_RFLWSS_alpha_diversity <- subset(RFLWSS_alpha_diversity, Estuary == "IRL")
#Test for total significance in SLE#
kruskal.test(Shannon ~ Estuary.Sampling, data = SLE_RFLWSS_alpha_diversity)
#Test for pairwise significance in SLE#
dunnTest(Shannon ~ Estuary.Sampling, data = SLE_RFLWSS_alpha_diversity, method="bh") 
#Test for total significance in IRL#
kruskal.test(Shannon ~ Estuary.Sampling, data = IRL_RFLWSS_alpha_diversity)
#Test for pairwise significance in IRL#
dunnTest(Shannon ~ Estuary.Sampling, data = IRL_RFLWSS_alpha_diversity, method="bh") 
#Determine summary stats#
ddply(RFLWSS_alpha_diversity, .(Estuary.Sampling), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

#Multiple test adjustment for Kruskal-Wallis tests#
alpha_div_p.value=c(0.01432,0.9244,0.0269,0.7267,0.005506,1.954e-07,2.20E-16, 1.002e-15, 0.0823, 3.493e-08, 2.2e-16, 2.2e-16, 1.149e-09, 7.842e-13, 0.001394)
alpha_div_p.value=data.frame(alpha_div_p.value)
alpha_div_p.value$padj <- p.adjust(alpha_div_p.value$alpha_div_p.value, method = "BH")
alpha_div_p.value

###MAKE FIGURE 5A###
#Make a boxplot#
ggplot(RFLWSS_alpha_diversity, aes(x=Estuary.Sampling, y=Shannon, color=Estuary.Sampling)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  xlab("Estuary by Sampling Period") +
  scale_color_manual(values=c("seagreen4", "seagreen3", "seagreen2", "seagreen1", "chocolate4", "chocolate3", "chocolate2", "chocolate1"), limits=c("IRL-W16", "IRL-D17", "IRL-W17", "IRL-D18", "SLE-W16", "SLE-D17", "SLE-W17", "SLE-D18"), labels=c("IRL Aug/Sept 2016", "IRL Mar/Apr 2017", "IRL Oct/Nov 2017", "IRL Apr 2018", "SLE Aug/Sept 2016", "SLE Mar/Apr 2017", "SLE Oct/Nov 2017", "SLE Apr 2018")) +
  scale_x_discrete(limits=c("IRL-W16", "IRL-D17", "IRL-W17", "IRL-D18", "SLE-W16", "SLE-D17", "SLE-W17", "SLE-D18"), labels=c("IRL Aug/Sept 2016", "IRL Mar/Apr 2017", "IRL Oct/Nov 2017", "IRL Apr 2018", "SLE Aug/Sept 2016", "SLE Mar/Apr 2017", "SLE Oct/Nov 2017", "SLE Apr 2018"))

#MAKE FIGURE 5B#
ggplot(RFLWSS_alpha_diversity, aes(x=Location, y=Shannon, color=Location)) + 
  scale_color_manual(values=c("seagreen4", "seagreen3", "seagreen2", "seagreen1", "chocolate4"), labels=c("North IRL", "North Central IRL", "South Central IRL", "South IRL", "SLE"), limits=c("North", "North Central", "South Central", "South", "SLE")) +
  scale_x_discrete(limits=c("North", "North Central", "South Central", "South", "SLE"), labels=c("North IRL", "North Central IRL", "South Central IRL", "South IRL", "SLE")) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))

###MAKE SETS OF MERGED BAR GRAPHS BY METADATA CATEGORY FOR LWSS SURVEY###

# get abundance in relative percentage#
LWSSphy <- transform_sample_counts(FLWSSDataA, function(x) 100*x/sum(x))


###CREATE FIGURES 6A AND S4A-D###
##CREATE ESTUARY FOCUSED TABLES FOR BAR GRAPHS##
##Create table ready for making stacked bar graph for Kingdoms >1% by Estuary##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Kingdom')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Kingdom to a character vector from a factor because R
dat$Kingdom <- as.character(dat$Kingdom)
# group dataframe by Kingdom, calculate mean rel. abundance
means <- ddply(dat, ~Kingdom, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Kingdom
# change their name to "Other Prokaryotes"
dat[dat$Kingdom %in% Other,]$Kingdom <- 'Other Prokaryotes'
#remove all Kingdomes labeled Other Prokaryotes
dat <-dat[!dat$Kingdom=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Kingdom))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Estuary","Kingdom"), na.rm = TRUE)
#remove unncessary columns
Kingdomdat2 <- subset(dat, select=c(Estuary, Abundance, Kingdom))


##Create table ready for making stacked bar graph for Phylums >1%  by Estuary##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Phylum')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)
# group dataframe by Phylum, calculate mean rel. abundance
means <- ddply(dat, ~Phylum, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Phylum
# change their name to "Other Prokaryotes"
dat[dat$Phylum %in% Other,]$Phylum <- 'Other Prokaryotes'
#remove all Phylums labeled Other Prokaryotes
dat <-dat[!dat$Phylum=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Phylum))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Estuary","Phylum"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Estuary, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Phylum<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Phylum))
#combine with original table
OnePhylumdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Classes >1% by Estuary##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Class')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Class to a character vector from a factor because R
dat$Class <- as.character(dat$Class)
# group dataframe by Class, calculate mean rel. abundance
means <- ddply(dat, ~Class, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Class
# change their name to "Other Prokaryotes"
dat[dat$Class %in% Other,]$Class <- 'Other Prokaryotes'
#remove all Classes labeled Other Prokaryotes
dat <-dat[!dat$Class=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Class))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Estuary","Class"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Estuary, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Class<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Class))
#combine with original table
OneClassdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Orders >1% by Estuary##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Order')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Order to a character vector from a factor because R
dat$Order <- as.character(dat$Order)
# group dataframe by Order, calculate mean rel. abundance
means <- ddply(dat, ~Order, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Order
# change their name to "Other Prokaryotes"
dat[dat$Order %in% Other,]$Order <- 'Other Prokaryotes'
#remove all Orders labeled Other Prokaryotes
dat <-dat[!dat$Order=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Order))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Estuary","Order"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Estuary, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Order<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Order))
#combine with original table
OneOrderdatAb23 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Families >1% by Estuary##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Family')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Family to a character vector from a factor because R
dat$Family <- as.character(dat$Family)
# group dataframe by Family, calculate mean rel. abundance
means <- ddply(dat, ~Family, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Family
# change their name to "Other Prokaryotes"
dat[dat$Family %in% Other,]$Family <- 'Other Prokaryotes'
#remove all Families labeled Other Prokaryotes
dat <-dat[!dat$Family=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Family))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Estuary","Family"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Estuary, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Family<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Family))
#combine with original table
OneFamilydatAb23 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Genera >1%by Estuary##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Estuary","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Estuary, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Genus))
#combine with original table
OneGenusdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Species >1% by Estuary##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family, Genus, Species levels#
dat<- unite(dat, Genus, Family:Genus, sep=';')
dat<- unite(dat, Species, Genus:Species, sep=';')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
# group dataframe by Species, calculate mean rel. abundance
means <- ddply(dat, ~Species, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Species
# change their name to "Other Prokaryotes"
dat[dat$Species %in% Other,]$Species <- 'Other Prokaryotes'
#remove all Species labeled Other Prokaryotes
dat <-dat[!dat$Species=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Estuary","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Estuary, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Species<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Estuary, Abundance, Species))
#combine with original table
OneSpeciesdatAb2 <- rbind(dat, Abundance)

##CREATE ESTUARY FOCUSED BAR GRAPHS##
##Make a series of graphs using ggplot based on above tables##

#Make IRL Phylum graph, FIGURE S2A#
spatial_plot_LWSS_Phylum2 <- ggplot(data=OnePhylumdatAb2, aes(x=Estuary, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple",  "rosybrown"))+
  xlab("Estuary") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))
spatial_plot_LWSS_Phylum2

#Make IRL Class graph, FIGURE S2B#
spatial_plot_LWSS_Class2 <- ggplot(OneClassdatAb2, aes(x=Estuary, y=Abundance, fill=Class)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "rosybrown"))+
  xlab("Estuary") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))
spatial_plot_LWSS_Class2

#Make IRL Order graph, FIGURE 6A#
spatial_plot_LWSS_Order2 <- ggplot(data=OneOrderdatAb2, aes(x=Estuary, y=Abundance, fill=Order)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "darkorchid3", "gold3", "darkolivegreen", "hotpink", "lawngreen", "rosybrown"))+
  xlab("Estuary") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))
spatial_plot_LWSS_Order2

#Make IRL Family graph, FIGURE S2C#
spatial_plot_LWSS_Family2 <- ggplot(data=OneFamilydatAb2, aes(x=Estuary, y=Abundance, fill=Family)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "rosybrown"))+
  xlab("Estuary") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))
spatial_plot_LWSS_Family2

#Make IRL Genus Graph, FIGURE S2D#
spatial_plot_LWSS_Genus2 <- ggplot(data=OneGenusdatAb2, aes(x=Estuary, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "rosybrown"))+
  xlab("Estuary") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_LWSS_Genus2

##CREATE ESTUARY BY SAMPLING PERIOD TABLE FOR BAR GRAPH, FIGURE S3##
##Create table ready for making stacked bar graph for Orders >1% by Estuary.Sampling##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Order')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Order to a character vector from a factor because R
dat$Order <- as.character(dat$Order)
# group dataframe by Order, calculate mean rel. abundance
means <- ddply(dat, ~Order, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Order
# change their name to "Other Prokaryotes"
dat[dat$Order %in% Other,]$Order <- 'Other Prokaryotes'
#remove all Orders labeled Other Prokaryotes
dat <-dat[!dat$Order=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Estuary.Sampling, Abundance, Order))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Estuary.Sampling","Order"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Estuary.Sampling, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Order<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Estuary.Sampling, Abundance, Order))
#combine with original table
Estuary.SamplingOrderdatAb2 <- rbind(dat, Abundance)

##CREATE ESTUARY BY SAMPLING PERIOD BAR GRAPH##
spatial_plot_LWSS_Estuary.Sampling_Order2 <- ggplot(data=Estuary.SamplingOrderdatAb2, aes(x=Estuary.Sampling, y=Abundance, fill=Order)) +
  geom_bar(aes(), stat="identity", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "darkorchid3", "gold3", "darkolivegreen", "hotpink", "lawngreen", "rosybrown"))+
  xlab("Estuary by Sampling Period") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("IRL-W16", "IRL-D17", "IRL-W17", "IRL-D18", "SLE-W16", "SLE-D17", "SLE-W17", "SLE-D18"), labels=c("IRL Aug/Sept 2016", "IRL Mar/Apr 2017", "IRL Oct/Nov 2017", "IRL Apr 2018", "SLE Aug/Sept 2016", "SLE Mar/Apr 2017", "SLE Oct/Nov 2017", "SLE Apr 2018"))
spatial_plot_LWSS_Estuary.Sampling_Order2

###CREATE FIGURES 6B AND S4A-D
##CREATE MUCK FOCUSED TABLES FOR BAR GRAPHS##
##Create table ready for making stacked bar graph for Phylums >1%  by Muck##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Phylum')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)
# group dataframe by Phylum, calculate mean rel. abundance
means <- ddply(dat, ~Phylum, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Phylum
# change their name to "Other Prokaryotes"
dat[dat$Phylum %in% Other,]$Phylum <- 'Other Prokaryotes'
#remove all Phylums labeled Other Prokaryotes
dat <-dat[!dat$Phylum=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Phylum))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Muck","Phylum"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Muck, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Phylum<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Phylum))
#combine with original table
MuckPhylumdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Classes >1% by Muck##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Class')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Class to a character vector from a factor because R
dat$Class <- as.character(dat$Class)
# group dataframe by Class, calculate mean rel. abundance
means <- ddply(dat, ~Class, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Class
# change their name to "Other Prokaryotes"
dat[dat$Class %in% Other,]$Class <- 'Other Prokaryotes'
#remove all Classes labeled Other Prokaryotes
dat <-dat[!dat$Class=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Class))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Muck","Class"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Muck, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Class<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Class))
#combine with original table
MuckClassdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Orders >1% by Muck##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Order')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Order to a character vector from a factor because R
dat$Order <- as.character(dat$Order)
# group dataframe by Order, calculate mean rel. abundance
means <- ddply(dat, ~Order, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Order
# change their name to "Other Prokaryotes"
dat[dat$Order %in% Other,]$Order <- 'Other Prokaryotes'
#remove all Orders labeled Other Prokaryotes
dat <-dat[!dat$Order=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Order))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Muck","Order"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Muck, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Order<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Order))
#combine with original table
MuckOrderdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Families >1% by Muck##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Family')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Family to a character vector from a factor because R
dat$Family <- as.character(dat$Family)
# group dataframe by Family, calculate mean rel. abundance
means <- ddply(dat, ~Family, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Family
# change their name to "Other Prokaryotes"
dat[dat$Family %in% Other,]$Family <- 'Other Prokaryotes'
#remove all Families labeled Other Prokaryotes
dat <-dat[!dat$Family=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Family))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Muck","Family"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Muck, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Family<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Family))
#combine with original table
MuckFamilydatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Genera >1%by Muck##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Muck","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Muck, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Genus))
#combine with original table
MuckGenusdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Species >1% by Muck##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family, Genus, Species levels#
dat<- unite(dat, Genus, Family:Genus, sep=';')
dat<- unite(dat, Species, Genus:Species, sep=';')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
# group dataframe by Species, calculate mean rel. abundance
means <- ddply(dat, ~Species, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Species
# change their name to "Other Prokaryotes"
dat[dat$Species %in% Other,]$Species <- 'Other Prokaryotes'
#remove all Species labeled Other Prokaryotes
dat <-dat[!dat$Species=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Muck","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Muck, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Species<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Muck, Abundance, Species))
#combine with original table
MuckSpeciesdatAb2 <- rbind(dat, Abundance)


##CREATE MUCK FOCUSED BAR GRAPHS##
#Make LWSS Muck Phylum graph, FIGURE S4A#
MuckPhylumdatAb2M <- filter(MuckPhylumdatAb2, Muck==c("Muck"))
MuckPhylumdatAb2N <- filter(MuckPhylumdatAb2, Muck==c("Not"))
MuckPhylumdatAb2MN <- rbind(MuckPhylumdatAb2M, MuckPhylumdatAb2N)
spatial_plot_LWSS_Muck_Phylum2 <- ggplot(data=MuckPhylumdatAb2MN, aes(x=Muck, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "darkorchid3", "gold3", "darkolivegreen", "hotpink", "lawngreen", "rosybrown"))+
  xlab("Muck") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))
spatial_plot_LWSS_Muck_Phylum2

#Make LWSS Muck Class graph, FIGURE S4B#
MuckClassdatAb2M <- filter(MuckClassdatAb2, Muck==c("Muck"))
MuckClassdatAb2N <- filter(MuckClassdatAb2, Muck==c("Not"))
MuckClassdatAb2MN <- rbind(MuckClassdatAb2M, MuckClassdatAb2N)
spatial_plot_LWSS_Muck_Class2 <- ggplot(MuckClassdatAb2MN, aes(x=Muck, y=Abundance, fill=Class)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "rosybrown"))+
  xlab("Muck") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))
spatial_plot_LWSS_Muck_Class2

#Make LWSS Muck Order graph, FIGURE 6B#
MuckOrderdatAb2M <- filter(MuckOrderdatAb2, Muck==c("Muck"))
MuckOrderdatAb2N <- filter(MuckOrderdatAb2, Muck==c("Not"))
MuckOrderdatAb2MN <- rbind(MuckOrderdatAb2M, MuckOrderdatAb2N)
spatial_plot_LWSS_Muck_Order2 <- ggplot(data=MuckOrderdatAb2MN, aes(x=Muck, y=Abundance, fill=Order)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "darkorchid3", "gold3", "darkolivegreen", "hotpink", "lawngreen", "rosybrown"))+
  xlab("Muck") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))
spatial_plot_LWSS_Muck_Order2

#Make LWSS Muck Family graph, FIGURE S4C#
MuckFamilydatAb2M <- filter(MuckFamilydatAb2, Muck==c("Muck"))
MuckFamilydatAb2N <- filter(MuckFamilydatAb2, Muck==c("Not"))
MuckFamilydatAb2MN <- rbind(MuckFamilydatAb2M, MuckFamilydatAb2N)
spatial_plot_LWSS_Muck_Family2 <- ggplot(data=MuckFamilydatAb2MN, aes(x=Muck, y=Abundance, fill=Family)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "rosybrown"))+
  xlab("Muck") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))
spatial_plot_LWSS_Muck_Family2

#Make LWSS Muck Genus Graph, FIGURE S4D#
MuckGenusdatAb2M <- filter(MuckGenusdatAb2, Muck==c("Muck"))
MuckGenusdatAb2N <- filter(MuckGenusdatAb2, Muck==c("Not"))
MuckGenusdatAb2MN <- rbind(MuckGenusdatAb2M, MuckGenusdatAb2N)
spatial_plot_LWSS_Muck_Genus2 <- ggplot(data=MuckGenusdatAb2MN, aes(x=Muck, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "rosybrown"))+
  xlab("Muck") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_LWSS_Muck_Genus2

###CREATE FIGURES 6C AND S5A-D
##CREATE LOI.Cu FOCUSED TABLES FOR BAR GRAPHS##
##Create table ready for making stacked bar graph for Phylums >1%  by LOI.Cu##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Phylum')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)
# group dataframe by Phylum, calculate mean rel. abundance
means <- ddply(dat, ~Phylum, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Phylum
# change their name to "Other Prokaryotes"
dat[dat$Phylum %in% Other,]$Phylum <- 'Other Prokaryotes'
#remove all Phylums labeled Other Prokaryotes
dat <-dat[!dat$Phylum=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Phylum))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("LOI.Cu","Phylum"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~LOI.Cu, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Phylum<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Phylum))
#combine with original table
LOI.CuPhylumdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Classes >1% by LOI.Cu##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Class')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Class to a character vector from a factor because R
dat$Class <- as.character(dat$Class)
# group dataframe by Class, calculate mean rel. abundance
means <- ddply(dat, ~Class, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Class
# change their name to "Other Prokaryotes"
dat[dat$Class %in% Other,]$Class <- 'Other Prokaryotes'
#remove all Classes labeled Other Prokaryotes
dat <-dat[!dat$Class=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Class))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("LOI.Cu","Class"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~LOI.Cu, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Class<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Class))
#combine with original table
LOI.CuClassdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Orders >1% by LOI.Cu##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Order')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Order to a character vector from a factor because R
dat$Order <- as.character(dat$Order)
# group dataframe by Order, calculate mean rel. abundance
means <- ddply(dat, ~Order, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Order
# change their name to "Other Prokaryotes"
dat[dat$Order %in% Other,]$Order <- 'Other Prokaryotes'
#remove all Orders labeled Other Prokaryotes
dat <-dat[!dat$Order=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Order))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("LOI.Cu","Order"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~LOI.Cu, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Order<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Order))
#combine with original table
LOI.CuOrderdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Families >1% by LOI.Cu##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Family')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Family to a character vector from a factor because R
dat$Family <- as.character(dat$Family)
# group dataframe by Family, calculate mean rel. abundance
means <- ddply(dat, ~Family, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Family
# change their name to "Other Prokaryotes"
dat[dat$Family %in% Other,]$Family <- 'Other Prokaryotes'
#remove all Families labeled Other Prokaryotes
dat <-dat[!dat$Family=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Family))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("LOI.Cu","Family"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~LOI.Cu, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Family<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Family))
#combine with original table
LOI.CuFamilydatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Genera >1%by LOI.Cu##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("LOI.Cu","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~LOI.Cu, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Genus))
#combine with original table
LOI.CuGenusdatAb2 <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Species >1% by LOI.Cu##
# agglomerate taxa
glom <- tax_glom(LWSSphy, taxrank = 'Species')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family, Genus, Species levels#
dat<- unite(dat, Genus, Family:Genus, sep=';')
dat<- unite(dat, Species, Genus:Species, sep=';')
# convert Species to a character vector from a factor because R
dat$Species <- as.character(dat$Species)
# group dataframe by Species, calculate mean rel. abundance
means <- ddply(dat, ~Species, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Species
# change their name to "Other Prokaryotes"
dat[dat$Species %in% Other,]$Species <- 'Other Prokaryotes'
#remove all Species labeled Other Prokaryotes
dat <-dat[!dat$Species=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Species))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("LOI.Cu","Species"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~LOI.Cu, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Species<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(LOI.Cu, Abundance, Species))
#combine with original table
LOI.CuSpeciesdatAb2 <- rbind(dat, Abundance)


##CREATE LOI.Cu FOCUSED BAR GRAPHS##
#Make LWSS LOI.Cu Phylum graph, FIGURE S5A#
LOI.CuPhylumdatAb2HH <- filter(LOI.CuPhylumdatAb2, LOI.Cu==c("High-High"))
LOI.CuPhylumdatAb2HL <- filter(LOI.CuPhylumdatAb2, LOI.Cu==c("High-Low"))
LOI.CuPhylumdatAb2HHHL <- rbind(LOI.CuPhylumdatAb2HH, LOI.CuPhylumdatAb2HL)
spatial_plot_LWSS_LOI.Cu_Phylum2 <- ggplot(data=LOI.CuPhylumdatAb2HHHL, aes(x=LOI.Cu, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "rosybrown"))+
  xlab("TOM/Cu") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("High/High", "High/Low"))

spatial_plot_LWSS_LOI.Cu_Phylum2

#Make LWSS LOI.Cu Class graph, FIGURE S5B#
LOI.CuClassdatAb2HH <- filter(LOI.CuClassdatAb2, LOI.Cu==c("High-High"))
LOI.CuClassdatAb2HL <- filter(LOI.CuClassdatAb2, LOI.Cu==c("High-Low"))
LOI.CuClassdatAb2HHHL <- rbind(LOI.CuClassdatAb2HH, LOI.CuClassdatAb2HL)
spatial_plot_LWSS_LOI.Cu_Class2 <- ggplot(LOI.CuClassdatAb2HHHL, aes(x=LOI.Cu, y=Abundance, fill=Class)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "rosybrown"))+
  xlab("TOM/Cu") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(labels=c("High/High", "High/Low"))
spatial_plot_LWSS_LOI.Cu_Class2

#Make LWSS LOI.Cu Order graph, FIGURE 6C#
LOI.CuOrderdatAb2HH <- filter(LOI.CuOrderdatAb2, LOI.Cu==c("High-High"))
LOI.CuOrderdatAb2HL <- filter(LOI.CuOrderdatAb2, LOI.Cu==c("High-Low"))
LOI.CuOrderdatAb2HHHL <- rbind(LOI.CuOrderdatAb2HH, LOI.CuOrderdatAb2HL)
spatial_plot_LWSS_LOI.Cu_Order2 <- ggplot(data=LOI.CuOrderdatAb2HHHL, aes(x=LOI.Cu, y=Abundance, fill=Order)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "darkorchid3", "gold3", "darkolivegreen", "hotpink", "lawngreen", "rosybrown"))+
  xlab("TOM/Cu") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(labels=c("High/High", "High/Low"))
spatial_plot_LWSS_LOI.Cu_Order2

#Make LWSS LOI.Cu Family graph, FIGURE S5C#
LOI.CuFamilydatAb2HH <- filter(LOI.CuFamilydatAb2, LOI.Cu==c("High-High"))
LOI.CuFamilydatAb2HL <- filter(LOI.CuFamilydatAb2, LOI.Cu==c("High-Low"))
LOI.CuFamilydatAb2HHHL <- rbind(LOI.CuFamilydatAb2HH, LOI.CuFamilydatAb2HL)
spatial_plot_LWSS_LOI.Cu_Family2 <- ggplot(data=LOI.CuFamilydatAb2HHHL, aes(x=LOI.Cu, y=Abundance, fill=Family)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "rosybrown"))+
  xlab("TOM/Cu") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(labels=c("High/High", "High/Low"))
spatial_plot_LWSS_LOI.Cu_Family2

#Make LWSS LOI.Cu Genus Graph, FIGURE S5D#
LOI.CuGenusdatAb2HH <- filter(LOI.CuGenusdatAb2, LOI.Cu==c("High-High"))
LOI.CuGenusdatAb2HL <- filter(LOI.CuGenusdatAb2, LOI.Cu==c("High-Low"))
LOI.CuGenusdatAb2HHHL <- rbind(LOI.CuGenusdatAb2HH, LOI.CuGenusdatAb2HL)
spatial_plot_LWSS_LOI.Cu_Genus2 <- ggplot(data=LOI.CuGenusdatAb2HHHL, aes(x=LOI.Cu, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("legend", values=c("salmon", "plum", "saddlebrown", "wheat", "deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "rosybrown"))+
  xlab("TOM/Cu") +
  ylab("Percentage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("High/High", "High/Low"))
spatial_plot_LWSS_LOI.Cu_Genus2

