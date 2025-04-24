# ALOUETTE RESERVOIR
# 2024 REPORT FIGURES AND STATS - 2024 Report
# phyto only
# Jennifer Sarchuk, BC Ministry of WLRS
# version used on work computer - R version 4.4.2
#####################################################

library("doBy")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("reshape2")
library("scales")

rm(list=ls(all=TRUE))

theme_set(theme_bw())

setwd("C:/Users/jsarchuk/OneDrive - Government of BC/ALU 2024/R files/alu.24")

### LOADING DATA FILES
phyto <- read.csv('C:/Users/jsarchuk/OneDrive - Government of BC/ALU 2024/R files/ALU_DAT_R_phyto.csv')

### FORMAT COLUMNS

phyto$date <- as.Date(phyto$date)
phyto$year <- factor(phyto$year)
phyto$month <- factor(phyto$month,levels=c('Jan','Feb','Mar','Apr','May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
phyto$basin <- factor(phyto$basin,levels=c('NB', 'SB'))
phyto$type <- factor(phyto$type)
phyto$class.alias <- factor(phyto$class.alias)
phyto$class.name <- factor(phyto$class.name)
phyto$spp <- factor(phyto$spp)
phyto$edibility <- factor(phyto$edibility,levels=c('B', 'I', 'E'))

### ADD FUNCTIONS

f.sum.stat <- function(x, ...){c(mean=mean(x, ...), sd=sd(x, ...), max=max(x, ...), min=min(x, ...))}
f.max.min <- function(x, ...){c(max=max(x, ...), min=min(x, ...))}
f.total <- function(x, ...){c(total=sum(x, ...))}

#####################################################
### PHYTOPLANKTON

# annual species richness
phyto.spp <- select(phyto, year, spp)
phyto.spp <- group_by(phyto.spp, year)
phyto.spp <- summarise(phyto.spp, annual.spp.count = n_distinct(spp))
write.csv(phyto.spp, "tables/ALU_DAT_R_phyto-annual-spp-richness.csv", row.names=FALSE)

phyto.spp.24 <- filter(phyto.spp, year=="2024")

#Figure species richness
f.phyto.spp <- ggplot(phyto.spp, aes(x=year, y=annual.spp.count))
f.phyto.spp <- f.phyto.spp + geom_point(shape=1, size=3.5, colour="black")
f.phyto.spp <- f.phyto.spp + scale_y_continuous(name="Phytoplankton Species Richness (no. species)", limits=c(0,70), breaks=seq(from=0, to=70, by=10))
f.phyto.spp <- f.phyto.spp + xlab("Year")
f.phyto.spp <- f.phyto.spp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                     axis.title.x = element_text(face = "bold"),
                                     axis.title.y = element_text(face = "bold"),
                                     legend.title = element_blank(),
                                     legend.position = "top")
f.phyto.spp

ggsave("plots/alu_phyto-spp-rich.png",width=8,height=5,units="in", dpi=600)

#2024 only
f.phyto.spp.24 <- ggplot(phyto.spp.24, aes(x=year, y=annual.spp.count))
f.phyto.spp.24 <- f.phyto.spp.24 + geom_point(shape=1, size=3.5, colour="black")
f.phyto.spp.24 <- f.phyto.spp.24 + scale_y_continuous(name="Phytoplankton Species Richness (no. species)", limits=c(0,70), breaks=seq(from=0, to=70, by=10))
f.phyto.spp.24 <- f.phyto.spp.24 + xlab("Year")
f.phyto.spp.24 <- f.phyto.spp.24 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                   axis.title.x = element_text(face = "bold"),
                                   axis.title.y = element_text(face = "bold"),
                                   legend.title = element_blank(),
                                   legend.position = "top")
f.phyto.spp.24

ggsave("plots/alu_phyto-spp-rich.24.png",width=6,height=5,units="in", dpi=600)


# annual species occurences
phyto.spp.id <- select(phyto, year, spp)
phyto.spp.id <- dcast(phyto.spp.id, spp~year, fun.aggregate=length)

#assign presence to 1 and not present to 0
phyto.spp.id[phyto.spp.id==0] <- NA
phyto.spp.id[phyto.spp.id>=1] <- 1
write.csv(phyto.spp.id, "tables/ALU_DAT_R_phyto-annual-spp-occurence.csv", row.names=FALSE)

#filter whole reservoir fert years
phyto.fert.yrs <- filter(phyto, date >= "2000-01-01")

#filter baseline year
phyto.1998 <- filter(phyto, date < "1999-01-01")

### GENERAL STATS

#get lt stats by basin and study era
phyto.fert.yrs.lt <- mutate(phyto.fert.yrs, era = "Whole Reservoir Nutrient Additions")
phyto.fert.yrs.lt$era <- factor(phyto.fert.yrs.lt$era)

phyto.1998.lt <- mutate(phyto.1998, era = "Baseline")
phyto.1998.lt$era <- factor(phyto.1998.lt$era)

phyto.lt.stats <- rbind(phyto.1998.lt, phyto.fert.yrs.lt)
phyto.lt.stats <- summaryBy(cells.ml + mm3.l ~ date + basin + era, 
                                  data=phyto.lt.stats, FUN=f.total, na.rm=TRUE)
phyto.lt.stats.basin <- summaryBy(cells.ml.total + mm3.l.total ~ basin + era, 
                                  data=phyto.lt.stats, FUN=f.sum.stat, na.rm=TRUE)
phyto.lt.stats <- summaryBy(cells.ml.total + mm3.l.total ~ era, 
                            data=phyto.lt.stats, FUN=f.sum.stat, na.rm=TRUE)
phyto.lt.stats <- mutate(phyto.lt.stats, basin = "Combined")
phyto.lt.stats$basin <- factor(phyto.lt.stats$basin)
phyto.lt.stats <- rbind(phyto.lt.stats.basin, phyto.lt.stats)
write.csv(phyto.lt.stats, "tables/ALU_DAT_R_phyto-lt-stats.csv", row.names=FALSE)

#get lt stats by class and study era
phyto.lt.class.stats <- rbind(phyto.1998.lt, phyto.fert.yrs.lt)
phyto.lt.class.stats <- summaryBy(cells.ml + mm3.l ~ date + basin + era + class.name + class.alias, data=phyto.lt.class.stats, FUN=f.total, na.rm=TRUE)
phyto.lt.class.stats <- summaryBy(cells.ml.total + mm3.l.total ~ era + class.name + class.alias, data=phyto.lt.class.stats, FUN=f.sum.stat, na.rm=TRUE)
phyto.lt.class.stats <- arrange(phyto.lt.class.stats, era, desc(cells.ml.total.mean))
write.csv(phyto.lt.class.stats, "tables/ALU_DAT_R_phyto-lt-class-stats.csv", row.names=FALSE)
  
#filter to standard sample months apr-nov
phyto.f <- filter(phyto, month == 'Apr' | month == 'May' | month == 'Jun'| month == 'Jul'| month == 'Aug'| month == 'Sep'| month == 'Oct'| month == 'Nov')

#get sample totals for boxplot
phyto.f.t1 <- summaryBy(cells.ml + mm3.l ~ date + basin + month + year, data = phyto.f, FUN=f.total, na.rm=TRUE)

#get total class values by basin, month and year for histogram
phyto.f.t <- summaryBy(cells.ml + mm3.l ~ date + basin + month + year + class.name, data = phyto.f, FUN=f.total, na.rm=TRUE)
phyto.f.t <- summaryBy(cells.ml.total + mm3.l.total ~ month + year + class.name + basin, data = phyto.f.t, FUN=f.total, na.rm=TRUE)

#get annual stats for NB, SB and combined values
phyto.f.stats.annual <- summaryBy(cells.ml + mm3.l ~ date + basin + year, data = phyto.f, FUN = f.total, na.rm=TRUE)
phyto.f.stats.annual <- summaryBy(cells.ml.total + mm3.l.total ~ year, data = phyto.f.stats.annual, FUN = f.sum.stat, na.rm=TRUE)
phyto.f.stats.annual <- mutate(phyto.f.stats.annual, basin = "Combined")
phyto.f.stats.annual$basin <- factor(phyto.f.stats.annual$basin)
phyto.f.stats <- summaryBy(cells.ml + mm3.l ~ date + +basin + year, data = phyto.f, FUN = f.total, na.rm=TRUE)
phyto.f.stats <- summaryBy(cells.ml.total + mm3.l.total ~ year + basin, data = phyto.f.stats, FUN = f.sum.stat, na.rm=TRUE)
phyto.f.stats <- rbind(phyto.f.stats, phyto.f.stats.annual)
write.csv(phyto.f.stats, "tables/ALU_DAT_R_phyto-annual-stats.csv", row.names=FALSE)


### EDIBILITY

### ANNUAL EDIBILITY STATS

#filter to standard sample months apr-nov, remove 2016 for 5-yr review
phyto.ei.stats <- filter(phyto, month == 'Apr' | month == 'May' | month == 'Jun'| month == 'Jul'| month == 'Aug'| month == 'Sep'| month == 'Oct'| month == 'Nov')
phyto.ei.stats <- summaryBy(cells.ml + mm3.l ~ date + basin + year + edibility, data = phyto.ei.stats, FUN = f.total, na.rm=TRUE)
phyto.ei.stats <- summaryBy(cells.ml.total + mm3.l.total ~ edibility + year, data = phyto.ei.stats, FUN = f.sum.stat, na.rm=TRUE)
phyto.ei.stats <- filter(phyto.ei.stats, edibility == "E" | edibility == "I"| edibility == "B")
write.csv(phyto.ei.stats, "tables/ALU_DAT_R_phyto-ed-annual-stats.csv", row.names=FALSE)


### MONTHLY EDIBILITY STATS

#get monthly totals for all years
phyto.ei.t <- summaryBy(cells.ml + mm3.l ~ date + basin + month + year + edibility, data = phyto.f, FUN = f.total, na.rm=TRUE)
phyto.ei.t <- summaryBy(cells.ml.total + mm3.l.total ~ edibility + month + year + basin, data = phyto.ei.t, FUN = f.total, na.rm=TRUE)
phyto.ei.t <- filter(phyto.ei.t, edibility == "E" | edibility == "I"| edibility == "B")

#get monthly stats for fert era
phyto.ei.fert.monthly <- summaryBy(cells.ml + mm3.l ~ date + basin + edibility + month, data = phyto.fert.yrs, FUN = f.total, na.rm=TRUE)
phyto.ei.fert.monthly <- summaryBy(cells.ml.total + mm3.l.total ~ edibility + month, data = phyto.ei.fert.monthly, FUN = f.sum.stat, na.rm=TRUE)
phyto.ei.fert.monthly <- filter(phyto.ei.fert.monthly, edibility == "E" | edibility == "I"| edibility == "B")
phyto.ei.fert.monthly <- filter(phyto.ei.fert.monthly, month == 'Apr' | month == 'May' | month == 'Jun'| month == 'Jul'| month == 'Aug'| month == 'Sep'| month == 'Oct'| month == 'Nov')
phyto.ei.fert.monthly <- mutate(phyto.ei.fert.monthly, era = "Whole Reservoir Nutrient Additions")
phyto.ei.fert.monthly$era <- factor(phyto.ei.fert.monthly$era)

#get monthly stats for baseline era
phyto.ei.1998.monthly <- summaryBy(cells.ml + mm3.l ~ date + basin + edibility + month, data = phyto.1998, FUN = f.total, na.rm=TRUE)
phyto.ei.1998.monthly <- summaryBy(cells.ml.total + mm3.l.total ~ edibility + month, data = phyto.ei.1998.monthly, FUN = f.sum.stat, na.rm=TRUE)
phyto.ei.1998.monthly <- mutate(phyto.ei.1998.monthly, era = "Baseline")
phyto.ei.1998.monthly$era <- factor(phyto.ei.1998.monthly$era)

#combine fert and baseline stats into single dataframe
phyto.ei.monthly <- rbind(phyto.ei.1998.monthly, phyto.ei.fert.monthly)


### EDIBILITY STATS BY BASIN
phyto.ei.fert.basin <- summaryBy(cells.ml + mm3.l ~ date + basin + edibility, data = phyto.fert.yrs, FUN = f.total, na.rm=TRUE)
phyto.ei.fert.basin <- summaryBy(cells.ml.total + mm3.l.total ~ edibility + basin, data = phyto.ei.fert.basin, FUN = f.sum.stat, na.rm=TRUE)
phyto.ei.1998.basin <- summaryBy(cells.ml + mm3.l ~ date + basin + edibility, data = phyto.1998, FUN = f.total, na.rm=TRUE)
phyto.ei.1998.basin <- summaryBy(cells.ml.total + mm3.l.total ~ edibility + basin, data = phyto.ei.1998.basin , FUN = f.sum.stat, na.rm=TRUE)


### FIGURE: SEASONAL PHYTO ABUNDANCE
#only 2024
phyto.f.t.24 <- filter(phyto.f.t, year == "2024")

f.phyto.abun.24 <- ggplot(phyto.f.t.24, aes(x=month, y=cells.ml.total.total, group=class.name, fill=class.name))
f.phyto.abun.24 <- f.phyto.abun.24 + geom_bar(stat="identity") 
f.phyto.abun.24 <- f.phyto.abun.24 + scale_fill_manual(values=c("#d73027","#fdae61","#ffffbf","#74add1","#4575b4"))
#f.phyto.abun.24 <- f.phyto.abun.24 + scale_fill_brewer(palette="Spectral")
f.phyto.abun.24 <- f.phyto.abun.24 + facet_grid(basin ~ year)
f.phyto.abun.24 <- f.phyto.abun.24 + scale_y_continuous(name=bquote(paste(bold("Abundance ("*cells%.%mL^-1*")"))))
f.phyto.abun.24 <- f.phyto.abun.24 + xlab("Month")
f.phyto.abun.24 <- f.phyto.abun.24 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                     axis.title.x = element_text(face = "bold"),
                                     axis.title.y = element_text(face = "bold"),
                                     legend.title = element_blank(),
                                     legend.position = "top")
f.phyto.abun.24

ggsave("plots/alu_phyto-abun-seasonal_2024.png",width=8,height=5,units="in", dpi=600)

#all years
f.phyto.abun <- ggplot(phyto.f.t, aes(x=month, y=cells.ml.total.total, group=class.name, fill=class.name))
f.phyto.abun <- f.phyto.abun + geom_bar(stat="identity")
f.phyto.abun <- f.phyto.abun + scale_fill_manual(values=c("#d73027","#fdae61","#ffffbf","#74add1","#4575b4"))
#f.phyto.abun <- f.phyto.abun + scale_fill_brewer(palette="Spectral")
f.phyto.abun <- f.phyto.abun + facet_grid(basin ~ year)
f.phyto.abun <- f.phyto.abun + scale_y_continuous(name=bquote(paste(bold("Abundance ("*cells%.%mL^-1*")"))), labels=comma, 
                                                  limits=c(0,210000), breaks=seq(from=0,to=210000,by=30000))
f.phyto.abun <- f.phyto.abun + xlab("Month")
f.phyto.abun <- f.phyto.abun + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                             axis.title.x = element_text(face = "bold"),
                             axis.title.y = element_text(face = "bold"),
                             legend.title = element_blank(),
                             legend.position = "top")
f.phyto.abun

ggsave("plots/alu_phyto-abun-seasonal.png",width=15.5,height=8,units="in", dpi=600)


### FIGURE: SEASONAL PHYTO BIOVOLUME
#only 2024
f.phyto.biov.24 <- ggplot(phyto.f.t.24, aes(x=month, y=mm3.l.total.total, group=class.name, fill=class.name))
f.phyto.biov.24 <- f.phyto.biov.24 + geom_bar(stat="identity") 
f.phyto.biov.24 <- f.phyto.biov.24 + scale_fill_manual(values=c("#d73027","#fdae61","#ffffbf","#74add1","#4575b4"))
#f.phyto.biov.24 <- f.phyto.biov.24 + scale_fill_brewer(palette="Spectral")
f.phyto.biov.24 <- f.phyto.biov.24 + facet_grid(basin ~ year)
f.phyto.biov.24 <- f.phyto.biov.24 + scale_y_continuous(bquote(paste(bold('Biovolume ('*mm^3%.%L^-1*')'))), limits=c(0,1.4), breaks=seq(from=0, to=1.4, by=.2))
f.phyto.biov.24 <- f.phyto.biov.24 + xlab("Month")
f.phyto.biov.24 <- f.phyto.biov.24 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                     axis.title.x = element_text(face = "bold"),
                                     axis.title.y = element_text(face = "bold"),
                                     legend.title = element_blank(),
                                     legend.position = "top")
f.phyto.biov.24

ggsave("plots/alu_phyto-biov-seasonal_2024.png",width=8,height=5,units="in", dpi=600)

#all years
f.phyto.biov <- ggplot(phyto.f.t, aes(x=month, y=mm3.l.total.total, group=class.name, fill=class.name))
f.phyto.biov <- f.phyto.biov + geom_bar(stat="identity") 
f.phyto.biov <- f.phyto.biov + scale_fill_manual(values=c("#d73027","#fdae61","#ffffbf","#74add1","#4575b4"))
#f.phyto.biov <- f.phyto.biov + scale_fill_brewer(palette="Spectral")
f.phyto.biov <- f.phyto.biov + facet_grid(basin ~ year)
f.phyto.biov <- f.phyto.biov + scale_y_continuous(bquote(paste(bold('Biovolume ( '*mm^3%.%L^-1*')'))), breaks=seq(from=0, to=7, by=1))
f.phyto.biov <- f.phyto.biov + xlab("Month")
f.phyto.biov <- f.phyto.biov + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                     axis.title.x = element_text(face = "bold"),
                                     axis.title.y = element_text(face = "bold"),
                                     legend.title = element_blank(),
                                     legend.position = "top")
f.phyto.biov

ggsave("plots/alu_phyto-biov-seasonal.png",width=15.5,height=8,units="in", dpi=600)


### FIGURE: BOXPLOT ANNUAL PHYTO ABUNDANCE
#all years
f.phyto.abun.yrx <- ggplot(phyto.f.t1, aes(x=year, y=cells.ml.total, colour=basin))
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + scale_colour_brewer(palette="Set1")
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), labels=comma, breaks=seq(from=0, to=150000, by=25000))
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + xlab("Year")
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                               axis.title.x = element_text(face = "bold"),
                                               axis.title.y = element_text(face = "bold"),
                                               legend.title = element_blank(),
                                               legend.position = "top")
f.phyto.abun.yrx 

ggsave("plots/alu_phyto-abun-annual-bx.png",width=8,height=5,units="in", dpi=600)

#zoomed in figure - annual phyto abundance

f.phyto.abun.yrx <- ggplot(phyto.f.t1, aes(x=year, y=cells.ml.total, colour=basin))
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + scale_colour_brewer(palette="Set1")
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), labels=comma)
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + coord_cartesian(ylim= c(0,25000))
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + xlab("Year")
f.phyto.abun.yrx  <- f.phyto.abun.yrx  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                               axis.title.x = element_text(face = "bold"),
                                               axis.title.y = element_text(face = "bold"),
                                               legend.title = element_blank(),
                                               legend.position = "top")
f.phyto.abun.yrx 

ggsave("plots/alu_phyto-abun-annual-bxzm.png",width=8,height=5,units="in", dpi=600)

#2024 only
phyto.f.t1.24 <- filter(phyto.f.t1, year=="2024")

f.phyto.abun.yrx.24 <- ggplot(phyto.f.t1.24, aes(x=year, y=cells.ml.total, colour=basin))
f.phyto.abun.yrx.24  <- f.phyto.abun.yrx.24  + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.abun.yrx.24  <- f.phyto.abun.yrx.24  + scale_colour_brewer(palette="Set1")
#f.phyto.abun.yrx.24  <- f.phyto.abun.yrx.24  + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), labels=comma, breaks=seq(from=0, to=4000, by=1000))
f.phyto.abun.yrx.24  <- f.phyto.abun.yrx.24  + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))))
f.phyto.abun.yrx.24  <- f.phyto.abun.yrx.24  + xlab("Year")
f.phyto.abun.yrx.24  <- f.phyto.abun.yrx.24  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                               axis.title.x = element_text(face = "bold"),
                                               axis.title.y = element_text(face = "bold"),
                                               legend.title = element_blank(),
                                               legend.position = "top")
f.phyto.abun.yrx.24 

ggsave("plots/alu_phyto-abun-annual-bx.24.png",width=8,height=5,units="in", dpi=600)


### FIGURE: BOXPLOT ANNUAL PHYTO BIOVOLUME

f.phyto.biov.yrx <- ggplot(phyto.f.t1, aes(x=year, y=mm3.l.total, colour=basin))
f.phyto.biov.yrx  <- f.phyto.biov.yrx  + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.biov.yrx  <- f.phyto.biov.yrx  + scale_colour_brewer(palette="Set1")
f.phyto.biov.yrx  <- f.phyto.biov.yrx  + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), 
                                                            labels=comma, breaks=seq(from=0, to=7, by=1))
f.phyto.biov.yrx  <- f.phyto.biov.yrx  + xlab("Year")
f.phyto.biov.yrx  <- f.phyto.biov.yrx  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                   axis.title.x = element_text(face = "bold"),
                                                   axis.title.y = element_text(face = "bold"),
                                                   legend.title = element_blank(),
                                                   legend.position = "top")
f.phyto.biov.yrx 

ggsave("plots/alu_phyto-biov-annual-bx.png",width=8,height=5,units="in", dpi=600)


#zoomed in figure - annual phyto biovolume

f.phyto.biov.yrxz <- ggplot(phyto.f.t1, aes(x=year, y=mm3.l.total, colour=basin))
f.phyto.biov.yrxz  <- f.phyto.biov.yrxz  + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.biov.yrxz  <- f.phyto.biov.yrxz  + scale_colour_brewer(palette="Set1")
f.phyto.biov.yrxz  <- f.phyto.biov.yrxz  + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), 
                                                              labels=comma, breaks=seq(from=0, to=2, by=0.5))
f.phyto.biov.yrxz  <- f.phyto.biov.yrxz  + coord_cartesian(ylim= c(0,2))
f.phyto.biov.yrxz  <- f.phyto.biov.yrxz  + xlab("Year")
f.phyto.biov.yrxz  <- f.phyto.biov.yrxz  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                               axis.title.x = element_text(face = "bold"),
                                               axis.title.y = element_text(face = "bold"),
                                               legend.title = element_blank(),
                                               legend.position = "top")
f.phyto.biov.yrxz 

ggsave("plots/alu_phyto-biov-annual-bxzm.png",width=8,height=5,units="in", dpi=600)

#2024 only
phyto.f.t1.24 <- filter(phyto.f.t1, year=="2024")

f.phyto.biov.yrx.24 <- ggplot(phyto.f.t1.24, aes(x=year, y=mm3.l.total, colour=basin))
f.phyto.biov.yrx.24  <- f.phyto.biov.yrx.24  + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.biov.yrx.24  <- f.phyto.biov.yrx.24  + scale_colour_brewer(palette="Set1")
f.phyto.biov.yrx.24  <- f.phyto.biov.yrx.24  + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), 
                                                            labels=comma)
f.phyto.biov.yrx.24  <- f.phyto.biov.yrx.24  + xlab("Year")
f.phyto.biov.yrx.24  <- f.phyto.biov.yrx.24  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                               axis.title.x = element_text(face = "bold"),
                                               axis.title.y = element_text(face = "bold"),
                                               legend.title = element_blank(),
                                               legend.position = "top")
f.phyto.biov.yrx.24 

ggsave("plots/alu_phyto-biov-annual-bx.24.png",width=8,height=5,units="in", dpi=600)


### FIGURE: BOXPLOT ANNUAL EDIBILITY - PHYTO ABUNDANCE
#all years
f.phyto.ei.abun.yrx <- ggplot(phyto.ei.t, aes(x=year, y=cells.ml.total.total, fill=edibility))
f.phyto.ei.abun.yrx <- f.phyto.ei.abun.yrx + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.ei.abun.yrx <- f.phyto.ei.abun.yrx + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.abun.yrx <- f.phyto.ei.abun.yrx + scale_fill_brewer(palette="Pastel1")
#f.phyto.ei.abun.yrx <- f.phyto.ei.abun.yrx + scale_fill_manual(values=c('#b3cde3', '#fbb4ae', '#ccebc5'))
f.phyto.ei.abun.yrx <- f.phyto.ei.abun.yrx + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), 
                                                                labels=comma, breaks=seq(from=0, to=150000, by=25000))
f.phyto.ei.abun.yrx <- f.phyto.ei.abun.yrx + xlab("Year")
f.phyto.ei.abun.yrx <- f.phyto.ei.abun.yrx + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                     axis.title.x = element_text(face = "bold"),
                                     axis.title.y = element_text(face = "bold"),
                                     legend.title = element_blank(),
                                     legend.position = "top")
f.phyto.ei.abun.yrx

ggsave("plots/alu_phyto-ed-abun-annual-bx.png",width=8,height=5,units="in", dpi=600)


#zoomed in figure

f.phyto.ei.abun.yrxz <- ggplot(phyto.ei.t, aes(x=year, y=cells.ml.total.total, fill=edibility))
f.phyto.ei.abun.yrxz <- f.phyto.ei.abun.yrxz + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.ei.abun.yrxz <- f.phyto.ei.abun.yrxz + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.abun.yrxz <- f.phyto.ei.abun.yrxz + scale_fill_brewer(palette="Pastel1")
#f.phyto.ei.abun.yrxz <- f.phyto.ei.abun.yrxz + scale_fill_manual(values=c('#b3cde3', '#fbb4ae', '#ccebc5'))
f.phyto.ei.abun.yrxz <- f.phyto.ei.abun.yrxz + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), 
                                                                  labels=comma)
f.phyto.ei.abun.yrxz <- f.phyto.ei.abun.yrxz + coord_cartesian(ylim= c(0,25000))
f.phyto.ei.abun.yrxz <- f.phyto.ei.abun.yrxz + xlab("Year")
f.phyto.ei.abun.yrxz <- f.phyto.ei.abun.yrxz + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                   axis.title.x = element_text(face = "bold"),
                                                   axis.title.y = element_text(face = "bold"),
                                                   legend.title = element_blank(),
                                                   legend.position = "top")
f.phyto.ei.abun.yrxz

ggsave("plots/alu_phyto-ed-abun-annual-bxzm.png",width=8,height=5,units="in", dpi=600)

#2024 only 
phyto.ei.t.24 <- filter(phyto.ei.t, year=="2024")

f.phyto.ei.abun.yrx.24 <- ggplot(phyto.ei.t.24, aes(x=year, y=cells.ml.total.total, fill=edibility))
f.phyto.ei.abun.yrx.24 <- f.phyto.ei.abun.yrx.24 + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.ei.abun.yrx.24 <- f.phyto.ei.abun.yrx.24 + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.abun.yrx.24 <- f.phyto.ei.abun.yrx.24 + scale_fill_brewer(palette="Pastel1")
#f.phyto.ei.abun.yrx.24 <- f.phyto.ei.abun.yrx.24 + scale_fill_manual(values=c('#b3cde3', '#fbb4ae', '#ccebc5'))
f.phyto.ei.abun.yrx.24 <- f.phyto.ei.abun.yrx.24 + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), 
                                                                labels=comma)
f.phyto.ei.abun.yrx.24 <- f.phyto.ei.abun.yrx.24 + xlab("Year")
f.phyto.ei.abun.yrx.24 <- f.phyto.ei.abun.yrx.24 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                   axis.title.x = element_text(face = "bold"),
                                                   axis.title.y = element_text(face = "bold"),
                                                   legend.title = element_blank(),
                                                   legend.position = "top")
f.phyto.ei.abun.yrx.24

ggsave("plots/alu_phyto-ed-abun-annual-bx.24.png",width=8,height=5,units="in", dpi=600)

### FIGURE: BOXPLOT ANNUAL EDIBILITY - PHYTO BIOVOLUME
#all years
f.phyto.ei.biov.yrx <- ggplot(phyto.ei.t, aes(x=year, y=mm3.l.total.total, fill=edibility))
f.phyto.ei.biov.yrx <- f.phyto.ei.biov.yrx + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.ei.biov.yrx <- f.phyto.ei.biov.yrx + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.biov.yrx <- f.phyto.ei.biov.yrx + scale_fill_brewer(palette="Pastel1")
#f.phyto.ei.biov.yrx <- f.phyto.ei.biov.yrx + scale_fill_manual(values=c('#b3cde3', '#fbb4ae', '#ccebc5'))
f.phyto.ei.biov.yrx <- f.phyto.ei.biov.yrx + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), 
                                                                labels=comma, breaks=seq(from=0, to=4, by=0.5))
f.phyto.ei.biov.yrx <- f.phyto.ei.biov.yrx + xlab("Year")
f.phyto.ei.biov.yrx <- f.phyto.ei.biov.yrx + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                   axis.title.x = element_text(face = "bold"),
                                                   axis.title.y = element_text(face = "bold"),
                                                   legend.title = element_blank(),
                                                   legend.position = "top")
f.phyto.ei.biov.yrx

ggsave("plots/alu_phyto-ed-biov-annual-bx.png",width=8,height=5,units="in", dpi=600)

#zoomed in figure

f.phyto.ei.biov.yrxz <- ggplot(phyto.ei.t, aes(x=year, y=mm3.l.total.total, fill=edibility))
f.phyto.ei.biov.yrxz <- f.phyto.ei.biov.yrxz + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.ei.biov.yrxz <- f.phyto.ei.biov.yrxz + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.biov.yrxz <- f.phyto.ei.biov.yrxz + scale_fill_brewer(palette="Pastel1")
#f.phyto.ei.biov.yrxz <- f.phyto.ei.biov.yrxz + scale_fill_manual(values=c('#b3cde3', '#fbb4ae', '#ccebc5'))
f.phyto.ei.biov.yrxz <- f.phyto.ei.biov.yrxz + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), 
                                                                  labels=comma, breaks=seq(from=0, to=4, by=0.5))
f.phyto.ei.biov.yrxz <- f.phyto.ei.biov.yrxz + coord_cartesian(ylim= c(0,2))
f.phyto.ei.biov.yrxz <- f.phyto.ei.biov.yrxz + xlab("Year")
f.phyto.ei.biov.yrxz <- f.phyto.ei.biov.yrxz + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                   axis.title.x = element_text(face = "bold"),
                                                   axis.title.y = element_text(face = "bold"),
                                                   legend.title = element_blank(),
                                                   legend.position = "top")
f.phyto.ei.biov.yrxz

ggsave("plots/alu_phyto-ed-biov-annual-bxzm.png",width=8,height=5,units="in", dpi=600)

#2024 only
f.phyto.ei.biov.yrx.24 <- ggplot(phyto.ei.t.24, aes(x=year, y=mm3.l.total.total, fill=edibility))
f.phyto.ei.biov.yrx.24 <- f.phyto.ei.biov.yrx.24 + geom_boxplot(position=position_dodge(width=1), width = 0.75)
f.phyto.ei.biov.yrx.24 <- f.phyto.ei.biov.yrx.24 + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.biov.yrx.24 <- f.phyto.ei.biov.yrx.24 + scale_fill_brewer(palette="Pastel1")
#f.phyto.ei.biov.yrx.24 <- f.phyto.ei.biov.yrx.24 + scale_fill_manual(values=c('#b3cde3', '#fbb4ae', '#ccebc5'))
f.phyto.ei.biov.yrx.24 <- f.phyto.ei.biov.yrx.24 + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), 
                                                                labels=comma, breaks=seq(from=0, to=6, by=0.5))
f.phyto.ei.biov.yrx.24 <- f.phyto.ei.biov.yrx.24 + xlab("Year")
f.phyto.ei.biov.yrx.24 <- f.phyto.ei.biov.yrx.24 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                   axis.title.x = element_text(face = "bold"),
                                                   axis.title.y = element_text(face = "bold"),
                                                   legend.title = element_blank(),
                                                   legend.position = "top")
f.phyto.ei.biov.yrx.24

ggsave("plots/alu_phyto-ed-biov-annual-bx.24.png",width=8,height=5,units="in", dpi=600)

### FIGURE: MONTHLY EDIBILITY - PHYTO ABUNDANCE

f.phyto.ei.abun <- ggplot(phyto.ei.monthly, aes(x=month, y=cells.ml.total.mean, group=edibility, fill=edibility))
f.phyto.ei.abun <- f.phyto.ei.abun + geom_bar(stat="identity") 
f.phyto.ei.abun <- f.phyto.ei.abun + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.abun <- f.phyto.ei.abun + scale_fill_brewer(palette="Set1")
#f.phyto.ei.abun <- f.phyto.ei.abun + scale_fill_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))
f.phyto.ei.abun <- f.phyto.ei.abun + facet_wrap(~ era, ncol=2)
f.phyto.ei.abun <- f.phyto.ei.abun + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), 
                                                        labels=comma, limits=c(0,25000), breaks=seq(from=0, to=25000, by=5000))
f.phyto.ei.abun <- f.phyto.ei.abun + scale_x_discrete(name="Month")
f.phyto.ei.abun <- f.phyto.ei.abun + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                           axis.title.x = element_text(face = "bold"),
                                           axis.title.y = element_text(face = "bold"),
                                           legend.title = element_blank(),
                                           legend.position = "top")
f.phyto.ei.abun

ggsave("plots/alu_phyto-ed-abun-month.png",width=8,height=5,units="in", dpi=600)


### FIGURE: MONTHLY EDIBILITY - PHYTO BIOVOLUME

f.phyto.ei.biov <- ggplot(phyto.ei.monthly, aes(x=month, y=mm3.l.total.mean, group=edibility, fill=edibility))
f.phyto.ei.biov <- f.phyto.ei.biov + geom_bar(stat="identity") 
f.phyto.ei.biov <- f.phyto.ei.biov + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.biov <- f.phyto.ei.biov + scale_fill_brewer(palette="Set1")
#f.phyto.ei.biov <- f.phyto.ei.biov + scale_fill_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))
f.phyto.ei.biov <- f.phyto.ei.biov + facet_wrap(~ era, ncol=2)
f.phyto.ei.biov <- f.phyto.ei.biov + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), 
                                                        limits=c(0,1), breaks=seq(from=0, to=1, by=0.2))
f.phyto.ei.biov <- f.phyto.ei.biov + xlab("Month")
f.phyto.ei.biov <- f.phyto.ei.biov + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                           axis.title.x = element_text(face = "bold"),
                                           axis.title.y = element_text(face = "bold"),
                                           legend.title = element_blank(),
                                           legend.position = "top")
f.phyto.ei.biov

ggsave("plots/alu_phyto-ed-biov-month.png",width=8,height=5,units="in", dpi=600)


### FIGURE: SEASONAL EDIBILITY ABUNDANCE
#2024 only 

phyto.ei.t.24 <- filter(phyto.ei.t, year == "2024")

f.phyto.ei.s.abun.24 <- ggplot(phyto.ei.t.24, aes(x=month, y=cells.ml.total.total, group=edibility, fill=edibility))
f.phyto.ei.s.abun.24 <- f.phyto.ei.s.abun.24 + geom_bar(stat="identity") 
f.phyto.ei.s.abun.24 <- f.phyto.ei.s.abun.24 + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.s.abun.24 <- f.phyto.ei.s.abun.24 + scale_fill_brewer(palette="Set1")
#f.phyto.ei.s.abun.24 <- f.phyto.ei.s.abun.24 + scale_fill_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))
f.phyto.ei.s.abun.24 <- f.phyto.ei.s.abun.24 + facet_grid(basin ~ year)
f.phyto.ei.s.abun.24 <- f.phyto.ei.s.abun.24 + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), 
                                                                  labels=comma)
f.phyto.ei.s.abun.24 <- f.phyto.ei.s.abun.24 + xlab("Month")
f.phyto.ei.s.abun.24 <- f.phyto.ei.s.abun.24 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                               axis.title.x = element_text(face = "bold"),
                                               axis.title.y = element_text(face = "bold"),
                                               legend.title = element_blank(),
                                               legend.position = "top")
f.phyto.ei.s.abun.24 

ggsave("plots/alu_phyto-ed-abun-seasonal_2024.png",width=8,height=5,units="in", dpi=600)

#all years
f.phyto.ei.s.abun <- ggplot(phyto.ei.t, aes(x=month, y=cells.ml.total.total, group=edibility, fill=edibility))
f.phyto.ei.s.abun <- f.phyto.ei.s.abun + geom_bar(stat="identity") 
f.phyto.ei.s.abun <- f.phyto.ei.s.abun + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.s.abun <- f.phyto.ei.s.abun + scale_fill_brewer(palette="Set1")
#f.phyto.ei.s.abun <- f.phyto.ei.s.abun + scale_fill_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))
f.phyto.ei.s.abun <- f.phyto.ei.s.abun + facet_grid(basin ~ year)
f.phyto.ei.s.abun <- f.phyto.ei.s.abun + scale_y_continuous(name=bquote(paste(bold("Phytoplankton Abundance ("*cells%.%mL^-1*")"))), 
                                                            labels=comma, breaks=seq(from=0, to=150000, by=25000))
f.phyto.ei.s.abun <- f.phyto.ei.s.abun + xlab("Month")
f.phyto.ei.s.abun <- f.phyto.ei.s.abun + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                     axis.title.x = element_text(face = "bold"),
                                     axis.title.y = element_text(face = "bold"),
                                     legend.title = element_blank(),
                                     legend.position = "top")
f.phyto.ei.s.abun 

ggsave("plots/alu_phyto-ed-abun-seasonal.png",width=15.5,height=8,units="in", dpi=600)


### FIGURE: SEASONAL EDIBILITY BIOVOLUME 
#2024 only
f.phyto.ei.s.biov.24 <- ggplot(phyto.ei.t.24, aes(x=month, y=mm3.l.total.total, group=edibility, fill=edibility))
f.phyto.ei.s.biov.24 <- f.phyto.ei.s.biov.24 + geom_bar(stat="identity") 
f.phyto.ei.s.biov.24 <- f.phyto.ei.s.biov.24 + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.s.biov.24 <- f.phyto.ei.s.biov.24 + scale_fill_brewer(palette="Set1")
#f.phyto.ei.s.biov.24 <- f.phyto.ei.s.biov.24 + scale_fill_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))
f.phyto.ei.s.biov.24 <- f.phyto.ei.s.biov.24 + facet_grid(basin ~ year)
f.phyto.ei.s.biov.24 <- f.phyto.ei.s.biov.24 + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), limits=c(0,1.4), breaks=seq(from=0, to=1.4, by=.2))
f.phyto.ei.s.biov.24 <- f.phyto.ei.s.biov.24 + xlab("Month")
f.phyto.ei.s.biov.24 <- f.phyto.ei.s.biov.24 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                               axis.title.x = element_text(face = "bold"),
                                               axis.title.y = element_text(face = "bold"),
                                               legend.title = element_blank(),
                                               legend.position = "top")
f.phyto.ei.s.biov.24

ggsave("plots/alu_phyto-ed-biov-seasonal_2024.png",width=8,height=5,units="in", dpi=600)

#all years
f.phyto.ei.s.biov <- ggplot(phyto.ei.t, aes(x=month, y=mm3.l.total.total, group=edibility, fill=edibility))
f.phyto.ei.s.biov <- f.phyto.ei.s.biov + geom_bar(stat="identity") 
f.phyto.ei.s.biov <- f.phyto.ei.s.biov + scale_fill_manual(values=c("#fdae61","#d73027","#74add1"))
#f.phyto.ei.s.biov <- f.phyto.ei.s.biov + scale_fill_brewer(palette="Set1")
#f.phyto.ei.s.biov <- f.phyto.ei.s.biov + scale_fill_manual(values=c('#377eb8', '#e41a1c', '#4daf4a'))
f.phyto.ei.s.biov <- f.phyto.ei.s.biov + facet_grid(basin ~ year)
f.phyto.ei.s.biov <- f.phyto.ei.s.biov + scale_y_continuous(name=bquote(paste(bold('Phytoplankton Biovolume ( '*mm^3%.%L^-1*')'))), 
                                                            breaks=seq(from=0, to=7, by=1))
f.phyto.ei.s.biov <- f.phyto.ei.s.biov + xlab("Month")
f.phyto.ei.s.biov <- f.phyto.ei.s.biov + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                     axis.title.x = element_text(face = "bold"),
                                     axis.title.y = element_text(face = "bold"),
                                     legend.title = element_blank(),
                                     legend.position = "top")
f.phyto.ei.s.biov

ggsave("plots/alu_phyto-ed-biov-seasonal.png",width=20,height=8,units="in", dpi=600)

#####################################################
