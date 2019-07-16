source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library("ape")
library("ade4")
source("https://install-github.me/r-lib/rlang")
library("ggplot2")
library("grid")
library("gridExtra")
library("lme4")
library("vegan")
library("dplyr")
library("ecodist")
library("MASS")
library("devtools")
library("dplyr")
library("stats")
library("emmeans")
library("stringi")
library("phyloseq")
library("cowplot")
source("https://bioconductor.org/biocLite.R")
library('ggtree')
library('ggtree')
library('digest')

######## Figure 1, Panel A #################

data<-read.csv("./R_file1.csv")
data<-data[(data$Bin_Name=="2kb+"),]
data$Tree_Groups <- factor(data$Tree_Groups,levels = c("nutrient.rich", "pseudo-oligotrophic", "oligotrophic"))

Fig.1.A<-ggplot(data,aes(Tree_Groups,Percent_Indels_plus_SNPs,colour=Tree_Groups,shape=Tree_Groups))+
  geom_jitter(aes(size =3,stroke=2))+
  annotate("text", size = 10, x=0.5, y=1.0, label= "LMER:", hjust=0.0)+
  annotate("text", size = 10, x=0.5, y=0.90, label= paste("list(Marginal~R^{2}==0.32)"),hjust=0.0,parse=TRUE,size=5)+
  annotate("text", size = 10, x=0.5, y=0.80, label= "p = 0.0005",hjust=0.0)+
  scale_color_manual(name="Phylogenetic Groups",labels=c("High P Lake/ High P Genotype","High P Lake/ Low P Genotype","Low P Lake/ Low P Genotype"),values = c("nutrient.rich"="green","oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9")) +
  scale_shape_manual(name="Phylogenetic Groups",labels=c("High P Lake/ High P Genotype","High P Lake/ Low P Genotype","Low P Lake/ Low P Genotype"),values = c("nutrient.rich"=16,"oligotrophic"=17,"pseudo-oligotrophic"=15)) +
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 30),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=30),legend.text=element_text(size=30))+
  guides(size=FALSE, colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))+
  ylab("\n Polymorphic Sites in the Genome (%)\n")+
  xlab("\nPhylogenetic Groups")+
  geom_segment(aes(x = 0.70, y = 0.25, xend = 2.3, yend = 0.25),size=1,color="black")+
  geom_segment(aes(x = 2.65, y = 1.0, xend = 3.5, yend = 1.0),size=1,color="black")+
  annotate("text", x = 1.5, y = 0.3, label = "A",size=10)+
  annotate("text", x = 3.05, y = 1.05, label = "B",size=10)+
  scale_x_discrete(labels=c("nutrient.rich" = "HL/HG", "oligotrophic" = "LL/LG", "pseudo-oligotrophic" = "HL/LG"))
Fig.1.A

######## Figure 1, Panel B #################

Fig.1.B<-ggplot(low_data,aes(Percent_Indels_plus_SNPs,Completeness,shape=Year,fill=Percent_Below_Cutoff))+
  geom_point(aes(size =3,stroke=2),color='blue')+
  stat_smooth(aes(group='none'),method='lm',col='black',fullrange = TRUE,se=FALSE)+
  scale_shape_manual(guide='none',name="Year",values=c("2011"=24,"2013"=24))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 30),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=30),legend.text=element_text(size=30))+
  guides(size=FALSE,shape=FALSE,color=FALSE,fill=FALSE)+
  scale_fill_gradient(name="% Contigs below\n Coverage Cutoff",low="white",high="black")+
  xlab("\nPolymorphic Sites in the Genome (%)\n")+
  ylab("Genome Completeness (%)\n")
Fig.1.B

######## Figure 1, Panel C #################

data<-read.csv("./R_file2.csv")
data$z.scores.calculated <- ave(data$Coverage, data$Name, FUN=scale)
mapping<-read.csv("./R_file3.csv")
data<-dplyr::full_join(data,mapping,by='Name')
test_data<-data[!is.na(data),]
test_data<-data[!is.na(data$Tree_Groups),]
levels(test_data$Tree_Groups) <- c("HL/HG", "HL/LG", "LL/LG")

Fig.1.C<-ggplot(test_data, aes(x=z.scores.calculated,fill=Tree_Groups)) +
  geom_histogram(binwidth=0.01) +
  scale_x_continuous(name="Coverage\n (standardized scores per strain)", limits=c(-1.25, 1.25))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 24),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=30),legend.text=element_text(size=30))+
  ylab("\n# of Contigs\n")+
  scale_fill_manual(name="Phylogenetic Group",values=c("LL/LG"="blue","HL/HG"="green","HL/LG"="#56B4E9"))+
  guides(fill='none')+
  facet_grid(Tree_Groups~.)
Fig.1.C

######## Figure 1, Panel D #################

data<-read.csv("./R_file4.csv")
lowcov<-read.csv("./R_file5.csv")
treegroups<-read.csv("./R_file1.csv")
data<-dplyr::full_join(data,lowcov,by="Name")
testing <- aggregate(Percent_Ref_Combined_FR ~  Name, data, median)
data<-dplyr::full_join(treegroups,testing,by="Name")
cutoff<-read.csv("./R_file6.csv")
data_all<-dplyr::full_join(data,cutoff,by="Name")
data_all<-data_all[!is.na(data_all$Tree_Groups),]
data_all<-data_all[(data_all$Bin_Name=="2kb+"),]

Fig.1.D<-ggplot(data_all, aes(x=Tree_Groups,y=Percent_Ref_Combined_FR*100,fill=Percent_Below_Cutoff.y,shape=Tree_Groups,color=Tree_Groups))+
  geom_jitter(size=4,stroke=2,width=0.2)+
  scale_shape_manual(guide='none',name="Phylogenetic Groups",values=c("pseudo-oligotrophic"=22,"nutrient.rich"=21,"oligotrophic"=24))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 30),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=30),legend.text=element_text(size=30))+
  guides(size=FALSE,shape=FALSE,color=FALSE,fill=FALSE)+
  scale_fill_gradient(name="% Contigs below\n Coverage Cutoff",low="white",high="black")+
  xlab("\nPhylogenetic Groups")+
  scale_color_manual(name="Phylogenetic Group",labels=c("oligotrophic"="LN","pseudo-oligotrophic"="HN2","nutrient.rich"="HN1"),values=c("oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","nutrient.rich"="green"))+
  ylab("Median % of Reference Variant\n")+
  scale_x_discrete(limits=c("nutrient.rich","pseudo-oligotrophic","oligotrophic"),labels=c("nutrient.rich" = "HL/HG", "oligotrophic" = "LL/LG", "pseudo-oligotrophic" = "HL/LG"))
Fig.1.D


####### Figure 2 ###############

metadata<-read.csv("./R_file7.csv")
dd <- data.frame(Strain  = metadata$Strain, TreeGroups = metadata$Tree_Groups)
streamlining<-read.csv("./R_file8.csv")
tree <- read.tree("./Vincent_tree_with_bootstrapping.txt")
tree<-read.csv("./R_file9.csv")
tree$node.label
original<-c("","96","69","97","93","54","81","49","100","100","100","89","100","99","100","100","100","41","18","21","59","68","91","100","20","25","32","38","21","21","48","47","94","94","94","41","68","64","70","96","60","78","58","42" ,"100")
tree$node.label<-c("","96","69","97","93","54","81","","100","100","100","89","100","99","100","100","100","","","","59","68","91","100","","","","","","","","","94","94","94","","68","64","70","96","60","78","58","","100")

p<-ggtree(tree) + geom_tree() + theme_tree() + geom_tiplab(size=5,hjust = -0.25)
my_tree <- p %<+% dd + geom_nodelab(nudge_y=1,nudge_x = -0.0015,size = 6, col= "black")+geom_tippoint(aes(color=TreeGroups,alpha=0.5),size=4)+ theme(legend.position="right") + scale_color_manual(values = c("nutrient-rich"="green", "oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9"))+guides(color=FALSE,alpha=FALSE)+xlim(NA,0.07)+expand_limits(y = c(0, 48))+geom_treescale(x=0.001,y=-0.5,offset=-1)

d2 <- data.frame(ID=streamlining$ID, Genome.Size=streamlining$GenomeSizeMb, TreeGroup=streamlining$Tree_Groups)
d3 <- data.frame(ID=streamlining$ID, Coding.DNA=streamlining$DNA_coding_number_of_bases_..ofTotal._2kb., TreeGroup=streamlining$Tree_Groups)
d4 <- data.frame(ID=streamlining$ID, Paralogs=streamlining$Paralog_., TreeGroup=streamlining$Tree_Groups)
d5 <- data.frame(ID=streamlining$ID, GC.Content=streamlining$GC_Content, TreeGroup=streamlining$Tree_Groups)
d6 <- data.frame(ID=streamlining$ID, Sigma.factors=streamlining$`Sigma_factors_percent`, TreeGroup=streamlining$Tree_Groups)
d7 <- data.frame(ID=streamlining$ID, Completness=streamlining$Completeness, TreeGroup=streamlining$Tree_Groups)

d2$ID <- factor(d2$ID,levels=unique(d2$ID))
scaleFUN <- function(x) sprintf("%.2f", x)
d2_mean<-aggregate(Genome.Size~TreeGroup,d2,mean)
d2_sd<-aggregate(Genome.Size~TreeGroup,d2,sd)
d2_length<-aggregate(Genome.Size~TreeGroup,d2,length)
d2_new<-cbind(d2_mean,d2_sd,d2_length)
d2_new<-d2_new[,c(1,2,4,6)]
d2_new$se<-d2_new$Genome.Size.1/sqrt(d2_new$Genome.Size.2)
d2_new<-d2_new[,c(1,2,5)]

p2<-ggplot(d2,aes(Genome.Size,ID,color=TreeGroup,alpha=0.5))+ 
  scale_color_manual(values=c("green","blue","#56B4E9"))+geom_point(aes(),size=4)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5,size=18),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 16),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=16),legend.text=element_text(size=16))+ 
  guides(size=FALSE,color=FALSE,alpha=FALSE)+
  ggtitle('Genome Size')+
  geom_segment(aes(x = 4.976299941, y = 0, xend = 4.976299941, yend = 46),colour="springgreen4",linetype='dashed',size=1)+
  geom_segment(aes(x=5.081535,xend=5.081535,y=0,yend=46), colour="blue",linetype='dashed',size=1)+
  expand_limits(y = c(0, 48))+
  geom_segment(aes(x=4.998920727,xend=4.998920727,y=0,yend=46), colour="#56B4E9",linetype='dashed',size=1)
d2_new$order<-c(8.5,37,23.5)

p2_new<-ggplot(d2_new,aes(Genome.Size,order,color=TreeGroup,alpha=0.5))+ 
  scale_color_manual(values=c("green","blue","#56B4E9"))+geom_point(aes(),size=4)+
  geom_errorbarh(aes(xmin=Genome.Size-se, xmax=Genome.Size+se),height=.0,size=3)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5,size=18),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 16),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=16),legend.text=element_text(size=16))+ 
  guides(size=FALSE,color=FALSE,alpha=FALSE)+
  ggtitle('Genome Size')+
  expand_limits(y = c(0, 48))

scaleFUN <- function(x) sprintf("%.1f", x)
d3$ID <- factor(d3$ID,levels=unique(d3$ID))
d3_mean<-aggregate(Coding.DNA~TreeGroup,d3,mean)
d3_sd<-aggregate(Coding.DNA~TreeGroup,d3,sd)
d3_length<-aggregate(Coding.DNA~TreeGroup,d3,length)

d3_new<-cbind(d3_mean,d3_sd,d3_length)
d3_new<-d3_new[,c(1,2,4,6)]
d3_new$se<-d3_new$Coding.DNA.1/sqrt(d3_new$Coding.DNA.2)
d3_new<-d3_new[,c(1,2,5)]
d3_new$order<-c(8.5,37,23.5)

scaleFUN <- function(x) sprintf("%.1f", x)
p3_new <- ggplot(d3_new,aes(Coding.DNA,order,color=TreeGroup,alpha=0.5))+ scale_color_manual(values=c("green","blue","#56B4E9"))+ geom_point(aes(),size=4)+theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5,size=18),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 16),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=16),legend.text=element_text(size=16))+ guides(size=FALSE,color=FALSE,alpha=FALSE)+ggtitle('Coding DNA*')+
  geom_errorbarh(aes(xmin=Coding.DNA-se, xmax=Coding.DNA+se),height=.0,size=3)+
  expand_limits(y = c(0, 48))

d4$ID <- factor(d4$ID,levels=unique(d4$ID))
scaleFUN <- function(x) sprintf("%.1f", x)

d4_mean<-aggregate(Paralogs~TreeGroup,d4,mean)
d4_sd<-aggregate(Paralogs~TreeGroup,d4,sd)
d4_length<-aggregate(Paralogs~TreeGroup,d4,length)
d4_new<-cbind(d4_mean,d4_sd,d4_length)
d4_new<-d4_new[,c(1,2,4,6)]
d4_new$se<-d4_new$Paralogs.1/sqrt(d4_new$Paralogs.2)
d4_new<-d4_new[,c(1,2,5)]
d4_new$order<-c(8.5,37,23.5)

scaleFUN <- function(x) sprintf("%.1f", x)
p4_new <- ggplot(d4_new,aes(Paralogs,order,color=TreeGroup,alpha=0.5))+ scale_color_manual(values=c("green","blue","#56B4E9"))+ geom_point(aes(),size=4)+theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5,size=18),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 16),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=16),legend.text=element_text(size=16))+ guides(size=FALSE,color=FALSE,alpha=FALSE)+ggtitle('Paralogs*')+scale_x_continuous(labels=scaleFUN)+
  geom_errorbarh(aes(xmin=Paralogs-se, xmax=Paralogs+se),height=.0,size=3)+
  expand_limits(y = c(0,48))

scaleFUN <- function(x) sprintf("%.1f", x)

d5$ID <- factor(d5$ID,levels=unique(d5$ID))
d5_mean<-aggregate(GC.Content~TreeGroup,d5,mean)
d5_sd<-aggregate(GC.Content~TreeGroup,d5,sd)
d5_length<-aggregate(GC.Content~TreeGroup,d5,length)
d5_new<-cbind(d5_mean,d5_sd,d5_length)
d5_new<-d5_new[,c(1,2,4,6)]
d5_new$se<-d5_new$GC.Content.1/sqrt(d5_new$GC.Content.2)
d5_new<-d5_new[,c(1,2,5)]
d5_new$order<-c(8.5,37,23.5)

scaleFUN <- function(x) sprintf("%.1f", x)
p5_new <- ggplot(d5_new,aes(GC.Content,order,color=TreeGroup,alpha=0.5))+ scale_color_manual(values=c("green","blue","#56B4E9"))+ geom_point(aes(),size=4)+theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5,size = 18),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 16),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=16),legend.text=element_text(size=16))+ guides(size=FALSE,color=FALSE,alpha=FALSE)+ggtitle("GC Content**")+scale_x_continuous(labels=scaleFUN)+
  geom_errorbarh(aes(xmin=GC.Content-se, xmax=GC.Content+se),height=.0,size=3)+
  annotate("text",label="A",color="blue", size=5,x=42.48333333, y= 39, hjust=0.5)+
  annotate("text",label="A",color="#56B4E9", size=5,x=42.50909091, y= 25.5, hjust=0.5)+
  annotate("text",label="B",color="springgreen4", size=5,x=42.70647059, y= 10.5, hjust=0.5)+
  expand_limits(y = c(0, 48))

d6$ID <- factor(d6$ID,levels=unique(d6$ID))
scaleFUN <- function(x) sprintf("%.2f", x)
d6_mean<-aggregate(Sigma.factors~TreeGroup,d6,mean)
d6_sd<-aggregate(Sigma.factors~TreeGroup,d6,sd)
d6_length<-aggregate(Sigma.factors~TreeGroup,d6,length)
d6_new<-cbind(d6_mean,d6_sd,d6_length)
d6_new<-d6_new[,c(1,2,4,6)]
d6_new$se<-d6_new$Sigma.factors.1/sqrt(d6_new$Sigma.factors.2)
d6_new<-d6_new[,c(1,2,5)]
d6_new$order<-c(8.5,37,23.5)

scaleFUN <- function(x) sprintf("%.2f", x)
p6_new <- ggplot(d6_new,aes(Sigma.factors,order,color=TreeGroup,alpha=0.5))+ scale_color_manual(values=c("green","blue","#56B4E9"))+ geom_point(aes(),size=4)+theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5,size = 18),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 16),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=16),legend.text=element_text(size=16))+ guides(size=FALSE,color=FALSE,alpha=FALSE)+ggtitle("Sigma Factors**")+scale_x_continuous(labels=scaleFUN)+
  geom_errorbarh(aes(xmin=Sigma.factors-se, xmax=Sigma.factors+se),height=.0,size=3)+
  annotate("text",label="A",color="springgreen4", size=5,x=0.2656341, y= 10.5, hjust=0.5)+
  annotate("text",label="A,B",color="blue", size=5,x=0.2809411, y= 39.0, hjust=0.5)+
  annotate("text",label="B",color="#56B4E9", size=5,x=0.2897346, y= 25.5, hjust=0.5)+
  expand_limits(y = c(0, 48))

d7$ID <- factor(d7$ID,levels=unique(d7$ID))
d7_mean<-aggregate(Completness~TreeGroup,d7,mean)
d7_sd<-aggregate(Completness~TreeGroup,d7,sd)
d7_length<-aggregate(Completness~TreeGroup,d7,length)
d7_new<-cbind(d7_mean,d7_sd,d7_length)
d7_new<-d7_new[,c(1,2,4,6)]
d7_new$se<-d7_new$Completness.1/sqrt(d7_new$Completness.2)
d7_new<-d7_new[,c(1,2,5)]
d7_new$order<-c(8.5,37,23.5)

scaleFUN <- function(x) sprintf("%.3f", x)
p7_new <- ggplot(d7_new,aes(Completness,order,color=TreeGroup,alpha=0.5))+ scale_color_manual(values=c("green","blue","#56B4E9"))+ geom_point(aes(),size=4)+theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5,size = 18),axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 16),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=16),legend.text=element_text(size=16))+ guides(size=FALSE,color=FALSE,alpha=FALSE)+ggtitle("Completeness**")+
  annotate("text",label="A",color="blue", size=5,x=96.77222222, y= 39, hjust=0.5)+
  annotate("text",label="B",color="#56B4E9", size=5,x=98.91181818, y= 25.5, hjust=0.5)+
  annotate("text",label="C",color="springgreen4", size=5,x=99.26588235, y= 10.5, hjust=0.5)+
  geom_errorbarh(aes(xmin=Completness-se, xmax=Completness+se),height=.0,size=3)+
  expand_limits(y = c(0, 48))+
  xlim(96.25,100)

grid.arrange(my_tree,p2_new,p7_new,p3_new,p5_new,p4_new,p6_new,ncol=7,widths=c(6,1,1,1,1,1,1))

############# Figure 3, Panel B ##############

data<-read.csv("./R_file10.csv")
data_table<-table(data)
Microcystis_Annotation<-read.csv("./R_file11.csv")
rownames(Microcystis_Annotation) <- Microcystis_Annotation[,1]
Microcystis_Annotation[,1]<-NULL
Microcystis_Annotation_dm<-vegan::vegdist(Microcystis_Annotation,method="bray",diag=TRUE,upper=TRUE)
library(ape)
result<-pcoa(Microcystis_Annotation_dm)
pcoa_vectors<-result$vectors[,1:5]
pcoa_vectors<-as.data.frame(pcoa_vectors)
pcoa_vectors$Name<-rownames(pcoa_vectors)
rownames(pcoa_vectors)<-NULL
treegroups<-read.csv("./R_file13.csv")
data<-dplyr::full_join(pcoa_vectors,treegroups,by="Name")
Microcystis_Annotation_dm<-as.matrix(Microcystis_Annotation_dm)

Fig.3.B<-ggplot(data,aes(Axis.1,Axis.2,colour=Tree_Groups,shape=Tree_Groups))+
  geom_point(aes(size =3,stroke=2))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 30),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=30),legend.text=element_text(size=30))+
  scale_color_manual(guide='none',name="Phylogenetic Group",labels=c("oligotrophic"="Oligotrophic","pseudo-oligotrophic"="Pseudo-oligotrophic","nutrient.rich"="Nutrient-rich"),values=c("oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","nutrient.rich"="green"))+
  scale_size(guide='none')+
  scale_shape(guide='none',name="Trophic Status")+
  xlab("\nAxis 1 (30.5% of variance)")+
  ylab("Axis 2 (11.0% of variance)\n")+
  annotate("text", x=-0.05, y=0.037,label=paste("list(Adonis:~F['2,45']==9.01", ", p==0.001", ", R^{2}==0.295)"),parse=TRUE,hjust=0.0,size=10)+
  annotate("text", x=-0.04, y=0.02,label="LL/LG", color="blue",size=10)+
  annotate("text", x=0.02, y=0.02,label="HL/LG", color="#56B4E9",size=10)+
  annotate("text", x=0.06, y=-0.02,label="HL/HG", color="green",size=10)
Fig.3.B

############# Figure 3, Panel A ##############

Fig.3.A<-ggplot(data,aes(Axis.1,Axis.2,colour=Lake,shape=Tree_Groups,size=Trophic_status))+
  geom_point(aes(stroke=2))+
  scale_size_manual(guide='none',name="Trophic Status",values=c("oligotrophic"=5,"mesotrophic"=9,"eutrophic"=13))+
  scale_color_manual(values=c("red","gray","#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442","blue","green","violet","deeppink4"))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 30),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=30),legend.text=element_text(size=30))+
  xlab("\n Axis 1 (30.5% of variance)")+
  ylab("Axis 2 (11.0% of variance)\n")+
  scale_shape_manual(guide='none',name="Phylogenetic Group",values=c("oligotrophic"=18,"pseudo-oligotrophic"=18,"nutrient.rich"=18),labels=c("oligotrophic"="Low-Nutrient","pseudo-oligotrophic"="High-Nutrient Group 2","nutrient.rich"="High-Nutrient Group 1"))+
  guides(color=guide_legend(override.aes=list(shape=18,size=10)))+
  annotate("text", x=-0.05, y=0.037,label=paste("list(Adonis:~F['9,40']==8.56", ", p==0.001", ", R^{2}==0.311)"),parse=TRUE,hjust=0.0,size=10)
Fig.3.A

######### Figure 4, Panel A #########################

phyloseq.mothur.scaled<-readRDS("./R_file12.rds")
otu_table<-otu_table(phyloseq.mothur.scaled)
otu_table<-as.matrix(otu_table)
otu_table<-t(otu_table)
taxonomic_dm<-vegan::vegdist(otu_table,method="bray",diag=TRUE,upper=TRUE)
result<-pcoa(taxonomic_dm)
result.vectors<-result$vectors[,1:6]
result.vectors<-as.data.frame(result.vectors)
result.vectors$Name<-row.names(result.vectors)
sampledf.lake.reps <- data.frame(sample_data(phyloseq.mothur.scaled))
metadata<-read.csv("./R_file13.csv")
result.vectors$Name<-substring(result.vectors$Name, 2)
results<-dplyr::full_join(metadata,result.vectors,by="Name")
results<-results[(!is.na(results$Tree_Groups)),]

Fig.4.A<-ggplot(results,aes(Axis.1,Axis.3,colour=Tree_Groups,shape=Tree_Groups))+
  geom_point(aes(size =3,stroke=2))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 30),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=30),legend.text=element_text(size=30))+
  scale_color_manual(guide='none',name="Phylogenetic Group",labels=c("oligotrophic"="Low-Nutrient","pseudo-oligotrophic"="High-Nutrient Group 2","nutrient.rich"="High-Nutrient Group 1"),values=c("oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","nutrient.rich"="green"))+
  scale_shape(guide='none',name="Trophic Status")+
  #guides(guide='none',colour = guide_legend(override.aes = list(size=6)))+
  scale_size(guide='none')+
  xlab("\n Axis 1 (15.5% of variance)")+
  ylab("Axis 3 (8.3% of variance)\n")+
  annotate("text", x=-0.5, y=0.45,label=paste("list(Adonis:~F['2,49']==2.4", ", p==0.001", ", R^{2}==0.092)"),parse=TRUE,hjust =0.0,size=10)
Fig.4.A

################ Figure 4, Panel B ##############

bacterial_data<-read.csv("./R_file14.csv")
metadata<-read.csv("./R_file15.csv")
rownames(bacterial_data)<-bacterial_data$Name
bacterial_data$Name<-NULL
Bacterial_data_trophic_bray_dm_of_pfams<-vegan::vegdist(bacterial_data,method="bray",diag=TRUE,upper=TRUE)
Bacterial_data_trophic_bray_dm_of_pfams<-stepacross(Bacterial_data_trophic_bray_dm_of_pfams, path="shortest",toolong=0.75) 
Bacteria_Trophic_pcoa<-pcoa(Bacterial_data_trophic_bray_dm_of_pfams)
Bacteria_Trophic_pcoa_vectors<-Bacteria_Trophic_pcoa$vectors[,1:5]
Bacteria_Trophic_pcoa_vectors<-as.data.frame(Bacteria_Trophic_pcoa_vectors)
Bacteria_Trophic_pcoa_vectors$Name<-rownames(Bacteria_Trophic_pcoa_vectors)
Bacteria_Trophic_data<-dplyr::full_join(metadata,Bacteria_Trophic_pcoa_vectors,by="Name")
Bacteria_Trophic_data<-Bacteria_Trophic_data[!is.na(Bacteria_Trophic_data$Axis.1),]
Bacteria_Trophic_pcoa_vectors$Name<-NULL

Fig.4.B<-ggplot(Bacteria_Trophic_data,aes(Axis.1,Axis.2,colour=Tree_Groups,shape=Tree_Groups))+
  geom_point(aes(size =3,stroke=2))+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 30),axis.title=element_text(size=30,face="plain"),legend.title=element_text(size=30),legend.text=element_text(size=30))+
  scale_color_manual(guide='none',name="Phylogenetic Group",labels=c("oligotrophic"="Oligotrophic","pseudo-oligotrophic"="Pseudo-oligotrophic","nutrient-rich"="Nutrient-rich"),values=c("oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","nutrient-rich"="green"))+
  scale_shape(guide='none',name="Trophic Status")+
  scale_size(guide='none')+
  annotate("text", x=-1.2, y=0.30,label=paste("list(Adonis:~F['2,45']==0.46", ", p==0.70", ", R^{2}==0.02)"),parse=TRUE,hjust=0.0,size=10)+
  xlab("\n Axis 1 (85.9% of variance)")+
  ylab("Axis 2 (6.51% of variance)\n")
Fig.4.B

################ Figure 4, Panel C ##############

data<-read.csv("./R_file16.csv")
full<-lmer(Bacteria_pfam_bray_distance~Microcystis_pfam_bray_distance+(1|Batch),data=data)
full3<-lmer(Bacteria_pfam_bray_distance~Microcystis_pfam_bray_distance3+(1|Batch),data=data)

scaleFUN <- function(x) sprintf("%.2f", x)
Fig.4.C<-ggplot(data=data,aes(Microcystis_pfam_bray_distance,Bacteria_pfam_bray_distance,colour=Batch))+
  geom_point(aes(size =1,stroke=2))+
  geom_line(aes(y=fitted.values(full3), group=Batch),size=3) +
  geom_line(aes(y=fitted.values(full), group=Batch),size=3) +
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 40),axis.title=element_text(size=40,face="plain"),legend.title=element_text(size=40),legend.text=element_text(size=40))+
  scale_color_manual(guide='none',name="Lake Environment",values=c("black", "gray", "#F0E442"))+
  scale_size(guide='none')+
  scale_y_continuous(labels=scaleFUN) +
  xlab("\nFunctional divergence: Hosts")+
  ylab("Functional divergence: Host Phycospheres\n")+
  annotate("text", x=0.01, y=0.80,label="LMER:",size=10, hjust =0.0)+
  annotate("text", x=0.01, y=0.77,label=paste("list(Linear~Model:~Marginal~R^{2}==0.093~p==0.060)"),parse=TRUE,hjust=0.0,size=10)+
  annotate("text", x=0.01, y=0.74,label=paste("list(Polynomial~Model:~Marginal~R^{2}==0.14~p==0.065)"),parse=TRUE,hjust=0.0,size=10)+
  annotate("text", x=0.01, y=0.18,label="Pairwise comparisons within\n W13-13, W13-15, W13-16",hjust=0.0,size=10,colour="darkgray")+
  annotate("text", x=0.08, y=0.32,label="W13-13, W13-15,\n W13-16 vs. W13-18",hjust=0.0,size=10,colour="gray")+
  annotate("text", x=0.077, y=0.58,label="W13-13, W13-15, W13-16,\n W13-18 vs. W13-11",hjust=0.0,size=10,colour="gray")+
  annotate("text", x=0.01, y=0.71,label="Lake Environment:",hjust=0.0,size=10,colour="black")+
  annotate("text", x=0.032, y=0.68,label="Gull 2013",hjust=0.0,size=10,colour="black")+
  annotate("text", x=0.032, y=0.65,label="Wintergreen 2013",hjust=0.0,size=10,colour="gray")
Fig.4.C

################ Figure 5 ##############

metadata<-read.csv("./R_file17.csv")
metadata<-metadata[!is.na(metadata$rmax..per.day.),]
Fig.5<-ggplot(metadata,aes(Tree_Groups,rmax..per.day.,color=Tree_Groups,shape=Tree_Groups))+
  geom_jitter(size=5,width=0.15)+
  theme(legend.text.align = 0,legend.key=element_rect(fill=NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size = 18),axis.title=element_text(size=18,face="plain"),legend.title=element_text(size=18),legend.text=element_text(size=18))+
  scale_color_manual(guide='none',name="Phylogenetic Group",labels=c("oligotrophic"="Oligotrophic","pseudo-oligotrophic"="Pseudo-oligotrophic","nutrient-rich"="Nutrient-rich"),values=c("oligotrophic"="blue","pseudo-oligotrophic"="#56B4E9","nutrient-rich"="green"))+
  scale_size(guide='none')+
  xlab("\nPhylogenetic Group")+
  ylab("Growth Rate\n(rmax/day)")+
  scale_x_discrete(labels=c("oligotrophic"="Low Nutrient","pseudo-oligotrophic"="High Nutrient-Low","nutrient-rich"="High Nutrient"),limits=c("oligotrophic","pseudo-oligotrophic","nutrient-rich"))+
  scale_shape(guide='none')+
  annotate("text", x=0.5, y=0.55,label=paste("list(LMER:~F['2,17']==2.8", ", p==0.053", ", R^{2}==0.27)"),parse=TRUE,hjust=0.0,size=6)+
  annotate("text", x=0.5, y=0.52,label=paste("pairwise comparisons, N.S."),hjust=0.0,size=6)+
  geom_segment(aes(y=0.326794475,yend=0.326794475,x=2.65,xend=3.35),colour="green",linetype='dotted',size=2)+
  geom_segment(aes(y=0.162463515,yend=0.162463515,x=0.65,xend=1.35), colour="blue",linetype='dotted',size=2)+
  geom_segment(aes(y=0.198410635,yend=0.198410635,x=1.65,xend=2.35), colour="#56B4E9",linetype='dotted',size=2)
Fig.5





