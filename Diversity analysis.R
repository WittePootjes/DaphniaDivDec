#Packages
library(dada2); packageVersion("dada2")
library(vegan)
library(igraph)
library(SpiecEasi)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(stringr)
theme_set(theme_bw())

options(scipen = 999)

setwd('C:/Users/Anna/Desktop/Daphnia_files/Complete analysis')
#Uploading ASV counts table
ASV_counts_non_norm <- read.table('ASV_table.csv', header = TRUE, sep =',', row.names = 1)

# Normalize counts by the sample sums
sums <- apply(ASV_counts_non_norm, 2, sum)
norm_otus <- ASV_counts_non_norm

for (i in 1:ncol(ASV_counts_non_norm)){
  norm_otus[,i] <- norm_otus[,i]/sums[i]
}

#Checking if it worked:
sum(norm_otus$SRR9091053) # the sum of AVSs in for the first samples gives me 1 
sum(norm_otus$SRR9091054) # the sum of AVSs in for the first samples gives me 1 
ASV_counts_table <- norm_otus

#Uploading ASV taxonomy table
Tax_table <- read.table('ASV_taxonomy.csv', header = TRUE, sep =',', row.names = 1)
class(Tax_table)

# Uploading meta-data
Meta_data <- read.table('meta_data.csv', header = TRUE, sep = ',', row.names =1)
Meta_data <- t(Meta_data)
class(Meta_data)

colnames(norm_otus) == colnames(Meta_data) #checking if the column names match before cutting some of them out
# Deleting samples containing COMBO media, DNA extraction kit etc.
Meta_data_table <- subset(Meta_data, select = -c(44:71, 13, 92))
norm_otus <- subset(norm_otus, select = -c(44:71, 13, 92))

colnames(norm_otus) == colnames(Meta_data_table) # Columns match so I didnt accidentaly lose some samples

#Need to do some modifications in the variable formats before feeding into ps
norm_otus <- as.matrix(norm_otus)
class(norm_otus)
Tax_table <- as.matrix(Tax_table)
class(Tax_table)
Meta_data_table <- t(Meta_data_table)
Meta_data_table <- as.data.frame(Meta_data_table)

ps <- phyloseq(otu_table((norm_otus), taxa_are_rows = TRUE),
               tax_table(Tax_table), sample_data(Meta_data_table))

ps


# Creating a long data with melt fuction 

tax <- data.frame(tax_table(ps))
otus <- data.frame(otu_table(ps))

long_data <- cbind(tax, otus)
long_data <- reshape2::melt(long_data)

#adding a treatment column to my long_data table

long_data$Sample <- as.character(long_data$variable)
for (i in 1:length(long_data$Sample)){
  sample_name <- long_data$variable[[i]]
  long_data$Sample[[i]] <- as.character(Meta_data_table$experiment_title[rownames(Meta_data_table) == sample_name])
}

# Plotting Phylum using ggplot2
ggplot(data=long_data, aes(x=variable, y=value, fill=Phylum)) + 
  geom_bar(position='stack', stat='identity') + 
  labs(x='Sample', y='Relative abundance') + theme_minimal()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  facet_wrap(vars(Sample), ncol=2, scales='free') 

# Plotting Class
ggplot(data=long_data, aes(x=variable, y=value, fill= Class)) + 
  geom_bar(position='stack', stat='identity') + 
  labs(x='Sample', y='Relative abundance') + theme_minimal()+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  facet_wrap(vars(Sample), ncol=2, scales='free') 

#It is not possible to plot relative abundances difference at a genus level since there is just too many of them
# they are also skewed towards a few of them and therefore I will only focus on plotting the most abundant ones
# Since phyloseq does not include any function in its package to filter out the rare taxa, I used a function created by Lisa

prevFilter =function(phy, prev){
  prev0 = apply(X = otu_table(phy),
                MARGIN = ifelse(taxa_are_rows(phy), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
  prevalenceThreshold = prev * nsamples(phy)
  nonrares = prune_taxa((prev0 > prevalenceThreshold), phy)
  rares = prune_taxa((prev0 < prevalenceThreshold), phy)
  rares = merge_taxa(rares, taxa_names(rares))
  otus = data.frame(otu_table(nonrares))
  otus = rbind(otus, data.frame(otu_table(rares)))
  tax = data.frame(tax_table(nonrares), stringsAsFactors = FALSE)
  tax = rbind(tax, rep("Other", 7))
  rownames(otus) <- c(taxa_names(nonrares), 'Bin')
  rownames(tax) <- c(taxa_names(nonrares), 'Bin')
  newphy = phyloseq(otu_table(otus, taxa_are_rows = TRUE), sample_data(phy), tax_table(as.matrix(tax)))
  return(newphy)
}

prev60 <-prevFilter(ps, 0.6) #Filtering at 60% prevalence 
tax <-data.frame(tax_table(prev60))
otus <-data.frame(otu_table(prev60))

#creating a new long_data containing filtered out genus at 60%
long_data <-cbind(tax, otus)
long_data <- reshape2::melt(long_data)
long_data
length(long_data)

#Adding treatment column
long_data$Sample <- as.character(long_data$variable)
for (i in 1:length(long_data$Sample)){
  sample_name <- long_data$variable[[i]]
  long_data$Sample[[i]] <- as.character(Meta_data_table$experiment_title[rownames(Meta_data_table) == sample_name])
}

#Creating a plot with filtered out genus at over 60% frequency
ggplot(data=long_data, aes(x=variable, y=value, fill=Genus)) + # base ggplot2 data
  geom_bar(position='stack', stat='identity') + # defines stacked bar plot
  labs(x='Sample', y='Relative abundance for genus at over 60% frequency') + theme_minimal() + # adds titles and minimal look
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # removes x-axis labels
  facet_wrap(vars(Sample), ncol=2, scales='free') # facets by Description


##Diversity analysis # Diversity indices assume equal depths therefore the data will need to be rarified

#First I need to re-load my non-normalised ASV table
no_norm_ASV_counts <- read.table('ASV_table.csv', header = TRUE, sep =',', row.names = 1)
no_norm_ASV_counts <- as.matrix(no_norm_ASV_counts)

#incorporating them into ps elements:
ps <- phyloseq(otu_table(no_norm_ASV_counts, taxa_are_rows = TRUE),
               tax_table(as.matrix(Tax_table)), sample_data(Meta_data_table))
ps #all good

#First I am going to plot a rarefraction curve to get a general idea of how my samples look like
#Creating a rarefraction curve with a vegan package function "rarecurve"
rarecurve(t(no_norm_ASV_counts), step=100, lwd=2, ylab="ASVs", label=TRUE)
abline(v=(min(rowSums(t(no_norm_ASV_counts)))))
max(sample_sums(ps)) #44802
min(sample_sums(ps)) #7

#Based on the rarefraction curve, if we rarefied to the lowest sample we would lose a lot of data, therefore 
# I will rarify to the sample size 10,000 losing some samples but instead preserving the abundance of others

#For further analysis, I decided to rarify to more than 10,000 sample size:
rarefied_data <-subset_samples(ps,sample_sums(ps)>10000)
rarefied_data <-rarefy_even_depth(rarefied_data) #rarify without replacement
sample_sums(rarefied_data) # I can see that all samples now have been rarified to 11833


#### Alpha Diversty

#Chao1 Richness
meta <- data.frame(sample_data(rarefied_data)) #Accessing my sample information from the ps object containing rarefied data
otus <- data.frame(otu_table(rarefied_data, taxa_are_rows = TRUE)) #Acessing an otu table from the ps object containing rarefied data
sum(otus$SRR9091053) #quick check if it worked, the sum gives back a number of 11833 so its good to go
rarecurve((otus), step=100, lwd=2, ylab="ASVs", label=TRUE) #visualising the rarefraction curve after rarefying 

# Getting singletons and doubletons for Chao richness
meta$Richness <- colSums(otus != 0) #Calculating species richness by summing all non-zero values of my data
singletons <- colSums(otus == 1) #Summing all values that appeared only once 
doubletons <- colSums(otus == 2) # Summing all values that apperad twice
rares = (singletons/(2*doubletons))
rares[doubletons==0] <- 0 
meta$Chao1 = meta$Richness + rares

#Shannon Diversity
shannon_div <- function(vector){
  vector <- vector*log(vector)
  # the vector has NA values if the species proportions are 0
  vectorsum <- sum(na.omit(vector))
  return((1)*vectorsum)
}

# Calculate diversities
meta$Shannon = apply(otus, 2, shannon_div)
meta$Pielou = meta$Shannon / log(meta$Richness)
meta$Simpson = colSums(otus * otus)
meta$InvSimpson = 1 / meta$Simpson

# Get long data with diversities by 'melting' wide data into a long data 
long_data <- reshape2::melt(meta)
long_data = long_data[long_data$variable %in% c('Chao1', 'Shannon', 'Pielou', 'InvSimpson'),] #selecting indices

#Ploting diveristy indices
ggplot(long_data,aes(x=experiment_title, y=value, colour=variable))+
  geom_boxplot()+
  facet_wrap(vars(variable), scales='free', ncol=1)+
  theme_minimal()+ 
  theme(axis.text.x =element_text(angle=90, vjust=0.5, hjust=1), legend.position="none")+
  coord_flip()+labs(colour='Index', x='Value')

#Phyloseq package also contains their own functions to map the diversity indices

plot_richness(rarefied_data, x='experiment_title', measures=c("Shannon", "Chao1", "InvSimpson"))+
  geom_boxplot()+
  theme(axis.text.x =element_text(angle=90, vjust=0.5, hjust=1), legend.position="none")


#Statistical test
# Are the differences observed between different treatment groups significant?
# I will need to use a non-parametric test that will work for more than two groups of unrelates samples
#Kruskal Wallis test is adequete in this situation

# Are the diversity indices different across sample types?
print(kruskal.test(Chao1~experiment_title, data=meta)) # p-value = 0.1049
print(kruskal.test(Shannon~experiment_title, data=meta)) # p-value = 0.06425
print(kruskal.test(Pielou~experiment_title, data=meta)) #p-value = 0.08147
print(kruskal.test(InvSimpson~experiment_title, data=meta)) #p-value = 0.2463
# The diversity indices between the samples are not significantly different

#####Beta diversity analysis or MICROBIAL COMPOSITION
#For beta diversity I will use the normalised values

norm_otus
norm_meta <- Meta_data_table
taxa_table <- Tax_table

bray <- vegdist(t(norm_otus)) #creating dissimilarity matrix with vegdist using bray curtis
norm_meta <- as.matrix(norm_meta)

pcoa.res.sim <- capscale(t(norm_otus)~norm_meta, distance='bray', na.action='na.omit') # Ordination
eigenvalues <- eigenvals(pcoa.res.sim) #xtracting eigen values

#Plotting PCoA samples 

plot(pcoa.res.sim$CA$u[,c(1,2)],
     xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))),
     pch = 19)


### Creating Envit centroids with treatment
bray <- vegdist(t(norm_otus)) #creating dissimilarity matrix with B.Curtis
bray
meta <- as.matrix(norm_meta)
sim_meta <- data.frame(meta)

pcoa.res.sim <- capscale(t(norm_otus)~meta, distance='bray', na.action='na.omit') # Ordination
ef.sim = envfit(pcoa.res.sim, sim_meta, perm=1000, choices=c(1,2), na.rm=TRUE) #Envit function finds vectors or factor averages of environmental variables. 
print(ef.sim) 
# r2 is (0.5518) which means that the model doesn't explain much of variation of the data but it is significant (p = 0.000999)

#How metadata relates to microbiome variation?

library(ape)
norm_otus
Y=t(norm_otus)  #transposing the normalized counts and get the sample number.
n <- nrow(Y) #number of rows = number of samples
ev.stand <- scale(pcoa.res.sim$CA$u[,c(1,2)]) # Standardizing the first and second eigenvectors
S <- cov(Y, ev.stand) #calculating covariance between taxon abundances and the standardized eigenvectors
U <- S %*% diag((pcoa.res.sim$CA$eig[c(1,2)]/(n-1))^(-0.5)) #scaling the covariance matrix by the eigenvalues
U <- U*0.8# scaling U vectors so they fit better 

# filteing taxa based on the length of their vector:

norm_func <- function(x){
  return(sqrt(sum(x^2)))
}

norms=apply(U,1,norm_func)
sorted=sort(norms,index.return=TRUE,decreasing=TRUE)
U.selected <- U[sorted$ix[1:5],]

## Add vectors for top taxa
x = U.selected[,1]*0.4
y = U.selected[,2]*0.4
x0 = rep(0, length(x))
y0 = rep(0, length(y))
labels = rownames(U.selected)
arrows(x1=x, y1=y, x0=x0, y0=y0)
text(x=x, y=y, labels=labels, pos=3)


## Change the U.selected rownames to taxonomic information
rownames(U.selected) <- taxa_table[rownames(taxa_table)%in% rownames(U.selected), 4]

##Making plots in ggplot2

#Creating data.frames containing: 
sample_df <- data.frame(PCoA1=pcoa.res.sim$CA$u[,1], PCoA2=pcoa.res.sim$CA$u[,2], group=sim_meta) # stores my Eigencevectors
centroid_df <- data.frame(PCoA1=ef.sim$factors$centroids[,1]*0.15, PCoA2=ef.sim$factors$centroids[,2]*0.15,
                          centroid=rownames(ef.sim$factors$centroids)) #stores my centroids
arrow_df <- data.frame(PCoA1=U.selected[,1]*0.50, PCoA2=U.selected[,2]*0.50,
                       taxon=rownames(U.selected), x0=rep(0, nrow(U.selected)), y0=rep(0, nrow(U.selected))) #stores the vectors


# I will change the name of the centroids from Treatmentcontrol to just control etc. so it looks nicer

centroid_df[centroid_df == 'experiment_titleaztreonam'] <- 'aztreonam'
centroid_df[centroid_df == 'experiment_titleerythromycin'] <- 'erythromycin'
centroid_df[centroid_df == 'experiment_titlesulfamethoxazole'] <- 'sulfamethoxazole'
centroid_df[centroid_df == 'mix'] <- 'aztreonam, erythromycin, sulfamethoxazole'
centroid_df[centroid_df == 'experiment_titleno antibiotics'] <- 'no antibiotics'


#Plotting the centroids and vectors:

gplot <- ggplot() + geom_point(data=sample_df, aes(x=PCoA1, y=PCoA2, color=experiment_title), cex = 1.5) +
  geom_text(data=centroid_df, aes(x=PCoA1, y=PCoA2, label=centroid)) + theme_minimal()+
  geom_segment(data=arrow_df, aes(x=x0, y=y0, xend=PCoA1, yend=PCoA2),arrow=arrow()) +
  geom_text(data=arrow_df, aes(x=PCoA1, y=PCoA2, label=taxon), nudge_x=-0.05, nudge_y=0.05)+
  xlim(-0.4, 0.35) + theme_minimal() + scale_color_brewer(palette='Dark2') + labs(color='Antibiotic')
  
gplot


######TAXA-TAXA CORRELATIONS
# Filtering taxa at 60% occurance with prev60
class(norm_otus)
Meta_data_table <- as.data.frame(t(Meta_data_table))
Tax_table <- as.matrix(Tax_table)


#incorporating of a normalised otu table into ps elements:
ps <- phyloseq(otu_table((norm_otus), taxa_are_rows = TRUE),
               tax_table(Tax_table), sample_data(Meta_data_table))



#Applying filter
prev60 <-prevFilter(ps, 0.6) #Filtering at 60% prevalence 

#Accessing filtered otu_table
norm_otus <- otu_table(prev60)
dim(norm_otus) #I am left with 74 ASVs left

# define a function to calculate the geometric mean of a vector
  geom.mean<-function(x){return(exp(mean(log(x))))} 
  
# initialize a new matrix where we will store transformed data 
otu_clr=matrix(0,nrow=nrow(norm_otus),ncol=ncol(norm_otus)) #empty matrix

# loop samples
for(sample.index in 1:ncol(norm_otus)){
  comp=norm_otus[,sample.index] # obtain sample values
  nonzero.indices=which(comp>0) # obtain non-zero indices of the current sample
  g=geom.mean(comp[nonzero.indices]) # compute the geometric mean
  otu_clr[nonzero.indices,sample.index]=log(comp[nonzero.indices]/g) #populate the empty matrix with transformed data
}
# transfer row and column names
rownames(otu_clr)=rownames(norm_otus)
colnames(otu_clr)=colnames(norm_otus)
otu_clr #transformed matrix

#defining a function that can compute the association between two vectors in three different ways:
get.assoc<-function(x,y,method="euclid"){
  if(method %in% c("pearson","spearman"))return(cor(x,y,method=method)) #Pearson and Spearman
  if(method=="euclid")return(sqrt(sum((x-y)^2))) #Euclid
}

#a function to shuffle abundances
shuffle<-function(abundances,by="row"){
  dim=ifelse(by=="row",1,2)
  return(apply(abundances,dim,sample))
}

#Combining both functions to collect a number of random association values for two species.

x.index=55 #index number 76
y.index=56 #index number 77
iter=100
random.values=c()
for(i in 1:iter){
  shuffled=shuffle(norm_otus)
  rand.val=get.assoc(shuffled[x.index,],shuffled[y.index,],method="pearson")
  random.values=c(random.values,rand.val)
}

#Now, we can visualise the random distribution in a histogram and add the observed distribution as a vertical line.

observed.val=cor(t(norm_otus[x.index,]),t(norm_otus[y.index,]))
hist(random.values,xlab="Pearson correlation")
abline(v=observed.val,col="red")

r=length(random.values[abs(rand.val)>abs(observed.val)[1]]) # there is an error here!!
pval=(r+1)/(iter+1) 


N=nrow(norm_otus) #extracting number if rows
adj.matrix=matrix(NA,nrow=N,ncol=N) #adjacency matrix
assoc.matrix=matrix(NA,nrow=N,ncol=N) #association matrix

for(i in 2:N){
  
  for(j in 1:(i-1)){
    
    
    out=cor.test(t(norm_otus[i,]), t(norm_otus[j,]), method="spearman", exact=FALSE)
    
    
    assoc.matrix[i,j]=out$estimate
    
# estimate holds the Spearman correlation
    
    adj.matrix[i,j]=out$p.value	
    
    
  }
  
  
}
rownames(adj.matrix)=rownames(norm_otus) # transfer ASVs names
colnames(adj.matrix)=rownames(norm_otus)


#Extracting p-values
p.values=adj.matrix[lower.tri(adj.matrix)]
p.values.adj=p.adjust(p.values,method="BH")
adj.matrix[lower.tri(adj.matrix)]=p.values.adj

#converting q-values into significances
pseudocount=0.0000000001 # maximal significance value will be -log10(pseudocount)=10
sig.threshold=1 # threshold on significance, here set to 1 adj.matrix[adj.matrix==0]=pseudocount
adj.matrix=-1*log10(adj.matrix)
adj.matrix[adj.matrix < sig.threshold]=0 

#graph_from_adjacency_matrix
  
microbe.network=graph_from_adjacency_matrix(adj.matrix, mode="lower",weighted=TRUE, diag=FALSE) 
plot(microbe.network)

microbe.network=delete.vertices(microbe.network,degree(microbe.network)==0) #Delete nodes with 0 edges
plot(microbe.network)

write.graph(graph=microbe.network, file='network_at_80', format='graphml') #saving before analysis in Cytoscape


#Running SPIEC-EASI

library(SpiecEasi)
library(phyloseq)

#incorporating of a normalised otu table into ps elements:
ps <- phyloseq(otu_table((norm_otus), taxa_are_rows = TRUE),
               tax_table(Tax_table), sample_data(Meta_data_table))


spiec.out=spiec.easi(ps, method="mb",icov.select.params=list(rep.num=20)) #Estimating the covariance matris
betaMat=as.matrix(symBeta(getOptBeta(spiec.out))) #Extracting adjacency matrix

#Filtering out correlation matrix
correls=cor(t(norm_otus), method="spearman")
correls[abs(correls)<0.55]=0 

#Creating a microbe network
microbe.network=graph_from_adjacency_matrix(correls,mode="lower",weighted=TRUE, diag=FALSE)

#Plotting it
plot(microbe.network)

# Removing nodes that are not connected to other nodes
microbe.network=delete.vertices(microbe.network,degree(microbe.network)==0)

write.graph(graph=microbe.network, file='network_spiec.graphml', format='graphml') #Saving the network



