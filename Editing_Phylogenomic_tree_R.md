This is some basic code for messing with phylogenies and making them look nice in R. It uses a bastard combination of multiple packages:

```
library(ape)             # v3.4
library(ggtree)          # v1.2.10
library(ggplot2)         # v2.0.0
library(wesanderson)     # v0.3.2 only needed for certain colors
library(phangorn)        # v2.0.1
```
Here’s how you can read in a generic tree file:

tree=read.tree("/Volumes/hoekstrafs1/taste/trees/01_15_16/phylo_tree")
If for whatever reason you have single nodes in the tree (nodes with a single descendant), here’s a quick ape function that takes care of that with magic:

tree<-collapse.singles(tree)
GGtree also has a specific function for reading in  phylogenies. I’ve found that while this makes initial interpretation easier, it seriously limits downstream actions because the tree object is harder to use with all the packages here. Regardless, here’s how to read in a phylogenomic tree
Let’s say you have ugly labels at the tips of your branches, like so:

head(tree$tip.label)
```
[1] "BW_NW_006501958.1_638666_641604_ORF_1"
[2] "CA_TRINITY_GG_101_c0_g1_i1_ORF_1"     
[3] "PO_TRINITY_GG_90_c3_g1_i1_ORF_1"      
[4] "LL_TRINITY_GG_87_c4_g1_i1_ORF_1"      
[5] "LL_TRINITY_GG_286_c0_g1_i1_ORF_1"     
[6] "BW_NW_006501958.1_646973_649914_ORF_1"
```
It’s easy to clean them up for easier presentation down the line:

```
seq = tree$tip.label  #Create vector of tip labels
dd = data.frame(seq)  #Make that the first column of a data frame
```
Below, we add a column to the data frame called 'taxa'. Using grep, we can easily populate this column with the names we want to use for our final
tree
```
dd$Species = 0
dd$Species[grep("PO", dd$seq)] = "Peromyscus polionotus"
dd$Species[grep("BW", dd$seq)] = "Peromyscus maniculatus"
dd$Species[grep("LL", dd$seq)] = "Peromyscus leucopus"
dd$Species[grep("CA", dd$seq)] = "Peromyscus californicus"
print(dd[1:4, ])
```
```  
                                   seq                 Species
1 BW_NW_006501958.1_638666_641604_ORF_1  Peromyscus maniculatus
2      CA_TRINITY_GG_101_c0_g1_i1_ORF_1 Peromyscus californicus
3       PO_TRINITY_GG_90_c3_g1_i1_ORF_1   Peromyscus polionotus
4       LL_TRINITY_GG_87_c4_g1_i1_ORF_1     Peromyscus leucopus
```
Now on to plotting. Here I’ve created a function that will assign colors to each taxon in the tree. Just input the dataframe you created earlier:

```
cls <- function(dd){
  ss<-sort(unique(dd$Species))
  colors<-setNames(palette()[1:length(ss)],ss)
  colors <- colors[dd$Species]
  return(colors)
}

colors = cls(dd)
```
And now we can plot a tree! Below is the most basic plotting function I use, from the ape package. You’ll notice that the tip labels are still the ugly ones from earlier; this is because I often use this method to very quickly plot a tree for my interpretation, rather than presenting to someone else. Plotting a nicely formatted tree takes a bit more tinkering, which I’ll show later.
```
plot(tree,tip.color=colors,cex=0.23,font=2,edge.width=c(1.5),no.margin=TRUE)
legend("topright",sort(unique(dd$Species)),fill=unique(colors))

```
