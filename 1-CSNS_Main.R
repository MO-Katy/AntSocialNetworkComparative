# Set directories, load packages & import data
MAINDIR  <- "/Volumes/Lacie/CSNS/data"
FIGDIR   <- "/Volumes/Lacie/CSNS/figures"
FACETDIR <- "/Volumes/Lacie/CSNS/FacetNet_Output"
library(igraph); library(viridis);library(jpeg); library(assertthat);library(devtools);
library(data.table); library(DescTools);library(diptest);library(ggfortify); library(lme4);
library(ggplot2); library(gridExtra); library(lmerTest); library(graphics); library(patchwork)

setwd(MAINDIR)
for (file in list.files()){
  assign(sub("\\..*", "", file), read.csv(file, check.names = FALSE))
}

# Create lists for each data type
nest_list <- list(nest_cfel1, nest_cfel2, nest_cfel3, nest_cfel4, nest_cfel5,
                  nest_pbar1, nest_pbar2, nest_pbar3, nest_pbar4, nest_pbar5,
                  nest_drug1, nest_drug2, nest_drug3, nest_drug4, nest_drug5,
                  nest_ipur1, nest_ipur2, nest_ipur3, nest_ipur4, nest_ipur5,
                  nest_rmet1, nest_rmet2, nest_rmet3, nest_rmet4, nest_rmet5)

forage_list <- list(forage_cfel1, forage_cfel2, forage_cfel3, forage_cfel4, forage_cfel5,
                    forage_pbar1, forage_pbar2, forage_pbar3, forage_pbar4, forage_pbar5,
                    forage_drug1, forage_drug2, forage_drug3, forage_drug4, forage_drug5,
                    forage_ipur1, forage_ipur2, forage_ipur3, forage_ipur4, forage_ipur5,
                    forage_rmet1, forage_rmet2, forage_rmet3, forage_rmet4, forage_rmet5)

md_list <- list(md_cfel1, md_cfel2, md_cfel3, md_cfel4, md_cfel5,
                md_pbar1, md_pbar2, md_pbar3, md_pbar4, md_pbar5,
                md_drug1, md_drug2, md_drug3, md_drug4, md_drug5,
                md_ipur1, md_ipur2, md_ipur3, md_ipur4, md_ipur5,
                md_rmet1, md_rmet2, md_rmet3, md_rmet4, md_rmet5)

interaction_list <- list(interaction_cfel1, interaction_cfel2, interaction_cfel3, interaction_cfel4, interaction_cfel5,
                         interaction_pbar1, interaction_pbar2, interaction_pbar3, interaction_pbar4, interaction_pbar5,
                         interaction_drug1, interaction_drug2, interaction_drug3, interaction_drug4, interaction_drug5,
                         interaction_ipur1, interaction_ipur2, interaction_ipur3, interaction_ipur4, interaction_ipur5,
                         interaction_rmet1, interaction_rmet2, interaction_rmet3, interaction_rmet4, interaction_rmet5)

# Uniformize row name format
for (i in 1:25){
  rownames(nest_list[[i]]) <- rownames(forage_list[[i]]) <- md_list[[i]]$ant_id
}

# Identify outliers based on total number of times seen
outliers_list <- list()
for (i in 1:25) {
  outliers <- md_list[[i]][(md_list[[i]]$counts > mean(md_list[[i]]$counts) + (2*sd(md_list[[i]]$counts))) | (md_list[[i]]$counts < mean(md_list[[i]]$counts) - (2*sd(md_list[[i]]$counts))),]$ant_id
  outliers_list[[i]] <- outliers
}

# In some replicates the queens are identified as outliers. Retain queens.
outliers_list[[15]] <- outliers_list[[15]][-which(outliers_list[[15]] == 68)] # drug5
outliers_list[[19]] <- outliers_list[[19]][-which(outliers_list[[19]] == 4)]  # ipur4
outliers_list[[21]] <- outliers_list[[21]][-which(outliers_list[[21]] == 38)]  # rmet1

# remove outliers from metadata, interaction data and spatial data
for (i in 1:25){
  md_list[[i]] <- md_list[[i]][!(md_list[[i]]$ant_id %in% outliers_list[[i]]),]
  interaction_list[[i]] <- interaction_list[[i]][!(interaction_list[[i]]$id1 %in% outliers_list[[i]] |
                                                     interaction_list[[i]]$id2 %in% outliers_list[[i]]),]
  forage_list[[i]] <- forage_list[[i]][!(forage_list[[i]][,1] %in% outliers_list[[i]]),]
  nest_list[[i]]   <- nest_list[[i]][!(nest_list[[i]][,1] %in% outliers_list[[i]]),]
}

# Restrict interactions to head-head contact only
interaction_hh_list <- list()
for (i in 1:25){
  if (i == 5){
    # In cfel5 head and bodies were inverted in FortStudio
    interaction_hh_list[[i]] <- interaction_list[[i]][interaction_list[[i]]$`2-2` == "True",c(1,2)]
  } else {
    interaction_hh_list[[i]] <- interaction_list[[i]][interaction_list[[i]]$`1-1` == "True",c(1,2)]
  }
}

# Order interactions first by id1 and then by id2
for (i in 1:25){
  interaction_hh_list[[i]] <- interaction_hh_list[[i]][order(interaction_hh_list[[i]]$id1,interaction_hh_list[[i]]$id2),]
}

# Collapse interaction data into edgelists
edgelist_list <- list()
for (i in 1:25){
  edgelist_list[[i]] <- as.data.frame(setDT(interaction_hh_list[[i]])[,list(Count=.N),names(interaction_hh_list[[i]])])
}

# Convert IDs to characters - otherwise iGraph acts up
for (i in 1:25){
  edgelist_list[[i]]$id1 <- as.character(edgelist_list[[i]]$id1)
  edgelist_list[[i]]$id2 <- as.character(edgelist_list[[i]]$id2)
}

# Save edgelists
edgelist_filesnames <- c("edgelist_cfel1.csv", "edgelist_cfel2.csv", "edgelist_cfel3.csv", "edgelist_cfel4.csv", "edgelist_cfel5.csv",
                         "edgelist_pbar1.csv", "edgelist_pbar2.csv", "edgelist_pbar3.csv", "edgelist_pbar4.csv", "edgelist_pbar5.csv",
                         "edgelist_drug1.csv", "edgelist_drug2.csv", "edgelist_drug3.csv", "edgelist_drug4.csv", "edgelist_drug5.csv",
                         "edgelist_ipur1.csv", "edgelist_ipur2.csv", "edgelist_ipur3.csv", "edgelist_ipur4.csv", "edgelist_ipur5.csv",
                         "edgelist_rmet1.csv", "edgelist_rmet2.csv", "edgelist_rmet3.csv", "edgelist_rmet4.csv", "edgelist_rmet5.csv")
for (i in 1:25){
  write.csv(edgelist_list[[i]], edgelist_filesnames[[i]], row.names = FALSE)
}

# Construct iGraph object
net_list <- list()
for (i in 1:25){
  net_list[[i]] <- graph_from_edgelist(as.matrix(edgelist_list[[i]][,1:2]), directed = FALSE)
}

# Assign edge weights
for (i in 1:25){
  E(net_list[[i]])$weight <- edgelist_list[[i]][,3]
}

# Calculate layouts with the force-directed Fruchterman-Reingold layout algorithm
layout_list <- list()
for (i in 1:25){
  layout_list[[i]] <- layout_with_fr(net_list[[i]])
}

# Color edges by their strength
for (i in 1:25){
  if(i %in% 1:5){
    E(net_list[[i]])$color <- gray((1-(E(net_list[[i]])$weight/max(E(net_list[[i]])$weight)))^3)
  } else if (i %in% 6:10){
    E(net_list[[i]])$color <- gray((1-(E(net_list[[i]])$weight/max(E(net_list[[i]])$weight)))^.5)
  } else if (i %in% 11:20){
    E(net_list[[i]])$color <- gray((1-(E(net_list[[i]])$weight/max(E(net_list[[i]])$weight))))
  } else if (i %in% 21:25){
    E(net_list[[i]])$color <- gray((1-(E(net_list[[i]])$weight/max(E(net_list[[i]])$weight)))^.2)
  }
}

# Define node shape
shape_list <- list()
for (i in 1:25){
  shape_list[[i]] <- rep("circle", length(V(net_list[[i]])))
}

# Make queens square
queen_list <- list(33, 4, 79, 66, 37,
                   38, 60, 32, 24, 49,
                   14, 28, 105, 17, 68,
                   1, 77, 68, 4, 67,
                   38, 68, 112, 34, 19)
for (i in 1:25){
  shape_list[[i]][which(names(V(net_list[[i]])) == queen_list[[i]])] <- "square"
}

# Apply shape vectors to network objects
for (i in 1:25){
  V(net_list[[i]])$shape <- shape_list[[i]]
}

# Define vertex colors
vertexCol_list <- list()
for (i in 1:25){
  col_ind <- c()
  for(j in names(V(net_list[[i]]))){
    col_ind <- c(col_ind, which(md_list[[i]]$ant_id == j))
  }
  vertexCol_list[[i]] <- gray(rowSums(nest_list[[i]][,-1]) / (rowSums(forage_list[[i]][,-1]) + rowSums(nest_list[[i]][,-1])))[col_ind]
}

# Recolor queens
for (i in 1:25){
  vertexCol_list[[i]][which(shape_list[[i]] == "square")] <- "magenta"
}

# Plot networks
setwd(FIGDIR)
jpeg('Cfel_Networks.jpg', width=6000, height=3000, unit='px')
par(mfrow = c(2,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
for (i in 1:5){
  plot(net_list[[i]], layout = layout_list[[i]],
       vertex.size=5, vertex.label=NA, vertex.color = vertexCol_list[[i]],
       edge.width = ((E(net_list[[i]])$weight/max(E(net_list[[i]])$weight))^1.2)*10)
}
dev.off()
jpeg('Pbar_Networks.jpg', width=6000, height=3000, unit='px')
par(mfrow = c(2,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
for (i in 6:10){
  plot(net_list[[i]], layout = layout_list[[i]],
       vertex.size=5, vertex.label=NA, vertex.color = vertexCol_list[[i]],
       edge.width = ((E(net_list[[i]])$weight/max(E(net_list[[i]])$weight))^.6)*30)
}
dev.off()
jpeg('Drug_Networks.jpg', width=6000, height=3000, unit='px')
par(mfrow = c(2,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
for (i in 11:15){
  plot(net_list[[i]], layout = layout_list[[i]],
       vertex.size=5, vertex.label=NA, vertex.color = vertexCol_list[[i]],
       edge.width = ((E(net_list[[i]])$weight/max(E(net_list[[i]])$weight))^1.2)*8)
}
dev.off()
jpeg('Ipur_Networks.jpg', width=6000, height=3000, unit='px')
par(mfrow = c(2,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
for (i in 16:20){
  plot(net_list[[i]], layout = layout_list[[i]],
       vertex.size=5, vertex.label=NA, vertex.color = vertexCol_list[[i]],
       edge.width = ((E(net_list[[i]])$weight/max(E(net_list[[i]])$weight))^.7)*20)
}
dev.off()
jpeg('Rmet_Networks.jpg', width=6000, height=3000, unit='px')
par(mfrow = c(2,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
for (i in 21:25){
  plot(net_list[[i]], layout = layout_list[[i]],
       vertex.size=5, vertex.label=NA, vertex.color = vertexCol_list[[i]],
       edge.width = ((E(net_list[[i]])$weight/max(E(net_list[[i]])$weight))^1.2)*20)
}
dev.off()

# Calculate and plot the strength distribution for each network/ species
strengthNorm_list <- list()
for (i in 1:25){
  strengthNorm_list[[i]] <- as.numeric(scale(strength(net_list[[i]])))[order(as.numeric(names(strength(net_list[[i]]))))]
}

jpeg('Strength_Distribution.jpg', width=6000, height=2000, unit='px')
par(mfrow = c(1,5), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
hist(c(strengthNorm_list[[1]], strengthNorm_list[[2]], strengthNorm_list[[3]], strengthNorm_list[[3]], strengthNorm_list[[5]]),
     breaks = 15, col = "black", main = substitute(paste(italic("C. fellah"))), xlab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12)
axis(1, lwd = 8, labels = FALSE)
axis(2, lwd = 8, labels = FALSE)
hist(c(strengthNorm_list[[11]], strengthNorm_list[[12]], strengthNorm_list[[13]], strengthNorm_list[[14]], strengthNorm_list[[15]]),
     breaks = 15, col = "black", main = substitute(paste(italic("D. rugosum"))), xlab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12)
axis(1, lwd = 8, labels = FALSE)
axis(2, lwd = 8, labels = FALSE)
hist(c(strengthNorm_list[[16]], strengthNorm_list[[17]], strengthNorm_list[[18]], strengthNorm_list[[19]], strengthNorm_list[[20]]),
     breaks = 15, col = "black", main = substitute(paste(italic("I. purpureus"))), xlab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12)
axis(1, lwd = 8, labels = FALSE)
axis(2, lwd = 8, labels = FALSE)
hist(c(strengthNorm_list[[6]], strengthNorm_list[[7]], strengthNorm_list[[8]], strengthNorm_list[[9]], strengthNorm_list[[10]]),
     breaks = 15, col = "black", main = substitute(paste(italic("P. rugosus"))), xlab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12)
axis(1, lwd = 8, labels = FALSE)
axis(2, lwd = 8, labels = FALSE)
hist(c(strengthNorm_list[[21]], strengthNorm_list[[22]], strengthNorm_list[[23]], strengthNorm_list[[24]], strengthNorm_list[[25]]),
     breaks = 15, col = "black", main = substitute(paste(italic("R. metallica"))), xlab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12)
axis(1, lwd = 8, labels = FALSE)
axis(2, lwd = 8, labels = FALSE)
dev.off()

all_normalized_strengths <- c()
for (i in 1:25){
  all_normalized_strengths <- c(all_normalized_strengths, strengthNorm_list[[i]])
}

jpeg('Strength_Distribution_Combined.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
hist(all_normalized_strengths,
     breaks = 30, col = "black", main = "Strength distribution", xlab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12)
axis(1, lwd = 8, labels = FALSE)
axis(2, lwd = 8, labels = FALSE)
dev.off() # This seems to follow a log-normal distribution

# Normalize foraging scores within colony
foraging_list <- list()
for(i in 1:25){
  foraging_list[[i]] <- rowSums(forage_list[[i]][,-1]) / (rowSums(forage_list[[i]][,-1]) + rowSums(nest_list[[i]][,-1]))
}

# Concatenate foraging scores per species
cfel_foraging <- c(foraging_list[[1]], foraging_list[[2]], foraging_list[[3]], foraging_list[[4]], foraging_list[[5]])
pbar_foraging <- c(foraging_list[[6]], foraging_list[[7]], foraging_list[[8]], foraging_list[[9]], foraging_list[[10]])
drug_foraging <- c(foraging_list[[11]], foraging_list[[12]], foraging_list[[13]], foraging_list[[14]], foraging_list[[15]])
ipur_foraging <- c(foraging_list[[16]], foraging_list[[17]], foraging_list[[18]], foraging_list[[19]], foraging_list[[20]])
rmet_foraging <- c(foraging_list[[21]], foraging_list[[22]], foraging_list[[23]], foraging_list[[24]], foraging_list[[25]])

# Concatenate normalized strength scores per species
cfel_strengthNorm <- c(strengthNorm_list[[1]], strengthNorm_list[[2]], strengthNorm_list[[3]], strengthNorm_list[[4]], strengthNorm_list[[5]])
pbar_strengthNorm <- c(strengthNorm_list[[6]], strengthNorm_list[[7]], strengthNorm_list[[8]], strengthNorm_list[[9]], strengthNorm_list[[10]])
drug_strengthNorm <- c(strengthNorm_list[[11]], strengthNorm_list[[12]], strengthNorm_list[[13]], strengthNorm_list[[14]], strengthNorm_list[[15]])
ipur_strengthNorm <- c(strengthNorm_list[[16]], strengthNorm_list[[17]], strengthNorm_list[[18]], strengthNorm_list[[19]], strengthNorm_list[[20]])
rmet_strengthNorm <- c(strengthNorm_list[[21]], strengthNorm_list[[22]], strengthNorm_list[[23]], strengthNorm_list[[24]], strengthNorm_list[[25]])

# Set colors for replicates
cfel_col      <- c(rep("midnightblue", length(foraging_list[[1]])), rep("sandybrown", length(foraging_list[[2]])), rep("cornflowerblue", length(foraging_list[[3]])), rep("#50C878", length(foraging_list[[4]])), rep("#FF7F50", length(foraging_list[[5]])))
pbar_col      <- c(rep("midnightblue", length(foraging_list[[6]])), rep("sandybrown", length(foraging_list[[7]])), rep("cornflowerblue", length(foraging_list[[8]])), rep("#50C878", length(foraging_list[[9]])), rep("#FF7F50", length(foraging_list[[10]])))
drug_col      <- c(rep("midnightblue", length(foraging_list[[11]])), rep("sandybrown", length(foraging_list[[12]])), rep("cornflowerblue", length(foraging_list[[13]])), rep("#50C878", length(foraging_list[[14]])), rep("#FF7F50", length(foraging_list[[15]])))
ipur_col      <- c(rep("midnightblue", length(foraging_list[[16]])), rep("sandybrown", length(foraging_list[[17]])), rep("cornflowerblue", length(foraging_list[[18]])), rep("#50C878", length(foraging_list[[19]])), rep("#FF7F50", length(foraging_list[[20]])))
rmet_col      <- c(rep("midnightblue", length(foraging_list[[21]])), rep("sandybrown", length(foraging_list[[22]])), rep("cornflowerblue", length(foraging_list[[23]])), rep("#50C878", length(foraging_list[[24]])), rep("#FF7F50", length(foraging_list[[25]])))

# Shuffle so that one color is not in the forefront of the plots
cfel_shuffle <- sample(1:length(cfel_col))
pbar_shuffle <- sample(1:length(pbar_col))
drug_shuffle <- sample(1:length(drug_col))
ipur_shuffle <- sample(1:length(ipur_col))
rmet_shuffle <- sample(1:length(rmet_col))


# Plot node strength versus foraging
jpeg('Strength_vs_Foraging.jpg', width=6000, height=2000, unit='px')
par(mfrow = c(1,5), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(cfel_foraging[cfel_shuffle] ~ cfel_strengthNorm[cfel_shuffle], pch = 16, main = substitute(paste(italic("C. fellah"))), xlab = "Strength", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylim = c(0,1), col = cfel_col)
axis(2, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
axis(1, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = cfel_strengthNorm, y = cfel_foraging), lwd = 10)
plot(drug_foraging[drug_shuffle] ~ drug_strengthNorm[drug_shuffle], pch = 16, main = substitute(paste(italic("D. rugosum"))), xlab = "Strength", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylim = c(0,1), col = drug_col)
axis(2, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
axis(1, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = drug_strengthNorm, y = drug_foraging), lwd = 10)
plot(pbar_foraging[pbar_shuffle] ~ pbar_strengthNorm[pbar_shuffle], pch = 16, main = substitute(paste(italic("P. rugosus"))), xlab = "Strength", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylim = c(0,1), col = pbar_col)
axis(2, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
axis(1, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = pbar_strengthNorm, y = pbar_foraging), lwd = 10)
plot(ipur_foraging[ipur_shuffle] ~ ipur_strengthNorm[ipur_shuffle], pch = 16, main = substitute(paste(italic("I. purpureus"))), xlab = "Strength", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylim = c(0,1), col = ipur_col)
axis(2, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
axis(1, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = ipur_strengthNorm, y = ipur_foraging), lwd = 10)
plot(rmet_foraging[rmet_shuffle] ~ rmet_strengthNorm[rmet_shuffle], pch = 16, main = substitute(paste(italic("R. metallica"))), xlab = "Strength", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylim = c(0,1), col = rmet_col)
axis(2, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
axis(1, labels = FALSE, tick = FALSE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = rmet_strengthNorm, y = rmet_foraging), lwd = 10)
dev.off()

# Import social maturity data
setwd(FACETDIR)
for (file in list.files()[which(grepl("_c", list.files(), fixed = TRUE))]){
  setwd(paste(FACETDIR, file, sep = "/"))
  assign(sub("\\..*", "", paste(file, "SoftMod", sep = "_")), read.csv(list.files()[which(grepl("soft_comm_step", list.files(), fixed = TRUE))]))
}

softMod_list <- list(cfel1_c_SoftMod, cfel2_c_SoftMod, cfel3_c_SoftMod, cfel4_c_SoftMod, cfel5_c_SoftMod,
                     pbar1_c_SoftMod, pbar2_c_SoftMod, pbar3_c_SoftMod, pbar4_c_SoftMod, pbar5_c_SoftMod,
                     drug1_c_SoftMod, drug2_c_SoftMod, drug3_c_SoftMod, drug4_c_SoftMod, drug5_c_SoftMod,
                     ipur1_c_SoftMod, ipur2_c_SoftMod, ipur3_c_SoftMod, ipur4_c_SoftMod, ipur5_c_SoftMod,
                     rmet1_c_SoftMod, rmet2_c_SoftMod, rmet3_c_SoftMod, rmet4_c_SoftMod, rmet5_c_SoftMod)

# Determine community identity
for (i in 1:25){
  foraging_list[[i]] <- foraging_list[[i]][forage_list[[i]][,1] %in% softMod_list[[i]]$id]
  strengthNorm_list[[i]] <- strengthNorm_list[[i]][forage_list[[i]][,1] %in% softMod_list[[i]]$id]
  if(cor.test(foraging_list[[i]], softMod_list[[i]]$cluster_0)$estimate < 0){
    colnames(softMod_list[[i]]) <- c("id", "nurse", "foraging")} else{colnames(softMod_list[[i]]) <- c("id", "foraging", "nurse")}
}

# Fpraging ~ social maturity dataframe
cfel_ForagingMaturity <- data.frame(Maturity = c(softMod_list[[1]]$foraging,
                                                 softMod_list[[2]]$foraging,
                                                 softMod_list[[3]]$foraging,
                                                 softMod_list[[4]]$foraging,
                                                 softMod_list[[5]]$foraging),
                                    Foraging = c(foraging_list[[1]],
                                                 foraging_list[[2]],
                                                 foraging_list[[3]],
                                                 foraging_list[[4]],
                                                 foraging_list[[5]]),
                                    Strength = c(strengthNorm_list[[1]],
                                                 strengthNorm_list[[2]],
                                                 strengthNorm_list[[3]],
                                                 strengthNorm_list[[4]],
                                                 strengthNorm_list[[5]]),
                                    col = c(rep("midnightblue", length(softMod_list[[1]]$foraging)),
                                            rep("sandybrown", length(softMod_list[[2]]$foraging)),
                                            rep("slategray", length(softMod_list[[3]]$foraging)),
                                            rep("red3", length(softMod_list[[4]]$foraging)),
                                            rep("darkcyan", length(softMod_list[[5]]$foraging))))
drug_ForagingMaturity <- data.frame(Maturity = c(softMod_list[[11]]$foraging,
                                                 softMod_list[[12]]$foraging,
                                                 softMod_list[[13]]$foraging,
                                                 softMod_list[[14]]$foraging,
                                                 softMod_list[[15]]$foraging),
                                    Foraging = c(foraging_list[[11]],
                                                 foraging_list[[12]],
                                                 foraging_list[[13]],
                                                 foraging_list[[14]],
                                                 foraging_list[[15]]),
                                    Strength = c(strengthNorm_list[[11]],
                                                 strengthNorm_list[[12]],
                                                 strengthNorm_list[[13]],
                                                 strengthNorm_list[[14]],
                                                 strengthNorm_list[[15]]),
                                    col = c(rep("midnightblue", length(softMod_list[[11]]$foraging)),
                                            rep("sandybrown", length(softMod_list[[12]]$foraging)),
                                            rep("slategray", length(softMod_list[[13]]$foraging)),
                                            rep("red3", length(softMod_list[[14]]$foraging)),
                                            rep("darkcyan", length(softMod_list[[15]]$foraging))))
ipur_ForagingMaturity <- data.frame(Maturity = c(softMod_list[[16]]$foraging,
                                                 softMod_list[[17]]$foraging,
                                                 softMod_list[[18]]$foraging,
                                                 softMod_list[[19]]$foraging,
                                                 softMod_list[[20]]$foraging),
                                    Foraging = c(foraging_list[[16]],
                                                 foraging_list[[17]],
                                                 foraging_list[[18]],
                                                 foraging_list[[19]],
                                                 foraging_list[[20]]),
                                    Strength = c(strengthNorm_list[[16]],
                                                 strengthNorm_list[[17]],
                                                 strengthNorm_list[[18]],
                                                 strengthNorm_list[[19]],
                                                 strengthNorm_list[[20]]),
                                    col = c(rep("midnightblue", length(softMod_list[[16]]$foraging)),
                                            rep("sandybrown", length(softMod_list[[17]]$foraging)),
                                            rep("slategray", length(softMod_list[[18]]$foraging)),
                                            rep("red3", length(softMod_list[[19]]$foraging)),
                                            rep("darkcyan", length(softMod_list[[20]]$foraging))))
pbar_ForagingMaturity <- data.frame(Maturity = c(softMod_list[[6]]$foraging,
                                                 softMod_list[[7]]$foraging,
                                                 softMod_list[[8]]$foraging,
                                                 softMod_list[[9]]$foraging,
                                                 softMod_list[[10]]$foraging),
                                    Foraging = c(foraging_list[[6]],
                                                 foraging_list[[7]],
                                                 foraging_list[[8]],
                                                 foraging_list[[9]],
                                                 foraging_list[[10]]),
                                    Strength = c(strengthNorm_list[[6]],
                                                 strengthNorm_list[[7]],
                                                 strengthNorm_list[[8]],
                                                 strengthNorm_list[[9]],
                                                 strengthNorm_list[[10]]),
                                    col = c(rep("midnightblue", length(softMod_list[[6]]$foraging)),
                                            rep("sandybrown", length(softMod_list[[7]]$foraging)),
                                            rep("slategray", length(softMod_list[[8]]$foraging)),
                                            rep("red3", length(softMod_list[[9]]$foraging)),
                                            rep("darkcyan", length(softMod_list[[10]]$foraging))))
rmet_ForagingMaturity <- data.frame(Maturity = c(softMod_list[[21]]$foraging,
                                                 softMod_list[[22]]$foraging,
                                                 softMod_list[[23]]$foraging,
                                                 softMod_list[[24]]$foraging,
                                                 softMod_list[[25]]$foraging),
                                    Foraging = c(foraging_list[[21]],
                                                 foraging_list[[22]],
                                                 foraging_list[[23]],
                                                 foraging_list[[24]],
                                                 foraging_list[[25]]),
                                    Strength = c(strengthNorm_list[[21]],
                                                 strengthNorm_list[[22]],
                                                 strengthNorm_list[[23]],
                                                 strengthNorm_list[[24]],
                                                 strengthNorm_list[[25]]),
                                    col = c(rep("midnightblue", length(softMod_list[[21]]$foraging)),
                                            rep("sandybrown", length(softMod_list[[22]]$foraging)),
                                            rep("slategray", length(softMod_list[[23]]$foraging)),
                                            rep("red3", length(softMod_list[[24]]$foraging)),
                                            rep("darkcyan", length(softMod_list[[25]]$foraging))))

# Shuffle for plotting
cfel_ForagingMaturity_SH <- cfel_ForagingMaturity[sample(1:nrow(cfel_ForagingMaturity)),]
drug_ForagingMaturity_SH <- drug_ForagingMaturity[sample(1:nrow(drug_ForagingMaturity)),]
ipur_ForagingMaturity_SH <- ipur_ForagingMaturity[sample(1:nrow(ipur_ForagingMaturity)),]
pbar_ForagingMaturity_SH <- pbar_ForagingMaturity[sample(1:nrow(pbar_ForagingMaturity)),]
rmet_ForagingMaturity_SH <- rmet_ForagingMaturity[sample(1:nrow(rmet_ForagingMaturity)),]

# Plot stacked histograms of social maturity distributions
common_theme <- theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.text = element_text(size = 35, family = "serif", margin = margin(t = -50)),
        axis.title = element_text(size = 35, family = "serif", margin = margin(t = -50)),
        plot.title = element_text(size = 50, family = "serif", hjust=0.5))
hist1 <- ggplot(data = cfel_ForagingMaturity_SH, aes(x = Maturity, fill = col)) +
  geom_histogram(binwidth = 0.1, position = "stack") +
  xlab("Maturity") +
  ylab("Frequency") +
  scale_fill_manual(values = c("midnightblue", "sandybrown", "slategray", "red3", "darkcyan")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,1),
                     labels = c(0,1)) +
  scale_y_continuous(breaks = c(0,140),
                     labels = c(0,140),
                     limits = c(0,140)) +
  labs(x = "Maturity", y = "Frequency", title = expression(italic("C. fellah"))) +
  common_theme
hist2 <- ggplot(data = drug_ForagingMaturity_SH, aes(x = Maturity, fill = col)) +
  geom_histogram(binwidth = 0.1, position = "stack") +
  xlab("Maturity") +
  ylab("Frequency") +
  scale_fill_manual(values = c("midnightblue", "sandybrown", "slategray", "red3", "darkcyan")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,1),
                     labels = c(0,1)) +
  scale_y_continuous(breaks = c(0,140),
                     labels = c(0,140),
                     limits = c(0,140)) +
  labs(x = "Maturity", y = "Frequency", title = expression(italic("D. rugosum")))  +
  common_theme
hist3 <- ggplot(data = ipur_ForagingMaturity_SH, aes(x = Maturity, fill = col)) +
  geom_histogram(binwidth = 0.1, position = "stack") +
  xlab("Maturity") +
  ylab("Frequency") +
  scale_fill_manual(values = c("midnightblue", "sandybrown", "slategray", "red3", "darkcyan")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,1),
                     labels = c(0,1)) +
  scale_y_continuous(breaks = c(0,140),
                     labels = c(0,140),
                     limits = c(0,140)) +
  labs(x = "Maturity", y = "Frequency", title = expression(italic("I. purpureus"))) +
  common_theme
hist4 <- ggplot(data = pbar_ForagingMaturity_SH, aes(x = Maturity, fill = col)) +
  geom_histogram(binwidth = 0.1, position = "stack") +
  xlab("Maturity") +
  ylab("Frequency") +
  ggtitle(expression(paste(italic("C. fellah")))) +
  scale_fill_manual(values = c("midnightblue", "sandybrown", "slategray", "red3", "darkcyan")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,1),
                     labels = c(0,1)) +
  scale_y_continuous(breaks = c(0,140),
                     labels = c(0,140),
                     limits = c(0,140)) +
  labs(x = "Maturity", y = "Frequency", title = expression(italic("P. barbatus"))) +
  common_theme
hist5 <- ggplot(data = rmet_ForagingMaturity_SH, aes(x = Maturity, fill = col)) +
  geom_histogram(binwidth = 0.1, position = "stack") +
  xlab("Maturity") +
  ylab("Frequency") +
  ggtitle(expression(paste(italic("C. fellah")))) +
  scale_fill_manual(values = c("midnightblue", "sandybrown", "slategray", "red3", "darkcyan")) +
  guides(fill = "none") +
  scale_x_continuous(breaks = c(0,1),
                     labels = c(0,1)) +
  scale_y_continuous(breaks = c(0,140),
                     labels = c(0,140),
                     limits = c(0,140)) +
  labs(x = "Maturity", y = "Frequency", title = expression(italic("R. metallica")))  +
  common_theme
setwd(FIGDIR)
jpeg('Maturity_Distribution.jpg', width=3000, height=1000, unit='px')
grid.arrange(hist1, hist2, hist3, hist4, hist5, ncol = 5)
dev.off()

setwd(FIGDIR)
jpeg('Maturity_vs_Foraging.jpg', width=6000, height=2000, unit='px')
par(mfrow = c(1,5), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(cfel_ForagingMaturity_SH$Foraging ~ cfel_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("C. fellah"))), xlab = "Maturity", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, ylim = c(0,1), col = cfel_ForagingMaturity_SH$col)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = cfel_ForagingMaturity_SH$Maturity, y = cfel_ForagingMaturity_SH$Foraging), lwd = 20)
plot(drug_ForagingMaturity_SH$Foraging ~ drug_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("D. rugosum"))), xlab = "Maturity", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, ylim = c(0,1), col = drug_ForagingMaturity_SH$col)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = drug_ForagingMaturity_SH$Maturity, y = drug_ForagingMaturity_SH$Foraging), lwd = 20)
plot(ipur_ForagingMaturity_SH$Foraging ~ ipur_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("I. purpureus"))), xlab = "Maturity", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, ylim = c(0,1), col = ipur_ForagingMaturity_SH$col)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = ipur_ForagingMaturity_SH$Maturity, y = ipur_ForagingMaturity_SH$Foraging), lwd = 20)
plot(pbar_ForagingMaturity_SH$Foraging ~ pbar_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("P. rugosus"))), xlab = "Maturity", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, ylim = c(0,1), col = pbar_ForagingMaturity_SH$col)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = pbar_ForagingMaturity_SH$Maturity, y = pbar_ForagingMaturity_SH$Foraging), lwd = 20)
plot(rmet_ForagingMaturity_SH$Foraging ~ rmet_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("R. metallica"))), xlab = "Maturity", ylab = "Foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, ylim = c(0,1), col = rmet_ForagingMaturity_SH$col)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(loess.smooth(x = rmet_ForagingMaturity_SH$Maturity, y = rmet_ForagingMaturity_SH$Foraging), lwd = 20)
dev.off()

jpeg('Maturity_vs_Strength.jpg', width=6000, height=2000, unit='px')
par(mfrow = c(1,5), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(cfel_ForagingMaturity_SH$Strength ~ cfel_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("C. fellah"))), xlab = "Maturity", ylab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, col = cfel_ForagingMaturity_SH$col, ylim = c(-2,8))
axis(2, at = c(-2,8), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(cfel_ForagingMaturity_SH$Strength~cfel_ForagingMaturity_SH$Maturity), lwd = 20)
plot(drug_ForagingMaturity_SH$Strength ~ drug_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("D. rugosum"))), xlab = "Maturity", ylab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, col = drug_ForagingMaturity_SH$col, ylim = c(-2,8))
axis(2, at = c(-2,8), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(drug_ForagingMaturity_SH$Strength~drug_ForagingMaturity_SH$Maturity), lwd = 20)
plot(ipur_ForagingMaturity_SH$Strength~ ipur_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("I. purpureus"))), xlab = "Maturity", ylab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, col = ipur_ForagingMaturity_SH$col, ylim = c(-2,8))
axis(2, at = c(-2,8), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(ipur_ForagingMaturity_SH$Strength~ipur_ForagingMaturity_SH$Maturity), lwd = 20)
plot(pbar_ForagingMaturity_SH$Strength~ pbar_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("P. rugosus"))), xlab = "Maturity", ylab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, col = pbar_ForagingMaturity_SH$col, ylim = c(-2,8))
axis(2, at = c(-2,8), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(pbar_ForagingMaturity_SH$Strength~pbar_ForagingMaturity_SH$Maturity), lwd = 20)
plot(rmet_ForagingMaturity_SH$Strength~ rmet_ForagingMaturity_SH$Maturity, pch = 16, main = substitute(paste(italic("R. metallica"))), xlab = "Maturity", ylab = "Normalized strength",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 8, col = rmet_ForagingMaturity_SH$col, ylim = c(-2,8))
axis(2, at = c(-2,8), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(rmet_ForagingMaturity_SH$Strength~rmet_ForagingMaturity_SH$Maturity), lwd = 20)
dev.off()

# Statistics for Foraging propensity and social maturity
Cfel_model <- lmer(Foraging ~ Maturity + (1 | col), data = cfel_ForagingMaturity)
Cfelmodel_summary <- coef(summary(Cfel_model))

Drug_model <- lmer(Foraging ~ Maturity + (1 | col), data = drug_ForagingMaturity)
Drugmodel_summary <- coef(summary(Drug_model))

Ipur_model <- lmer(Foraging ~ Maturity + (1 | col), data = ipur_ForagingMaturity)
Ipurmodel_summary <- coef(summary(Ipur_model))

Pbar_model <- lmer(Foraging ~ Maturity + (1 | col), data = pbar_ForagingMaturity)
Pbarmodel_summary <- coef(summary(Pbar_model))

Rmet_model <- lmer(Foraging ~ Maturity + (1 | col), data = rmet_ForagingMaturity)
Rmetmodel_summary <- coef(summary(Rmet_model))

# Statistics for interaction frequency and social maturity
Cfel_modeli <- lmer(Strength ~ Maturity + (1 | col), data = cfel_ForagingMaturity)
Cfelmodeli_summary <- coef(summary(Cfel_modeli))

Drug_modeli <- lmer(Strength ~ Maturity + (1 | col), data = drug_ForagingMaturity)
Drugmodeli_summary <- coef(summary(Drug_modeli))

Ipur_modeli <- lmer(Strength ~ Maturity + (1 | col), data = ipur_ForagingMaturity)
Ipurmodeli_summary <- coef(summary(Ipur_modeli))

Pbar_modeli <- lmer(Strength ~ Maturity + (1 | col), data = pbar_ForagingMaturity)
Pbarmodeli_summary <- coef(summary(Pbar_modeli))

Rmet_modeli <- lmer(Strength ~ Maturity + (1 | col), data = rmet_ForagingMaturity)
Rmetmodeli_summary <- coef(summary(Rmet_modeli))

# Quantify the extent of division of labor
DOL_DF <- data.frame(
  Camponotus = c(sd(softMod_list[[1]]$foraging), sd(softMod_list[[2]]$foraging), sd(softMod_list[[3]]$foraging), sd(softMod_list[[4]]$foraging), sd(softMod_list[[5]]$foraging)),
  Diacamma = c(sd(softMod_list[[11]]$foraging), sd(softMod_list[[12]]$foraging), sd(softMod_list[[13]]$foraging), sd(softMod_list[[14]]$foraging), sd(softMod_list[[15]]$foraging)),
  Iridomyrmex = c(sd(softMod_list[[16]]$foraging), sd(softMod_list[[17]]$foraging), sd(softMod_list[[18]]$foraging), sd(softMod_list[[19]]$foraging), sd(softMod_list[[20]]$foraging)),
  Pogonomyrmex = c(sd(softMod_list[[6]]$foraging), sd(softMod_list[[7]]$foraging), sd(softMod_list[[8]]$foraging), sd(softMod_list[[9]]$foraging), sd(softMod_list[[10]]$foraging)),
  Rhytidoponera = c(sd(softMod_list[[21]]$foraging), sd(softMod_list[[22]]$foraging), sd(softMod_list[[23]]$foraging), sd(softMod_list[[24]]$foraging), sd(softMod_list[[25]]$foraging)))

setwd(FIGDIR)
jpeg('Extent_of_DOL.jpg', width=3000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
boxplot(DOL_DF, pch = 16, ylab = "D. O. L.",
        yaxt="n", xaxt = "n", cex.lab = 10, cex.main = 12, cex = 5, ylim = c(0.25,0.5), xlab.cex = 5, main = "",
        lwd = 8, col = c("#FF7F50", "#808080", "#6A0DAD", "#50C878", "cornflowerblue"))
axis(2, at = c(0.25,0.5), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(1,2,3,4,5), labels = c(expression(italic("C. fellah")), expression(italic("D. rugosum")), expression(italic("I. purpureus")), expression(italic("P. rugosus")), expression(italic("R. metallica"))),
     tick = TRUE, lwd = 8, cex.axis = 6)
# Add a horizontal line
segments(1, 0.47, 3, 0.47, lwd = 5, col = "black")

# Add asterisk to indicate significance
text(2, 0.49, "*", cex = 8, col = "black")
dev.off()

# Construct a table of network features for PCA
MetaDF     <- as.data.frame(rep(colnames(DOL_DF), each = 5))
MetaDF$DOL <-  as.numeric(as.vector(as.matrix(DOL_DF)))
colnames(MetaDF) <- c("Species", "DOL")
rownames(MetaDF) <- c("C.fel_1", "C.fel_2", "C.fel_3", "C.fel_4", "C.fel_5",
                      "D.rug_1", "D.rug_2", "D.rug_3", "D.rug_4", "D.rug_5",
                      "I.pur_1", "I.pur_2", "I.pur_3", "I.pur_4", "I.pur_5",
                      "P.rug_1", "P.rug_2", "P.rug_3", "P.rug_4", "P.rug_5",
                      "R.met_1", "R.met_2", "R.met_3", "R.met_4", "R.met_5")

# Add soft modularity vals per species
cfelMod <- SoftModularity[1:5,2]
drugMod <- SoftModularity[6:10,2]
ipurMod <- SoftModularity[11:15,2]
pbarMod <- SoftModularity[16:20,2]
rmetMod <- SoftModularity[21:25,2]
MetaDF$Modularity <- c(cfelMod, drugMod, ipurMod,pbarMod, rmetMod)

t.test(cfelMod, drugMod)

# The proportion of individuals with extremal maturity scores
Extremal <- c()
for (i in 1:25){
  Extremal <- c(Extremal, sum(softMod_list[[i]]$foraging > 0.9 | softMod_list[[i]]$foraging < 0.1) / length(softMod_list[[i]]$foraging))
}
MetaDF$Extremal <- c(Extremal)

prop.test(sum(softMod_list[[1]]$foraging > 0.75 | softMod_list[[1]]$foraging < 0.25),
          sum(softMod_list[[1]]$foraging < 0.75 & softMod_list[[1]]$foraging > 0.25) + sum(softMod_list[[1]]$foraging > 0.75 | softMod_list[[1]]$foraging < 0.25),
          p=0.5, alternative="two.sided")

# The balance of individuals between the two social communities 
Balance <- c()
for (i in 1:25){
  Balance <- c(Balance, sum(softMod_list[[i]]$foraging > 0.9) / (sum(softMod_list[[i]]$foraging < 0.1) + sum(softMod_list[[i]]$foraging > 0.9)))
}
MetaDF$Balance <- c(Balance)

# The correlation between social maturity score and foraging
MatFor <- c()
for (i in 1:25){
  MatFor <- c(MatFor, cor(softMod_list[[i]]$foraging, foraging_list[[i]]))
}
MetaDF$MatFor <- c(MatFor)^2

# The correlation between social maturity score and node strength
MatStr <- c()
for (i in 1:25){
  MatStr <- c(MatStr, cor(softMod_list[[i]]$foraging, strengthNorm_list[[i]]))
}
MetaDF$MatStr <- c(MatStr)^2

# Total foraging effort for each colony
Foraging_MDF <- c()
for (i in 1:25){
  Foraging_MDF <- c(Foraging_MDF, mean(foraging_list[[i]]))
}
MetaDF$Foraging <- c(Foraging_MDF)

# Proportion of time FORAGERS spend foraging
foraging_list_foragers <- list()
foragingFreqForagers <- c()
for(i in 1:25){
  foraging_list_foragers[[i]] <- foraging_list[[i]][foraging_list[[i]]>0]
  foragingFreqForagers <- c(foragingFreqForagers, median(foraging_list[[i]][foraging_list[[i]]>0]))
  }

# Plotting the boxplots
setwd(FIGDIR)
jpeg('ForagersForagingFreq.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
boxplot(foraging_list_foragers, yaxt="n", xaxt="n",
        col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5),  # color of the boxes
        border = "black", # color of the box borders
        ylab = "Proportion of time foraging",
        xlab = "", cex.lab = 10, cex.main = 12, lwd = 5)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(1,25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(3, 8, 13, 18, 23), labels = c(expression(italic("C. fellah")), expression(italic("P. rugosus")), expression(italic("D. rugosum")), expression(italic("I. purpureus")), expression(italic("R. metallica"))),
     tick = TRUE, lwd = 8, cex.axis = 7)
dev.off()

# Total foragers per colony
Foragers_MDF <- c()
for (i in 1:25){
  Foragers_MDF <- c(Foragers_MDF, sum(foraging_list[[i]] > 0))
}
MetaDF$Foragers <- c(Foragers_MDF)

# Construct and plot PCA
MetaPCA <- prcomp(MetaDF[,-c(1,9)], center = TRUE,scale. = TRUE)
colours <- MetaDF[,1]
colours[colours == "Camponotus"] <- "#FF7F50"
colours[colours == "Iridomyrmex"] <- "#6A0DAD"
colours[colours == "Diacamma"] <- "#808080"
colours[colours == "Pogonomyrmex"] <- "#50C878"
colours[colours == "Rhytidoponera"] <- "cornflowerblue"

setwd(FIGDIR)
jpeg('PCA.jpg', width=2000, height=1000, unit='px')
autoplot(MetaPCA, loadings = TRUE, loadings.label = TRUE, loadings.size = 5, loadings.label.size = 15,
         loadings.colour = 'black', loadings.label.colour = 'black') +
  geom_point(size = 20, col = colours) +
  theme(text = element_text(size=50, family = "serif"))
dev.off()

jpeg('DOL_vs_Modularity.jpg', width=1200, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(MetaDF$DOL~MetaDF$Modularity, yaxt="n", xaxt="n", 
     col = colours, pch = 16, cex = 10, ylim = c(0.25,0.45), xlim = c(0.02,0.25),
     xlab = "Soft Modularity",
     ylab = "D.O.L", cex.lab = 7, cex.main = 12,)
axis(2, at = c(0.25,0.45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 8)
axis(1, at = c(0.02,0.25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 8)
abline(lm(MetaDF$DOL~MetaDF$Modularity), lwd=20, col = "darkgray")
dev.off()

# Statistics - relation between division of labor and modularity
cor.test(MetaDF$DOL, MetaDF$Modularity)

# Statistics - relation species differences in proportion of extremal
summary(aov(Extremal ~ Species, data = MetaDF))
summary(aov(Balance ~ Species, data = MetaDF))

# How does foraging output affect the correlation between strength and maturity
jpeg('MatrStr_vs_Foraging.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(MetaDF$MatStr~foragingFreqForagers, yaxt="n", xaxt="n", 
     col = colours, pch = 16, cex = 10, ylim = c(0,1), xlim = c(0,0.6),
     xlab = "Proportion of time foragers spend foraging",
     ylab = "strength ~ modularity", cex.lab = 10, cex.main = 12,)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,0.6), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(MetaDF$MatStr~foragingFreqForagers), lwd=20, col = "darkgray")
dev.off()

# Statistics on foraging output between species
summary(aov(Foraging ~ Species, data = MetaDF))

# Calculate and plot queen entropy and strength
Qrole_list <- list()
for(i in 1:25){
  Qrole_list[[i]] <- data.frame(
    ID = rep(NA,length(unique(c(edgelist_list[[i]]$id1, edgelist_list[[i]]$id2)))),
    Strength = rep(NA,length(unique(c(edgelist_list[[i]]$id1, edgelist_list[[i]]$id2)))),
    Entropy = rep(NA,length(unique(c(edgelist_list[[i]]$id1, edgelist_list[[i]]$id2)))),
    queen = rep("darkgray",length(unique(c(edgelist_list[[i]]$id1, edgelist_list[[i]]$id2))))
  )
  for(j in 1:length(unique(c(edgelist_list[[i]]$id1, edgelist_list[[i]]$id2)))){
    id <- unique(c(edgelist_list[[i]]$id1, edgelist_list[[i]]$id2))[j]
    weights <- edgelist_list[[i]][edgelist_list[[i]]$id1 == id | edgelist_list[[i]]$id2 == id,]$Count
    missing <- length(unique(c(edgelist_list[[i]]$id1, edgelist_list[[i]]$id2))) - length(weights)
    weights <- c(weights, rep(0, missing))
    Qrole_list[[i]]$ID[j] <- id
    Qrole_list[[i]]$Strength[j] <- sum(weights)
    Qrole_list[[i]]$Entropy[j] <- Entropy(weights)
    Qrole_list[[i]]$Specialization[j] <- sd(weights)/sum(weights)
  }
  Qrole_list[[i]]$queen[which(Qrole_list[[i]]$ID == queen_list[[i]])] <- "magenta"
  Qrole_list[[i]]$Strength <- scale(Qrole_list[[i]]$Strength) # Normalise such that all means are 0 and SDs are equal to 1
  Qrole_list[[i]]$Entropy <- scale(Qrole_list[[i]]$Entropy) # Normalise such that all means are 0 and SDs are equal to 1
}

jpeg('QueenPosition.jpg', width=5000, height=1000, unit='px')
par(mfrow = c(1,5), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(c(Qrole_list[[1]]$Strength, Qrole_list[[2]]$Strength, Qrole_list[[3]]$Strength, Qrole_list[[4]]$Strength, Qrole_list[[5]]$Strength), 
     c(Qrole_list[[1]]$Entropy, Qrole_list[[2]]$Entropy, Qrole_list[[3]]$Entropy, Qrole_list[[4]]$Entropy, Qrole_list[[5]]$Entropy),
     yaxt="n", xaxt="n", cex.lab = 10, pch = 16,
     col = "darkgray", cex = 13,
     ylab = "Entropy", xlab = "Strength", main = expression(italic("C. fellah")), cex.main = 10)
points(c(Qrole_list[[1]]$Strength[Qrole_list[[1]]$queen=="magenta"], Qrole_list[[2]]$Strength[Qrole_list[[2]]$queen=="magenta"], Qrole_list[[3]]$Strength[Qrole_list[[3]]$queen=="magenta"], Qrole_list[[4]]$Strength[Qrole_list[[4]]$queen=="magenta"], Qrole_list[[5]]$Strength[Qrole_list[[5]]$queen=="magenta"]), 
       c(Qrole_list[[1]]$Entropy[Qrole_list[[1]]$queen=="magenta"], Qrole_list[[2]]$Entropy[Qrole_list[[2]]$queen=="magenta"], Qrole_list[[3]]$Entropy[Qrole_list[[3]]$queen=="magenta"], Qrole_list[[4]]$Entropy[Qrole_list[[4]]$queen=="magenta"], Qrole_list[[5]]$Entropy[Qrole_list[[5]]$queen=="magenta"]),
       col = "magenta", cex = 13, pch = 16)
axis(1, at = c(-1,3), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(2, at = c(-3,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
plot(c(Qrole_list[[11]]$Strength, Qrole_list[[12]]$Strength, Qrole_list[[13]]$Strength, Qrole_list[[14]]$Strength, Qrole_list[[15]]$Strength), 
     c(Qrole_list[[11]]$Entropy, Qrole_list[[12]]$Entropy, Qrole_list[[13]]$Entropy, Qrole_list[[14]]$Entropy, Qrole_list[[15]]$Entropy),
     yaxt="n", xaxt="n", cex.lab = 10, pch = 16,
     col = "darkgray", cex = 13,
     ylab = "Entropy", xlab = "Strength", main = expression(italic("D. rugosum")), cex.main = 10)
points(c(Qrole_list[[11]]$Strength[Qrole_list[[11]]$queen=="magenta"], Qrole_list[[12]]$Strength[Qrole_list[[12]]$queen=="magenta"], Qrole_list[[13]]$Strength[Qrole_list[[13]]$queen=="magenta"], Qrole_list[[14]]$Strength[Qrole_list[[14]]$queen=="magenta"], Qrole_list[[15]]$Strength[Qrole_list[[15]]$queen=="magenta"]), 
       c(Qrole_list[[11]]$Entropy[Qrole_list[[11]]$queen=="magenta"], Qrole_list[[12]]$Entropy[Qrole_list[[12]]$queen=="magenta"], Qrole_list[[13]]$Entropy[Qrole_list[[13]]$queen=="magenta"], Qrole_list[[14]]$Entropy[Qrole_list[[14]]$queen=="magenta"], Qrole_list[[15]]$Entropy[Qrole_list[[15]]$queen=="magenta"]),
       col = "magenta", cex = 13, pch = 16)
axis(1, at = c(-1,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(2, at = c(-6,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
plot(c(Qrole_list[[16]]$Strength, Qrole_list[[17]]$Strength, Qrole_list[[18]]$Strength, Qrole_list[[19]]$Strength, Qrole_list[[20]]$Strength), 
     c(Qrole_list[[16]]$Entropy, Qrole_list[[17]]$Entropy, Qrole_list[[18]]$Entropy, Qrole_list[[19]]$Entropy, Qrole_list[[20]]$Entropy),
     yaxt="n", xaxt="n", cex.lab = 10, pch = 16,
     col = "darkgray", cex = 13,
     ylab = "Entropy", xlab = "Strength", main = expression(italic("I. purpureus")), cex.main = 10)
axis(1, at = c(-1,3), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(2, at = c(-7,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
points(c(Qrole_list[[16]]$Strength[Qrole_list[[16]]$queen=="magenta"], Qrole_list[[17]]$Strength[Qrole_list[[17]]$queen=="magenta"], Qrole_list[[18]]$Strength[Qrole_list[[18]]$queen=="magenta"], Qrole_list[[19]]$Strength[Qrole_list[[19]]$queen=="magenta"], Qrole_list[[20]]$Strength[Qrole_list[[20]]$queen=="magenta"]), 
       c(Qrole_list[[16]]$Entropy[Qrole_list[[16]]$queen=="magenta"], Qrole_list[[17]]$Entropy[Qrole_list[[17]]$queen=="magenta"], Qrole_list[[18]]$Entropy[Qrole_list[[18]]$queen=="magenta"], Qrole_list[[19]]$Entropy[Qrole_list[[19]]$queen=="magenta"], Qrole_list[[20]]$Entropy[Qrole_list[[20]]$queen=="magenta"]),
       col = "magenta", cex = 13, pch = 16)
plot(c(Qrole_list[[6]]$Strength, Qrole_list[[7]]$Strength, Qrole_list[[8]]$Strength, Qrole_list[[9]]$Strength, Qrole_list[[10]]$Strength), 
     c(Qrole_list[[6]]$Entropy, Qrole_list[[7]]$Entropy, Qrole_list[[8]]$Entropy, Qrole_list[[9]]$Entropy, Qrole_list[[10]]$Entropy),
     yaxt="n", xaxt="n", cex.lab = 10, pch = 16,
     col = "darkgray", cex = 13,
     ylab = "Entropy", xlab = "Strength", main = expression(italic("P. rugosus")), cex.main = 10)
points(c(Qrole_list[[6]]$Strength[Qrole_list[[6]]$queen=="magenta"], Qrole_list[[7]]$Strength[Qrole_list[[7]]$queen=="magenta"], Qrole_list[[8]]$Strength[Qrole_list[[8]]$queen=="magenta"], Qrole_list[[9]]$Strength[Qrole_list[[9]]$queen=="magenta"], Qrole_list[[10]]$Strength[Qrole_list[[10]]$queen=="magenta"]), 
       c(Qrole_list[[6]]$Entropy[Qrole_list[[6]]$queen=="magenta"], Qrole_list[[7]]$Entropy[Qrole_list[[7]]$queen=="magenta"], Qrole_list[[8]]$Entropy[Qrole_list[[8]]$queen=="magenta"], Qrole_list[[9]]$Entropy[Qrole_list[[9]]$queen=="magenta"], Qrole_list[[10]]$Entropy[Qrole_list[[10]]$queen=="magenta"]),
       col = "magenta", cex = 13, pch = 16)
axis(1, at = c(-1,7), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(2, at = c(-5,2), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
plot(c(Qrole_list[[21]]$Strength, Qrole_list[[22]]$Strength, Qrole_list[[23]]$Strength, Qrole_list[[24]]$Strength, Qrole_list[[25]]$Strength), 
     c(Qrole_list[[21]]$Entropy, Qrole_list[[22]]$Entropy, Qrole_list[[23]]$Entropy, Qrole_list[[24]]$Entropy, Qrole_list[[25]]$Entropy),
     yaxt="n", xaxt="n", cex.lab = 10, pch = 16,
     col = "darkgray", cex = 13,
     ylab = "Entropy", xlab = "Strength", main = expression(italic("R. metallica")), cex.main = 10)
points(c(Qrole_list[[21]]$Strength[Qrole_list[[21]]$queen=="magenta"], Qrole_list[[22]]$Strength[Qrole_list[[22]]$queen=="magenta"], Qrole_list[[23]]$Strength[Qrole_list[[23]]$queen=="magenta"], Qrole_list[[24]]$Strength[Qrole_list[[24]]$queen=="magenta"], Qrole_list[[25]]$Strength[Qrole_list[[25]]$queen=="magenta"]), 
       c(Qrole_list[[21]]$Entropy[Qrole_list[[21]]$queen=="magenta"], Qrole_list[[22]]$Entropy[Qrole_list[[22]]$queen=="magenta"], Qrole_list[[23]]$Entropy[Qrole_list[[23]]$queen=="magenta"], Qrole_list[[24]]$Entropy[Qrole_list[[24]]$queen=="magenta"], Qrole_list[[25]]$Entropy[Qrole_list[[25]]$queen=="magenta"]),
       col = "magenta", cex = 13, pch = 16)
axis(1, at = c(-1,4), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(2, at = c(-6,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
dev.off()

# Combine the data from the list into a single data frame
CFQ_data <- rbind(rbind(rbind(rbind(Qrole_list[[1]], Qrole_list[[2]]), Qrole_list[[3]]), Qrole_list[[4]]), Qrole_list[[5]]) 
PFQ_data <- rbind(rbind(rbind(rbind(Qrole_list[[6]], Qrole_list[[7]]), Qrole_list[[8]]), Qrole_list[[9]]), Qrole_list[[10]]) 
DFQ_data <- rbind(rbind(rbind(rbind(Qrole_list[[11]], Qrole_list[[12]]), Qrole_list[[13]]), Qrole_list[[14]]), Qrole_list[[15]]) 
IFQ_data <- rbind(rbind(rbind(rbind(Qrole_list[[16]], Qrole_list[[17]]), Qrole_list[[18]]), Qrole_list[[19]]), Qrole_list[[20]]) 
RFQ_data <- rbind(rbind(rbind(rbind(Qrole_list[[21]], Qrole_list[[22]]), Qrole_list[[23]]), Qrole_list[[24]]), Qrole_list[[25]]) 


# Create blank plot for spacing
blank_plot <- ggplot() +
  geom_blank() +
  theme_void()

# Create the entorpy boxplot
entropy_theme <- theme(legend.position = "none",
      axis.title.x = element_blank(),    # Set x-axis label appearance
      axis.text.x = element_blank(),     # Set x-axis text appearance
      axis.ticks.x = element_blank(),    # Set x-axis ticks appearance
      axis.line.x = element_blank(),
      axis.line.y = element_line(),
      axis.title.y = element_text(size = 30, angle = 90, family = "serif"),
      axis.text.y = element_text(size = 20, family = "serif"))
Cboxplot1 <- ggplot(CFQ_data, aes(x = factor(queen), y = Entropy, fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "", y = "Entropy") +
  scale_x_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  entropy_theme +
  annotate("segment", x = 1, xend = 2, y = 2, yend = 2, color = "black", size = 1.5) +
  annotate("text", x = 1.5, y = 2.3, label = "*", family = "serif", size = 10, color = "black", fontface = "bold")
Dboxplot1 <- ggplot(DFQ_data, aes(x = factor(queen), y = Entropy, fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "", y = "Entropy") +
  scale_x_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  entropy_theme +
  annotate("segment", x = 1, xend = 2, y = 2, yend = 2, color = "black", size = 1.5) +
  annotate("text", x = 1.5, y = 2.3, label = "*", family = "serif", size = 10, color = "black", fontface = "bold")
Iboxplot1 <- ggplot(IFQ_data, aes(x = factor(queen), y = Entropy, fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "", y = "Entropy") +
  scale_x_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  entropy_theme +
  annotate("segment", x = 1, xend = 2, y = 2, yend = 2, color = "black", size = 1.5) +
  annotate("text", x = 1.5, y = 2.3, label = "*", family = "serif", size = 10, color = "black", fontface = "bold")
Pboxplot1 <- ggplot(PFQ_data, aes(x = factor(queen), y = Entropy, fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "", y = "Entropy") +
  scale_x_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  entropy_theme
Rboxplot1 <- ggplot(RFQ_data, aes(x = factor(queen), y = Entropy, fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "", y = "Entropy") +
  scale_x_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  entropy_theme


# Create scatter plot
Cscatter <- ggplot(CFQ_data, aes(x = Strength, y = Entropy, color = queen)) +
  geom_point(size = 6) +
  labs(x = "Strength", y = "") +
  scale_color_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  theme(legend.position = "none") + 
  geom_point(data = subset(CFQ_data, queen == "magenta"), color = "magenta", size = 6) + # Add magenta points
  geom_blank(data = CFQ_data, aes(x = -Inf), show.legend = FALSE)  # Add blank layer for x-axis
Dscatter <- ggplot(DFQ_data, aes(x = Strength, y = Entropy, color = queen)) +
  geom_point(size = 6) +
  labs(x = "Strength", y = "") +
  scale_color_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  theme(legend.position = "none") + 
  geom_point(data = subset(DFQ_data, queen == "magenta"), color = "magenta", size = 6) + # Add magenta points
  geom_blank(data = CFQ_data, aes(x = -Inf), show.legend = FALSE)  # Add blank layer for x-axis
Iscatter <- ggplot(IFQ_data, aes(x = Strength, y = Entropy, color = queen)) +
  geom_point(size = 6) +
  labs(x = "Strength", y = "") +
  scale_color_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  theme(legend.position = "none") + 
  geom_point(data = subset(IFQ_data, queen == "magenta"), color = "magenta", size = 6) + # Add magenta points
  geom_blank(data = CFQ_data, aes(x = -Inf), show.legend = FALSE)  # Add blank layer for x-axis
Pscatter <- ggplot(PFQ_data, aes(x = Strength, y = Entropy, color = queen)) +
  geom_point(size = 6) +
  labs(x = "Strength", y = "") +
  scale_color_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  theme(legend.position = "none") + 
  geom_point(data = subset(PFQ_data, queen == "magenta"), color = "magenta", size = 6) + # Add magenta points
  geom_blank(data = CFQ_data, aes(x = -Inf), show.legend = FALSE)
Rscatter <- ggplot(RFQ_data, aes(x = Strength, y = Entropy, color = queen)) +
  geom_point(size = 6) +
  labs(x = "Strength", y = "") +
  scale_color_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  theme(legend.position = "none") + 
  geom_point(data = subset(RFQ_data, queen == "magenta"), color = "magenta", size = 6) + # Add magenta points
  geom_blank(data = CFQ_data, aes(x = -Inf), show.legend = FALSE)


#Create strength boxplots
strength_theme <- theme(legend.position = "none",
                          axis.title.y = element_blank(),    # Set x-axis label appearance
                          axis.text.y = element_blank(),     # Set x-axis text appearance
                          axis.ticks.y = element_blank(),    # Set x-axis ticks appearance
                          axis.line.y = element_blank(),
                          axis.line.x = element_line(),
                          axis.title.x = element_text(size = 30, family = "serif"),
                          axis.text.x = element_text(size = 20, family = "serif")) # Set x-axis line appearance
Cboxplot2<- ggplot(CFQ_data, aes(x = Strength, y = factor(queen), fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Strength", y = "") +
  scale_y_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  strength_theme
Dboxplot2<- ggplot(DFQ_data, aes(x = Strength, y = factor(queen), fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Strength", y = "") +
  scale_y_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  strength_theme
Iboxplot2<- ggplot(IFQ_data, aes(x = Strength, y = factor(queen), fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Strength", y = "") +
  scale_y_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  strength_theme
Pboxplot2<- ggplot(PFQ_data, aes(x = Strength, y = factor(queen), fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Strength", y = "") +
  scale_y_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  strength_theme +
  annotate("segment", x = 8, xend = 8, y = 1, yend = 2, color = "black", size = 1.5) +
  annotate("text", y = 1.5, x = 8.6, label = "**", family = "serif", size = 10, color = "black", fontface = "bold", angle = 90)
Rboxplot2<- ggplot(RFQ_data, aes(x = Strength, y = factor(queen), fill = queen)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Strength", y = "") +
  scale_y_discrete(labels = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  scale_fill_manual(values = c("darkgray" = "darkgray", "magenta" = "magenta")) +
  theme_void() +
  strength_theme

# Combine the plots
Cfinalplot <- (Cboxplot1 / Cscatter) + blank_plot + Cboxplot2 +
  plot_layout(ncol = 2, widths = c(0.2, 1), heights = c(1,0.2)) +
  plot_annotation(title = expression(italic("C. fellah")), theme = theme(plot.title = element_text(hjust = 0.5, size = 35, family = "serif")))
Dfinalplot <- (Dboxplot1 / Dscatter) + blank_plot + Dboxplot2 +
  plot_layout(ncol = 2, widths = c(0.2, 1), heights = c(1,0.2)) +
  plot_annotation(title = expression(italic("D. rugosum")), theme = theme(plot.title = element_text(hjust = 0.5, size = 35, family = "serif")))
Ifinalplot <- (Iboxplot1 / Iscatter) + blank_plot + Iboxplot2 +
  plot_layout(ncol = 2, widths = c(0.2, 1), heights = c(1,0.2)) +
  plot_annotation(title = expression(italic("I. purpureus")), theme = theme(plot.title = element_text(hjust = 0.5, size = 35, family = "serif")))
Pfinalplot <- (Pboxplot1 / Pscatter) + blank_plot + Pboxplot2 +
  plot_layout(ncol = 2, widths = c(0.2, 1), heights = c(1,0.2)) +
  plot_annotation(title = expression(italic("P. rugosus")), theme = theme(plot.title = element_text(hjust = 0.5, size = 35, family = "serif")))
Rfinalplot <- (Rboxplot1 / Rscatter) + blank_plot + Rboxplot2 +
  plot_layout(ncol = 2, widths = c(0.2, 1), heights = c(1,0.2)) +
  plot_annotation(title = expression(italic("R. metallica")), theme = theme(plot.title = element_text(hjust = 0.5, size = 35, family = "serif")))

AllQueenPlots <- (Cfinalplot / Dfinalplot / Ifinalplot / Pfinalplot / Rfinalplot) + plot_layout(ncol = 5) +
  plot_annotation(title = expression(italic(
    "C. fellah                        D. rugosum                        I. purpureus                        P. rugosus                        R. metallica")),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 35, family = "serif")))
ggsave("AllQueenPlots.jpg", AllQueenPlots, width = 25, height = 5, units = "in")



# Define the additional annotation
additional_annotation <- annotate("text", x = 0.5, y = 1.2, label = "Additional Annotation", size = 10, family = "serif")

# Combine plots and add the additional annotation
AllQueenPlots <- (Cfinalplot / Dfinalplot / Ifinalplot / Pfinalplot / Rfinalplot) +
  plot_layout(
    ncol = 5) +
  plot_annotation(title = expression(italic("C. fellah")), theme = theme(plot.title = element_text(hjust = 0.1, size = 35, family = "serif")))

# Add the additional annotation
AllQueenPlots <- AllQueenPlots + additional_annotation



# Compare queen and worker strengths and entropies for each species
t.test(c(Qrole_list[[1]]$Entropy[Qrole_list[[1]]$queen == "darkgray",],
         Qrole_list[[2]]$Entropy[Qrole_list[[2]]$queen == "darkgray",],
         Qrole_list[[3]]$Entropy[Qrole_list[[3]]$queen == "darkgray",],
         Qrole_list[[4]]$Entropy[Qrole_list[[4]]$queen == "darkgray",],
         Qrole_list[[5]]$Entropy[Qrole_list[[5]]$queen == "darkgray",]),
       c(Qrole_list[[1]]$Entropy[Qrole_list[[1]]$queen == "magenta",],
         Qrole_list[[2]]$Entropy[Qrole_list[[2]]$queen == "magenta",],
         Qrole_list[[3]]$Entropy[Qrole_list[[3]]$queen == "magenta",],
         Qrole_list[[4]]$Entropy[Qrole_list[[4]]$queen == "magenta",],
         Qrole_list[[5]]$Entropy[Qrole_list[[5]]$queen == "magenta",]))
t.test(c(Qrole_list[[1]]$Strength[Qrole_list[[1]]$queen == "darkgray",],
         Qrole_list[[2]]$Strength[Qrole_list[[2]]$queen == "darkgray",],
         Qrole_list[[3]]$Strength[Qrole_list[[3]]$queen == "darkgray",],
         Qrole_list[[4]]$Strength[Qrole_list[[4]]$queen == "darkgray",],
         Qrole_list[[5]]$Strength[Qrole_list[[5]]$queen == "darkgray",]),
       c(Qrole_list[[1]]$Strength[Qrole_list[[1]]$queen == "magenta",],
         Qrole_list[[2]]$Strength[Qrole_list[[2]]$queen == "magenta",],
         Qrole_list[[3]]$Strength[Qrole_list[[3]]$queen == "magenta",],
         Qrole_list[[4]]$Strength[Qrole_list[[4]]$queen == "magenta",],
         Qrole_list[[5]]$Strength[Qrole_list[[5]]$queen == "magenta",]))
t.test(c(Qrole_list[[11]]$Entropy[Qrole_list[[11]]$queen == "darkgray",],
         Qrole_list[[12]]$Entropy[Qrole_list[[12]]$queen == "darkgray",],
         Qrole_list[[13]]$Entropy[Qrole_list[[13]]$queen == "darkgray",],
         Qrole_list[[14]]$Entropy[Qrole_list[[14]]$queen == "darkgray",],
         Qrole_list[[15]]$Entropy[Qrole_list[[15]]$queen == "darkgray",]),
       c(Qrole_list[[11]]$Entropy[Qrole_list[[11]]$queen == "magenta",],
         Qrole_list[[12]]$Entropy[Qrole_list[[12]]$queen == "magenta",],
         Qrole_list[[13]]$Entropy[Qrole_list[[13]]$queen == "magenta",],
         Qrole_list[[14]]$Entropy[Qrole_list[[14]]$queen == "magenta",],
         Qrole_list[[15]]$Entropy[Qrole_list[[15]]$queen == "magenta",]))
t.test(c(Qrole_list[[11]]$Strength[Qrole_list[[11]]$queen == "darkgray",],
         Qrole_list[[12]]$Strength[Qrole_list[[12]]$queen == "darkgray",],
         Qrole_list[[13]]$Strength[Qrole_list[[13]]$queen == "darkgray",],
         Qrole_list[[14]]$Strength[Qrole_list[[14]]$queen == "darkgray",],
         Qrole_list[[15]]$Strength[Qrole_list[[15]]$queen == "darkgray",]),
       c(Qrole_list[[11]]$Strength[Qrole_list[[11]]$queen == "magenta",],
         Qrole_list[[12]]$Strength[Qrole_list[[12]]$queen == "magenta",],
         Qrole_list[[13]]$Strength[Qrole_list[[13]]$queen == "magenta",],
         Qrole_list[[14]]$Strength[Qrole_list[[14]]$queen == "magenta",],
         Qrole_list[[15]]$Strength[Qrole_list[[15]]$queen == "magenta",]))
t.test(c(Qrole_list[[16]]$Entropy[Qrole_list[[16]]$queen == "darkgray",],
         Qrole_list[[17]]$Entropy[Qrole_list[[17]]$queen == "darkgray",],
         Qrole_list[[18]]$Entropy[Qrole_list[[18]]$queen == "darkgray",],
         Qrole_list[[19]]$Entropy[Qrole_list[[19]]$queen == "darkgray",],
         Qrole_list[[20]]$Entropy[Qrole_list[[20]]$queen == "darkgray",]),
       c(Qrole_list[[16]]$Entropy[Qrole_list[[16]]$queen == "magenta",],
         Qrole_list[[17]]$Entropy[Qrole_list[[17]]$queen == "magenta",],
         Qrole_list[[18]]$Entropy[Qrole_list[[18]]$queen == "magenta",],
         Qrole_list[[19]]$Entropy[Qrole_list[[19]]$queen == "magenta",],
         Qrole_list[[20]]$Entropy[Qrole_list[[20]]$queen == "magenta",]))
t.test(c(Qrole_list[[16]]$Strength[Qrole_list[[16]]$queen == "darkgray",],
         Qrole_list[[17]]$Strength[Qrole_list[[17]]$queen == "darkgray",],
         Qrole_list[[18]]$Strength[Qrole_list[[18]]$queen == "darkgray",],
         Qrole_list[[19]]$Strength[Qrole_list[[19]]$queen == "darkgray",],
         Qrole_list[[20]]$Strength[Qrole_list[[20]]$queen == "darkgray",]),
       c(Qrole_list[[16]]$Strength[Qrole_list[[16]]$queen == "magenta",],
         Qrole_list[[17]]$Strength[Qrole_list[[17]]$queen == "magenta",],
         Qrole_list[[18]]$Strength[Qrole_list[[18]]$queen == "magenta",],
         Qrole_list[[19]]$Strength[Qrole_list[[19]]$queen == "magenta",],
         Qrole_list[[20]]$Strength[Qrole_list[[20]]$queen == "magenta",]))
t.test(c(Qrole_list[[6]]$Entropy[Qrole_list[[6]]$queen == "darkgray",],
         Qrole_list[[7]]$Entropy[Qrole_list[[7]]$queen == "darkgray",],
         Qrole_list[[8]]$Entropy[Qrole_list[[8]]$queen == "darkgray",],
         Qrole_list[[9]]$Entropy[Qrole_list[[9]]$queen == "darkgray",],
         Qrole_list[[10]]$Entropy[Qrole_list[[10]]$queen == "darkgray",]),
       c(Qrole_list[[6]]$Entropy[Qrole_list[[6]]$queen == "magenta",],
         Qrole_list[[7]]$Entropy[Qrole_list[[7]]$queen == "magenta",],
         Qrole_list[[8]]$Entropy[Qrole_list[[8]]$queen == "magenta",],
         Qrole_list[[9]]$Entropy[Qrole_list[[9]]$queen == "magenta",],
         Qrole_list[[10]]$Entropy[Qrole_list[[10]]$queen == "magenta",]))
t.test(c(Qrole_list[[6]]$Strength[Qrole_list[[6]]$queen == "darkgray",],
         Qrole_list[[7]]$Strength[Qrole_list[[7]]$queen == "darkgray",],
         Qrole_list[[8]]$Strength[Qrole_list[[8]]$queen == "darkgray",],
         Qrole_list[[9]]$Strength[Qrole_list[[9]]$queen == "darkgray",],
         Qrole_list[[10]]$Strength[Qrole_list[[10]]$queen == "darkgray",]),
       c(Qrole_list[[6]]$Strength[Qrole_list[[6]]$queen == "magenta",],
         Qrole_list[[7]]$Strength[Qrole_list[[7]]$queen == "magenta",],
         Qrole_list[[8]]$Strength[Qrole_list[[8]]$queen == "magenta",],
         Qrole_list[[9]]$Strength[Qrole_list[[9]]$queen == "magenta",],
         Qrole_list[[10]]$Strength[Qrole_list[[10]]$queen == "magenta",]))
t.test(c(Qrole_list[[21]]$Entropy[Qrole_list[[21]]$queen == "darkgray",],
         Qrole_list[[22]]$Entropy[Qrole_list[[22]]$queen == "darkgray",],
         Qrole_list[[23]]$Entropy[Qrole_list[[23]]$queen == "darkgray",],
         Qrole_list[[24]]$Entropy[Qrole_list[[24]]$queen == "darkgray",],
         Qrole_list[[25]]$Entropy[Qrole_list[[25]]$queen == "darkgray",]),
       c(Qrole_list[[21]]$Entropy[Qrole_list[[21]]$queen == "magenta",],
         Qrole_list[[22]]$Entropy[Qrole_list[[22]]$queen == "magenta",],
         Qrole_list[[23]]$Entropy[Qrole_list[[23]]$queen == "magenta",],
         Qrole_list[[24]]$Entropy[Qrole_list[[24]]$queen == "magenta",],
         Qrole_list[[25]]$Entropy[Qrole_list[[25]]$queen == "magenta",]))
t.test(c(Qrole_list[[21]]$Strength[Qrole_list[[21]]$queen == "darkgray",],
         Qrole_list[[22]]$Strength[Qrole_list[[22]]$queen == "darkgray",],
         Qrole_list[[23]]$Strength[Qrole_list[[23]]$queen == "darkgray",],
         Qrole_list[[24]]$Strength[Qrole_list[[24]]$queen == "darkgray",],
         Qrole_list[[25]]$Strength[Qrole_list[[25]]$queen == "darkgray",]),
       c(Qrole_list[[21]]$Strength[Qrole_list[[21]]$queen == "magenta",],
         Qrole_list[[22]]$Strength[Qrole_list[[22]]$queen == "magenta",],
         Qrole_list[[23]]$Strength[Qrole_list[[23]]$queen == "magenta",],
         Qrole_list[[24]]$Strength[Qrole_list[[24]]$queen == "magenta",],
         Qrole_list[[25]]$Strength[Qrole_list[[25]]$queen == "magenta",]))

# Quantify queen space-use - nothing interesting here
spatialEntropy <- c()
HomeRange <- c()
for (i in 1:25){
  spatialEntropy <- c(spatialEntropy, Entropy(as.numeric(nest_list[[i]][rownames(nest_list[[i]]) == queen_list[[i]],])))
  mat <- rev(sort(as.numeric(nest_list[[i]][rownames(nest_list[[i]]) == queen_list[[i]],])))
  range <- c()
  for(j in 1:length(mat)){
    if(sum(range) < sum(mat)*0.9){
      range <- c(range, mat[j]) 
    }}
  HomeRange <- c(HomeRange, length(range)/length(mat))
}

socialEntropy <- c()
for (i in 1:25){
  socialEntropy <- c(socialEntropy, Qrole_list[[i]]$Entropy[Qrole_list[[i]]$queen == "magenta",])
}
socialStrength <- c()
for (i in 1:25){
  socialStrength <- c(socialStrength, Qrole_list[[i]]$Strength[Qrole_list[[i]]$queen == "magenta",])
}

# Quantify interaction profile over time
interaction_Q_list <- list()
for (i in 1:25){
interaction_Q_list[[i]] <- interaction_list[[i]][interaction_list[[i]]$id1 == queen_list[[i]] | interaction_list[[i]]$id2 == queen_list[[i]],]
}

# Bin interactions by hour. IMPORTANTLY USING ALL CONTACTS; NOT ONLY HH HERE
for (i in 1:25){
  interaction_Q_list[[i]]$hour <- round.POSIXt(as.POSIXct(interaction_Q_list[[i]]$start, format = "%Y-%m-%dT%H:%M:%OS"), units = "hours")
  interaction_Q_list[[i]]$hour <- as.numeric((interaction_Q_list[[i]]$hour - interaction_Q_list[[i]]$hour[1])/3600)
  interaction_Q_list[[i]]$min <- round.POSIXt(as.POSIXct(interaction_Q_list[[i]]$start, format = "%Y-%m-%dT%H:%M:%OS"), units = "mins")
  interaction_Q_list[[i]]$min <- as.numeric((interaction_Q_list[[i]]$min - interaction_Q_list[[i]]$min[1])/60)
  print(i)
}

hourSplits <- list()
minSplits <- list()
interactionTotHour_list <- list()
interactionTotIDHour_list <- list()
IDconservationHour_list <- list()
interactionTotMin_list <- list()
interactionTotIDMin_list <- list()
IDconservationMin_list <- list()
WeightedMeanMod_list <- list()
for (i in 1:25){
  hourSplits[[i]] <- split(interaction_Q_list[[i]], as.character(interaction_Q_list[[i]]$hour))
  minSplits[[i]] <- split(interaction_Q_list[[i]], as.character(interaction_Q_list[[i]]$min))
  interactionTotHour <- c()
  interactionTotIDHour <- c()
  IDconservationHour <- c()
  interactionTotMin <- c()
  interactionTotIDMin <- c()
  IDconservationMin <- c()
  WeightedMeanMod <- c()
  for (j in 1:length(hourSplits[[i]])){
    interactionTotHour <- c(interactionTotHour, nrow(hourSplits[[i]][[j]]))
    interactionTotIDHour <- c(interactionTotIDHour, length(unique(c(hourSplits[[i]][[j]]$id1, hourSplits[[i]][[j]]$id2))) - 1)
    table <- table(c(hourSplits[[i]][[j]]$id1, hourSplits[[i]][[j]]$id2))
    socmodtable <- data.frame(counts = table,
                              maturity = NA)
    for(k in names(table)){
      if(k %in% softMod_list[[i]]$id){
        socmodtable[socmodtable$counts.Var1 == k, 3] <- softMod_list[[i]][softMod_list[[i]]$id == k,]$foraging
      }
    }
    weightedmeanmod <- sum(socmodtable$counts.Freq*socmodtable$maturity)/sum(socmodtable$counts.Freq)
    WeightedMeanMod <- c(WeightedMeanMod, weightedmeanmod)
    if (j < length(hourSplits[[i]])){
      # Calculate the conservation from one hour to the next as a proportion of the total number of individuals across both hours
      IDconservationHour <- c(IDconservationHour, length(intersect(unique(c(hourSplits[[i]][[j]]$id1, hourSplits[[i]][[j]]$id2)), unique(c(hourSplits[[i]][[j+1]]$id1, hourSplits[[i]][[j+1]]$id2)))) /
        length(unique(c(unique(c(hourSplits[[i]][[j]]$id1, hourSplits[[i]][[j]]$id2)), unique(c(hourSplits[[i]][[j+1]]$id1, hourSplits[[i]][[j+1]]$id2))))))
    }
  }
  for (j in 1:length(minSplits[[i]])){
    interactionTotMin <- c(interactionTotMin, nrow(minSplits[[i]][[j]]))
    interactionTotIDMin <- c(interactionTotIDMin, length(unique(c(minSplits[[i]][[j]]$id1, minSplits[[i]][[j]]$id2))) - 1)
    if (j < length(minSplits[[i]])){
      # Calculate the conservation from one hour to the next as a proportion of the total number of individuals across both hours
      IDconservationMin <- c(IDconservationMin, length(intersect(unique(c(minSplits[[i]][[j]]$id1, minSplits[[i]][[j]]$id2)), unique(c(minSplits[[i]][[j+1]]$id1, minSplits[[i]][[j+1]]$id2)))) /
                              length(unique(c(unique(c(minSplits[[i]][[j]]$id1, minSplits[[i]][[j]]$id2)), unique(c(minSplits[[i]][[j+1]]$id1, minSplits[[i]][[j+1]]$id2))))))
    }
  }
  interactionTotHour_list[[i]] <- interactionTotHour
  interactionTotIDHour_list[[i]] <- interactionTotIDHour
  IDconservationHour_list[[i]] <- IDconservationHour
  interactionTotMin_list[[i]] <- interactionTotMin
  interactionTotIDMin_list[[i]] <- interactionTotIDMin
  IDconservationMin_list[[i]] <- IDconservationMin
  WeightedMeanMod_list[[i]] <- WeightedMeanMod
  print(i)
}

contactsPerWorker <- c() 
for (i in 1:25){
  contactsPerWorker <- c(contactsPerWorker, mean(interactionTotHour_list[[i]] / interactionTotIDHour_list[[i]]))
}

# Calculate the average number of individuals queens of each species interact with per hour
hourlyQinteractionsVec <- c()
for(i in 1:25){
  hourlyQinteractionsVec <- c(hourlyQinteractionsVec, mean(interactionTotIDHour_list[[i]]/(nrow(md_list[[i]])-1)))
}

minlyQinteractionsVec <- c()
for(i in 1:25){
  minlyQinteractionsVec <- c(minlyQinteractionsVec, mean(interactionTotIDMin_list[[i]]/(nrow(md_list[[i]])-1)))
}

# Calculate average conservation per hour
hourlyQConservationVec <- c()
for(i in 1:25){
  hourlyQConservationVec <- c(hourlyQConservationVec, mean(IDconservationHour_list[[i]]))
}

#Calculate average weighted mean partner maturity
WeightedMeanModVec <- c()
for(i in 1:25){
  WeightedMeanModVec <- c(WeightedMeanModVec, mean(na.omit(WeightedMeanMod_list[[i]])))
}

jpeg('QueenDynamics.jpg', width = 1200, height = 1200, units = 'px')
par(mfrow = c(4, 2), mar = c(5, 5, 4, 2) + 0.1, oma = c(2, 2, 0, 0), bty = "n", mgp = c(5, 1, 0), family = "serif")
plot(HomeRange ~ socialStrength, pch = 16, cex = 3,
     col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5), ylab = "",
     cex.axis = 1.5, cex.lab = 2, xlab = "")
mtext("90% home range size", side = 2, line = 3, cex = 1.5)
mtext("Social strength", side = 1, line = 3, cex = 1.5)
abline(lm(HomeRange ~ socialStrength))

plot(HomeRange ~ socialEntropy, pch = 16, cex = 3,
     col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5), ylab = "",
     cex.axis = 1.5, cex.lab = 2, xlab = "")
mtext("90% home range size", side = 2, line = 3, cex = 1.5)
mtext("Social entropy", side = 1, line = 3, cex = 1.5)
abline(lm(HomeRange ~ socialEntropy))

plot(hourlyQinteractionsVec ~ socialStrength, pch = 16, cex = 3,
     col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5), ylab = "",
     cex.axis = 1.5, cex.lab = 2, xlab = "")
mtext("Workers contacted per hour", side = 2, line = 3, cex = 1.5)
mtext("Social strength", side = 1, line = 3, cex = 1.5)
abline(lm(hourlyQinteractionsVec ~ socialStrength), col = "black")

plot(hourlyQinteractionsVec ~ socialEntropy, pch = 16, cex = 3,
     col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5), ylab = "",
     cex.axis = 1.5, cex.lab = 2, xlab = "")
mtext("Workers contacted per hour", side = 2, line = 3, cex = 1.5)
mtext("Social entropy", side = 1, line = 3, cex = 1.5)
abline(lm(hourlyQinteractionsVec ~ socialEntropy), col = "black")

plot(contactsPerWorker ~ socialStrength, pch = 16, cex = 3,
     col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5), ylab = "",
     cex.axis = 1.5, cex.lab = 2, xlab = "")
mtext("Contacts per worker", side = 2, line = 3, cex = 1.5)
mtext("Social strength", side = 1, line = 3, cex = 1.5)
abline(lm(contactsPerWorker ~ socialStrength), col = "black")

plot(contactsPerWorker ~ socialEntropy, pch = 16, cex = 3,
     col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5), ylab = "",
     cex.axis = 1.5, cex.lab = 2, xlab = "")
mtext("Contacts per worker", side = 2, line = 3, cex = 1.5)
mtext("Social entropy", side = 1, line = 3, cex = 1.5)
abline(lm(contactsPerWorker ~ socialEntropy))

plot(hourlyQConservationVec ~ socialStrength, pch = 16, cex = 3,
     col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5), ylab = "",
     cex.axis = 1.5, cex.lab = 2, xlab = "")
mtext("Worker identities conserved between hours", side = 2, line = 3, cex = 1.5)
mtext("Social strength", side = 1, line = 3, cex = 1.5)
abline(lm(hourlyQConservationVec ~ socialStrength), col = "black")

plot(hourlyQConservationVec ~ socialEntropy, pch = 16, cex = 3,
     col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5), ylab = "",
     cex.axis = 1.5, cex.lab = 2, xlab = "")
mtext("Worker identities conserved between hours", side = 2, line = 3, cex = 1.5)
mtext("Social entropy", side = 1, line = 3, cex = 1.5)
abline(lm(hourlyQConservationVec ~ socialEntropy))
dev.off()

# How large are the arenas, and does area size correlate with spatial entropy
ArenaSize <- c()
for(i in 1:25){
  ArenaSize <- c(ArenaSize, ncol(nest_list[[i]]))
}
cor.test(HomeRange, ArenaSize)

# Statistics for understanding drivers of queen strength and entropy
QueenAllDF <- data.frame(HR = HomeRange,
                         socialEnt = socialEntropy,
                         socialStr = socialStrength,
                         species = MetaDF$Species,
                         ContPerHour = hourlyQinteractionsVec,
                         Conserv = hourlyQConservationVec,
                         PerWorker = contactsPerWorker,
                         MeanMod = WeightedMeanModVec)

# How does home range size relate to entropy and strength
QueenHR_ent_model <- lmer(socialEnt ~ HR + (1|species), data = QueenAllDF)
QueenHR_ent_summary <- coef(summary(QueenHR_ent_model))

QueenHR_str_model <- lmer(socialStr ~ HR + (1|species), data = QueenAllDF)
QueenHR_str_summary <- coef(summary(QueenHR_str_model))

# How does proportion of workers contacted per hour relate to entropy and strength
QueenCPH_ent_model <- lmer(socialEnt ~ ContPerHour + (1|species), data = QueenAllDF)
QueenCPH_ent_summary <- coef(summary(QueenCPH_ent_model))

QueenCPH_str_model <- lmer(socialStr ~ ContPerHour + (1|species), data = QueenAllDF)
QueenCPH_str_summary <- coef(summary(QueenCPH_str_model))

# How does hourly turn-over rate relate to entropy and strength
QueenHTO_ent_model <- lmer(socialEnt ~ Conserv + (1|species), data = QueenAllDF)
QueenHTO_ent_summary <- coef(summary(QueenHTO_ent_model))

QueenHTO_str_model <- lmer(socialStr ~ ContPerHour + (1|species), data = QueenAllDF)
QueenHTO_str_summary <- coef(summary(QueenHTO_str_model))

# How does contacts per worker relate to entropy and strength
QueenPW_ent_model <- lmer(socialEnt ~ contactsPerWorker + (1|species), data = QueenAllDF)
QueenPW_ent_summary <- coef(summary(QueenPW_ent_model))

QueenPW_str_model <- lmer(socialStr ~ contactsPerWorker + (1|species), data = QueenAllDF)
QueenPW_str_summary <- coef(summary(QueenPW_str_model))


plot(QueenAllDF$socialEnt, QueenAllDF$MeanMod,
     pch = 16, xlab = "Social entropy", ylab = "Weighted mean maturity of interaction partners",
    col = rep(c("#FF7F50", "#50C878", "#808080", "#6A0DAD", "cornflowerblue"), each = 5))

summary(lm(QueenAllDF$PerWorker ~ QueenAllDF$MeanMod))
EntMod_model <- lmer(socialEnt ~ MeanMod + (1 | species), data = QueenAllDF)
EntMod_summary <- coef(summary(EntMod_model))

range(QueenAllDF$HR)

# Calculate the homerange size of each worker
WorkerHomeRange_list <- list()
for (i in 1:25){
  WorkerHomeRangeVec <- c()
  for(j in 1:nrow(nest_list[[i]])){
    mat <- rev(sort(as.numeric(nest_list[[i]][j,])))
    range <- c()
    for(k in 1:length(mat)){
      if(sum(range) < sum(mat)*0.9){
        range <- c(range, mat[k]) 
      }}
    WorkerHomeRangeVec <- c(WorkerHomeRangeVec, length(range)/length(mat))
    }
    WorkerHomeRange_list[[i]] <- WorkerHomeRangeVec
    print(i)
}

# Compare homerange size with social maturity
HR_Mod_list <- list()
for(i in 1:25){
  df <- data.frame(id = rownames(nest_list[[i]]),
             HomeRange = WorkerHomeRange_list[[i]],
             rep = i)
  HR_Mod_list[[i]] <- merge(softMod_list[[i]], df, by = "id")
}

HR_Mod_cfel <- rbind(HR_Mod_list[[1]], HR_Mod_list[[2]], HR_Mod_list[[3]], HR_Mod_list[[4]], HR_Mod_list[[5]])
HR_Mod_cfel$rep[HR_Mod_cfel$rep == 1] <- "midnightblue"
HR_Mod_cfel$rep[HR_Mod_cfel$rep == 2] <- "sandybrown"
HR_Mod_cfel$rep[HR_Mod_cfel$rep == 3] <- "slategray"
HR_Mod_cfel$rep[HR_Mod_cfel$rep == 4] <- "red3"
HR_Mod_cfel$rep[HR_Mod_cfel$rep == 5] <- "darkcyan"

HR_Mod_drug <- rbind(HR_Mod_list[[6]], HR_Mod_list[[7]], HR_Mod_list[[8]], HR_Mod_list[[9]], HR_Mod_list[[10]])
HR_Mod_drug$rep[HR_Mod_drug$rep == 6] <- "midnightblue"
HR_Mod_drug$rep[HR_Mod_drug$rep == 7] <- "sandybrown"
HR_Mod_drug$rep[HR_Mod_drug$rep == 8] <- "slategray"
HR_Mod_drug$rep[HR_Mod_drug$rep == 9] <- "red3"
HR_Mod_drug$rep[HR_Mod_drug$rep == 10] <- "darkcyan"

HR_Mod_ipur <- rbind(HR_Mod_list[[11]], HR_Mod_list[[12]], HR_Mod_list[[13]], HR_Mod_list[[14]], HR_Mod_list[[15]])
HR_Mod_ipur$rep[HR_Mod_ipur$rep == 11] <- "midnightblue"
HR_Mod_ipur$rep[HR_Mod_ipur$rep == 12] <- "sandybrown"
HR_Mod_ipur$rep[HR_Mod_ipur$rep == 13] <- "slategray"
HR_Mod_ipur$rep[HR_Mod_ipur$rep == 14] <- "red3"
HR_Mod_ipur$rep[HR_Mod_ipur$rep == 15] <- "darkcyan"

HR_Mod_prug <- rbind(HR_Mod_list[[16]], HR_Mod_list[[17]], HR_Mod_list[[18]], HR_Mod_list[[19]], HR_Mod_list[[20]])
HR_Mod_prug$rep[HR_Mod_prug$rep == 16] <- "midnightblue"
HR_Mod_prug$rep[HR_Mod_prug$rep == 17] <- "sandybrown"
HR_Mod_prug$rep[HR_Mod_prug$rep == 18] <- "slategray"
HR_Mod_prug$rep[HR_Mod_prug$rep == 19] <- "red3"
HR_Mod_prug$rep[HR_Mod_prug$rep == 20] <- "darkcyan"

HR_Mod_rmet <- rbind(HR_Mod_list[[21]], HR_Mod_list[[22]], HR_Mod_list[[23]], HR_Mod_list[[24]], HR_Mod_list[[25]])
HR_Mod_rmet$rep[HR_Mod_rmet$rep == 21] <- "midnightblue"
HR_Mod_rmet$rep[HR_Mod_rmet$rep == 22] <- "sandybrown"
HR_Mod_rmet$rep[HR_Mod_rmet$rep == 23] <- "slategray"
HR_Mod_rmet$rep[HR_Mod_rmet$rep == 24] <- "red3"
HR_Mod_rmet$rep[HR_Mod_rmet$rep == 25] <- "darkcyan"

# Shuffle dataframes for plotting
HR_Mod_cfel_SH <- HR_Mod_cfel[sample(1:nrow(HR_Mod_cfel)),]
HR_Mod_drug_SH <- HR_Mod_drug[sample(1:nrow(HR_Mod_drug)),]
HR_Mod_ipur_SH <- HR_Mod_ipur[sample(1:nrow(HR_Mod_ipur)),]
HR_Mod_prug_SH <- HR_Mod_prug[sample(1:nrow(HR_Mod_prug)),]
HR_Mod_rmet_SH <- HR_Mod_rmet[sample(1:nrow(HR_Mod_rmet)),]

setwd(FIGDIR)
jpeg('Maturity_vs_HomeRange.jpg', width=6000, height=2000, unit='px')
par(mfrow = c(1,5), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(HR_Mod_cfel_SH$HomeRange ~ HR_Mod_cfel_SH$foraging,
     pch = 16, ylab = "90% home range size", xlab = "Maturity", col = HR_Mod_cfel_SH$rep,
     yaxt="n", xaxt="n", cex.lab = 10, cex = 8, ylim = c(0, 0.6))
axis(2, at = c(0,0.6), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(HR_Mod_cfel_SH$HomeRange ~ HR_Mod_cfel_SH$foraging), lwd = 20)
plot(HR_Mod_drug_SH$HomeRange ~ HR_Mod_drug_SH$foraging,
     pch = 16, ylab = "90% home range size", xlab = "Maturity", col = HR_Mod_drug_SH$rep,
     yaxt="n", xaxt="n", cex.lab = 10, cex = 8, ylim = c(0, 0.6))
axis(2, at = c(0,0.6), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(HR_Mod_drug_SH$HomeRange ~ HR_Mod_drug_SH$foraging), lwd = 20)
plot(HR_Mod_ipur_SH$HomeRange ~ HR_Mod_ipur_SH$foraging,
     pch = 16, ylab = "90% home range size", xlab = "Maturity", col = HR_Mod_ipur_SH$rep,
     yaxt="n", xaxt="n", cex.lab = 10, cex = 8, ylim = c(0, 0.6))
axis(2, at = c(0,0.6), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(HR_Mod_ipur_SH$HomeRange ~ HR_Mod_ipur_SH$foraging), lwd = 20)
plot(HR_Mod_prug_SH$HomeRange ~ HR_Mod_prug_SH$foraging,
     pch = 16, ylab = "90% home range size", xlab = "Maturity", col = HR_Mod_prug_SH$rep,
     yaxt="n", xaxt="n", cex.lab = 10, cex = 8, ylim = c(0, 0.6))
axis(2, at = c(0,0.6), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(HR_Mod_prug_SH$HomeRange ~ HR_Mod_prug_SH$foraging), lwd = 20)
plot(HR_Mod_rmet_SH$HomeRange ~ HR_Mod_rmet_SH$foraging,
     pch = 16, ylab = "90% home range size", xlab = "Maturity", col = HR_Mod_rmet_SH$rep,
     yaxt="n", xaxt="n", cex.lab = 10, cex = 8, ylim = c(0, 0.6))
axis(2, at = c(0,0.6), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
abline(lm(HR_Mod_rmet_SH$HomeRange ~ HR_Mod_rmet_SH$foraging), lwd = 20)
dev.off()

HR_Mod_cfel$species = "cfel"
HR_Mod_drug$species = "drug"
HR_Mod_ipur$species = "ipur"
HR_Mod_prug$species = "prug"
HR_Mod_rmet$species = "rmet"

HR_Mod_All <- rbind(rbind(rbind(rbind(HR_Mod_cfel, HR_Mod_drug),HR_Mod_ipur), HR_Mod_prug), HR_Mod_rmet)

# Statistics of home range versus foraging
HR_Mod_model <- lmer(foraging ~ HomeRange + (1 | species), data = HR_Mod_All)
HR_Mod_summary <- coef(summary(HR_Mod_model))

