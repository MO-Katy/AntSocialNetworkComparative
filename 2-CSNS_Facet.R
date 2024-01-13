# Set directories and import data
FACETDIR <- "/Volumes/Lacie/CSNS/FacetNet_Output"
MAINDIR <- "/Volumes/Lacie/CSNS/data"
FIGDIR   <- "/Volumes/Lacie/CSNS/figures"
setwd(MAINDIR)
AllSoftMod <- read.csv("SoftModularity.csv", row.names = 1)
library(tidyr); library(lme4)

for (file in list.files()[which(grepl("randMod", list.files(), fixed = TRUE))]){
  assign(sub("\\..*", "", file), read.csv(file, check.names = FALSE))
}

randMod_list <- list(randMod_cfel1, randMod_cfel2, randMod_cfel3, randMod_cfel4, randMod_cfel5,
                     randMod_pbar1, randMod_pbar2, randMod_pbar3, randMod_pbar4, randMod_pbar5,
                     randMod_drug1, randMod_drug2, randMod_drug3, randMod_drug4, randMod_drug5,
                     randMod_ipur1, randMod_ipur2, randMod_ipur3, randMod_ipur4, randMod_ipur5,
                     randMod_rmet1, randMod_rmet2, randMod_rmet3, randMod_rmet4, randMod_rmet5)

# Remove iteration number
for (i in 1:length(randMod_list)){
  randMod_list[[i]] <- randMod_list[[i]]$`0`[seq(1, nrow(randMod_list[[i]]), 2)]
}

colnames(AllSoftMod) <- c("2", "3", "4", "5")
SpeciesSoftMod <- data.frame(cfel = colMeans(AllSoftMod[1:5,], na.rm = TRUE),
                             drug = colMeans(AllSoftMod[6:10,], na.rm = TRUE),
                             ipur = colMeans(AllSoftMod[11:15,], na.rm = TRUE),
                             pbar = colMeans(AllSoftMod[16:20,], na.rm = TRUE),
                             rmet = colMeans(AllSoftMod[21:25,], na.rm = TRUE))

# Plot soft modularity for no. communities = 2 -> 5. All species together.
setwd(FIGDIR)
jpeg('Soft_Modularity.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(as.numeric(SpeciesSoftMod[,1])~ c(2:5), type = "l", ylim = c(0,0.25), col = "#FF7F50", lwd = 15,
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylab = "Soft Modularity", xlab = "Community Number")
axis(2, at = c(0,0.25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(2,5), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
lines(as.numeric(SpeciesSoftMod[,2])~ c(2:5), type = "l", col = "#808080", lwd = 15)
lines(as.numeric(SpeciesSoftMod[,3])~ c(2:5), type = "l", col = "#6A0DAD", lwd = 15)
lines(as.numeric(SpeciesSoftMod[,4])~ c(2:5), type = "l", col = "#50C878", lwd = 15)
lines(as.numeric(SpeciesSoftMod[,5])~ c(2:5), type = "l", col = "cornflowerblue", lwd = 15)
points(as.numeric(AllSoftMod[,1])~ rep(2,25), cex = 4, pch = 16, col = rep(c("#FF7F50", "#808080", "#6A0DAD", "#50C878", "cornflowerblue"), each = 5))
points(as.numeric(AllSoftMod[,2])~ rep(3,25), cex = 4, pch = 16, col = rep(c("#FF7F50", "#808080", "#6A0DAD", "#50C878", "cornflowerblue"), each = 5))
points(as.numeric(AllSoftMod[,3])~ rep(4,25), cex = 4, pch = 16, col = rep(c("#FF7F50", "#808080", "#6A0DAD", "#50C878", "cornflowerblue"), each = 5))
points(as.numeric(AllSoftMod[,4])~ rep(5,25), cex = 4, pch = 16, col = rep(c("#FF7F50", "#808080", "#6A0DAD", "#50C878", "cornflowerblue"), each = 5))
legend(x=4, y = 0.25, col = c("#FF7F50", "#808080", "#6A0DAD", "#50C878", "cornflowerblue"), legend = c("C. fellah", "D. rugosum", "I. purpureus", "P. rugosus", "R. metallica"),
       fill = c("#FF7F50", "#808080", "#6A0DAD", "#50C878", "cornflowerblue"), cex = 6, bty = "n", text.font = 3)
dev.off()

# Plot soft modularity for no. communities = 2 -> 5. Separatgin species.
jpeg('Soft_Modularity_Colony.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,5), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(as.numeric(AllSoftMod[1,])~ c(2:5), type = "l", ylim = c(0,0.25), col = "#FF7F50", lwd = 15,
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylab = "Soft Modularity", xlab = "Community Number", main = substitute(paste(italic("C. fellah"))))
axis(2, at = c(0,0.25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 8)
axis(1, at = c(2,5), labels = FALSE, tick = TRUE, lwd = 8, cex.axis = 8)
lines(as.numeric(AllSoftMod[2,])~ c(2:5), type = "l", col = "#FF7F50", lwd = 15)
lines(as.numeric(AllSoftMod[3,])~ c(2:5), type = "l", col = "#FF7F50", lwd = 15)
lines(as.numeric(AllSoftMod[4,])~ c(2:5), type = "l", col = "#FF7F50", lwd = 15)
lines(as.numeric(AllSoftMod[5,])~ c(2:5), type = "l", col = "#FF7F50", lwd = 15)
plot(as.numeric(AllSoftMod[6,])~ c(2:5), type = "l", ylim = c(0,0.25), col = "#808080", lwd = 15,
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylab = "Soft Modularity", xlab = "Community Number", main = substitute(paste(italic("D. rugosum"))))
axis(2, at = c(0,0.25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 8)
axis(1, at = c(2,5), labels = FALSE, tick = TRUE, lwd = 8, cex.axis = 8)
lines(as.numeric(AllSoftMod[7,])~ c(2:5), type = "l", col = "#808080", lwd = 15)
lines(as.numeric(AllSoftMod[8,])~ c(2:5), type = "l", col = "#808080", lwd = 15)
lines(as.numeric(AllSoftMod[9,])~ c(2:5), type = "l", col = "#808080", lwd = 15)
lines(as.numeric(AllSoftMod[10,])~ c(2:5), type = "l", col = "#808080", lwd = 15)
plot(as.numeric(AllSoftMod[11,])~ c(2:5), type = "l", ylim = c(0,0.25), col = "#6A0DAD", lwd = 15,
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylab = "Soft Modularity", xlab = "Community Number", main = substitute(paste(italic("I. purpureus"))))
axis(2, at = c(0,0.25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 8)
axis(1, at = c(2,5), labels = FALSE, tick = TRUE, lwd = 8, cex.axis = 8)
lines(as.numeric(AllSoftMod[12,])~ c(2:5), type = "l", col = "#6A0DAD", lwd = 15)
lines(as.numeric(AllSoftMod[13,])~ c(2:5), type = "l", col = "#6A0DAD", lwd = 15)
lines(as.numeric(AllSoftMod[14,])~ c(2:5), type = "l", col = "#6A0DAD", lwd = 15)
lines(as.numeric(AllSoftMod[15,])~ c(2:5), type = "l", col = "#6A0DAD", lwd = 15)
plot(as.numeric(AllSoftMod[16,])~ c(2:5), type = "l", ylim = c(0,0.25), col = "#50C878", lwd = 15,
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylab = "Soft Modularity", xlab = "Community Number", main = substitute(paste(italic("P. rugosus"))))
axis(2, at = c(0,0.25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 8)
axis(1, at = c(2,5), labels = FALSE, tick = TRUE, lwd = 8, cex.axis = 8)
lines(as.numeric(AllSoftMod[17,])~ c(2:5), type = "l", col = "#50C878", lwd = 15)
lines(as.numeric(AllSoftMod[18,])~ c(2:5), type = "l", col = "#50C878", lwd = 15)
lines(as.numeric(AllSoftMod[19,])~ c(2:5), type = "l", col = "#50C878", lwd = 15)
lines(as.numeric(AllSoftMod[20,])~ c(2:5), type = "l", col = "#50C878", lwd = 15)
plot(as.numeric(AllSoftMod[21,])~ c(2:5), type = "l", ylim = c(0,0.25), col = "cornflowerblue", lwd = 15,
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 12, cex = 5, ylab = "Soft Modularity", xlab = "Community Number", main = substitute(paste(italic("R. metallica"))))
axis(2, at = c(0,0.25), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 8)
axis(1, at = c(2,5), labels = FALSE, tick = TRUE, lwd = 8, cex.axis = 8)
lines(as.numeric(AllSoftMod[22,])~ c(2:5), type = "l", col = "cornflowerblue", lwd = 15)
lines(as.numeric(AllSoftMod[23,])~ c(2:5), type = "l", col = "cornflowerblue", lwd = 15)
lines(as.numeric(AllSoftMod[24,])~ c(2:5), type = "l", col = "cornflowerblue", lwd = 15)
lines(as.numeric(AllSoftMod[25,])~ c(2:5), type = "l", col = "cornflowerblue", lwd = 15)
dev.off()

# Plot the expected distribution of modularity scores for randomly rewired networks compared to observed values
jpeg('Modularity_Significance.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(5,5), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(5,5,0), family = "serif")
for(i in 1:25){
  hist(randMod_list[[i]], breaks = 5, xlim = c(0,0.25),
       xlab = "Soft modularity", main = "", yaxt="n", xaxt="n", cex.lab = 8)
  abline(v=AllSoftMod[i,1], col = "#FF7F50", lwd = 8)
  axis(2, at = c(0,60), labels = FALSE, tick = TRUE, lwd = 2, cex.axis = 8)
  axis(1, at = c(0,0.25), labels = FALSE, tick = TRUE, lwd = 2, cex.axis = 8)
}
dev.off()

# Is there a significant effect of species on modularity?
AllSoftMod$species <- as.factor(rep(c("cfel", "drug", "ipur", "prug", "rmet"), each = 5))
summary(aov(`2` ~ species, data = AllSoftMod))

# Statistics comparing 2 communities to 3, 4 & 5:
t.test(AllSoftMod[1:5,1]- AllSoftMod[1:5,2])
t.test(AllSoftMod[1:5,1]- AllSoftMod[1:5,3])
t.test(AllSoftMod[1:5,1]- AllSoftMod[1:5,4])

t.test(AllSoftMod[6:10,1]- AllSoftMod[6:10,2])
t.test(AllSoftMod[6:10,1]- AllSoftMod[6:10,3])
t.test(AllSoftMod[6:10,1]- AllSoftMod[6:10,4])

t.test(AllSoftMod[11:15,1]- AllSoftMod[11:15,2])
t.test(AllSoftMod[11:15,1]- AllSoftMod[11:15,3])
t.test(AllSoftMod[11:15,1]- AllSoftMod[11:15,4])

t.test(AllSoftMod[16:20,1]- AllSoftMod[16:20,2])
t.test(AllSoftMod[16:20,1]- AllSoftMod[16:20,3])
t.test(AllSoftMod[16:20,1]- AllSoftMod[16:20,4])

t.test(AllSoftMod[21:25,1]- AllSoftMod[21:25,2])
t.test(AllSoftMod[21:25,1]- AllSoftMod[21:25,3])
t.test(AllSoftMod[21:25,1]- AllSoftMod[21:25,4])
