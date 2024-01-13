# Set directories and import data
FACETNETDIR <- "/Users/tkay/Desktop/Work/facet_unil"
INPUTDIR <- "/Volumes/Lacie/CSNS/FacetNet_Input"
OUTPUTDIR <- "/Volumes/Lacie/CSNS/FacetNet_Output"
MAINDIR <- "/Volumes/Lacie/CSNS/data"
FIGDIR <- "/Volumes/Lacie/CSNS/figures"

setwd(MAINDIR)
for (file in list.files()[which(grepl("edgelist", list.files(), fixed = TRUE))]){
  assign(sub("\\..*", "", paste(file, "Contacts", sep = "_")),read.csv(file))
}

edgelist_list <- list(edgelist_cfel1, edgelist_cfel2, edgelist_cfel3, edgelist_cfel4, edgelist_cfel5,
                      edgelist_pbar1, edgelist_pbar2, edgelist_pbar3, edgelist_pbar4, edgelist_pbar5,
                      edgelist_drug1, edgelist_drug2, edgelist_drug3, edgelist_drug4, edgelist_drug5,
                      edgelist_ipur1, edgelist_ipur2, edgelist_ipur3, edgelist_ipur4, edgelist_ipur5,
                      edgelist_rmet1, edgelist_rmet2, edgelist_rmet3, edgelist_rmet4, edgelist_rmet5)

names_list <- c("cfel1", "cfel2", "cfel3", "cfel4", "cfel5",
                "pbar1", "pbar2", "pbar3", "pbar4", "pbar5",
                "drug1", "drug2", "drug3", "drug4", "drug5",
                "ipur1", "ipur2", "ipur3", "ipur4", "ipur5",
                "rmet1", "rmet2", "rmet3", "rmet4", "rmet5")

names_list_c <- paste0(names_list, "_c")
names_list_r <- paste0(names_list, "_r")

# Relabel dataframes for FacetNet
for (i in 1:25){
  setwd(paste(INPUTDIR, names_list[[i]], sep = "/"))
  write.table(edgelist_list[[i]], "0.edgelist", row.names = FALSE, col.names=FALSE)
}

# Set FacetNet commands & run FacetNet
alpha <- 1 
t_steps <- 1
for (m in c(5,4,3,2)){
  for (i in 1:25){
    command <- paste("python3", paste(FACETNETDIR, "facetnet_step.py", sep = "/"), paste(INPUTDIR,names_list[[i]], "0.edgelist", sep = "/"), alpha, m, paste(OUTPUTDIR, names_list_c[[i]], sep = "/"), t_steps, sep = " ")
    system(command)
  }
}

for (i in 1:25){
  for(ITERATION in 0:99){
    setwd(paste(INPUTDIR, names_list_r[[i]], sep = "/"))
    edgelist_list[[i]]$Count <- sample(edgelist_list[[i]]$Count)
    write.table(edgelist_list[[i]], paste(ITERATION, ".edgelist", sep = ""), row.names = FALSE, col.names=FALSE)
  }
}

alpha <- 1
t_steps <- 100
m <- 2
for (i in 1:25){
  command <- paste("python3", paste(FACETNETDIR, "facetnet_evol.py", sep = "/"), t_steps, alpha, m, paste(INPUTDIR,names_list_r[[i]], sep = "/"), paste(OUTPUTDIR,names_list_r[[i]], sep = "/"), sep = " ")
  system(command)
}