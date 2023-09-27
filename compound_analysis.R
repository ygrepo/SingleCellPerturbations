install.packages("rcdk", verbose = FALSE)

install.packages("arrow")
install.packages("data.table")
install.packages('dbscan')
library(rcdk)
library(arrow)
library(data.table)
library(fingerprint)
library(dendextend)
library(dbscan)

setwd('~/github/SingleCellPerturbations')
de_train <- read_parquet('data/de_train.parquet')

feats = de_train[,6:length(colnames(de_train))]
chems = de_train$sm_name

mols <- parse.smiles(unique(de_train$SMILES))
fps <- lapply(mols, get.fingerprint, type='circular')
names(fps) <- unique(de_train$sm_name)
head(fps)


de_train <- setDT(de_train)
head(de_train)
finger_features <- cbind(data.table("sm_name" = names(fps)),
                         fp.to.matrix(fps))
head(finger_feature)
ff_tab <- setDT(finger_features)
fwrite(ff_tab, "data/finger_features.csv")

options(repr.plot.width = 25, repr.plot.height = 10)
fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
fp.dist <- 1 - fp.sim
cls <- hclust(as.dist(fp.dist))
cls$labels <- names(fps)
plot(cls, main = 'Clustering of the compounds based on their molecular 
     fingerprints and the Tanimoto similarity metric')



heights <- cls$height
n_clusters <- length(heights) + 1  # Number of clusters including individual data points
wss <- numeric(n_clusters)

# Calculate within-cluster sum of squares (WSS) for different levels
for (k in 1:n_clusters) {
  cut_tree <- cutree(cls, k = k)
  cluster_heights <- heights[cut_tree]  # Heights of clusters at level k
  wss[k] <- sum((cluster_heights - mean(cluster_heights))^2)  # WSS calculation
}

# Find the "Elbow Point" using a heuristic
elbow_point <- 1
for (k in 2:(n_clusters - 1)) {
  if (wss[k] < wss[k - 1] && wss[k] < wss[k + 1]) {
    elbow_point <- k
    break
  }
}

# Visualize the "Elbow Point" on the dendrogram
dend <- as.dendrogram(cls)
dend %>% 
  set("labels_col", rainbow(elbow_point), k = elbow_point) %>% # change color
  set("branches_k_color", k = elbow_point) %>% 
  plot(main = 'Clustering of the compounds based on Elbow point (k = 11)')
abline(h = heights[elbow_point], col = "red", lty = 2)

clusters <- cutree(cls, k = elbow_point)
finger_features[, cluster:= clusters]

for (cl in unique(finger_features$cluster)) {
  
  cat("\ncluster №", cl, "\n")
  print(finger_features[cluster == cl, sm_name])
}
