# install.packages("rcdk", verbose = FALSE)
#
# install.packages("arrow")
# install.packages("plotly")
rm(list = ls())
library(rcdk)
library(arrow)
library(data.table)
library(fingerprint)
library(dendextend)
library(dbscan)
library(umap)
library(RColorBrewer)
library(plotly)


setwd("~/github/SingleCellPerturbations")
de_train <- read_parquet("data/de_train.parquet")

# Retrieve compound names ----
feats <- de_train[, 6:length(colnames(de_train))]
chems <- de_train$sm_name

# Get fingerprints from drug smiles ----
mols <- parse.smiles(unique(de_train$SMILES))
fps <- lapply(mols, get.fingerprint, type = "circular")
names(fps) <- unique(de_train$sm_name)
head(fps)


de_train <- setDT(de_train)
head(de_train)

# Create a table with rows the compounds and columns their fingerprints ----
finger_features <- cbind(
  data.table("sm_name" = names(fps)),
  fp.to.matrix(fps)
)
head(finger_features)
ff_tab <- setDT(finger_features)
fwrite(ff_tab, "data/finger_features.csv")

options(repr.plot.width = 25, repr.plot.height = 10)
# Compute Tanimoto similarity between fps ----
fp.sim <- fingerprint::fp.sim.matrix(fps, method = "tanimoto")
fp.dist <- 1 - fp.sim

# Clustering of the compounds based on their fps and Tanimoto similarity metrics ----
cls <- hclust(as.dist(fp.dist))
cls$labels <- names(fps)
plot(cls, main = "Clustering of the compounds based on their molecular
     fingerprints and the Tanimoto similarity metric")



heights <- cls$height
n_clusters <- length(heights) + 1 # Number of clusters including individual data points
wss <- numeric(n_clusters)

# Calculate within-cluster sum of squares (WSS) for different levels
for (k in 1:n_clusters) {
  cut_tree <- cutree(cls, k = k)
  cluster_heights <- heights[cut_tree] # Heights of clusters at level k
  wss[k] <- sum((cluster_heights - mean(cluster_heights))^2) # WSS calculation
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
  plot(main = "Clustering of the compounds based on Elbow point (k = 11)")
abline(h = heights[elbow_point], col = "red", lty = 2)

clusters <- cutree(cls, k = elbow_point)
finger_features[, cluster := clusters]

for (cl in unique(finger_features$cluster)) {
  cat("\ncluster №", cl, "\n")
  print(finger_features[cluster == cl, sm_name])
}

# Umap of the drugs using Tanimoto similarity metric ----
fp.umap <- umap(fp.sim, n_components = 3)
layout <- fp.umap[["layout"]]
layout <- data.frame(layout)
clusters <- data.frame(clusters)
layout <- cbind(layout, clusters$clusters)
layout <- cbind(layout, rownames(clusters))
colnames(layout)[colnames(layout) == "clusters$clusters"] <- "cluster"
colnames(layout)[colnames(layout) == "rownames(clusters)"] <- "molecule"
layout$cluster_name <- as.factor(layout$cluster)
head(layout)


display.brewer.all(type = "qual")
pal1 <- brewer.pal(12, "Set3")
pal2 <- brewer.pal(12, "Paired")
bigpal <- c(pal1, pal2)

plot_ly(layout,
  x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers",
  color = ~cluster_name, text = ~cluster_name, colors = bigpal
)

# Umap using fingerprints ----
ff.umap <- umap(finger_features[, 2:length(colnames(finger_features))], n_components = 3)
layout <- ff.umap[["layout"]]
layout <- data.frame(layout)
clusters <- data.frame(clusters)
layout <- cbind(layout, clusters$clusters)
layout <- cbind(layout, rownames(clusters))
colnames(layout)[colnames(layout) == "clusters$clusters"] <- "cluster"
colnames(layout)[colnames(layout) == "rownames(clusters)"] <- "molecule"
layout$cluster_name <- as.factor(layout$cluster)
head(layout)
fig3 <- plot_ly(layout,
  x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers",
  colors = bigpal, color = ~cluster_name
)
fig3

# Clustering using fingerprints ----
cl <- hdbscan(layout[, 1:3], minPts = 3)
layout$hbscan_cluster <- as.factor(cl$cluster)
plot_ly(layout,
  x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "text",
  color = ~hbscan_cluster, text = ~molecule, colors = bigpal
)
just_cls <- layout[, c("molecule", "hbscan_cluster", "cluster_name")]
head(just_cls)

sub <- just_cls[just_cls$hbscan_cluster == 15, ]
sub

# Umap using gene expressions ----
feat.umap <- umap(feats, n_components = 3)
fproj <- feat.umap[["layout"]]
fproj <- data.frame(fproj)
head(fproj)
fproj$chems <- chems
plot_ly(fproj,
  x = ~X1, y = ~X2, z = ~X3, type = "scatter3d",
  mode = "markers", text = ~chems
)
fproj$molecule <- fproj$chems
merged <- merge(fproj, just_cls, by = "molecule")
head(merged)
plot_ly(merged,
  x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", colors = bigpal,
  mode = "markers", text = ~molecule, color = ~hbscan_cluster
)

new_sub <- merged[merged$molecule %in% c("Scriptaid", "Vorinostat"), ]
new_sub
plot_ly(new_sub,
  x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", colors = bigpal,
  mode = "markers", text = ~molecule, color = "cluster_name"
)

to_merge <- merged[, c("cluster_name", "molecule")]
names(to_merge) <- c("cluster_name", "molecule")
head(to_merge)

# Plot of average of Umap coordinates ----
avg_proj <- aggregate(merged[, c("X1", "X2", "X3")], by = list(merged$molecule), FUN = mean)
names(avg_proj) <- c("molecule", "X1", "X2", "X3")
avg_proj <- merge(avg_proj, unique(to_merge), by = "molecule")
plot_ly(avg_proj,
  x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers",
  colors = bigpal, text = ~molecule, color = ~cluster_name
)
