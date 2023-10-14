rm(list = ls())

library(data.table)
library(qs) # read qs files
library(arrow) # read .parquet files
library(dendextend) # add info in dendrograms
library(pheatmap)
library(rcdk)
library(fingerprint) # fp.sim.matrix()
library(ggplot2)


setwd("~/github/SingleCellPerturbations")

source("R/utils.R")
# Load data ----

filename <- "data/de_train.parquet"
de_train <- read_parquet(filename)
de_train <- setDT(de_train)

head(de_train)
dataSummary(de_train)

# Molecular Descriptor ----

descriptors <- get.desc.names(type = "all")

desc_names <- gsub("org.openscience.cdk.qsar.descriptors.molecular.", "", descriptors)
desc_names <- gsub("Descriptor", "", desc_names)
desc_names <- gsub("org.openscience.cdk.qsar.descriptors.protein.", "", desc_names)

cat("Molecular descriptors available in the rcdk package:", length(desc_names), "\n")
cat(desc_names, sep = "\n")


drugs <- unique(de_train$sm_name)
mols <- parse.smiles(unique(de_train$SMILES))

for (j in 1:length(descriptors)) {
  tryCatch(
    expr = {
      # print(paste(j, ":", descriptors[j], "\n"))
      # print(eval.desc(mols[[1]], descriptors[j]))
      eval.desc(mols[[1]], descriptors[j])
    },
    error = function(e) {
      print(descriptors[j])
    },
    warning = function(w) {
      print(descriptors[j])
      print(w)
    }
  )
}

excl_desc <- c(
  "org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor", # 3D
  "org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor", # 3D
  "org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor", # 3D
  "org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor", # 3D
  "org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor", # 3D
  "org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor", # 3D
  "org.openscience.cdk.qsar.descriptors.protein.TaeAminoAcidDescriptor", # The molecule should be of type IBioPolymer
  "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor", # R would hang on mol[[2]]
  "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor" # geomShape is NA for most molecules
)
descriptors <- descriptors[!descriptors %in% excl_desc]
cat("Molecular descriptors remained:", length(descriptors), "\n")


desc_names <- gsub("org.openscience.cdk.qsar.descriptors.molecular.", "", descriptors)
desc_names <- gsub("Descriptor", "", desc_names)
desc_names <- gsub("org.openscience.cdk.qsar.descriptors.protein.", "", desc_names)

desc.list <- list()
col_names <- c()
start_time <- Sys.time()
for (i in 1:length(mols)) {
  tmp <- c()
  for (j in 1:length(descriptors)) {
    # calculate descriptor - each one returns a vector
    values <- eval.desc(mols[[i]], descriptors[j])
    if (any(is.na(values))) {
      cat(
        "For", drugs[[i]],
        "the descriptor", desc_names[j], "has NA"
      )
    }
    tmp <- c(tmp, values)
    # only on the first iteration, come up with column names
    if (i == 1) {
      col_names <- c(col_names, paste(desc_names[j], names(values), sep = "."))
    }
  }
  desc.list[[i]] <- tmp
}

mol_descriptors_features <- as.data.table(do.call("rbind", desc.list))
names(mol_descriptors_features) <- col_names
mol_descriptors_features[, (names(mol_descriptors_features)) :=
  lapply(.SD, unlist),
.SDcols = names(mol_descriptors_features)
]
mol_descriptors_features <- cbind(
  "sm_name" = drugs,
  mol_descriptors_features
)

Sys.time() - start_time
head(mol_descriptors_features)

# Remove constant features and cols with NA ----
# drugs with NAs (save as above)
cols <- names(mol_descriptors_features)
na_cols <- sapply(mol_descriptors_features[, ..cols], function(x) any(is.na(x)))
na_cols <- names(na_cols[na_cols])
# check_na <- cbind("sm_name" = drugs,
#                    mol_descriptors_features[, ..na_cols])
cat(
  "Removed", length(na_cols),
  "columns with NA for 2 compounds (CEP-18770 and MLN 2238, see above)"
)
mol_descriptors_features <- mol_descriptors_features[, !na_cols, with = FALSE]

# constant features
cols <- names(mol_descriptors_features)
const_cols <- sapply(mol_descriptors_features[, ..cols], function(x) uniqueN(x) == 1)
const_cols <- names(const_cols[const_cols])
cat(
  "\nRemoved", length(const_cols),
  "columns with the same value for all compounds"
)

mol_descriptors_features <- mol_descriptors_features[, !const_cols, with = FALSE]
cat("\n The resulting number of features:", ncol(mol_descriptors_features) - 1)

# binary columns
# cols <- names(mol_descriptors_features)[-1]
# binary_cols <- sapply(mol_descriptors_features[, ..cols], function(x) uniqueN(x) == 2)
# binary_cols[binary_cols]

# Number of parameters per descriptor ----

desc_summary <- data.table(
  Descriptor = gsub("\\..*", "", names(mol_descriptors_features))[-1],
  Parameters = gsub(".*\\.", "", names(mol_descriptors_features))[-1]
)
# desc_summary

# Groups are my educated guess!!
desc_map <- list(
  "Molecular complexity" = c(
    "FragmentComplexity", "KappaShapeIndices", "ZagrebIndex",
    "PetitjeanNumber", "AtomCount", "FractionalCSP3"
  ),
  "Physicochemical properties" = c(
    "AutocorrelationMass", "AutocorrelationPolarizability",
    "LargestPiSystem", "HybridizationRatio", "MannholdLogP",
    "ALOGP", "XLogP", "BPol", "APol", "FractionalPSA",
    "TPSA", "RuleOfFive", "Weight"
  ),
  "Structural features" = c(
    "AromaticBondsCount", "AromaticAtomsCount", "LargestChain",
    "SmallRing", "MDE", "EccentricConnectivityIndex",
    "ChiPath", "ChiPathCluster", "ChiCluster", "ChiChain",
    "BondCount", "RotatableBondsCount", "KierHallSmarts",
    "FMF", "BasicGroupCount", "AcidicGroupCount",
    "AminoAcidCount", "CarbonTypes", "HBondDonorCount",
    "HBondAcceptorCount", "WienerNumbers",
    "WeightedPath", "VAdjMa"
  )
)

desc_summary <- desc_summary[, .N, by = "Descriptor"]
desc_summary[, Group := sapply(Descriptor, function(d) {
  l_num <- which(sapply(desc_map, function(group) d %in% group))
  names(desc_map[l_num])
})]
setcolorder(desc_summary, neworder = c("Group"))

desc_summary[order(Group)]

# check <- "Weight."
# check <- c("sm_name", names(mol_descriptors_features)[
#     grepl(check, names(mol_descriptors_features)) ])
# head(mol_descriptors_features[, ..check])
# range(mol_descriptors_features[, ..check][,2])


getDig(15, 7)
group_sum <- copy(desc_summary)
group_sum <- group_sum[, .(count = sum(N)), by = "Group"][, prop := count / sum(count)]

ggplot(group_sum, aes(x = "", y = count, fill = factor(Group))) +
  geom_bar(stat = "identity", width = 1) +
  # coord_polar("y", start = 0) +
  labs(
    fill = "Descriptor group", x = NULL, y = NULL,
    title = "Number of features by the Molecular Descriptor groups"
  ) +
  # scale_fill_manual(values = c("#7BA23F", "#65a1c7", "#a674b5")) +
  theme_bw(base_size = 18) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    legend.position = "bottom"
  ) +
  geom_text(aes(label = paste0(count, " (", scales::percent(prop), ")")),
    position = position_stack(vjust = 0.5), size = 6
  )

# lot correlations between the molecular descriptors ----


filename <- "data/mol_descriptors_features.csv"
mol_descriptors_features <- fread(filename)
getwd()


getFig(25, 25)

cols <- names(mol_descriptors_features)[-1]
desc_corr <- cor(mol_descriptors_features[, ..cols], method = "spearman")
myColors <- c(
  colorRampPalette(c("darkblue", "white"))(100),
  colorRampPalette(c("white", "red"))(100)
)
breaksList <- seq(-1, 1, by = 0.01)

plt <- pheatmap(desc_corr,
  fontsize = 12,
  color = myColors,
  breaks = breaksList, # col = colorRampPalette(c("darkblue","white","red"))(200),
  main = "Molecular descriptors correlations"
)
plt


# PCA ----

pca_mol_data <- prcomp(mol_descriptors_features[, -1, with = FALSE],
  scale. = TRUE,
  center = TRUE
)

pca_targets_encoded <- de_train[, .(sm_name, cell_type)]
pca_targets_encoded <- cbind(
  pca_targets_encoded,
  pca_mol_data$x
)
head(pca_targets_encoded)
data.summary(pca_targets_encoded)

getFig(20, 10)

cumulative_var <- cumsum(pca_mol_data$sdev^2 / sum(pca_mol_data$sdev^2))

plot(cumulative_var,
  xlab = "Principal Component",
  ylab = "Cumulative Proportion of Variance Explained",
  ylim = c(0, 1), type = "b", cex.axis = 1.5, cex.lab = 1.5, xaxt = "n"
)
breaks <- seq(0, ncol(pca_mol_data$x), by = 10)
axis(1, at = breaks, labels = breaks)


abline(v = c(which.max(cumulative_var >= 0.80)), col = "red", lty = 2)
abline(h = 0.80, col = "red", lty = 2)

filename <- "data/molecular_data_pca.csv"
fwrite(pca_targets_encoded, filename)
