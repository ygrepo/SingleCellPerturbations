rm(list = ls())
library(data.table)
library(arrow)
library(dplyr)
library(tidyr)
library(stringdist)
library(qs)
library(rcdk)
library(ggplot2)
library(visNetwork)
library(dbparser)


setwd("~/github/SingleCellPerturbations")

source("R/utils.R")


# Load data ----

filename <- "data/de_train.parquet"
de_train <- read_parquet(filename)
de_train <- setDT(de_train)

op_drugs <- unique(de_train[, .(sm_name, SMILES)])
head(op_drugs)
unique(op_drugs$sm_name) %>% sort()


# Load Drugs ----

# dvobj <- parseDrugBank(db_path            = "drugs/full database.xml",
#                        drug_options       = drug_node_options(),
#                        parse_salts        = TRUE,
#                        parse_products     = TRUE,
#                        references_options = references_node_options(),
#                        cett_options       = cett_nodes_options())

# salts <- dvobj$salts
# qsave(salts, "drugs/drug_bank/salts.qs")

# products <- dvobj$products
# qsave(products, "drugs/drug_bank/products.qs")

# references <- dvobj$references
# qsave(references, "drugs/drug_bank/references.qs")

# cett <- dvobj$cett
# qsave(cett, "drugs/drug_bank/cett.qs")

# drugs <- dvobj$drugs
# qsave(drugs, "drugs/drug_bank/drugs.qs")


drugs <- qread("data/drugs.qs")
vocab <- fread("data/drugbank vocabulary.csv")

cat("DrugBank Vocabulary (identifiers, names, and synonyms to permit
    easy linking and integration into any type of project)")
head(vocab)
cat("\nShape: ", ncol(vocab), " columns x ",
  nrow(vocab), " rows",
  sep = ""
)

# Check smiles ----
b_smiles <- left_join(
  select(drugs$general_information, c("primary_key", "name")) %>%
    rename(sm_name = name),
  drugs$calculated_properties %>%
    filter(kind == "SMILES") %>%
    select(value, parent_key) %>%
    rename(Drug_bank_SMILES = value),
  by = c("primary_key" = "parent_key")
) %>%
  as.data.table()

head(db_smiles[op_drugs, on = "sm_name"])

# Find the most similar Drug Bank names to that in the de_train ----

# Calculate min distances
op_drugs[, min_dist := sapply(sm_name, function(x) {
  distances <- stringdist(tolower(x), tolower(vocab$`Common name`))
  min(distances)
})]

# Get the closest name
op_drugs[, DrugBank_name := sapply(sm_name, function(x) {
  # the threshold is set based on manual check)
  if (op_drugs[sm_name == x]$min_dist < 2) {
    distances <- stringdist(tolower(x), tolower(vocab$`Common name`))
    vocab$`Common name`[which.min(distances)]
  } else {
    NA
  }
})]

head(op_drugs[!is.na(DrugBank_name)][order(-min_dist)], 10)
cat("Not found drugs:", op_drugs[is.na(DrugBank_name), .N])

# Match compound names using synonyms ----

check <- op_drugs[min_dist > 1]$sm_name
vocab[grepl(paste0(check, collapse = "|"), vocab$Synonyms), ]
# only one found

op_drugs[sm_name == "PD-0325901", DrugBank_name := "Mirdametinib"]

cat("Not found drugs:", op_drugs[is.na(DrugBank_name), .N])


# Match compounds by chemical formula ----
# Get the formula of competition drugs
op_drugs[, formula := sapply(SMILES, function(x) {
  mol <- parse.smiles(x)
  f <- get.mol2formula(mol[[1]])
  f@string
})]

# Get the formula of Drug Bank drugs
db_form <- drugs$calculated_properties %>%
  filter(kind == "Molecular Formula") %>%
  select(value, parent_key) %>%
  rename(formula = value) %>%
  as.data.table()

# Subset not yet found drugs and merge with formula data
vocab <- left_join(vocab, db_form, by = c("DrugBank ID" = "parent_key"))
vocab <- vocab[, .(`DrugBank ID`, DrugBank_name = `Common name`, CAS, formula)]

vocab <- vocab[
  op_drugs[
    is.na(DrugBank_name),
    .(
      sm_name,
      formula, SMILES
    )
  ],
  on = "formula"
]
cat("Drug Bank vocabulary merged with competition drugs by formula")
head(vocab)

# Manual Check ----
# Pick one if several compounds have the same formula (searched manually in pubchem)
dupl_formula <- vocab[, .N, by = formula][N > 1]
cat("Several compounds have the same formula:")
vocab[formula %in% dupl_formula$formula]

op_drugs[
  sm_name == "5-(9-Isopropyl-8-methyl-2-morpholino-9H-purin-6-yl)pyrimidin-2-amine",
  DrugBank_name := "VS-5584"
]
op_drugs[
  sm_name == "RG7090",
  DrugBank_name := "Basimglurant"
]
op_drugs[
  sm_name == "RVX-208",
  DrugBank_name := "Apabetalone"
]
op_drugs[
  sm_name == "Ganetespib (STA-9090)",
  DrugBank_name := "Ganetespib"
]

# manuall check
# one wrong match https://pubchem.ncbi.nlm.nih.gov/compound/46173038#section=Synonyms
vocab[sm_name == "K-02288"]$DrugBank_name <- NA
vocab[sm_name == "K-02288"]$`DrugBank ID` <- NA
vocab[sm_name == "K-02288"]$CAS <- NA

# note Pitavastatin Calcium - is a wrong name in the de_train, should be Pitavastatin:
# https://pubchem.ncbi.nlm.nih.gov/compound/5282452
# https://pubchem.ncbi.nlm.nih.gov/compound/5282451

not_found <- vocab[is.na(DrugBank_name)]
vocab <- vocab[!is.na(DrugBank_name), .(DrugBank_name, sm_name)]

add_names <- op_drugs[is.na(DrugBank_name), !"DrugBank_name"]
op_drugs_map <- rbind(
  op_drugs[!is.na(DrugBank_name)],
  add_names[vocab, on = c("sm_name"), nomatch = 0]
)
op_drugs_map <- op_drugs_map[, .(sm_name, DrugBank_name)]

cat(
  "The resulting mapping of the competition's compound names and Drug Bank names",
  "\n(As an example I show rows with the same compounds but different names)"
)
head(op_drugs_map[sm_name != DrugBank_name])
cat("\nShape: ", ncol(op_drugs_map), " columns x ",
  nrow(op_drugs_map), " rows",
  sep = ""
)
cat("\n\nNot found (", nrow(not_found), "compounds):")
not_found$sm_name

# EDA for the drugs found in Drug Bank ----

names(drugs$general_information)

# Select relevant columns for further analysis
drugs_data <- op_drugs_map %>%
  left_join(rename(drugs$general_information, drugbank_id = primary_key),
    by = c("DrugBank_name" = "name")
  ) %>%
  select(
    sm_name, DrugBank_name, drugbank_id, cas_number,
    type, description,
    average_mass, monoisotopic_mass, state
  )

cat("General information data for the competition's drugs")
head(drugs_data, 3)

cat("The counts of drugs by their type (can be one of Small molecule or Biotech )\n")
drugs_data[, .N, by = "type"]

cat("The counts of drugs by their state\n")
drugs_data[, .N, by = "state"]

# Drug groups ----

getFig(20, 5)
groups <- left_join(drugs_data[, .(sm_name, drugbank_id)],
  drugs$groups,
  by = c("drugbank_id" = "drugbank-id")
)
groups %>%
  ggplot(aes(x = group, fill = "coral")) +
  geom_bar() +
  labs(
    x = "Drug Group",
    y = "Quantity",
    title = "Drug Group Distribution"
  ) +
  coord_flip() +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")


# Drug classification ----

drugs_class <- left_join(drugs_data[, .(sm_name, drugbank_id)],
  rename(drugs$drug_classification,
    description_class = description
  ),
  by = "drugbank_id"
)
cat("Drug classification data")
head(drugs_class, 3)

cat("The distribution of drugs by kingdom\n")
drugs_class[, .N, by = "kingdom"][order(-N)]
cat("\nThe distribution of drugs by superclass\n")
drugs_class[, .N, by = "superclass"][order(-N)]
cat("\nThe distribution of drugs by class\n")
drugs_class[, .N, by = "class"][order(-N)]

# Graph representation ----

graph <- drugs_class[
  !is.na(kingdom),
  .(
    kingdom, superclass, class,
    subclass, sm_name
  )
]

unique_nodes <- unique(na.omit(c(
  graph$kingdom, graph$superclass, graph$class, graph$subclass
)))
unique_nodes <- unique_nodes[unique_nodes != ""]

# Create a data frame of nodes
nodes <- data.frame(id = unique_nodes, label = unique_nodes)
nodes$group <- ifelse(nodes$id %in% graph$kingdom, "kingdom",
  ifelse(nodes$id %in% graph$superclass, "superclass",
    ifelse(nodes$id %in% graph$class, "class", "subclass")
  )
)
nodes$value <- as.numeric(factor(
  nodes$group,
  levels = c("subclass", "class", "superclass", "kingdom")
))

# Create a data frame of edges
edges <- data.frame(
  from = c(graph$kingdom, graph$superclass, graph$class),
  to = c(graph$superclass, graph$class, graph$subclass)
) %>%
  mutate(to = ifelse(to == "", NA, to)) %>%
  unique()

# Plot graph
visNetwork(nodes, edges, main = "The graph of compounds classification (Zoom in)", width = "100%") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend()

# Pharmacology ----

drugs_pharm <- left_join(drugs_data[, .(sm_name, drugbank_id)],
  drugs$pharmacology,
  by = "drugbank_id"
)

cat("Pharmacology of Midostaurin:\n\n")
unlist(drugs_pharm[sm_name == "Midostaurin",
  c(names(drugs$pharmacology)),
  with = FALSE
])

# ATC Classification ----

atc <- left_join(drugs_data[, .(sm_name, drugbank_id)],
  drugs$atc_codes,
  by = c("drugbank_id" = "drugbank-id")
) %>%
  select(sm_name, `drugbank_id`, level_1, level_2, level_3, level_4) %>%
  unique()
head(atc)

cat("Drug counts at top 10 1st level pharmacological groups")
atc[, .N, by = "level_1"][order(-N)][1:10]
cat("Drug counts at top 10 2nd level pharmacological groups")
atc[, .N, by = "level_2"][order(-N)][1:10]
cat("Drug counts at top 10 3rd level pharmacological groups")
atc[, .N, by = "level_3"][order(-N)][1:10]
cat("Drug counts at top 10 4th level pharmacological groups")
atc[, .N, by = "level_4"][order(-N)][1:10]

# Graph representation ----

graph <- atc[!is.na(level_1)]
unique_nodes <- unique(na.omit(c(
  graph$sm_name, graph$level_4, graph$level_3, graph$level_2, graph$level_1
)))

# Create a data frame of nodes
nodes <- data.frame(id = unique_nodes, label = unique_nodes)
nodes$group <- ifelse(nodes$id %in% graph$level_1, "level_1",
  ifelse(nodes$id %in% graph$level_2, "level_2",
    ifelse(nodes$id %in% graph$level_3, "level_3",
      ifelse(nodes$id %in% graph$level_4, "level_4", "sm_name")
    )
  )
)
nodes$value <- as.numeric(factor(
  nodes$group,
  levels = c("sm_name", "level_1", "level_2", "level_3", "level_4")
))

# Create a data frame of edges
edges <- data.frame(
  from = c(graph$sm_name, graph$level_1, graph$level_2, graph$level_3),
  to = c(graph$level_1, graph$level_2, graph$level_3, graph$level_4)
) %>% unique()

# Plot graph
visNetwork(nodes, edges,
  main = "Graph representation of ATC Classification of
  the competition compounds (Zoom in)", width = "100%"
) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend()

#  Calculated properties ----

# https://docs.drugbank.com/xml/#calculated-properties


cal_prop <- left_join(drugs_data[, .(sm_name, drugbank_id)],
  drugs$calculated_properties,
  by = c("drugbank_id" = "parent_key")
)

# There are two sources of data:
# unique(cal_prop$source)
# cal_prop[source == "ALOGPS", unique(kind)]
# cal_prop[source == "ChemAxon", unique(kind)]

# There are 24-26 observations per drug
# cal_prop[, .N, by = "sm_name"][, range(N)]

# As logP sometimes provided from two sources, we calculate mean logP
logp <- cal_prop[kind == "logP", !"source"]
logp[, value := mean(as.numeric(value)),
  by = "sm_name"
]
cal_prop <- rbind(
  unique(logp),
  cal_prop[kind != "logP", !"source"]
)
# Correct water solubility data as numeric
cal_prop[kind == "Water Solubility", value := gsub(" g/l.*", "", value)]

# cal_prop[sm_name == "HYDROXYUREA"]

getFig(20, 20)
non_num_cols <- c(
  "IUPAC Name", "Traditional IUPAC Name", "SMILES",
  "Molecular Formula", "InChI", "InChIKey", "Bioavailability",
  "Rule of Five", "Ghose Filter", "MDDR-Like Rule"
)
cal_prop %>%
  filter(!kind %in% non_num_cols) %>%
  mutate(value = as.numeric(value)) %>%
  ggplot(aes(value)) +
  geom_histogram(fill = "coral") +
  facet_wrap(~kind, scales = "free") +
  theme_bw(base_size = 22)

# Plot binary properties

getFig(10, 10)
cal_prop %>%
  filter(kind %in% c(
    "Bioavailability", "Rule of Five",
    "Ghose Filter", "MDDR-Like Rule"
  )) %>%
  ggplot(aes(value)) +
  geom_bar(fill = "coral") +
  facet_wrap(~kind, scales = "free") +
  theme_bw(base_size = 22)

# All properties in wide format ----


cal_prop <- pivot_wider(cal_prop,
  names_from = "kind",
  values_from = "value",
  values_fill = NA
) %>%
  as.data.table()

head(cal_pro)

# Pathways ----

pathways <- left_join(drugs_data[, .(sm_name, drugbank_id)],
  drugs$pathway$general_information,
  by = c("drugbank_id" = "parent_key")
)

pathways[!is.na(smpdb_id)]
