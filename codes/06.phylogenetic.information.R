########################################################################
####### PHYLOGENETIC INFORMATION  ######################################
########################################################################

# Filter species data -----
load("rdata/temperate.all.rates.RData")
temperate_rates

gimnosperm <- c("Tsuga canadensis", "Abies balsamea","Larix laricina","Picea abies","Picea glauca",
                "Picea mariana", "Picea rubens", "Pinus banksiana", "Pinus echinata",
                "Pinus ponderosa","Pinus resinosa","Pinus rigida","Pinus strobus","Juniperus virginiana",
                "Pinus taeda", "Pinus virginiana", "Taxodium distichum","Thuja occidentalis")
temperate_rates$group[temperate_rates$sp %in% gimnosperm]= "Gymnosperm"
temperate_rates$group[is.na(temperate_rates$group)] = "Angiosperm"
temperate_rates <- temperate_rates %>%
  mutate (group = as.factor(group)) 


# species phylogenetic information ------
# Brownian Motion
df.analy <- temperate_rates
df.analy$species <- sub(" ", "_", df.analy$sp)
df.analy<- df.analy %>% remove_rownames %>% column_to_rownames(var="species")

load("rdata/phylotree.RData")
df.analy
BM <- corBrownian(1, phy = phylotree.r)
df.analy <- df.analy[phylotree.r$tip.label,]
rownames(df.analy)==phylotree.r$tip.label #match species and tip labels
str(df.analy)

# save(df.analy, file = "rdata/matrix.trait.preds.RData")
