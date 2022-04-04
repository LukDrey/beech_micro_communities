big_alg <- prune_samples(sample_data(phy_algae)$size_cat %in% c("big"), phy_algae)

medium_alg <- prune_samples(sample_data(phy_algae)$size_cat %in% c("medium"), phy_algae)

small_alg <- prune_samples(sample_data(phy_algae)$size_cat %in% c("small"), phy_algae)

big_genus <- aggregate_taxa(big_alg, "Genus")

medium_genus <- aggregate_taxa(medium_alg, "Genus")

small_genus <- aggregate_taxa(small_alg, "Genus")

big_taxa <- subset_taxa(big_genus, Genus %in% c("Apatococcus", "Dictyochloropsis", "Trebouxia", "Symbiochloris"))

medium_taxa <- subset_taxa(medium_genus, Genus %in% c("Apatococcus", "Dictyochloropsis", "Trebouxia", "Symbiochloris"))

small_taxa <- subset_taxa(small_genus, Genus %in% c("Apatococcus", "Dictyochloropsis", "Trebouxia", "Symbiochloris"))

taxa_sums(big_taxa)

taxa_sums(medium_taxa)

taxa_sums(small_taxa)
