

# Custom Traits

# PHA:
- K01895
- K00925
- K00625
- K00626
- K00023
- K03821
```{r}

traits <- list("PhosphatePhoU" = c(RNAseq.data$features$annotation.db$module.dict$`Phosphate transporter`, "K02039"))


custom_trait_attributes <- Add_Custom_Trait(RNAseq.data, traits, pairwise.distances, distance.metrics, bkgd.traits)

#Incorporate the custom traits into the existing library

trait.attributes.pruned <- c(trait.attributes.pruned, custom_trait_attributes)


#Incorporate the custom trait definition into the annotation.db
RNAseq.data$features$annotation.db$module.dict <- c(RNAseq.data$features$annotation.db$module.dict, traits)

```



```{r}
Plot_Trait_Attribute_Expression("PhosphatePhoU", trait.attributes.pruned, RNAseq.data)
```


