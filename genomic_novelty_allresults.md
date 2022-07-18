genomic\_novelty\_analysis
================

``` r
library(tidyverse)
library(ggplot2)
theme_set(
  theme_light() + theme(legend.position = "right"))
```

**Figure1: Density distribution of Alignment length, bitscores and
percent identity in MAGs and isolates**

``` r
# loading the file the mags and isolate samples
mags_iso_samp <- read.csv("~/R/samp_all_density_new.csv")
## Alignment_length

align_length_dist <- ggplot(mags_iso_samp) + 
  geom_density(aes(x=log10(Alignment_Length), colour=cat),show_guide=FALSE)+
  stat_density(aes(x=log10(Alignment_Length), colour=cat),
               geom="line",position="identity")+
  scale_color_brewer(palette = "Dark2")+
  ylab("Density")+
  xlab("Alignment length")+
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"))

print(align_length_dist)
```

![](genomic_novelty_allresults_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
## Bitscore density distribution

bitscore_dist <- ggplot(mags_iso_samp) + 
  geom_density(aes(x=norm_bitscore, colour=cat),show_guide=FALSE)+
  stat_density(aes(x=norm_bitscore, colour=cat),
               geom="line",position="identity")+
  scale_color_brewer(palette = "Dark2")+
  ylab("Density")+
  xlab("Bit scores distribution")+
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"))

print(bitscore_dist)
```

![](genomic_novelty_allresults_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
## Percentage identity distribution
percent_ident_dist <- ggplot(mags_iso_samp) + 
  geom_density(aes(x=Pident_matches, colour=cat),show_guide=FALSE)+
  stat_density(aes(x=Pident_matches, colour=cat),
               geom="line",position="identity")+
  scale_color_brewer(palette = "Dark2")+
  ylab("Density")+
  xlab("Percent Identity")+
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"))


print(percent_ident_dist)
```

![](genomic_novelty_allresults_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

**figure 2: Genomic novelty Culture level distribution in mags and
isolates** Since we are not sure of ideal bitscore category we decided
to use several cutoffs for our bitscores categories

``` r
## loading the big diamond files with unique cultured level(cl) classification
diamond_unique_taxa <- read.csv("~/R/diamond_unique_taxa_corrected.csv")
# we want to filter out the domain level
diamond_unique_taxa_d <- filter(diamond_unique_taxa, !cultured_level == "domain")
#loading the files that have the isolates archaea and bacteria
isolate_archae_new <- read.csv("~/R/isolate_archae_new.csv")
isolate_bacteria_new <- read.csv("~/R/isolate_bacteria_new.csv")
diamond_arc_mag <- filter(diamond_unique_taxa_d, domain == "Archaea")
diamond_bac_mag <- filter(diamond_unique_taxa_d, domain == "Bacteria")

#Archaea mag category, calculating the norm_bitscore, grouping by the cultured level and sampling
arc_cl_samp <- diamond_arc_mag %>%
  #filter(Bit_score >= 50) %>%
  select(Pident_matches, Alignment_Length, cultured_level, Bit_score)%>%
  mutate(norm_bitscore = Bit_score/Alignment_Length)%>%
  group_by(cultured_level) %>%
  sample_n(1000)

#  Archaea iso category

iso_arc_samp <- isolate_archae_new %>%
  #filter(Bit_score >= 50) %>%
  select(Pident_matches, Alignment_Length, Bit_score) %>%
  mutate(norm_bitscore = Bit_score/Alignment_Length)%>%
  sample_n(1000)

iso_arc_samp$cultured_level <-"isolates"
#binding the archaea and mags files 
arc_cl_samp <- rbind(iso_arc_samp, arc_cl_samp)
arc_cl_samp$cultured_level<-factor(arc_cl_samp$cultured_level, 
                                   levels=c("isolates","species","genus","family","order","class","phylum"))
#plot of normalized bit_score
arc_cl_bitsplot <- ggplot(arc_cl_samp, aes(x= cultured_level, y= norm_bitscore)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.3, height = 0), alpha = 0.03) +
  ylab("Normalized Bit scores") +
  xlab("Taxonomic cultured level") +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=12),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"))

print(arc_cl_bitsplot)
```

![](genomic_novelty_allresults_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#plot of percentage identity 
arc_cl_identplot <- ggplot(arc_cl_samp, aes(x= cultured_level, y= Pident_matches)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.3, height = 0), alpha = 0.03) +
  ylab("Identity to the SwissProt") +
  xlab("Taxonomic cultured level") +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=12),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"))

print(arc_cl_identplot)
```

![](genomic_novelty_allresults_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
#Comments
 
#we filter our bitscores from all bitscores, >50, >200, >500, >1000, >1500 
#we replicate the same analysis for bacteria
```

\*\*Figure3: correlation between the taxonomic “unculturedness” of mag
and isolate phyla (i.e., how far it is from a cultured category(refseq))
and the similarity of its genes to genes in Swissprot

``` r
# We filtered the dominant selected mag phyla 
selected_mag_phyla <- diamond_unique_taxa_d %>%
  filter(phylum %in% c("Proteobacteria", "Firmicutes", "Actinobacteriota",
                      "Desulfobacterota", "Patescibacteria", "Planctomycetota",
                       "Bacteroidota",  "Acidobacteriota", "Spirochaetota",
                      "Halobacterota", "Thermoplasmatota", "Crenaracheota", "Synergistota"))
# calculating the mean identity of each cultured level 
mag_phyla_all <- selected_mag_phyla %>%
  #filter(Bit_score >= 50) %>%
  select(Pident_matches, Alignment_Length, Bit_score, cultured_level, phylum, taxonomic_dist) %>%
  mutate(norm_bitscore = Bit_score/Alignment_Length)%>%
  group_by(taxonomic_dist, phylum)%>%
  summarize(mean_identity = mean(Pident_matches), na.rm = TRUE)

## Repeating the same analysis for bit score >= 50, bit score >= 200

#plotting taxonomic distance and mean identity
ggplot(na.omit(mag_phyla_all), aes(x = taxonomic_dist, y = mean_identity, color = phylum)) +
  geom_line() + 
  geom_point()+
  xlab("Taxonomic distance to RefSeq relative")+
  ylab("Mean identity to SwissProt")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
```

![](genomic_novelty_allresults_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
  #facet_wrap(~phylum) + 
  theme(
    strip.text.x = element_text(
      size = 9, color = "black", face = "bold.italic"
    ))
```

    ## List of 1
    ##  $ strip.text.x:List of 11
    ##   ..$ family       : NULL
    ##   ..$ face         : chr "bold.italic"
    ##   ..$ colour       : chr "black"
    ##   ..$ size         : num 9
    ##   ..$ hjust        : NULL
    ##   ..$ vjust        : NULL
    ##   ..$ angle        : NULL
    ##   ..$ lineheight   : NULL
    ##   ..$ margin       : NULL
    ##   ..$ debug        : NULL
    ##   ..$ inherit.blank: logi FALSE
    ##   ..- attr(*, "class")= chr [1:2] "element_text" "element"
    ##  - attr(*, "class")= chr [1:2] "theme" "gg"
    ##  - attr(*, "complete")= logi FALSE
    ##  - attr(*, "validate")= logi TRUE

**Environment group classification**

``` r
#loading files  containg the ecosystem category
env_mag <- read.csv("~/R/tax_env.csv")
#merge the files with cultured level and taxonomic distance 
select_mag_env <-diamond_unique_taxa_d %>%
  dplyr::select(genome_id, Bit_score, Pident_matches, Alignment_Length, cultured_level, taxonomic_dist) %>%
  merge(env_mag, by="genome_id")
# Terrestrial subset
select_mag_env_ter <- select_mag_env %>%
  filter(ecosystem_category %in% "Terrestrial")

select_mag_env_ter$habitat <- "Terrestrial"
# Aquatic subset
select_mag_env_aqu <- select_mag_env %>%
  filter(ecosystem_category %in% "Aquatic")

select_mag_env_aqu$habitat <- "Aquatic"
# Host associated subsets
select_mag_env_hostass <- select_mag_env %>%
  filter(ecosystem_category %in% c("Human", "Mammals", "Plants", "Anthropoda", "Fungi"))

select_mag_env_hostass$habitat <- "Host_associated"
# Engineered subsets
select_mag_env_eng <- select_mag_env %>%
  filter(ecosystem_category %in% c("Solid waste", "Wastewater", "Lab enrichment", "Built environment", "Biotransformation"))


select_mag_env_eng$habitat <- "Engineered"

#Bind all the subsets of the files together 
mag_env_groups <- rbind(select_mag_env_aqu, select_mag_env_eng, select_mag_env_hostass, select_mag_env_ter)

mag_env_all <- mag_env_groups %>%
  #filter(Bit_score >= 200) %>%
  select(Pident_matches, Alignment_Length, Bit_score, cultured_level, phylum, taxonomic_dist, habitat) %>%
  mutate(norm_bitscore = Bit_score/Alignment_Length)%>%
  group_by(taxonomic_dist, habitat)%>%
  summarize(mean_identity = mean(Pident_matches), na.rm = TRUE)

#Repeat the same analysis for bitscore >= 50, bitscore >= 200 
ggplot(na.omit(mag_env_all), aes(x = taxonomic_dist, y = mean_identity, color = habitat)) +
  geom_line() + 
  geom_point()+
  xlab("Taxonomic distance to RefSeq relative")+
  ylab("Mean identity to SwissProt")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
```

![](genomic_novelty_allresults_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
  #facet_wrap(~habitat) + 
  theme(
    strip.text.x = element_text(
      size = 9, color = "black", face = "bold.italic"
    ))
```

    ## List of 1
    ##  $ strip.text.x:List of 11
    ##   ..$ family       : NULL
    ##   ..$ face         : chr "bold.italic"
    ##   ..$ colour       : chr "black"
    ##   ..$ size         : num 9
    ##   ..$ hjust        : NULL
    ##   ..$ vjust        : NULL
    ##   ..$ angle        : NULL
    ##   ..$ lineheight   : NULL
    ##   ..$ margin       : NULL
    ##   ..$ debug        : NULL
    ##   ..$ inherit.blank: logi FALSE
    ##   ..- attr(*, "class")= chr [1:2] "element_text" "element"
    ##  - attr(*, "class")= chr [1:2] "theme" "gg"
    ##  - attr(*, "complete")= logi FALSE
    ##  - attr(*, "validate")= logi TRUE

``` r
#Repeat the same plot for the bitscore >=50 and bitscore >=200
```

**Normalized bitscores distribution of the environment**

``` r
mag_env_samp <- mag_env_groups %>%
  #filter(Bit_score >= 1000) %>%
  select(Pident_matches, Alignment_Length, Bit_score, cultured_level, phylum, taxonomic_dist, habitat) %>%
  mutate(norm_bitscore = Bit_score/Alignment_Length)%>%
  group_by(cultured_level) %>%
  sample_n(1000)

magenv_bitsplot <- ggplot(mag_env_samp, aes(x= habitat, y= norm_bitscore)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.3, height = 0), alpha = 0.03) +
  ylab("Normalized Bit scores") +
  xlab("Dominant MAGs phyla") +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=12),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) +
  facet_grid(~cultured_level)

print(magenv_bitsplot)
```

![](genomic_novelty_allresults_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
