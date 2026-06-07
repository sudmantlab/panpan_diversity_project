setwd("/Users/joanocha/Google Drive/My Drive/POSTDOC/PANPAN/analysis/Figure3_VEP_snvs_svs/snvs_vep/")
library(vcfR)
library(tidyverse)
library(janitor)
library(data.table)
library(ggupset)
library(ggh4x)
library(knitr)
library(cowplot)
library(dplyr)
library(vcfR)

############## CROSS_SPECIES ANALYSIS ###################

# --- FILE PATHS ---
vep_human_file <- "ht2t/hprc_mapped2ht2t.phased_snvs.vep.txt.gz"
vep_allhuman_file <- "ht2t/hprc_hgsvc_mapped2ht2t.phased_snvs.vep.txt.gz"
vep_chimp_file <- "mPanTro3/lifted_pantros_mapped2mPanTro3.phased_snvs.vep.txt" ### nicolas will make lifted vcf for svs
vep_bonobo_file <- "mPanPan1/panpa_mapped2mPanTro3.phased_snvs.vep.txt"

gff_human_file <- "ht2t/ht2t.gff.gz"
gtf_chimp_file <- "mPanTro3/mPanTro3.gtf.gz"
gtf_bonobo_file <- "mPanPan1/mPanPan1.gtf.gz"

# --- WRANGLE HUMAN ---
## human SNV VEP results - Read and process human VEP data
vep_human <- fread(cmd = paste("zgrep -v '^##'", vep_allhuman_file)) %>%
  clean_names() %>%
  separate_rows(extra, sep = ";") %>%
  separate(extra, into = c("key", "value"), sep = "=", fill = "right") %>%
  pivot_wider(names_from = key, values_from = value) %>%
  clean_names() %>%
  mutate(impact = fct_relevel(impact, c("HIGH", "MODERATE", "LOW", "MODIFIER"))) %>%
  mutate(species = "human")
gff_human <- read_tsv(gff_human_file, comment = "#", col_names = FALSE) %>%
  filter(X3 %in% c("gene", "pseudogene")) %>%
  mutate(index = row_number()) %>%
  dplyr::select(index, X2, X3, X9) %>%
  separate_rows(X9, sep = ";") %>%
  separate(X9, into = c("key", "value"), sep = "=") %>%
  filter(key %in% c("Dbxref", "ID", "Name", "description")) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  clean_names() %>%
  dplyr::select(-index) %>%
  mutate(gene = str_remove(dbxref, ".*GeneID:") %>% str_remove(",.*")) %>%
  rename(annotation_tool = x2, type = x3)
vep_human_by_gene <- vep_human %>%
  filter(gene != "-") %>%
  separate_rows(consequence, sep = ",") %>%
  distinct(gene, impact, consequence, species, number_uploaded_variation) %>%
  filter(impact == sort(unique(impact))[1], .by = c(species, gene)) %>%
  summarise(
    consequence = str_c(unique(consequence), collapse = ","),
    variant_id = str_c(unique(number_uploaded_variation), collapse = ","),
    impact = sort(impact)[1],
    .by = c(species, gene)
  ) %>%
  left_join(gff_human, by = "gene")
vep_human_by_gene_filtered <- vep_human_by_gene %>%
  dplyr::select(species, type, name, description, impact, consequence, variant_id) %>%
  filter(
    type == "gene",
    !(str_detect(name, "^LOC") & is.na(description)),
    !(str_detect(name, "^LOC") & str_detect(description, "^uncharacterized")),
    !(str_detect(name, "^LOC") & str_detect(description, "^putatively uncharacterized")),
    !(str_detect(name, "^LOC") & str_detect(description, "^putative uncharacterized")),
    !str_detect(name, "^TRNA"),
    !str_detect(name, "^LINC"),
    !str_detect(name, "^SNORA"),
    !str_detect(name, "^RNA5S"),
    !str_detect(description, "ribosomal RNA$"),
    !str_detect(description, "^microRNA")
  ) %>%
  distinct()
head(vep_human_by_gene_filtered)


# --- WRANGLE CHIMP ---
# Read and process chimp VEP data
vep_chimp <- fread(cmd = paste("zgrep -v '^##'", vep_chimp_file)) %>%
  clean_names() %>%
  separate_rows(extra, sep = ";") %>%
  separate(extra, into = c("key", "value"), sep = "=", fill = "right") %>%
  pivot_wider(names_from = key, values_from = value) %>%
  clean_names() %>%
  mutate(impact=fct_relevel(impact, c("HIGH", "MODERATE", "LOW", "MODIFIER"))) %>%
  mutate(species="chimpanzee")
gtf_chimp <- read_tsv(gtf_chimp_file, comment = "#", col_names = FALSE) %>%
  filter(X3 %in% c("gene", "pseudogene")) %>%
  mutate(index=row_number()) %>%
  dplyr::select(index, X2, X3, X9) %>%
  separate_rows(X9, sep = "; ") %>%
  separate(X9, into = c("key", "value"), sep = " ", extra = "merge") %>%
  mutate(value=str_remove_all(value, "\"")) %>%
  filter(key %in% c("gene_id", "description")) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  clean_names() %>%
  dplyr::select(-index) %>%
  #mutate(gene=str_remove(dbxref, ".*GeneID:") %>% str_remove(",.*")) %>%
  rename(gene=gene_id, annotation_tool=x2, type=x3)
vep_chimp_by_gene <- vep_chimp %>%
  filter(gene != "-") %>%
  separate_rows(consequence, sep = ",") %>%
  distinct(gene, impact, consequence, species, number_uploaded_variation) %>%
  filter(impact==sort(unique(impact))[1], .by=c(species, gene)) %>%
  summarise(consequence=str_c(unique(consequence), collapse = ","), variant_id=str_c(unique(number_uploaded_variation), collapse = ","), impact=sort(impact)[1], .by=c(species, gene)) %>%
  left_join(gtf_chimp, by="gene") 
# Filter out uncharacterized/non-coding genes
vep_chimp_by_gene_filtered <- vep_chimp_by_gene %>%
  select(species, type, gene, description, impact, consequence, variant_id) %>%
  rename(name = gene) %>%
  filter(
    type == "gene",
    !(str_detect(name, "^LOC") & is.na(description)),
    !(str_detect(name, "^LOC") & str_detect(description, "^uncharacterized")),
    !(str_detect(name, "^LOC") & str_detect(description, "^putatively uncharacterized")),
    !(str_detect(name, "^LOC") & str_detect(description, "^putative uncharacterized")),
    !str_detect(name, "^TRNA"),
    !str_detect(name, "^LINC"),
    !str_detect(name, "^SNORA"),
    !str_detect(name, "^RNA5S"),
    !str_detect(description, "ribosomal RNA$"),
    !str_detect(description, "^microRNA")
  ) %>%
  distinct()


# --- WRANGLE BONOBO --- #
vep_bonobo <- fread(cmd = paste("zgrep -v '^##'", vep_bonobo_file))%>% 
  clean_names() %>%
  separate_rows(extra, sep = ";") %>%
  separate(extra, into = c("key", "value"), sep = "=", fill = "right") %>%
  pivot_wider(names_from = key, values_from = value) %>%
  clean_names() %>%
  mutate(impact=fct_relevel(impact, c("HIGH", "MODERATE", "LOW", "MODIFIER"))) %>%
  mutate(species="bonobo")
gtf_bonobo <- read_tsv(gtf_bonobo_file, comment = "#", col_names = FALSE) %>%
  filter(X3 %in% c("gene", "pseudogene")) %>%
  mutate(index=row_number()) %>%
  dplyr::select(index, X2, X3, X9) %>%
  separate_rows(X9, sep = "; ") %>%
  separate(X9, into = c("key", "value"), sep = " ", extra = "merge") %>%
  mutate(value=str_remove_all(value, "\"")) %>%
  filter(key %in% c("gene_id", "description")) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  clean_names() %>%
  dplyr::select(-index) %>%
  #mutate(gene=str_remove(dbxref, ".*GeneID:") %>% str_remove(",.*")) %>%
  rename(gene=gene_id, annotation_tool=x2, type=x3)
vep_bonobo_by_gene <- vep_bonobo %>%
  filter(gene != "-") %>%
  separate_rows(consequence, sep = ",") %>%
  distinct(gene, impact, consequence, species, number_uploaded_variation) %>%
  filter(impact==sort(unique(impact))[1], .by=c(species, gene)) %>%
  summarise(consequence=str_c(unique(consequence), collapse = ","), variant_id=str_c(unique(number_uploaded_variation), collapse = ","), impact=sort(impact)[1], .by=c(species, gene)) %>%
  left_join(gtf_bonobo, by="gene")
vep_bonobo_by_gene_filtered <- vep_bonobo_by_gene %>%
  #filter(impact=="HIGH") %>%
  dplyr::select(species, type, gene, description, impact, consequence, variant_id) %>%
  rename(name=gene) %>%
  filter(
    type == "gene",
    !(str_detect(name, "^LOC") & is.na(description)),
    !(str_detect(name, "^LOC") & str_detect(description, "^uncharacterized")),
    !(str_detect(name, "^LOC") & str_detect(description, "^putatively uncharacterized")),
    !(str_detect(name, "^LOC") & str_detect(description, "^putative uncharacterized")),
    !str_detect(name, "^TRNA"),
    !str_detect(name, "^LINC"),
    !str_detect(name, "^SNORA"),
    !str_detect(name, "^RNA5S"),
    !str_detect(description, "ribosomal RNA$"),
    !str_detect(description, "^microRNA")
  ) %>%
  distinct()

# --- WRANGLE ALL --- #
vep_all_species_by_gene_wide <- vep_human_by_gene_filtered %>% 
  dplyr::select(species, name, description, impact) %>%
  bind_rows(vep_chimp_by_gene_filtered %>% dplyr::select(species, name, description, impact)) %>%
  bind_rows(vep_bonobo_by_gene_filtered %>% dplyr::select(species, name, description, impact)) %>%
  dplyr::select(-description) %>%
  mutate(impact=(5-as.numeric(impact))) %>%
  pivot_wider(names_from = species, values_from = impact, values_fill = 0) %>%
  left_join(gff_human %>% filter(!is.na(description)) %>% dplyr::select(name, description) %>% distinct()) %>%
  left_join(gtf_chimp %>% filter(!is.na(description)) %>% transmute(name=gene, description_chimp=description) %>% distinct()) %>%
  left_join(gtf_bonobo %>% filter(!is.na(description)) %>% transmute(name=gene, description_bonobo=description) %>% distinct()) %>%
  mutate(description=ifelse(is.na(description), description_chimp, description),
         description=ifelse(is.na(description), description_bonobo, description),) %>%
  dplyr::select(-description_chimp, -description_bonobo)
vep_all_species_by_gene_wide %>%
  mutate(human=human>0, chimpanzee=chimpanzee>0, bonobo=bonobo>0) %>%
  count(human, chimpanzee, bonobo)
vep_all_species_by_gene_wide %>%
  mutate(human=human>1, chimpanzee=chimpanzee>1, bonobo=bonobo>1) %>%
  count(human, chimpanzee, bonobo)
vep_all_species_by_gene_wide %>%
  filter((human>3)+(chimpanzee>3)+(bonobo>3)>1)%>%
  arrange(name)
vep_all_species_by_gene_wide %>%
  filter((human>1)+(chimpanzee>1)+(bonobo>1)>2) 
head(vep_all_species_by_gene_wide)

# --- OUTPUT GENE TABLE --- #
vep_all_species_by_gene_wide %>%
  write_tsv("gene_impact_wide_allhumans.tsv")

vep_all_species_by_gene_wide %>%
  filter((human>3)+(chimpanzee>3)+(bonobo>3)>1) %>%
  arrange(name) %>%
  write_csv("genes_impacted_by_svns_shared_across_species_allhumans.csv")


############## SUBSPECIES ANALYSIS ###################

### needs vep_chimp (defined above)
#vep_chimp
### needs vcf_for_vep_join 
#vcf_for_vep_join<-fread("vcf_for_vep_join.tsv")

### needs vcf_long_chimp_by_species 
#vcf_long_chimp_by_specie<-fread("mPanTro3/vcf_long_chimp_by_species.tsv")

#### output the SNV table and define snv_impact_long
snv_impact_long <- bind_rows(vep_bonobo, vep_chimp, vep_human %>% left_join(gff_human %>% distinct(gene, name)) %>% mutate(gene=name) %>% dplyr::select(-name)) %>% 
  transmute(species, gene, impact, id=number_uploaded_variation, feature, consequence) %>%
  group_by(species, id) %>%
  filter(impact==unique(impact) %>% sort() %>% .[1]) %>%
  summarise(impact=unique(impact),
            gene=unique(gene) %>% str_c(collapse = ","), 
            feature=unique(feature) %>% str_c(collapse = ","), 
            consequence=str_c(consequence, collapse = ",") %>% str_split(pattern=",") %>% lapply(unique) %>% map_chr(function(x){str_c(x, collapse = ",")})) %>%
  ungroup() %>%
  left_join(vcf_for_vep_join, by=c("id", "species"))
write_tsv(snv_impact_long, "snv_impact_long_allhuman.tsv")

#### sharing of variant effect on the variant level
variant_level_sharing <- snv_impact_long %>%
  filter(species=="chimpanzee") %>%
  dplyr::select(-species) %>%
  left_join(vcf_long_chimp_by_species %>% dplyr::select(id, common_name, allele_count, sample_size), ., by="id") %>%
  filter(impact=="HIGH") 
write_tsv(variant_level_sharing, "mPanTro3/high_impact_variant_by_subspecies.tsv")

#### sharing of variant effect on the gene level
gene_level_sharing <- vep_chimp %>%
  filter(gene != "-") %>%
  separate_rows(consequence, sep = ",") %>%
  distinct(gene, impact, consequence, number_uploaded_variation) %>%
  inner_join(vcf_long_chimp_by_species %>% filter(allele_count>0, allele_count<(sample_size*2)) %>% transmute(number_uploaded_variation=id, common_name), ., relationship = "many-to-many") %>%
  filter(impact==sort(unique(impact))[1], .by=c(common_name, gene)) %>%
  summarise(consequence=str_c(unique(consequence), collapse = ","), variant_id=str_c(unique(number_uploaded_variation), collapse = ","), impact=sort(impact)[1], .by=c(common_name, gene)) %>%
  left_join(gtf_chimp, by="gene") %>%
  dplyr::select(-annotation_tool, -type)
gene_level_sharing_wide <- gene_level_sharing %>%
  dplyr::select(gene, description, common_name, impact) %>%
  mutate(impact=(5-as.numeric(impact))) %>%
  pivot_wider(names_from = common_name, values_from = impact, values_fill = 0) 
gene_level_sharing_wide %>%
  write_tsv("mPanTro3/gene_impact_wide_by_subspecies.tsv")




