library(tidyverse)
library(hrbrthemes)
library(ggtext)

theme_set(
  theme_minimal() +
  theme(
    plot.title.position = "plot",
    plot.margin = margin(25,25,25,25),
    axis.title.x = element_markdown(hjust = .5, size = 12),
    axis.title.y = element_markdown(hjust = .5, size = 12),
		legend.position = "top",
		panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, colour = 'grey20')
  )
)

options(readr.show_col_types = FALSE)


clean_alleles <- function(s) {
  s |>
    str_replace_all('"','') |>
    str_replace('\\[','') |>
    str_replace('\\]','') 
}

get_allele_frequency <- function(s){
  s |>
    str_replace('.*"AF":\\[','') |>
    str_replace('\\].*','') 
}

process_variant_qc <- function(my_tibble, prs_snps){
  my_tibble |>
    mutate(alleles = clean_alleles(alleles)) |>
    mutate(af = get_allele_frequency(variant_qc)) |>
    dplyr::select(locus, alleles, af) |>
    filter(str_count(af,',') == str_count(alleles,',')) |>
    mutate(allele = str_split(alleles, ",")) |>
    mutate(af = str_split(af, ",")) |>
    unnest(c(allele, af)) |>
    mutate(id = str_c(locus, allele, sep = ':')) |>
    distinct(id, af) |>
    group_by(id) |> mutate(n = n()) |>
    filter(n == 1) |> 
    filter(id %in% prs_snps$id_hg38) |>
    dplyr::select(-n)
}

prs_snps <-
  read_tsv('mavaddat2019-313snps-liftover.tsv', col_types = 'ccccccc') |>
    mutate(id_hg19 = str_c("chr",old_locus,":",allele2)) |>
    mutate(id_hg38 = str_c(locus,':',allele2))

my_data_grar <-
  read_tsv('grar313_variantqc.tsv') |>
  process_variant_qc(prs_snps) 

my_data_sabe <-
  read_tsv('sabe313_variantqc.tsv') |>
  process_variant_qc(prs_snps) 

my_data_afr <-
  read_tsv('afr1kgp3_variantqc.tsv') |>
  mutate(variant_qc = info) |>
  process_variant_qc(prs_snps) 

my_data_eur <-
  read_tsv('eur1kgp3_variantqc.tsv') |>
  mutate(variant_qc = info) |>
  process_variant_qc(prs_snps) 

my_data_eas <-
  read_tsv('eas1kgp3_variantqc.tsv') |>
  mutate(variant_qc = info) |>
  process_variant_qc(prs_snps) 

snps_in_common <-
  tibble(
    id = c(
      my_data_grar$id, 
      my_data_sabe$id, 
      my_data_eur$id, 
      my_data_afr$id, 
      my_data_eas$id
    )
  ) |>
  count(id) |>
  filter(n == 5) |>
  pull(id)


# SABE/GRAR vs 1KGP3
data2plot <- bind_rows(
  my_data_grar |>
    filter(id %in% snps_in_common) |>
    left_join(my_data_eur |> rename(af_ref_pop = 'af'), by = 'id') |>
    mutate(brazil_pop = 'GRAR', ref_pop = 'EUR'),
  my_data_grar |>
    filter(id %in% snps_in_common) |>
    left_join(my_data_afr |> rename(af_ref_pop = 'af'), by = 'id') |>
    mutate(brazil_pop = 'GRAR', ref_pop = 'AFR'),
  my_data_grar |>
    filter(id %in% snps_in_common) |>
    left_join(my_data_eas |> rename(af_ref_pop = 'af'), by = 'id') |>
    mutate(brazil_pop = 'GRAR', ref_pop = 'EAS'),
  my_data_sabe |>
    filter(id %in% snps_in_common) |>
    left_join(my_data_eur |> rename(af_ref_pop = 'af'), by = 'id') |>
    mutate(brazil_pop = 'SABE', ref_pop = 'EUR'),
  my_data_sabe |>
    filter(id %in% snps_in_common) |>
    left_join(my_data_afr |> rename(af_ref_pop = 'af'), by = 'id') |>
    mutate(brazil_pop = 'SABE', ref_pop = 'AFR'),
  my_data_sabe |>
    filter(id %in% snps_in_common) |>
    left_join(my_data_eas |> rename(af_ref_pop = 'af'), by = 'id') |>
    mutate(brazil_pop = 'SABE', ref_pop = 'EAS')) |>
  mutate(af = as.double(af), af_ref_pop = as.double(af_ref_pop)) |>
  mutate(delta = abs(af_ref_pop - af))

pearson_r <- 
  data2plot |>
  group_by(ref_pop, brazil_pop) |>
  summarise(r = cor(af, af_ref_pop, method = 'pearson'))

p <-
  data2plot |>
  ggplot(aes(x = af, y = af_ref_pop, )) +
  geom_point(size = .6, aes(color = delta)) +
  geom_abline(color = 'red') +
  geom_text(data = pearson_r, 
            x = .8, 
            y = .1, 
            aes(label = scales::comma(r, accuracy = .01)), 
            size = 3.5) +
  viridis::scale_color_viridis(direction = -1) +
  facet_grid(brazil_pop~ref_pop) +
  coord_fixed() +
  theme(axis.ticks = element_line(color = 'black')) +
  scale_x_continuous(labels = function(x) x * 100, breaks = c(0,.50,1)) +
  scale_y_continuous(labels = function(x) x * 100, breaks = c(0,.50,1)) +
  labs(x = 'RAF (%)', y = 'RAF (%)')

ggsave(p, filename = 'freq.pdf', device = 'pdf', width = 4, height = 4)

# SABE vs GRAR
data2plot2 <-
  my_data_grar |> mutate(af = as.numeric(af)) |> rename(af_grar = 'af') |>
  inner_join(my_data_sabe |> mutate(af = as.numeric(af)) |> rename(af_sabe = 'af'), by = 'id')

my_r <- cor(data2plot2$af_grar, data2plot2$af_sabe) 

p2 <-
  data2plot2 |>
  mutate(delta = abs(af_sabe - af_grar)) |>
  ggplot(aes(x = af_sabe, y = af_grar)) +
  geom_point(aes(color = delta)) +
  geom_abline(color = 'red') +
  annotate(geom = 'text',
            x = .8,
            y = .1,
            label = scales::comma(my_r, accuracy = .01),
            size = 3.5) +
  viridis::scale_color_viridis(direction = -1, limits = c(0,.5)) +
  coord_fixed() +
  theme(axis.ticks = element_line(color = 'black')) +
  scale_x_continuous(labels = function(x) x * 100) +
  scale_y_continuous(labels = function(x) x * 100) +
  labs(x = 'RAF SABE (%)', y = 'RAF GRAR (%)')

ggsave(p2, filename = 'freq_sabegrar.pdf', device = 'pdf', width = 4, height = 4)