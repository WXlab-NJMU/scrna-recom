library(readxl)
data <- readxl::read_excel("data-raw/Cell_marker_Seq.xlsx")
cellmarkers <- data %>% 
  filter(Symbol != "")  %>% 
  filter(PMID != "")  %>% 
  filter(cellontology_id != "")  %>% 
  mutate(across('cell_name', ~ gsub(' ', '_', .x))) %>%
  mutate(across('cell_name', ~ gsub('\\.', '_', .x))) %>%
  mutate(across('cell_name', ~ gsub('\\+', 'p_', .x))) %>%
  mutate(across('cell_name', ~ gsub('\\-', 'n_', .x))) %>%
  mutate(across('cell_name', ~ gsub(r"{\s*\([^\)]+\)}","", .x))) %>% 
  group_by(tissue_class, cell_name) %>%
  summarise(marker_genes = paste(Symbol, collapse = "|"), .groups = "keep")
usethis::use_data(cellmarkers, overwrite = T)
# set internal to TRUE if data is too large
##usethis::use_data(cellmarkers, internal = TRUE)

