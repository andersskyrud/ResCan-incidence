##################################################
######  -- Data formatting OUS lab data -- #######
##################################################

## -- Setup -- ##

# Load necessary libraries
pacman::p_load(
  tidyverse,
  markdown,
  haven,
  janitor,
  gtsummary,
  broom,
  ggplot2,
  viridis,
  RColorBrewer,
  lubridate,
  rio,
  here
)

# Setting the working directory
setwd("~")

# Import lab data

miclis <- import("ResCanMiclis_PID.xlsx")
unilab <- import_list(c("ResCanUnilab_PID_1.xlsx", "ResCanUnilab_PID_2.xlsx"), 
                      rbind = T, readxl = FALSE)
swisslab <- import("ResCanSwisslab_PID.xlsx", readxl = FALSE)

# Import cancer data

study1_data <- import("~/study1_data.rds")

# Creating a diagnosed date variable, linking it to lab data, 
# and keeping only samples taken from one year before to three years after cancer diagnosis

# Fetch first cancer occurrence

diagnosedato <- study1_data %>% 
  select(PID, DIAGNOSEDATO, diagnose_3) 

# Process for Unilab

unilab <- unilab %>%
  mutate(PID = as.character(PID)) %>%
  left_join(diagnosedato, by = "PID", keep = T) %>%
  rename(PID = PID.x) %>%
  # This variable is directly extracted from the SQL database as UNIX time, 
  # but it obviously originates from 1st January 1900
  mutate(OrderDate_date = as_date(OrderDate, origin = "1900-01-01")) %>%
  filter(OrderDate_date >= DIAGNOSEDATO - months(6) & OrderDate_date <= DIAGNOSEDATO + years(3))

# Process for Swisslab

swisslab <- swisslab %>%
  mutate(PID = as.character(PID)) %>%
  left_join(diagnosedato, by = "PID", keep = T) %>%
  rename(PID = PID.x) %>%
  # Same here
  mutate(OrderDate_date = as_date(Ordredato, origin = "1900-01-01")) %>%
  filter(OrderDate_date >= DIAGNOSEDATO - months(6) & OrderDate_date <= DIAGNOSEDATO + years(3))

# Process for Miclis

miclis <- miclis %>%
  mutate(PID = as.character(PID)) %>%
  left_join(diagnosedato, by = "PID", keep = T) %>%
  rename(PID = PID.x) %>%
  filter(SAMPLETIME >= DIAGNOSEDATO - months(6) & SAMPLETIME <= DIAGNOSEDATO + years(3))

# Selecting only blood and urine samples

# Unilab 

unilab <- unilab %>%
  group_by(CDOrderID) %>%
  filter(any(str_detect(Name.1, regex("blod(?!sk)|urin", ignore_case = TRUE)))) %>%
  ungroup()

# Swisslab 

swisslab <- swisslab %>%
  group_by(Ordrenummer) %>%
  filter(any(str_detect(Materialnavn, regex("blod|urin", ignore_case = TRUE)))) %>%
  ungroup()

# Miclis 

miclis <- miclis %>%
  group_by(REQID) %>%
  filter(any(str_detect(NAME, regex("blod|urin", ignore_case = TRUE)))) %>%
  ungroup()

# Unilab

# If 'ResultCode' is NA and the following 'ResultCode' is 'UT', set it to 'SP'
unilab$ResultCode[is.na(unilab$ResultCode) & lead(unilab$ResultCode)=="UT"] <- "SP"

# Selecting codes that often appear together with (positive results for) species
unilab_species <- unilab %>% 
  # Filtering on result codes indicating various levels of bacterial growth
  filter(ResultCode == "RV" |  # Rich growth
           ResultCode == "MV" |  # Moderate growth
           ResultCode == "VF" |  # Usual flora
           ResultCode == "VE" |  # Growth
           (ResultCode == "SV" & Word == "Sparsom vekst") |  # Sparse growth
           ResultCode == "EN" |  # Some colonies
           ResultCode == "SP" |  # Species from malditof etc.
           grepl("[0-9]", ResultCode)) %>%  # All result codes with numbers, includes cultivation in broth 
  # Removing all rate codes and incorrectly categorized NA rows
  filter(Name.1 != "Aerob dyrkning") %>% 
  filter(Name.1 != "Bakteriologisk dyrkning urin") %>%
  filter(Name.1 != "Bakterier dyrkning") %>%
  filter(Name.1 != "Buljonganrikning") %>%
  filter(Name.1 != "Anaerob dyrkning") 

unilab_NA <- unilab %>% 
  filter(is.na(ResultCode)) # Some species are found under NA, species names do not appear without positive test

# Assign species names based on pattern in the name.1 variable for unilab_species and unilab_NA datasets
unilab_species <- unilab_species %>%
  mutate(Species = case_when(
    str_detect(string = Name.1, pattern = regex("escherichia coli", ignore_case = T)) ~ "Escherichia coli",
    str_detect(string = Name.1, pattern = regex("staphylococcus aureus|(Staph. aureus)", ignore_case = T)) ~ "Staphylococcus aureus",
    str_detect(string = Name.1, pattern = regex("blandings|gram positive (.*) gram negative|gram negative (.*) gram positive", ignore_case = T)) ~ "Blandingsvekst",
    str_detect(string = Name.1, pattern = regex("enterococcus faeca", ignore_case = T)) ~ "Enterococcus faecalis",
    str_detect(string = Name.1, pattern = regex("enterococcus faeci", ignore_case = T)) ~ "Enterococcus faecium",
    str_detect(string = Name.1, pattern = regex("candida albicans", ignore_case = T)) ~ "Candida albicans",
    str_detect(string = Name.1, pattern = regex("normal|vanlig (.*)flora", ignore_case = T)) ~ "Normal flora",
    str_detect(string = Name.1, pattern = regex("(?!.gram positive*)gram negative(?!.*gram positive)", ignore_case = T)) ~ "Gram-negative bakterier",
    str_detect(string = Name.1, pattern = regex("(?!.gram negative*)gram positive(?!.*gram negative)|enterococcus(?! faec)", ignore_case = T)) ~ "Gram-positive bakterier",
    str_detect(string = Name.1, pattern = regex("klebsiella pneum", ignore_case = T)) ~ "Klebsiella pneumoniae",
    str_detect(string = Name.1, pattern = regex("klebsiella(?! pneum)", ignore_case = T)) ~ "Klebsiella spp.",
    str_detect(string = Name.1, pattern = regex("Pseudomonas aeruginosa", ignore_case = T)) ~ "Pseudomonas aeruginosa",
    str_detect(string = Name.1, pattern = regex("haemophilus influenzae", ignore_case = T)) ~ "Haemophilus influenzae",
    str_detect(string = Name.1, pattern = regex("streptococcus pneum|pneumokokk", ignore_case = T)) ~ "Streptococcus pneumoniae",
    str_detect(string = Name.1, pattern = regex("Proteus", ignore_case = T)) ~ "Proteus spp.",
    str_detect(string = Name.1, pattern = regex("(beta-hemolytiske gruppe A)|( streptokokker gr. A)|(Streptococcus pyogenes)", ignore_case = T)) ~ "Streptokokker gruppe A",
    str_detect(string = Name.1, pattern = regex("(beta-hemolytiske gruppe C)|( streptokokker gr. C)", ignore_case = T)) ~ "Streptokokker gruppe C",
    str_detect(string = Name.1, pattern = regex("(beta-hemolytiske gruppe B)|(Streptococcus agal)|( streptokokker gr. B)", ignore_case = T)) ~ "Streptokokker gruppe B",
    str_detect(string = Name.1, pattern = regex("(beta-hemolytiske gruppe G)|(streptokokker gr. G)", ignore_case = T)) ~ "Streptokokker gruppe G",
    str_detect(string = Name.1, pattern = regex("^Streptococcus(?! pneum)(?! pyog)(?! agal)|(Streptokokker, alfa-hemolytiske)|alfahemolystiske", ignore_case = T)) ~ "Streptococcus spp.",
    str_detect(string = Name.1, pattern = regex("Moraxella catarrhalis", ignore_case = T)) ~ "Moraxella catarrhalis",
    str_detect(string = Name.1, pattern = regex("Enterobacter cloacae", ignore_case = T)) ~ "Enterobacter cloacae",
    str_detect(string = Name.1, pattern = regex("Candida(?! albi)", ignore_case = T)) ~ "Candida spp.",
    str_detect(string = Name.1, pattern = regex("Enterobacter(?! cloa)", ignore_case = T)) ~ "Enterobacter spp.",
    str_detect(string = Name.1, pattern = regex("Staphylococcus(?! aur)|stafylokokker, hvite|koagulase-negative", ignore_case = T)) ~ "Staphyloccocus spp.", 
    str_detect(string = Name.1, pattern = regex("aerococcus", ignore_case = T)) ~ "Aerococcus spp.",
    str_detect(string = Name.1, pattern = regex("serratia", ignore_case = T)) ~ "Serratia spp.",
    str_detect(string = Name.1, pattern = regex("pseudomonas(?! aeru)", ignore_case = T)) ~ "Pseudomonas spp.",
    str_detect(string = Name.1, pattern = regex("(Clostridium difficile)|(Clostridioides difficile)", ignore_case = T)) ~ "Clostridioides difficile",
    str_detect(string = Name.1, pattern = regex("Stenotrophomonas maltophilia", ignore_case = T)) ~ "Stenotrophomonas maltophilia",
    str_detect(string = Name.1, pattern = regex("Propionibacterium", ignore_case = T)) ~ "Propionibacterium spp.",
    str_detect(string = Name.1, pattern = regex("citrobacter", ignore_case = T)) ~ "Citrobacter spp.",
    str_detect(string = Name.1, pattern = regex("Morganella", ignore_case = T)) ~ "Morganella spp.",
    str_detect(string = Name.1, pattern = regex("syrefaste|mycobacterium|mykobakt", ignore_case = T)) ~ "Mycobacterium spp. (el. syrefaste staver)",
    str_detect(string = Name.1, pattern = regex("bacteroides", ignore_case = T)) ~ "Bacteroides spp.",
    str_detect(string = Name.1, pattern = regex("acinetobacter", ignore_case = T)) ~ "Acinetobacter spp.",
    str_detect(string = Name.1, pattern = regex("Corynebacterium", ignore_case = T)) ~ "Corynebacterium spp.",
    str_detect(string = Name.1, pattern = regex("campylobacter", ignore_case = T)) ~ "Campylobacter spp.",
    str_detect(string = Name.1, pattern = regex("(aspergillus fumigatus)|fumigatus", ignore_case = T)) ~ "Aspergillus fumigatus",
    str_detect(string = Name.1, pattern = regex("Lactobacillus", ignore_case = T)) ~ "Lactobacillus spp.",
    str_detect(string = Name.1, pattern = regex("escherichia(?! coli)", ignore_case = T)) ~ "Andre Escherichia spp.",
    str_detect(string = Name.1, pattern = regex("Lactobacillus", ignore_case = T)) ~ "Lactobacillus spp.",
    str_detect(string = Name.1, pattern = regex("salmonella", ignore_case = T)) ~ "Salmonella spp.",
    TRUE ~ "Andre funn"
  )) 

unilab_NA <- unilab_NA %>% 
  mutate(Species = case_when(
    str_detect(string = Name.1, pattern = regex("escherichia coli", ignore_case = T)) ~ "Escherichia coli",
    str_detect(string = Name.1, pattern = regex("staphylococcus aureus|(Staph. aureus)", ignore_case = T)) ~ "Staphylococcus aureus",
    str_detect(string = Name.1, pattern = regex("blandings|gram positive (.*) gram negative|gram negative (.*) gram positive", ignore_case = T)) ~ "Blandingsvekst",
    str_detect(string = Name.1, pattern = regex("enterococcus faeca", ignore_case = T)) ~ "Enterococcus faecalis",
    str_detect(string = Name.1, pattern = regex("enterococcus faeci", ignore_case = T)) ~ "Enterococcus faecium",
    str_detect(string = Name.1, pattern = regex("candida albicans", ignore_case = T)) ~ "Candida albicans",
    str_detect(string = Name.1, pattern = regex("normal|vanlig (.*)flora", ignore_case = T)) ~ "Normal flora",
    str_detect(string = Name.1, pattern = regex("(?!.gram positive*)gram negative(?!.*gram positive)", ignore_case = T)) ~ "Gram-negative bakterier",
    str_detect(string = Name.1, pattern = regex("(?!.gram negative*)gram positive(?!.*gram negative)|enterococcus(?! faec)", ignore_case = T)) ~ "Gram-positive bakterier",
    str_detect(string = Name.1, pattern = regex("klebsiella pneum", ignore_case = T)) ~ "Klebsiella pneumoniae",
    str_detect(string = Name.1, pattern = regex("klebsiella(?! pneum)", ignore_case = T)) ~ "Klebsiella spp.",
    str_detect(string = Name.1, pattern = regex("Pseudomonas aeruginosa", ignore_case = T)) ~ "Pseudomonas aeruginosa",
    str_detect(string = Name.1, pattern = regex("haemophilus influenzae", ignore_case = T)) ~ "Haemophilus influenzae",
    str_detect(string = Name.1, pattern = regex("streptococcus pneum|pneumokokk", ignore_case = T)) ~ "Streptococcus pneumoniae",
    str_detect(string = Name.1, pattern = regex("Proteus", ignore_case = T)) ~ "Proteus spp.",
    str_detect(string = Name.1, pattern = regex("(beta-hemolytiske gruppe A)|( streptokokker gr. A)|(Streptococcus pyogenes)", ignore_case = T)) ~ "Streptokokker gruppe A",
    str_detect(string = Name.1, pattern = regex("(beta-hemolytiske gruppe C)|( streptokokker gr. C)", ignore_case = T)) ~ "Streptokokker gruppe C",
    str_detect(string = Name.1, pattern = regex("(beta-hemolytiske gruppe B)|(Streptococcus agal)|( streptokokker gr. B)", ignore_case = T)) ~ "Streptokokker gruppe B",
    str_detect(string = Name.1, pattern = regex("(beta-hemolytiske gruppe G)|(streptokokker gr. G)", ignore_case = T)) ~ "Streptokokker gruppe G",
    str_detect(string = Name.1, pattern = regex("^Streptococcus(?! pneum)(?! pyog)(?! agal)|(Streptokokker, alfa-hemolytiske)|alfahemolystiske", ignore_case = T)) ~ "Streptococcus spp.",
    str_detect(string = Name.1, pattern = regex("Moraxella catarrhalis", ignore_case = T)) ~ "Moraxella catarrhalis",
    str_detect(string = Name.1, pattern = regex("Enterobacter cloacae", ignore_case = T)) ~ "Enterobacter cloacae",
    str_detect(string = Name.1, pattern = regex("Candida(?! albi)", ignore_case = T)) ~ "Candida spp.",
    str_detect(string = Name.1, pattern = regex("Enterobacter(?! cloa)", ignore_case = T)) ~ "Enterobacter spp.",
    str_detect(string = Name.1, pattern = regex("Staphylococcus(?! aur)|stafylokokker, hvite|koagulase-negative", ignore_case = T)) ~ "Staphyloccocus spp.", 
    str_detect(string = Name.1, pattern = regex("aerococcus", ignore_case = T)) ~ "Aerococcus spp.",
    str_detect(string = Name.1, pattern = regex("serratia", ignore_case = T)) ~ "Serratia spp.",
    str_detect(string = Name.1, pattern = regex("pseudomonas(?! aeru)", ignore_case = T)) ~ "Pseudomonas spp.",
    str_detect(string = Name.1, pattern = regex("(Clostridium difficile)|(Clostridioides difficile)", ignore_case = T)) ~ "Clostridioides difficile",
    str_detect(string = Name.1, pattern = regex("Stenotrophomonas maltophilia", ignore_case = T)) ~ "Stenotrophomonas maltophilia",
    str_detect(string = Name.1, pattern = regex("Propionibacterium", ignore_case = T)) ~ "Propionibacterium spp.",
    str_detect(string = Name.1, pattern = regex("citrobacter", ignore_case = T)) ~ "Citrobacter spp.",
    str_detect(string = Name.1, pattern = regex("Morganella", ignore_case = T)) ~ "Morganella spp.",
    str_detect(string = Name.1, pattern = regex("syrefaste|mycobacterium|mykobakt", ignore_case = T)) ~ "Mycobacterium spp. (el. syrefaste staver)",
    str_detect(string = Name.1, pattern = regex("bacteroides", ignore_case = T)) ~ "Bacteroides spp.",
    str_detect(string = Name.1, pattern = regex("acinetobacter", ignore_case = T)) ~ "Acinetobacter spp.",
    str_detect(string = Name.1, pattern = regex("Corynebacterium", ignore_case = T)) ~ "Corynebacterium spp.",
    str_detect(string = Name.1, pattern = regex("campylobacter", ignore_case = T)) ~ "Campylobacter spp.",
    str_detect(string = Name.1, pattern = regex("(aspergillus fumigatus)|fumigatus", ignore_case = T)) ~ "Aspergillus fumigatus",
    str_detect(string = Name.1, pattern = regex("Lactobacillus", ignore_case = T)) ~ "Lactobacillus spp.",
    str_detect(string = Name.1, pattern = regex("escherichia(?! coli)", ignore_case = T)) ~ "Andre Escherichia spp.",
    str_detect(string = Name.1, pattern = regex("Lactobacillus", ignore_case = T)) ~ "Lactobacillus spp.",
    str_detect(string = Name.1, pattern = regex("salmonella", ignore_case = T)) ~ "Salmonella spp."
  )) %>% 
  filter(!is.na(Species))

unilab_all_species <- 
  bind_rows(unilab_species, unilab_NA) %>% 
  rename(SampleDate = OrderDate_date) 
  

# Swisslab

# Filter the swisslab dataframe to exclude non-positive results
swisslab_positiv <- swisslab %>%
  filter(!str_detect(string = Funn, pattern = regex("ikke pevist|ingen|se kommentar|^Negativ$|^Negativ.$|mislykket|uavklart|inkonklusiv|usikkert", ignore_case = T)))

# Assign the species names to each record based on the patterns in 'Funn' field using regex
swisslab_positiv <- swisslab_positiv %>%
  mutate(Species = case_when(
    str_detect(string = Funn, pattern = regex("escherichia coli", ignore_case = T)) ~ "Escherichia coli",
    str_detect(string = Funn, pattern = regex("staphylococcus aureus", ignore_case = T)) ~ "Staphylococcus aureus",
    str_detect(string = Funn, pattern = regex("blandingsvekst|gram positive (.*) gram negative|gram negative (.*) gram positive", ignore_case = T)) ~ "Blandingsvekst",
    str_detect(string = Funn, pattern = regex("enterococcus faeca", ignore_case = T)) ~ "Enterococcus faecalis",
    str_detect(string = Funn, pattern = regex("enterococcus faeci", ignore_case = T)) ~ "Enterococcus faecium",    
    str_detect(string = Funn, pattern = regex("candida albicans", ignore_case = T)) ~ "Candida albicans",
    str_detect(string = Funn, pattern = regex("normal|vanlig (.*)flora", ignore_case = T)) ~ "Normal flora",
    str_detect(string = Funn, pattern = regex("(?!.gram positive*)gram negative(?!.*gram positive)", ignore_case = T)) ~ "Gram-negative bakterier",
    str_detect(string = Funn, pattern = regex("(?!.gram negative*)gram positive(?!.*gram negative)|enterococcus(?! faec)", ignore_case = T)) ~ "Gram-positive bakterier",
    str_detect(string = Funn, pattern = regex("klebsiella pneum", ignore_case = T)) ~ "Klebsiella pneumoniae",
    str_detect(string = Funn, pattern = regex("klebsiella(?! pneum)", ignore_case = T)) ~ "Klebsiella spp.",
    str_detect(string = Funn, pattern = regex("Pseudomonas aeruginosa", ignore_case = T)) ~ "Pseudomonas aeruginosa",
    str_detect(string = Funn, pattern = regex("haemophilus influenzae", ignore_case = T)) ~ "Haemophilus influenzae",
    str_detect(string = Funn, pattern = regex("streptococcus pneum|pneumokokk", ignore_case = T)) ~ "Streptococcus pneumoniae",
    str_detect(string = Funn, pattern = regex("Proteus", ignore_case = T)) ~ "Proteus spp.",
    str_detect(string = Funn, pattern = regex("(beta-hemolytiske gruppe A)|( streptokokker gr. A)|(Streptococcus pyogenes)", ignore_case = T)) ~ "Streptokokker gruppe A",
    str_detect(string = Funn, pattern = regex("(beta-hemolytiske gruppe C)|( streptokokker gr. C)", ignore_case = T)) ~ "Streptokokker gruppe C",
    str_detect(string = Funn, pattern = regex("(beta-hemolytiske gruppe B)|(Streptococcus agal)|( streptokokker gr. B)", ignore_case = T)) ~ "Streptokokker gruppe B",
    str_detect(string = Funn, pattern = regex("(beta-hemolytiske gruppe G)|(streptokokker gr. G)", ignore_case = T)) ~ "Streptokokker gruppe G",
    str_detect(string = Funn, pattern = regex("^Streptococcus(?! pneum)(?! pyog)(?! agal)|(Streptokokker, alfa-hemolytiske)|alfahemolystiske|(Streptokokker, non-hemolytiske)", ignore_case = T)) ~ "Streptococcus spp.",
    str_detect(string = Funn, pattern = regex("Moraxella catarrhalis", ignore_case = T)) ~ "Moraxella catarrhalis",
    str_detect(string = Funn, pattern = regex("Enterobacter cloacae", ignore_case = T)) ~ "Enterobacter cloacae",
    str_detect(string = Funn, pattern = regex("Candida(?! albi)", ignore_case = T)) ~ "Candida spp.",
    str_detect(string = Funn, pattern = regex("Enterobacter(?! cloa)", ignore_case = T)) ~ "Enterobacter spp.",
    str_detect(string = Funn, pattern = regex("Staphylococcus(?! aur)|stafylokokker, hvite", ignore_case = T)) ~ "Staphyloccocus spp.", # klarer ikke skille ut Staphylococcus epidermidis fordi "stafylokokker, hvite" er egen kategori
    str_detect(string = Funn, pattern = regex("aerococcus", ignore_case = T)) ~ "Aerococcus spp.",
    str_detect(string = Funn, pattern = regex("serratia", ignore_case = T)) ~ "Serratia spp.",
    str_detect(string = Funn, pattern = regex("pseudomonas(?! aeru)", ignore_case = T)) ~ "Pseudomonas spp.",
    str_detect(string = Funn, pattern = regex("Clostridium difficile", ignore_case = T)) ~ "Clostridioides difficile",
    str_detect(string = Funn, pattern = regex("Stenotrophomonas maltophilia", ignore_case = T)) ~ "Stenotrophomonas maltophilia",
    str_detect(string = Funn, pattern = regex("Propionibacterium", ignore_case = T)) ~ "Propionibacterium spp.",
    str_detect(string = Funn, pattern = regex("citrobacter", ignore_case = T)) ~ "Citrobacter spp.",
    str_detect(string = Funn, pattern = regex("Morganella", ignore_case = T)) ~ "Morganella spp.",
    str_detect(string = Funn, pattern = regex("syrefaste|mycobacterium|mykobakt", ignore_case = T)) ~ "Mycobacterium spp. (el. syrefaste staver)",
    str_detect(string = Funn, pattern = regex("bacteroides", ignore_case = T)) ~ "Bacteroides spp.",
    str_detect(string = Funn, pattern = regex("acinetobacter", ignore_case = T)) ~ "Acinetobacter spp.",
    str_detect(string = Funn, pattern = regex("Corynebacterium", ignore_case = T)) ~ "Corynebacterium spp.",
    str_detect(string = Funn, pattern = regex("campylobacter", ignore_case = T)) ~ "Campylobacter spp.",
    str_detect(string = Funn, pattern = regex("(aspergillus fumigatus)|fumigatus", ignore_case = T)) ~ "Aspergillus fumigatus",
    str_detect(string = Funn, pattern = regex("providencia", ignore_case = T)) ~ "Providencia spp.",
    str_detect(string = Funn, pattern = regex("escherichia(?! coli)", ignore_case = T)) ~ "Andre Escherichia spp.",
    str_detect(string = Funn, pattern = regex("Lactobacillus", ignore_case = T)) ~ "Lactobacillus spp.",
    str_detect(string = Funn, pattern = regex("salmonella", ignore_case = T)) ~ "Salmonella spp.",
    TRUE ~ "Andre funn"
  ))

swisslab_all_species <- swisslab_positiv %>%
  rename(SampleDate = OrderDate_date)

# Miclis

# Assign the species names to each record in the miclis dataset based on the patterns in 'MIC_NAME' field using regex
miclis_all_species <- miclis %>% 
  mutate(Species = case_when(
    str_detect(string = MIC_NAME, pattern = regex("escherichia coli", ignore_case = T)) ~ "Escherichia coli",
    str_detect(string = MIC_NAME, pattern = regex("staphylococcus aureus|(Staph. aureus)", ignore_case = T)) ~ "Staphylococcus aureus",
    str_detect(string = MIC_NAME, pattern = regex("blandings|gram positive (.*) gram negative|gram negative (.*) gram positive", ignore_case = T)) ~ "Blandingsvekst",
    str_detect(string = MIC_NAME, pattern = regex("enterococcus faeca|ENTEROCOCCUS FC???CALIS", ignore_case = T)) ~ "Enterococcus faecalis",
    str_detect(string = MIC_NAME, pattern = regex("enterococcus faeci", ignore_case = T)) ~ "Enterococcus faecium",     str_detect(string = MIC_NAME, pattern = regex("candida albicans", ignore_case = T)) ~ "Candida albicans",
    str_detect(string = MIC_NAME, pattern = regex("normal|vanlig (.*)flora", ignore_case = T)) ~ "Normal flora",
    str_detect(string = MIC_NAME, pattern = regex("(?!.gram positive*)gram negative(?!.*gram positive)", ignore_case = T)) ~ "Gram-negative bakterier",
    str_detect(string = MIC_NAME, pattern = regex("(?!.gram negative*)gram positive(?!.*gram negative)|enterococcus(?! faec)", ignore_case = T)) ~ "Gram-positive bakterier",
    str_detect(string = MIC_NAME, pattern = regex("klebsiella pneum", ignore_case = T)) ~ "Klebsiella pneumoniae",
    str_detect(string = MIC_NAME, pattern = regex("klebsiella(?! pneum)", ignore_case = T)) ~ "Klebsiella spp.",
    str_detect(string = MIC_NAME, pattern = regex("Pseudomonas aeruginosa", ignore_case = T)) ~ "Pseudomonas aeruginosa",
    str_detect(string = MIC_NAME, pattern = regex("haemophilus influenzae", ignore_case = T)) ~ "Haemophilus influenzae",
    str_detect(string = MIC_NAME, pattern = regex("streptococcus pneum|pneumokokk", ignore_case = T)) ~ "Streptococcus pneumoniae",
    str_detect(string = MIC_NAME, pattern = regex("Proteus", ignore_case = T)) ~ "Proteus spp.",
    str_detect(string = MIC_NAME, pattern = regex("(beta-hemolytiske gruppe A)|( streptokokker gr. A)|(Streptococcus pyogenes)", ignore_case = T)) ~ "Streptokokker gruppe A",
    str_detect(string = MIC_NAME, pattern = regex("(beta-hemolytiske gruppe C)|( streptokokker gr. C)|GR.C", ignore_case = T)) ~ "Streptokokker gruppe C",
    str_detect(string = MIC_NAME, pattern = regex("(beta-hemolytiske gruppe B)|(Streptococcus agal)|( streptokokker gr. B)|GR.B", ignore_case = T)) ~ "Streptokokker gruppe B",
    str_detect(string = MIC_NAME, pattern = regex("(beta-hemolytiske gruppe G)|(streptokokker gr. G)|GR.G", ignore_case = T)) ~ "Streptokokker gruppe G",
    str_detect(string = MIC_NAME, pattern = regex("^Streptococcus(?! pneum)(?! pyog)(?! agal)|(Streptokokker, alfa-hemolytiske)|alfahemolystiske|STREPTOKOKKER$", ignore_case = T)) ~ "Streptococcus spp.",
    str_detect(string = MIC_NAME, pattern = regex("Moraxella catarrhalis", ignore_case = T)) ~ "Moraxella catarrhalis",
    str_detect(string = MIC_NAME, pattern = regex("Enterobacter cloacae", ignore_case = T)) ~ "Enterobacter cloacae",
    str_detect(string = MIC_NAME, pattern = regex("Candida(?! albi)", ignore_case = T)) ~ "Candida spp.",
    str_detect(string = MIC_NAME, pattern = regex("Enterobacter(?! cloa)", ignore_case = T)) ~ "Enterobacter spp.",
    str_detect(string = MIC_NAME, pattern = regex("Staphylococcus(?! aur)|stafylokokker, hvite|koagulase-negative|KOAGULASE", ignore_case = T)) ~ "Staphyloccocus spp.", 
    str_detect(string = MIC_NAME, pattern = regex("aerococcus", ignore_case = T)) ~ "Aerococcus spp.",
    str_detect(string = MIC_NAME, pattern = regex("serratia", ignore_case = T)) ~ "Serratia spp.",
    str_detect(string = MIC_NAME, pattern = regex("pseudomonas(?! aeru)", ignore_case = T)) ~ "Pseudomonas spp.",
    str_detect(string = MIC_NAME, pattern = regex("(Clostridium difficile)|(Clostridioides difficile)", ignore_case = T)) ~ "Clostridioides difficile",
    str_detect(string = MIC_NAME, pattern = regex("Stenotrophomonas maltophilia", ignore_case = T)) ~ "Stenotrophomonas maltophilia",
    str_detect(string = MIC_NAME, pattern = regex("Propionibacterium", ignore_case = T)) ~ "Propionibacterium spp.",
    str_detect(string = MIC_NAME, pattern = regex("citrobacter", ignore_case = T)) ~ "Citrobacter spp.",
    str_detect(string = MIC_NAME, pattern = regex("Morganella", ignore_case = T)) ~ "Morganella spp.",
    str_detect(string = MIC_NAME, pattern = regex("syrefaste|mycobacterium|mykobakt", ignore_case = T)) ~ "Mycobacterium spp. (el. syrefaste staver)",
    str_detect(string = MIC_NAME, pattern = regex("bacteroides", ignore_case = T)) ~ "Bacteroides spp.",
    str_detect(string = MIC_NAME, pattern = regex("acinetobacter", ignore_case = T)) ~ "Acinetobacter spp.",
    str_detect(string = MIC_NAME, pattern = regex("Corynebacterium", ignore_case = T)) ~ "Corynebacterium spp.",
    str_detect(string = MIC_NAME, pattern = regex("campylobacter", ignore_case = T)) ~ "Campylobacter spp.",
    str_detect(string = MIC_NAME, pattern = regex("(aspergillus fumigatus)|fumigatus", ignore_case = T)) ~ "Aspergillus fumigatus",
    str_detect(string = MIC_NAME, pattern = regex("Lactobacillus|LAKTOBASILLER", ignore_case = T)) ~ "Lactobacillus spp.",
    str_detect(string = MIC_NAME, pattern = regex("escherichia(?! coli)", ignore_case = T)) ~ "Andre Escherichia spp.",
    str_detect(string = MIC_NAME, pattern = regex("Lactobacillus", ignore_case = T)) ~ "Lactobacillus spp.",
    str_detect(string = MIC_NAME, pattern = regex("salmonella", ignore_case = T)) ~ "Salmonella spp.",
    TRUE ~ "Andre funn"
  )) %>%
  # Convert 'SAMPLETIME' to Date type
  mutate(SAMPLETIME = as.Date(SAMPLETIME)) %>%
  # Rename the 'SAMPLETIME' column to 'SampleDate'
  rename(SampleDate = SAMPLETIME)

# Combine all species from swisslab, unilab, and miclis into a single dataframe
OUS_all_species <- bind_rows(
  list(swisslab_all_species,
       unilab_all_species, 
       miclis_all_species)) %>%
  # Select required columns
  select(PID, SampleDate, Species) %>% 
  # Add a new column 'SampleYear' extracted from 'SampleDate'
  mutate(SampleYear = year(SampleDate)) %>% 
  # Recategorize species based on certain conditions
  mutate(Species = case_when(
    str_detect(string = Species, pattern = regex("Aerococcus spp.|Gram-positive bakterier|Corynebacterium spp.|Lactobacillus spp.|Mycobacterium spp.|Propionibacterium spp.")) ~ "Other Gram-positives",
    str_detect(string = Species, pattern = regex("Gram-negative bakterier|Andre Escherichia spp.|Haemophilus influenzae|Moraxella catarrhalis|Morganella spp.|Providencia spp.|Serratia spp.|Stenotrophomonas maltophilia")) ~ "Other Gram-negatives",
    str_detect(string = Species, pattern = regex("Andre funn|Aspergillus fumigatus")) ~ "Other findings",
    str_detect(string = Species, pattern = regex("Candida albicans|Candida spp.")) ~ "Yeasts",
    str_detect(string = Species, pattern = regex("Citrobacter spp.|Salmonella spp.")) ~ "Other Enterobacteriaceae",
    str_detect(string = Species, pattern = regex("Blandingsvekst")) ~ "Mixed growth",
    str_detect(string = Species, pattern = regex("Enterobacter spp.|Enterobacter cloacae")) ~ "Enterobacter spp.",
    str_detect(string = Species, pattern = regex("Klebsiella pneumoniae|Klebsiella spp.")) ~ "Klebsiella spp.",
    str_detect(string = Species, pattern = regex("Pseudomonas aeruginosa|Pseudomonas spp.")) ~ "Pseudomonas spp.",
    str_detect(string = Species, pattern = regex("Staphyloccocus spp.")) ~ "Staphylococcus spp. (mostly coagulase-negative)",
    str_detect(string = Species, pattern = regex("Streptococcus spp.")) ~ "Viridans and non-haemolytic streptococci",
    str_detect(string = Species, pattern = regex("Streptokokker gruppe A")) ~ "Streptococcus pyogenes",
    str_detect(string = Species, pattern = regex("Streptokokker gruppe B")) ~ "Streptococcus agalactiae",
    str_detect(string = Species, pattern = regex("Streptokokker gruppe C|Streptokokker gruppe G")) ~ "Beta-haemolytic streptococci group C and G",
    TRUE ~ Species
  ))

# Define custom order for species categories
custom_order <- c("Staphylococcus aureus", "Staphylococcus spp. (mostly coagulase-negative)", 
                  "Streptococcus pneumoniae", "Streptococcus pyogenes", "Streptococcus agalactiae",
                  "Beta-haemolytic streptococci group C and G",
                  "Viridans and non-haemolytic streptococci",
                  "Enterococcus faecalis", "Enterococcus faecium", "Other Gram-positives", "Escherichia coli", 
                  "Klebsiella spp.", "Enterobacter spp.", "Proteus spp.", "Other Enterobacteriaceae", 
                  "Pseudomonas spp.", "Acinetobacter spp.", "Other Gram-negatives", "Bacteroides spp.",
                  "Yeasts", "Mixed growth", "Other findings")
OUS_all_species$Species <- factor(OUS_all_species$Species, levels = custom_order)  

# Calculate proportions for each species within each year
OUS_all_species_percent <- OUS_all_species %>%
  # Group by 'SampleYear' and 'Species' and count records in each group
  group_by(SampleYear, Species) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  # Calculate proportion for each species in each year
  group_by(SampleYear) %>%
  mutate(percent = Count / sum(Count)) %>%
  ungroup()

# Define color palette for the plot
color_palette <- brewer.pal(12, "Set3")

# Create a function to generate repeating colors for the plot
color_generator <- colorRampPalette(colors = color_palette)

# Number of unique species categories in the data
n_categories <- length(unique(OUS_all_species_percent$Species))

# Create a stacked bar plot showing distribution of species by sample year
figureS3 <- ggplot(OUS_all_species_percent, aes(x = factor(SampleYear), y = percent, fill = Species)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color_generator(n_categories), guide = guide_legend(title = "Species", ncol = 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Sample Year", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Distribution of Species by Sample Year")

# Save the plot to a PNG file
ggsave("G:/Prosjekter/PDB 3169 - The ResCan study - P_/Artikkel 1/figureS3.png", plot = figureS3, dpi = 320, width = 25, height = 25, units = "cm")
