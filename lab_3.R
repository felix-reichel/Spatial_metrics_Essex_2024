# 3 practicals ----
## 3.1 use contdird_COW_v3.2.csv, perhaps with the help of contdirs_COW_v3.2.csv, 
#to create a row-standardized contiguity W for all sovereign countries in the world as of 2016.

## 3.2 use ideal_point_UNGA_Mar2021.csv
### 3.2.1 to see ten nearest countries to China and US respectively in terms of IdealPointDistance (that's the variable name),
### 3.2.3 to create a W based on the inverse IdealPointDistance (no row-standardization needed),
### *should time permit* 3.2.2 to visualize IdealPointDistance among G20 countries in a lower-triangular way.

## 3.3 use the dichotomous democracy-autocracy classification data (democracy-v2.0_1.csv) to
###### to create a regime similarity W as of 2010 (ignoring the missing Kosovo) of which
###### w_ij = 1 if i and j have the same regime type (i != j).


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

library(cshapes)
library(countrycode)
library(spdep)
library(tidyverse)
library(vroom)


## 3.1 Create row-standardized contiguity matrix for 2016 ----
contiguity_files <- list.files(path = "data/", pattern = "contdir", full.names = TRUE)

contiguity_data <- map(contiguity_files, read_csv) %>%
  map(~ filter(.x, year == max(.x$year)))

names(contiguity_data) <- c("dyadic", "monadic")

hist(contiguity_data$monadic$total, xlab = "Number of neighboring countries", ylab = "", main = "")

missing_countries_state1 <- anti_join(distinct(contiguity_data$monadic, stateab), 
                                      distinct(contiguity_data$dyadic, state1ab), 
                                      by = c("stateab" = "state1ab"))

missing_countries_state2 <- anti_join(distinct(contiguity_data$monadic, stateab), 
                                      distinct(contiguity_data$dyadic, state2ab), 
                                      by = c("stateab" = "state2ab"))

contiguity_matrix <- select(contiguity_data$dyadic, state1ab, state2ab, conttype) %>%
  add_row(state1ab = c("ICE", "NEW"), state2ab = c("ICE", "NEW")) %>%
  complete(state1ab, state2ab) %>%
  mutate(contig = if_else(is.na(conttype), 0, 1), .keep = "unused") %>%
  pivot_wider(values_from = contig, names_from = state2ab)

contiguity_matrix <- rowwise(contiguity_matrix) %>%
  mutate(rowsum = sum(c_across(where(is.double))), .before = AAB) %>%
  mutate(across(where(is.double) & !c(rowsum), ~ if_else(rowsum != 0, .x / rowsum, 0))) %>%
  select(AAB:ncol(.)) %>%
  as.matrix()

rownames(contiguity_matrix) <- colnames(contiguity_matrix)
rowSums(contiguity_matrix)

# 3.2 National ideal point distance based on UNGA voting ----
# 3.3 Create regime similarity matrix for 2010 ----
regime_data <- read_csv("data/democracy-v2.0_1.csv") %>%
  filter(year == 2010) %>%
  drop_na(democracy) %>%
  select(abbreviation, democracy)

table(regime_data$democracy)

regime_dyadic <- uncount(regime_data, weights = nrow(regime_data)) %>%
  rename(cntry1 = abbreviation) %>%
  group_by(cntry1) %>%
  mutate(cntry2 = regime_data$abbreviation, .before = democracy) %>%
  ungroup()

regime_dyadic <- full_join(regime_dyadic, regime_dyadic, by = c("cntry1" = "cntry2", "cntry2" = "cntry1"))

regime_similarity_matrix <- mutate(regime_dyadic, same_regime = if_else(democracy.x == democracy.y & cntry1 != cntry2, 1, 0)) %>%
  select(cntry1, cntry2, same_regime) %>%
  arrange(cntry1, cntry2) %>%
  pivot_wider(names_from = cntry2, values_from = same_regime)

print(regime_similarity_matrix)


# essex summer school: spatial econometrics
# solution 3
# Muzhou Zhang & Martin Steinwand
# 24 July 2023
rm(list=ls())

library(cshapes)
library(countrycode)
library(spdep)
library(tidyverse)
library(vroom)


# 3.2 national ideal point distance, based on the UNGA voting ----
ipd <- read_csv("data/ideal_point_UNGA_Mar2021.csv")

### 3.2.1 to see ten nearest countries to China and US respectively 
#in terms of IdealPointDistance (that's the variable name),
ipd_10farthest <- filter(ipd, iso3c.x %in% c("CHN", "USA")) %>%
  group_by(iso3c.x) %>% # for each country
  slice_min(order_by = IdealPointDistance, n = 10) %>%
  ungroup(.) %>%
  select(iso3c.x, Countryname.y, IdealPointDistance)

ipd_10farthest <- map(c("CHN", "USA"), ~ filter(ipd_10farthest, iso3c.x == .x))

names(ipd_10farthest) <- c("China", "United States")

ipd_10farthest <- map(ipd_10farthest, select, -iso3c.x)

ipd_10china <- ipd_10farthest$China
ipd_10USA <- ipd_10farthest$'United States'

### 3.2.3 to create a W based on the inverse IdealPointDistance (no row-standardization needed),
W_ipd_tibble <- complete(ipd, Countryname.x, Countryname.y) %>%
  mutate(IdealPointDistance = if_else(Countryname.x == Countryname.y, 0, IdealPointDistance)) %>% 
  select(Countryname.x, Countryname.y, IdealPointDistance) %>%
  pivot_wider(names_from = Countryname.y, values_from = IdealPointDistance)

W_ipd <- mutate(W_ipd_tibble, across(where(is.double), ~ if_else(.x != 0, 1 / .x, .x)), .keep = "used") %>%
  as.matrix()

rownames(W_ipd) <- colnames(W_ipd)


#3.2.2 to visualize IdealPointDistance among G20 countries in a lower-triangular way.
# W for G20 countries
G20 <- c("Argentina", "Australia", "Brazil", "Canada", "China", "France", "Germany", "India", 
         "Indonesia", "Italy", "Japan", "Mexico", "Russia", "Saudi Arabia", "South Africa", 
         "South Korea", "Turkey", "United Kingdom", "United States")

W_ipd_G20 <- filter(W_ipd_tibble, Countryname.x %in% G20) %>%
  select(Countryname.x, all_of(G20))

W_ipd_G20 <- as.matrix(select(W_ipd_G20, Argentina:`United States`))

W_ipd_G20[lower.tri(W_ipd_G20, diag = TRUE)] <- NA

dyadic_ipd_G2O <- as_tibble(W_ipd_G20) %>%
  mutate(cname1 = G20, .before = Argentina) %>%
  pivot_longer(cols = Argentina:`United States`, names_to = "cname2", values_to = "ipd")

ggplot(dyadic_ipd_G2O) + 
  geom_tile(aes(x = cname1, y = reorder(cname2, desc(cname2)), fill = ipd)) +
  scale_fill_gradient(low = "snow2", high = "black", na.value = "white", name = "Ideal Point Distance") +
  labs(x = "", y = "", title = "Revealed National Preference at the UNGA, Dayadic Comparison of G20 Countries", 
       caption = "Source: Bailey, Strezhnev, and Voeten (2017)") +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=45, vjust=.8),
    legend.position = "bottom",
    legend.justification = c(0, 0),
    panel.grid = element_blank(),
    text = element_text(family = "Georgia")
  )

rm(list = ls())
dev.off()
