###### **********************
## GROUP 6 - TRICHINELLA SP. T8
##
## Group 6 - Ruhi Sajabi, Satchit Bhattacharya, Saniya Chadha, Phuc Pham, and Grace Enike.
##
## 2026-11-02
##
##***************************

## _ Packages used -------
rm(list = ls())
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(stringr)
library(viridis)

#trich_main is the trichinella sp. we will be comparing our assigned trich. group to
#trich_t8 is the trichinella sp. our group was assigned

setwd("./R/")

trich_main = read_csv("../data/GCA_001447595.1.csv")
trich_main |> view()

trich_t8 = read_csv("../data/GCA_001447745.1.csv")
trich_t8

## Identifying number of unique classes in the trichinella t8 sp.

trich_t8$matching_class |> unique() |> length()

## Identifying unique classes in the trichinella spiralis

trich_main$matching_class |> unique() |> length()

## Unique classes by matching class

trich_t8 |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n)) 

trich_main |> 
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n))

## Identify the retrotransposons in the trich_t8 class. From Copilot, everything starting with LTR

trich_t8_retro = trich_t8 |> 
  filter(str_detect(matching_class, "^(LTR|LINE|SINE|PLE)"))

str_detect(trich_t8$matching_class, "^(LTR|LINE|SINE|PLE)")

## Identify the retrotransposons in the trich_main class. From Copilot, everything starting with LTR

trich_main_retro = trich_main |> 
  filter(str_detect(matching_class, "^(LTR|PLE|SINE|LINE)"))

str_detect(trich_main$matching_class, "^(LTR|PLE|SINE|LINE)")

unique(trich_t8_retro$matching_class)


## Percentage deletions for trich_t8_retro

trich_t8_retro |> 
  ggplot(aes(x = matching_class, y = perc_del)) +
  geom_boxplot() +
  labs(title = "Trich T8 Retrotransposons", y ="Percent Deletions (%)" , x = "Retrotransposons class")

## Percentage deletions for trich_main_retro

trich_main_retro |> 
  ggplot(aes(x = matching_class, y = perc_del)) +
  geom_boxplot() +
  labs(title = "T. spiralis Retrotransposons", y ="Percent Deletions (%)" , x = "Retrotransposons class")

## Percentage insertions for trich_t8_retro

trich_t8_retro |> 
  ggplot(aes(x = matching_class, y = perc_insert)) +
  geom_boxplot() +
  labs(title = "Trich T8 Retrotransposons", y ="Percent Insertions (%)" , x = "Retrotransposons class")

## Percentage insertions for trich_main_retro

trich_main_retro |> 
  ggplot(aes(x = matching_class, y = perc_insert)) +
  geom_boxplot() +
  labs(title = "T. spiarlis Retrotransposons", y ="Percent Insertions (%)" , x = "Retrotransposons class")

## Percentage deletion and insertion in trich_te_retro

trich_t8_retro |> 
  ggplot(aes(x = perc_insert, y = perc_del))  +
  geom_point(aes(color = matching_class)) +
  geom_smooth()

## 2nd Attempt at easier to interpret visual via gpplot

trich_t8_retro |> 
  ggplot(aes(x = perc_insert, y = perc_del))  +
  geom_boxplot(aes(color = matching_class)) +
  geom_smooth()

## Percentage deletion and insertion in trich_main

trich_main_retro |> 
  ggplot(aes(x = perc_insert, y = perc_del))  +
  geom_point(aes(color = matching_class)) +
  geom_smooth()

## Alternative visuals 
trich_main_retro |> 
  ggplot(aes(x = perc_insert, y = perc_del))  +
  geom_boxplot(aes(color = matching_class)) +
  geom_smooth()


## Comparing retrotransposons insertion and deletion patterns between both species

retro_combined = bind_rows(trich_main_retro |> mutate(species = "T. spiralis"),
                           trich_t8_retro |> mutate(species = "T8"))

retro_combined |>
  ggplot(aes(x = perc_del, y = perc_insert, color = species)) +
  geom_point() +
  geom_smooth() +
  labs(title = "Insertion and Deletion Patterns in Retrotransposons", x = "Percentage deletions (%)", y = "Percentage insertions (%)", color = "Organism")

summary(trich_t8_retro$perc_insert)
summary(trich_t8_retro$perc_del)

############################################### Show individual classes
retro_combined |>
  ggplot(aes(x = matching_class,
             y = perc_insert,
             fill = species)) + #fill species colours the box plot
  labs(
    title = "Insertion Percentages Across Retrotransposon Classes",
    x = "Retrotransposon Class",
    y = "Percentage Insertions (%)"
  ) +
  geom_boxplot() +
  facet_wrap(~species) + #make separate panels for each species instead of one graph 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# angle 45 = rotate the x axis to fit 
#h_just = aligns the labels 
ggsave("../output/figure1_insertions.png",
       width = 9, height = 9, dpi = 600)

###
retro_combined |>
  ggplot(aes(x = matching_class,
             y = perc_del,
             fill = species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  labs(
    title = "Deletion Percentages Across Retrotransposon Classes",
    x = "Retrotransposon Class",
    y = "Percentage Deletions (%)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../output/figure2_deletions.png",
       width = 9, height = 9, dpi = 600)

## Grace's code template stops here ----
## Novelty - Retrotransposon Composition - Stacked Bar Chart ----
retro_combined |>
  count(species, matching_class) |>
  group_by(species) |> #Normalize counts into percentages for each genome
  mutate(percent = n / sum(n) * 100) |>
  ggplot(aes(x = species,
             y = percent,
             fill = matching_class)) +
  geom_col() +
  labs(
    title = "Retrotransposon Composition",
    x = "Genome",
    y = "Percentage of Annotations",
    fill = "Retrotransposon Class"
  ) +
  theme_minimal()
ggsave("../output/figure3_retrotransposon_composition.png",
       width = 9, height = 9, dpi = 600)

## Novelty - Retrotransposon Composition - Percentage Differences ----
retro_combined |>
  count(species, matching_class) |>
  group_by(species) |>
  mutate(percent = n / sum(n) * 100) |>
  select(species, matching_class, percent) |> #Place genome percentages side by side
  pivot_wider(
    names_from = species,
    values_from = percent,
    values_fill = 0
  ) |>
  mutate(difference = T8 - `T. spiralis`) |> #Calculate percentage-point differences (T8 minus T. spiralis)
  arrange(desc(abs(difference)))

## Novelty - Statistical Analysis - Chi-square Test ----
class_counts <- retro_combined |>
  count(matching_class, species) |>
  pivot_wider(
    names_from = species,
    values_from = n,
    values_fill = 0
  )# Count annotations in each class for both genomes

count_matrix <- class_counts |>
  select(`T. spiralis`, T8) |>
  as.matrix()
# Creates a table for the statistical test

chi_result <- chisq.test(count_matrix)
chi_result
# Run the chi-square test and view results




########## Minaal Edits to code ----

### Edit 1 ----
# I used two-way ANOVA test to check if there was as significance between species in insertion/deletion percentages, retrotransposon classes, or depending on the interaction between species and retrotransposon class. Results show no significant difference between species for  insertion patterns, deletion patterns, but significant p value between retrotransposon classes. I used a ANOVA test to see another statistical tests result since your chi-square results p value was 0.340, which is not significant. 


# Statistical Analysis - Two-way ANOVA for Deletions

anova_del <- aov(perc_del ~ species * matching_class,
                 data = retro_combined)

summary(anova_del)

par(mfrow = c(2, 2))
plot(anova_del)

## Statistical Analysis - Two-way ANOVA for Insertions 

anova_insert <- aov(perc_insert ~ species * matching_class,
                    data = retro_combined)

summary(anova_insert)

par(mfrow = c(2, 2))
plot(anova_insert)

### Edit 2 ----
# I expanded the dataset by adding my group's dataset from Trichinella Nelsoni. Using similar coding of your groups novelty approach for Percentage Differences, I compared the TE composition between T8, T.spiralis, and T.Nelsoni. I also added a ggplot for visualization. Similar to your groups finding, LTR/GYPSY was the highest in T8 compared to T.spiralis, and T.Nelsoni.

### Comparing Retrotransposon groups between T.spiralis, T8, T.Nelsoni

trich_nelsoni = read_csv("../data/GCA_001447455.1 (1).csv")
trich_nelsoni


## Unique classes by matching class

trich_nelsoni |>
  group_by(matching_class) |> 
  count() |> 
  arrange(desc(n))

## Identify the retrotransposons in the trich_nelsoni class. From Copilot, everything starting with LTR

trich_nelsoni_retro = trich_nelsoni |> 
  filter(str_detect(matching_class, "^(LTR|PLE|SINE|LINE)"))

str_detect(trich_nelsoni$matching_class, "^(LTR|PLE|SINE|LINE)")


##### binding between, T8, main, and nelsoni 
retro_three = bind_rows(trich_main_retro |> mutate(species = "T. spiralis"),
                        trich_t8_retro |> mutate(species = "T8"),
                        trich_nelsoni_retro |> mutate(species = "T. nelsoni"))



##  Retrotransposon Composition - Percentage Differences Three Species 

retro_three |>
  count(species, matching_class) |>
  group_by(species) |>
  mutate(percent = n / sum(n) * 100) |>
  select(species, matching_class, percent) |>
  pivot_wider(
    names_from = species,
    values_from = percent,
    values_fill = 0
  ) |>
  mutate(
    difference_T8_vs_spiralis = T8 - `T. spiralis`,
    difference_nelsoni_vs_spiralis = `T. nelsoni` - `T. spiralis`,
    difference_nelsoni_vs_T8 = `T. nelsoni` - T8
  ) |> view()


## ggplot Comparing Retrotransposon Counts Across Three Species 

retro_three |>
  count(species, matching_class) |>
  ggplot(aes(x = matching_class, y = n, fill = species)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Retrotransposon Counts Across Three Trichinella Species",
    x = "Retrotransposon Class",
    y = "Number of Retrotransposons",
    fill = "Species"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))

ggsave("../output/2026-11-22_Fig_retro_three_species_comparison.png",
       width = 8, height = 6, dpi = 600)


### Edit 3 ----
# I shortened your code for creating the boxplots of insertion and deletion percentages. instead of two seperate boxplots I combined them into one plot with the four boxplots.  

###Combine Insertion and Deletion Boxplots

retro_combined |>
  pivot_longer(
    cols = c(perc_insert, perc_del),
    names_to = "measurement",
    values_to = "percentage"
  ) |>
  ggplot(aes(x = matching_class,
             y = percentage,
             fill = species)) +
  geom_boxplot() +
  facet_grid(measurement ~ species, scales = "free_y") +
  labs(
    title = "Insertion and Deletion Percentages Across Retrotransposon Classes",
    x = "Retrotransposon Class",
    y = "Percentage (%)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../output/Boxplot_ins_del_combined.png",
       width = 9, height = 9, dpi = 600)
