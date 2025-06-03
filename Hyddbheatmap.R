library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)


#(usually I use the import function 
# in the top right corner, as it is less error prone
# and is easier to assign the E-value column as numeric)
# Import data 
HydrogenasegenesThesis <- read_excel("C:/Users/sipe3/Downloads/HydrogenasegenesThesis.xlsx", sheet = "collected", col_types = c("text", "text", "text", "numeric"))

#View data
View(HydrogenasegenesThesis)

#assign as a dataframe
df <- data.frame(HydrogenasegenesThesis)

#group by sequence name
df1 <- group_by(df, Sequence_name)

#Assign the lowest E-value for each gene
df2 <- summarise(df1, `E.value_full` = min(`E.value_full`))

#join HydDB class name to the lowest E-value of each sequence
good <- semi_join(df1, df2, by = c("Sequence_name","E.value_full"))

#import text file (desired_cultureid.txt) with list of desired order for culture ids 
good$Culture <- factor(good$Culture, levels = c(desired_cultureid$V1))

#Assigning desired order of HydDB classifications, also via text file with prioritized list
desired_order1 <- c(desired_order$V1) 
good$target_name <- factor(good$target_name, levels = desired_order1)

#E-value above 1 is more convenient for ggplot, edited out in final product
"Evalue+one" <- good$E.value_full + 1 

#Plot the heatmap
ggplot(good, aes(x = target_name, y = Culture, fill = `Evalue+one`)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("green", "yellow", "white"),
                       trans = "log",
                       limits = c((1+1.0E-50), (1+1E-4))) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_line(color = "darkred", size = 0.25, linetype = "dotted"),
    panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = "solid")
  ) +
  theme(
    axis.text.x = element_text(color = "darkred", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black")
  ) +
  labs(title = "HydDB HMM scan", x = "Hydrogenase Class", y = "Strain")


