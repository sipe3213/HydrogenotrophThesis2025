library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(readxl)

#This protocol is to sort only the lowest E-value hydDB classification (target_name) per annotated sequence (query_name)


#Using import function in top right corner is recommended, to make sure E-values are numeric.
hmmscan_output <- read_excel("C:/Users/sipe3/Downloads/hmmscan_output.xlsx", 
                         sheet = "PC2DesulfHydDBgroups_sum", col_types = c("text", 
                                                                          "text", "text", "text", "numeric", 
                                                                          "text", "text", "text", "text", "text", 
                                                                          "text", "text", "text", "text", "text", 
                                                                          "text", "text", "text", "numeric"))
View(hmmscan_output)

#Assign data to dataframe
df <- data.frame(hmmscan_output)

#group by sequence name (query_name) and find lowest E-value for each unique sequence name
df1 <- group_by(df, query_name)
df2 <- summarise(df1, `E.value_full` = min(`E.value_full`))

#join together
good <- semi_join(df1, df2, by = c("query_name","E.value_full"))

#assign them together with hyddb classifications (target_name)
dfdone <- data.frame(good$query_name, good$target_name, good$`E.value_full`)

#Convert dataframe to csv table and downloaded
write.csv(dfdone, "C:/Users/sipe3/Downloads/PC2.csv", row.names = T)

# All csv files are mended together in excel, and continued in Hyddbheatmap




