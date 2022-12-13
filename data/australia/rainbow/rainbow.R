install.packages("rainbow")
library("rainbow")

Australiafertility$x
Australiafertility$y

write.csv(Australiafertility$y, file = "data/rainbow/Australia_age_specific_fertility.csv")

