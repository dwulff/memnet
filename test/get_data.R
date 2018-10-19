#
# options(stringsAsFactors = F)
#
# a = "/Users/dwulff/Dropbox (2.0)/Work/Projects/Memory/-- AgingLexicon/1 Raw Data/0_StanMidus/1 Data/0 Stanford/2_Corrected2.0/readydata_wave1.txt"
# a = read.table(a)
# names(a)[1] = 'PID'
#
# b = "/Users/dwulff/Dropbox (2.0)/Work/Projects/Memory/-- AgingLexicon/1 Raw Data/0_StanMidus/1 Data/0 Stanford/2_Corrected2.0/participant_info.txt"
# b = read.table(b, header = T)
#
# require(dplyr)
# a = a %>% inner_join(b[,c('PID','age')])
#
# anis = split(a$V3, a$PID)
# ages = split(a$age, a$PID)
#
# names(anis) = sapply(ages, function(x) x[1])
#
# animal_fluency = anis
#
# devtools::use_data(animal_fluency, overwrite = T)
