library(tidyverse)

load('censutractlevelcommuting/w_dc/2010 _odMatrix_US.RData')

## Converting census tract level commuting data into county level data and restrict to CONUS

odmatrix_US_tbl <- odmatrix_US %>% as.tibble() %>%
    mutate(w_cntyid = substr(w_tractid, 1, 5),
           h_cntyid = substr(h_tractid, 1, 5)) %>%
    filter( !(substr(w_cntyid,1,2) %in% OutsideCONUS | substr(h_cntyid,1,2) %in% OutsideCONUS)) %>%
    group_by(w_cntyid, h_cntyid) %>%
    summarise(
        commuters = sum(S000, na.rm=T) ### S000 is the total number of workers
    )

#rm(odmatrix_US, odmatrix_tr, year, years)


G2 <- graph.data.frame(odmatrix_US_tbl) %>%
    as.undirected(mode = "collapse", edge.attr.comb = "sum")

C2010 <- as_long_data_frame(G2)[, c('ver[el[, 1], ]', 'ver2[el[, 2], ]', 'commuters')]
names(C2010) <- c('node1', 'node2', 'commuters')
write.csv(C2010, file="C2010_withDC.csv")