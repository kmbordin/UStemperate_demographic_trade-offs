########################################################################
####### GROWTH ESTIMATES ACROSS SPECIES ################################
########################################################################

# New dataframe to be used -----
load("rdata/census1to3.data.RData")
df <- data.frame(tag = census1$tmt.tree.id, #treeid
                 sp = census1$species.cor,  #species
                 dbh1 = census1$d, #diameter in mm
                 dbh2 = census2$d, 
                 dbh3 = census3$d,
                 ba1 = census1$ba, #basal area
                 ba2 = census2$ba,
                 ba3 = census3$ba,
                 date1 = census1$census.date, #census date
                 date2 = census2$census.date, 
                 date3 = census3$census.date, 
                 pom1 = census1$pom, #point of measurement
                 pom2 = census2$pom, 
                 pom3 = census3$pom, 
                 DFstatus1 = census1$tree.status, #0= alive; 1= dead
                 DFstatus2 = census2$tree.status, 
                 DFstatus3 = census3$tree.status, 
                 database.code = census1$database.code, # database code
                 plot.id = census1$tmt.plot.id) #plot.id

df <- df %>%
  mutate(pom1 = as.numeric(pom1), 
         pom2 = as.numeric(pom2),
         pom3 = as.numeric(pom3)) %>%
  mutate(pom2 = coalesce(pom2,pom1)) %>% # fill missing info
  mutate(pom3 = coalesce(pom3,pom1)) %>%
  mutate(dbh1 = dbh1/10, # convert dbh to cm if the data are in mm
         dbh2 = dbh2/10,
         dbh3 = dbh3/10) 

range(df$date3-df$date1)  

# Growth calculations -----
# subset data to get cohort to be used for calculating growth in census 1-2 and exclude NA's
live_data <- df$DFstatus1=="0" & df$DFstatus2=="0" & df$dbh1 >= 10 #0 means alive
live <- subset(df, live_data); print(head(live))

# check good data for calculating growth
gooddate <- live$date1>0 & live$date2>0 #only positive census = two censuses
samepom <- live$pom1 == live$pom2 #same pom
hasdbh  <- !is.na(live$dbh1) & !is.na(live$dbh2) #dbh measurement
good2 <- gooddate & samepom & hasdbh

# subset live data by above conditions
live <- subset(live, good2); print(dim(live))

# calculate growth parameters and unite in dataframe 
# note that dbh is in cm and time in the data is in fractional years

time <- (live$date2-live$date1) #time interval
growth.dbh <- (live$dbh2-live$dbh1)/time #growth in dbh

growth <-data.frame(live, time=time, growth.dbh = growth.dbh)

# Growth transformation ----- 
skewness(growth$growth.dbh[growth$growth.dbh >= 0]^0.5)
skewness(growth$growth.dbh[growth$growth.dbh >= 0]^0.47) # BEST ONE
skewness(growth$growth.dbh[growth$growth.dbh >= 0]^0.465)
