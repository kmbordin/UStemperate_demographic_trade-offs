########################################################################
####### GROWTH AND MORTALITY ACROSS SPECIES ############################
########################################################################

# load dataset ----
load("rdata/live.RData")
load("rdata/mortality.probability.at.0gr_cov.RData")

# filtering species 
gr <- live %>%  filter (sp %in% pred.all$sp)

gr <- gr %>% 
  mutate(time.interv = date3-date2,
         growth.dbh = (dbh3-dbh2)/time.interv) %>% 
  filter(growth.dbh >= 0) %>% 
  group_by(sp) %>% 
  summarise(mean.gr = mean(growth.dbh), 
            gr.95th = quantile(growth.dbh, probs= 0.95), 
            group = NA) %>% 
  mutate(std.q.95th = gr.95th/max(gr.95th))

temperate_rates <- merge(pred.all, gr, by.x = "sp", by.y = "sp")
# save(temperate_rates, file = "rdata/temperate.all.rates.RData")

# bootstrap growth rates -----
gr <- live %>% 
  filter (sp %in% pred.all$sp)

names(gr)
gr <- gr %>% 
  mutate(time.interv = date3-date2,
         growth.dbh = (dbh3-dbh2)/time.interv) %>% 
  filter(growth.dbh >= 0) 
sp.gr <- split(gr, gr$sp)

sp.gr.all <- lapply(sp.gr, obtain_sd_growth95, nboot = 999)

#summarise the results per species 
sp.gr.all.boot <- base::do.call(rbind, sp.gr.all)
#save(sp.gr.all.boot, file = "rdata/sp_growth_boot.RData")
