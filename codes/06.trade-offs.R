########################################################################
####### TRADE-OFFS ESTIMATES  ##########################################
########################################################################

stand <- readRDS("rdata/stand_dev_us.RDS") # stand development charaterisation

# load dataset for all species ----
load("rdata/temperate.all.rates.RData")
load("rdata/sp_growth_boot.RData")
df.analy <- temperate_rates %>% 
  mutate(mort.prob = pred_sp_m_invl,# rename variables but keep the original name
         q95th = gr.95th)

sp.gr.all.boot = as_tibble (sp.gr.all.boot) %>% 
  mutate(sp = rownames(sp.gr.all.boot),
         gr.upper = `97.5%`,
         gr.lower = `2.5%`)
df.analy <- merge(df.analy,sp.gr.all.boot, by.x = "sp", by.y = "sp")
angio <- filter (df.analy, group=="Angiosperm")
gimno <- filter (df.analy, group=="Gymnosperm")

# SMA regressions for all species ----
model <- sma(mort.prob ~ gr.95th, data = df.analy)
summary(model)
model.a <- sma(mort.prob ~ gr.95th, data = angio)
summary(model.a)
model.g <- sma(mort.prob ~ gr.95th, data = gimno)
summary(model.g)

# early development trade-offs ----
load("rdata/data_early_succ_25.RData")
load("rdata/mortality.probability.at.0gr_EARLY_cov_25.RData")
load("rdata/live.RData")

early <- filter (stand, for.dev <= "0.25")
early.data = live %>% filter(plot.id %in% early$tmt.plot.id)
gr_early <- early.data %>% 
  mutate(time.interv = date3-date2,
         growth.dbh = (dbh3-dbh2)/time.interv) %>% 
  filter(growth.dbh >= 0) %>% 
  group_by(sp) %>% 
  summarise(mean.gr = mean(growth.dbh), 
            gr.95th = quantile(growth.dbh, probs= 0.95),
            group = NA) %>% 
  mutate(std.q.95th = gr.95th/max(gr.95th))

temp_early_rates <- merge(pred.early, gr_early, by.x = "sp", by.y = "sp")

gr_early <- early.data %>% 
  mutate(time.interv = date3-date2,
         growth.dbh = (dbh3-dbh2)/time.interv) %>% 
  filter(growth.dbh >= 0) 
sp.gr.early <- split(gr_early, gr_early$sp)

# this might be lazy, use the RData below
# sp.gr.early <- lapply(sp.gr.early, obtain_sd_growth95, nboot = 999)
# sp.gr.early.boot <- base::do.call(rbind, sp.gr.early)
# sp.gr.early.boot = as_tibble (sp.gr.early.boot) %>% 
#   mutate(sp = rownames(sp.gr.early.boot),
#          gr.upper = `97.5%`,
#          gr.lower = `2.5%`)

load("rdata/sp.gr.early.boot.rdata")
temp_early_rates <- merge(temp_early_rates, sp.gr.early.boot, by.x = "sp", by.y = "sp")
gimnosperm <- c("Tsuga canadensis", "Abies balsamea","Larix laricina","Picea abies","Picea glauca",
                "Picea mariana", "Picea rubens", "Pinus banksiana", "Pinus echinata",
                "Pinus ponderosa","Pinus resinosa","Pinus rigida","Pinus strobus","Juniperus virginiana",
                "Pinus taeda", "Pinus virginiana", "Taxodium distichum","Thuja occidentalis")
temp_early_rates$group[temp_early_rates$sp %in% gimnosperm]= "Gymnosperm"
temp_early_rates$group[is.na(temp_early_rates$group)] = "Angiosperm"
temp_early_rates <- temp_early_rates %>%  mutate (group = as.factor(group)) 

# SMA regressions for early development ----
model_early <- sma(pred_sp_m_invl ~ gr.95th, data = temp_early_rates)
summary(model_early)

# late development trade-offs ----
load("rdata/data_late_succ_75.RData")
load("rdata/mortality.probability.at.0gr_LATE_cov_75.RData")

late <- filter (stand, for.dev >= "0.75")
late.data = live %>% filter(plot.id %in% late$tmt.plot.id)
gr_late <- late.data %>% 
  mutate(time.interv = date3-date2,
         growth.dbh = (dbh3-dbh2)/time.interv) %>% 
  filter(growth.dbh >= 0) %>% 
  group_by(sp) %>% 
  summarise(mean.gr = mean(growth.dbh), 
            gr.95th = quantile(growth.dbh, probs= 0.95),
            group = NA) %>% 
  mutate(std.q.95th = gr.95th/max(gr.95th))

temp_late_rates <- merge(pred.late, gr_late, by.x = "sp", by.y = "sp")
gr_late <- late.data %>% 
  mutate(time.interv = date3-date2,
         growth.dbh = (dbh3-dbh2)/time.interv) %>% 
  filter(growth.dbh >= 0) 
sp.gr.late <- split(gr_late, gr_late$sp)

# this might be lazy, use the RData below
# sp.gr.late <- lapply(sp.gr.late, obtain_sd_growth95, nboot = 999)
# sp.gr.late.boot <- base::do.call(rbind, sp.gr.late)
# sp.gr.late.boot = as_tibble (sp.gr.late.boot) %>% 
#   mutate(sp = rownames(sp.gr.late.boot),
#          gr.upper = `97.5%`,
#          gr.lower = `2.5%`)

load("rdata/sp.gr.late.boot.rdata")
temp_late_rates <- merge(temp_late_rates, sp.gr.late.boot, by.x = "sp", by.y = "sp")

gimnosperm <- c("Tsuga canadensis", "Abies balsamea","Larix laricina","Picea abies","Picea glauca",
                "Picea mariana", "Picea rubens", "Pinus banksiana", "Pinus echinata",
                "Pinus ponderosa","Pinus resinosa","Pinus rigida","Pinus strobus","Juniperus virginiana",
                "Pinus taeda", "Pinus virginiana", "Taxodium distichum","Thuja occidentalis")
temp_late_rates$group[temp_late_rates$sp %in% gimnosperm]= "Gymnosperm"
temp_late_rates$group[is.na(temp_late_rates$group)] = "Angiosperm"
temp_late_rates <- temp_late_rates %>%
  mutate (group = as.factor(group)) 

# SMA regressions for late development ----
model_late <- sma(pred_sp_m_invl ~ gr.95th, data = temp_late_rates)
summary(model_late)

# commom species ----
commom_sp <- intersect(unique(data_early$sp), unique(data_late$sp)) # select the commom species 
late2 = temp_late_rates %>% filter(sp %in% commom_sp)
early2 = temp_early_rates %>% filter(sp %in% commom_sp)
test_stage <- data.frame(mort.early = early2$pred_sp_m_invl,
                         mort.late = late2$pred_sp_m_invl,
                         sp.late = late2$sp,
                         sp.early = early2$sp,
                         gr.early = early2$gr.95th,
                         gr.late = late2$gr.95th, 
                         group = early2$group)
mort.test <- melt(test_stage[,c(1:2,3,7)])
mort.test <- mort.test %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = replace(variable, variable == "mort.early", "Early")) %>%
  mutate(variable = replace(variable, variable == "mort.late", "Late"))

m1 <- t.test(value ~ variable, data=mort.test, paired=TRUE)
m1
gr.test <- melt(test_stage[,4:7])
gr.test <- gr.test %>%
  mutate(variable=as.character(variable)) %>%
  mutate(variable = replace(variable, variable == "gr.early", "Early")) %>%
  mutate(variable = replace(variable, variable == "gr.late", "Late"))

m2 <- t.test(value ~ variable, data=gr.test, paired=TRUE)
m2
