########################################################################
####### FIGURES WITHIN THE MANUSCRIPT ##################################
########################################################################
# careful! some plots can take some time to be generated. 

# Figure 1 -----
Figure1 <- ggplot() +
  geom_ellipse(aes(x0 = 1.95, y0 = 2.1, a = 1, b = 0.75, angle = 0), fill = "grey70", color = "grey70", alpha=1) +  
  geom_ellipse(aes(x0 = 1.45, y0 = 1.45, a = 1.7, b = 0.6, angle = 340), fill = "grey80", color = "grey80", alpha=.9) +   
  annotate("text", x = 1, y = 1, label = "No disturbance", size = 6, angle = 30) +
  annotate("text", x = 1.9, y = 2.7, label = "Disturbance", size = 6) +
  geom_polygon(aes(x = c(0, 1.5, 0), y = c(3, 3, 1.5)), fill = "lightblue", alpha = 0.5) +  # Top-left
  geom_polygon(aes(x = c(1.5, 3,3), y = c(0, 0, 1.5)), fill = "lightblue", alpha = 0.5) +  # Bottom-right
  annotate("text", x = 0.45, y = 2.9, label = "Selected against", size = 4, color = "black") +
  annotate("text", x = 2.4, y = 0.3, label = "Biophysical and \nevolutionary constraints", size = 4, color = "black") + 
  annotate("text", x = 0, y = 0.3, label = "Low", size = 4.5, color = "black", angle = 90) +
  annotate("text", x = 0.1, y = 0, label = "Slow", size = 4.5, color = "black", angle = 0) +
  annotate("text", x = 2.8, y = 3, label = "Fast", size = 4.5, color = "black", angle = 0) +
  annotate("text", x = 2.9, y = 2.7, label = "High", size = 4.5, color = "black", angle = -90) +
  scale_x_continuous(name = "Species maximum growth rate", limits = c(0, 3), breaks = seq(0, 5, 1)) +
  scale_y_continuous(name = "Mortality probability at zero growth", limits = c(0, 3), breaks = seq(0, 5, 1)) +
  theme_light() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())


# Figure 2 -----
load("rdata/filtered_fit.stan.rdata")
load("rdata/post.mod.all_cov_matrix.rdata")
data_std <- data2
data_std$dbh <- as.numeric(scale(data_std$dbh))
data_std$gr <- as.numeric(scale(data_std$gr))

## Predictions based on 0-0.5 cm range
# Build a grid of covariate value combinations for which we want posterior predictions:
data = data2
dbh_val <- (log(12.7) - mean(data$dbh)) / sd(data$dbh) # corresponds to dbh of 12.7 cm
gr_val <- c((0^0.47 - mean(data$gr)) / sd(data$gr),   # corresponds to gr of 0 cm yr-1
            (0.5^0.47 - mean(data$gr)) / sd(data$gr)) # corresponds to gr of 0.5 cm yr-1

grid <- modelr::data_grid(data = data_std,
                          gr = seq(gr_val[1], gr_val[2], length.out = 100), 
                          dbh = dbh_val)

# Predictions from the fitted model only:
predict <- post %>% 
  # add grid of covariate values to the posterior draws:
  crossing(grid) %>% 
  # Prediction across species (using hyperparameters only):
  mutate(pred = a + b_dbh * dbh + b_gr * gr) %>% 
  # IF we work from the MORTALITY model (otherwise silence the line below):
  rename(pred_m = pred) %>% 
  mutate(pred_m_invl = inv.logit(pred_m)) %>% 
  # Prediction per species:
  mutate(pred_sp = spcoef_a + spcoef_b_dbh * dbh + spcoef_b_gr * gr) %>% 
  rename(pred_sp_m = pred_sp) %>% 
  mutate(pred_sp_m_invl = inv.logit(pred_sp_m)) %>% 
  # back-transform covariate, if they were standardised prior to running the model:
  mutate(dbh_raw = (dbh * sd(data$dbh)) + mean(data$dbh),
         dbh_raw = exp(dbh_raw),
         gr_raw = (gr * sd(data$gr)) + mean(data$gr),
         gr_raw = gr_raw^(1/0.47), 
         Group = as.character(NA))

gimnosperm <- c("Tsuga canadensis", "Abies balsamea","Larix laricina","Picea abies","Picea glauca",
                "Picea mariana", "Picea rubens", "Pinus banksiana", "Pinus echinata",
                "Pinus ponderosa","Pinus resinosa","Pinus rigida","Pinus strobus","Juniperus virginiana",
                "Pinus taeda", "Pinus virginiana", "Taxodium distichum","Thuja occidentalis")
predict$Group[predict$sp %in% gimnosperm]= "Gymnosperm"
predict$Group[is.na(predict$Group)] = "Angiosperm"
predict <- predict %>%
  mutate (group = as.factor(Group)) 


# note that gr_raw should still be back-transformed with dbh_raw = dbh_raw^(1/0.47)
# to be on the raw scale.

# Summarise and plot central estimate and interval(s):
a = predict %>% 
  group_by(dbh_raw, gr_raw) %>% 
  point_interval(pred_m_invl, 
                 .point = mean, .interval = qi, .width = 0.9) 
# Figure:
b= ggplot(a,aes(gr_raw, pred_m_invl)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5, fill = "#4393C3") +
  geom_line() +
  labs(x = "Previous growth rate\n(0 to 0.5 cm yr-¹)",
       y = "Annual mortality \nprobability")+
  mytheme+
  theme_light()+
  theme(axis.text=element_text(size=8),
    axis.title=element_text(size=8),
    text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

predict_summ <- predict %>% 
  group_by(sp, dbh_raw, gr_raw, Group) %>% 
  point_interval(pred_sp_m_invl, 
                 .point = mean, .interval = qi, .width = 0.9) %>% 
  ungroup

Figure2 <- predict_summ %>% 
  ggplot(aes(x=gr_raw, y=pred_sp_m_invl, group = sp)) + #predict and real growth
  # Silence the geom_ribbon to hide the sp-level uncertainty:
  # geom_ribbon(aes(ymin = .lower, ymax = .upper,
  #                 fill = sp), alpha = 0.5) +
  geom_line(aes(colour = Group), alpha = 0.5) +
  scale_colour_manual(values = c("#009E73", "#D55E00"))+
  labs(x = "Previous growth rate\n(0 to 0.5 cm yr-¹)",
       y = "Annual mortality probability")+
  mytheme +
  theme_light()+
  coord_cartesian(ylim = c(0, 0.4)) + # to focus on specific range of the y-axis
  theme(#legend.position = "bottom",
    axis.text=element_text(size=14),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    legend.position = c(0.2,0.85),
    legend.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())+
  annotation_custom(ggplotGrob(b), xmin = 0.25, xmax = 0.5, ymin = 0.15, ymax = 0.4)

# Figure 3 
load("rdata/matrix.trait.preds.RData")
load("rdata/sp_growth_boot.rdata")

df.analy <- df.analy %>% 
  mutate(mort.prob = pred_sp_m_invl) %>% 
  bind_cols(sp.gr.all.boot) %>% 
  rename(gr.all.low = "2.5%",
         gr.all.high = "97.5%")

# overall figure
figure <- ggplot()+
  geom_linerange(df.analy, mapping= aes(x= gr.95th, y= mort.prob,ymin = .lower, ymax = .upper), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_linerange(df.analy, mapping= aes(x= gr.95th, y= mort.prob,xmin = gr.all.low, xmax = gr.all.high), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_point(df.analy, mapping = aes(x= gr.95th, y= mort.prob), size=2) +
  stat_ma_line(df.analy, mapping = aes(x= gr.95th, y= mort.prob),method = "SMA", se=T, colour="black")+
  ylim(0,0.5)+
  annotate("text", x = 1.6, y = 0.5, label = "SMA slope = 0.35", hjust = 1)+
  annotate("text", x = 1.6, y = 0.48, label = "R² = 0.07",hjust = 1)+
  annotate("text", x = 1.6, y = 0.46, label = "p-value = 0.02",hjust = 1)+
  my_theme

# add Angiosperm data
angio <- df.analy %>% filter(group=="Angiosperm") 

f2 <- figure + 
  geom_point(angio, mapping = aes(x= gr.95th, y= pred_sp_m_invl, color="#009E73"), size=2) +
  theme(legend.position = "bottom",
        text = element_text(size=10))+
  labs (y= "Annual mortality probability \nat zero growth",
        x="Maximum growth (cm yr-¹)")+
  ylim(0,0.5)+ 
  scale_colour_manual(values = c("#009E73"))

# add Gymnosperm data
gimno <- df.analy %>% filter(group=="Gymnosperm")
f3 <- f2 +
  geom_point(gimno,mapping = aes(x= gr.95th, y= mort.prob, color="#D55E00"),size=2)+
  my_theme +
  scale_colour_manual(values = c("#009E73","#D55E00"))+
  theme(text = element_text(size=15))
f3

# separate Angiosperm plot
angio.1 <- ggplot() +
  geom_linerange(angio, mapping= aes(x= gr.95th, y= mort.prob,ymin = .lower, ymax = .upper), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_linerange(angio, mapping= aes(x= gr.95th, y= mort.prob,xmin = gr.all.low, xmax = gr.all.high), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_point(angio, mapping = aes(x= gr.95th, y= pred_sp_m_invl, color="#009E73"), size=2) +
  my_theme+
  scale_colour_manual(values = c( "#009E73"))+
  labs(x="Maximum growth (cm yr-¹)", y="Annual mortality \nprobability at zero growth",title = "Angiosperms")+
  ylim(0,0.5)+
  annotate("text", x = 1.6, y = 0.5, label = "SMA slope = 0.38", hjust = 1)+
  annotate("text", x = 1.6, y = 0.46, label = "R² = 0.06",hjust = 1)+
  annotate("text", x = 1.6, y = 0.42, label = "p-value = 0.07",hjust = 1)+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(0.2,1.6)

# separate Gymnosperm plot
gimno.1 <- ggplot()+
  geom_linerange(gimno, mapping= aes(x= gr.95th, y= mort.prob,ymin = .lower, ymax = .upper), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_linerange(gimno, mapping= aes(x= gr.95th, y= mort.prob,xmin = gr.all.low, xmax = gr.all.high), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_point(gimno,mapping = aes(x= gr.95th, y= pred_sp_m_invl, color="#D55E00"),size=2)+
  my_theme+
  #ylim(0,0.32)+
  scale_colour_manual(values = c( "#D55E00"))+
  labs(x="Maximum growth (cm yr-¹)", y="Annual mortality \nprobability at zero growth",title = "Gymnosperms")+
  ylim(0,0.5)+
  annotate("text", x = 1.6, y = 0.5, label = "SMA slope = 0.24", hjust = 1)+
  annotate("text", x = 1.6, y = 0.46, label = "R² = 0.09",hjust = 1)+
  annotate("text", x = 1.6, y = 0.42, label = "p-value = 0.22",hjust = 1)+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(0.2,1.6)

# combine plots
plots <- (f3|(angio.1/gimno.1)) +plot_annotation(tag_levels = c("A"))+ plot_layout(widths = c(2, 1))


# Figure 4 ----
load("rdata/development.rates_4bins.rdata")
load("rdata/mort.paired_4bins.rdata")
load("rdata/gr.paired_4bins.rData")

temp_early_rates <- development.rates %>% filter(rate=="early")
temp_late_rates <- development.rates %>% filter(rate=="late")

p1 <-  ggplot()+
  geom_linerange(temp_early_rates, mapping= aes(x= gr.95th, y= pred_sp_m_invl,ymin = .lower, ymax = .upper), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_linerange(temp_early_rates, mapping= aes(x= gr.95th, y= pred_sp_m_invl,xmin = gr.lower, xmax = gr.upper), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_point(temp_early_rates, mapping = aes(x= gr.95th, y= pred_sp_m_invl, colour=group), size=1.5) +
  my_theme+
  scale_colour_manual(values = c("#009E73", "#D55E00"))+
  labs(x="", y="Annual mortality \nprobability at zero growth", title="Early development")+
  ylim(0,0.5)+xlim(0.3,1.6)+
  annotate("text", x = 1.6, y = 0.5, label = "SMA slope = 0.33", hjust = 1)+
  annotate("text", x = 1.6, y = 0.46, label = "R² = 0.02",hjust = 1)+
  annotate("text", x = 1.6, y = 0.42, label = "p-value = 0.49",hjust = 1)+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p2 <- ggplot()+
  geom_linerange(temp_late_rates, mapping= aes(x= gr.95th, y= pred_sp_m_invl,ymin = .lower, ymax = .upper), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_linerange(temp_late_rates, mapping= aes(x= gr.95th, y= pred_sp_m_invl,xmin = gr.lower, xmax = gr.upper), linewidth = 0.5, colour = "grey80") + # colour = sign_gr
  geom_point(temp_late_rates, mapping = aes(x= gr.95th, y= pred_sp_m_invl, colour=group), size=1.5) + 
  stat_ma_line(temp_late_rates, mapping = aes(x= gr.95th, y= pred_sp_m_invl),method = "SMA", se=TRUE, colour="black")+
  my_theme+  ylim(0,0.8)+xlim(0.2,1.3)+
  stat_ma_line(method = "SMA", se=TRUE, colour="black")+
  scale_colour_manual(values = c("#009E73", "#D55E00"))+
  labs(x="Maximum growth (cm yr-¹)", y="Annual mortality \nprobability at zero growth", title = "Late development")+
  annotate("text", x = 1.3, y = 0.8, label = "SMA slope = 0.53", hjust = 1)+
  annotate("text", x = 1.3, y = 0.74, label = "R² = 0.17",hjust = 1)+
  annotate("text", x = 1.3, y = 0.68, label = "p-value = 0.01",hjust = 1)+
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p3 <- ggplot(mort.test, aes(x = variable, y = value, group=sp.late)) + 
  geom_line(aes(color = group)) + 
  geom_point(size = 3, aes(color = group)) + 
  my_theme+ 
  labs(x="", y="Annual mortality probability \nat zero growth")+
  scale_colour_manual(values = c("#009E73", "#D55E00"))#+

p4 <- ggplot(gr.test, aes(x = variable, y = value, group=sp.early)) + 
  geom_line(aes(color = group)) + 
  geom_point(size = 3, aes(color = group)) + 
  my_theme+ 
  labs(x="Stand development", y="Maximum growth (cm yr-¹)")+
  annotate("text", x = 1, y = 1.45, label = "*",
           colour = "black", size=6, hjust = 0)+
  scale_colour_manual(values = c("#009E73", "#D55E00"))

plots_stage <- ((p1|p3)/(p2|p4)) +plot_annotation(tag_levels = c("A"))



