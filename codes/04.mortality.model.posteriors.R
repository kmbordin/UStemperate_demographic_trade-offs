########################################################################
####### EXTRACT MORTALITY PROBABILITIES ################################
########################################################################
# r datas ----
load("rdata/post.mod.all_cov_matrix.rdata")
load("rdata/post.mod.late_cov_matrix_75.rdata")
load("rdata/post.mod.early_cov_matrix_25.rdata")

#load("rdata/filtered_fit.stan.rdata")
data_all <- data2
#load("rdata/data_late_succ_75.RData")
data_late
#load("rdata/data_early_succ_25.RData")
data_early

# extract model predictions -----
pred.all = preds(data = data_all,post = post)
pred.late = preds(data = data_late,post = post_late)
pred.early = preds(data = data_early,post = post_early)

# plot posterior mean and CI per species ----
all = plot_post(pred = pred.all, var = "All species")
early = plot_post(pred = pred.early, var = "Early successional species")
late = plot_post(pred = pred.late, var = "Late successional species")

# save(pred.all, file="rdata/mortality.probability.at.0gr_cov.RData")
# save(pred.early, file="rdata/mortality.probability.at.0gr_EARLY_cov_25.RData")
# save(pred.late, file="rdata/mortality.probability.at.0gr_LATE_cov_75.RData")