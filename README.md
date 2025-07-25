[![DOI](https://zenodo.org/badge/892592324.svg)](https://doi.org/10.5281/zenodo.14402211)

------------------------------------------------------------------------

### This repository contains all files related to the analyses performed for the manuscript "Growth-survival trade-off in temperate trees is weak and restricted to late-successional stages", led by Kauane M. Bordin

Important information: the raw data is openly available at <https://research.fs.usda.gov/products/dataandtools/tools/fia-datamart> For this manuscript's purposes, the processed data was standardised and harmonised in the context of the TreeMort project (<https://more.bham.ac.uk/treemort/>).

For specific requests and questions, please send an email to kauanembordin\@gmail.com

#### *Content of the repository:*

1.  `codes`: contains the R scripts to run the analyses and generate the graphical visualization

    1.  `01.loading.packages`: contains all required packages and functions created for the data analysis

    2.  `02.growth.rates`: create the files to estimate maximum tree growth (i.e., 95th quantile of species growth rates)

    3.  `03.mortality.model`: run the Bayesian model to estimate species mortality probability at zero growth for all species, as well as species from early and late developmental stands

    4.  `04.mortality.model.posteriors`: extract the posterior distribution of the species mortality probability at zero growth for all species and species from early and late developmental stands

    5.  `05.growth.and.mortality`: generates the 95% confidence intervals of maximum growth and combine the estimated mortality probability in a single file

    6.  `06.trade-offs`: test the relationship between growth and survival across the temperate tree species

    7.  `07.figures`: generate the four figures presented in the main manuscript file

        `stan_model_for_mortality_probability.stan`: Bayesian model to estimate the mortality probability at zero growth, using the time census interval *t* as an exponential, and alpha priors as (-1, 2) on logit scale: -3 corresponds to a mort. prob. of 0.04 and 1 to 0.73; see main files for details

2.  `rdata`: due to the nature of the analyses (sometimes very have and lazy to run), this folder contains the files .RData to allow model tests

    1.  `census1to3.RData`: stem-level three census data

    2.  `data_early_succ_25.RData`: stem-level data from stands of early development

    3.  `data_late_succ_75.RData`: stem-level data from stands of late development

    4.  `development.rates_4bins.RData`: all estimated demographic parameters of species in stands from early and late development

    5.  `filtered_fit.stan.RData`: filtered data of all stems to fir stan model

    6.  `fit.stan.total.RData`: filtered data of all stems to fir stan model

    7.  `gr.paired_4bins.RData`: paired data of species to test for differences between growth rates in early and late developmental stands

    8.  `live.RData`: first filtered data - alive and dead information

    9.  `mort.paired_4bins.RData`: paired data of species to test for differences between mortality probability in early and late developmental stands

    10. `mortality.probability.at.0gr_cov.RData`: estimated demographic parameters of species across all sample

    11. `mortality.probability.at.0gr_EARLY_cov_25.RData`: estimated demographic parameters of species across early developmental stands

    12. `mortality.probability.at.0gr_LATE_cov_75.RData`: estimated demographic parameters of species across late developmental stands

    13. `post.mod.all_cov_matrix.RData`: posterior distributions of the mortality probability model across all stands

    14. `post.mod.early_cov_matrix_25.RData`: posterior distributions of the mortality probability model across early stands

    15. `post.mod.late_cov_matrix_75.RData`: posterior distributions of the mortality probability model across late stands

    16. `sp.gr.early.boot.RData`: confidence itervals of 95th percetile of growth rates of species across early developmental stands

    17. `sp.gr.late.boot.RData`: confidence itervals of 95th percetile of growth rates of species across late developmental stands

    18. `sp_growth.boot.RData`: confidence itervals of 95th percetile of growth rates of species across all stands

    19. `stan_mort_mod_output_nc_mod_early_cov_matrix.RData`: stan model output for species across early developmental stands (careful - very heavy)

    20. `stan_mort_mod_output_nc_mod_late_cov_matrix.RData`: stan model output for species across latedevelopmental stands (careful - very heavy)

    21. `stan_mort_mod_output_nc_mod2_cov_matrix.RData`: stan model output for species across all developmental stands (careful - very heavy)

    22. `stand_dev_us.rds`: plot-level information of development status across all stands

    23. `temperate.all.rates.RData`: all estimated demographic parameters of species in all stands

3.  `results`: main results from the analyses shown within the manuscript

    1.  `conceptual.fig.png:` Figure 1

    2.  `sp.level.range.all.0-0.5cm_groups.png`: Figure 2

    3.  `figure 3.png`: Figure 3

    4.  `early_and_late_developm.png`: Figure 4

    5.  `species_summary.rtf`: Appendix 1

    6.  `sp.level.predicts.zerogr_all.png`, `sp.level.predicts.zerogr_EARLY.png`, `sp.level.predicts.zerogr_LATE.png`: Appendix 2 to 4

        1.  `data.total.csv`: data frame with all the results dataset (see below)

            | Variable              | README                                                                                |
            |-----------------------|-------------------------------------------------|
            | group                 | phylogenetic clade (Angiosperm or Gymnosperm)                                         |
            | family                | Family                                                                                |
            | genus                 | Genus                                                                                 |
            | sp                    | Species                                                                               |
            | all.mort.prob         | Mortality probability at zero growth for all species                                  |
            | all.lower.mort.prob   | 2.5% of mortality probability at zero growth for all species                          |
            | all.upper.mort.prob   | 97.5% of mortality probability at zero growth for all species                         |
            | all.95gr              | 95th quantile of growth for all species                                               |
            | all.sd.gr             | standard deviation of 95th quantile of growth                                         |
            | all.lower.gr          | 2.5% of 95th quantile of growth for all species                                       |
            | all.upper.gr          | 97.5% of 95th quantile of growth for all species                                      |
            | late.mort.prob        | Mortality probability at zero growth for species in late development stands           |
            | late.lower.mort.prob  | 2.5% of mortality probability at zero growth for species in late development stands   |
            | late.upper.mort.prob  | 97.5% of mortality probability at zero growth for species in late development stands  |
            | late.95gr             | 95th quantile of growth for species in late development stands                        |
            | late.sd.gr            | Standard deviation of 95th quantile of growth for species in late development stands  |
            | late.upper.gr         | 97.5% of 95th quantile of growth for species in late development stands               |
            | late.lower.gr         | 2.5% of 95th quantile of growth for species in late development stands                |
            | early.mort.prob       | Mortality probability at zero growth for species in early development stands          |
            | early.lower.mort.prob | 2.5% of mortality probability at zero growth for species in early development stands  |
            | early.upper.mort.prob | 97.5% of mortality probability at zero growth for species in early development stands |
            | early.95gr            | 95th quantile of growth for species in early development stands                       |
            | early.sd.gr           | Standard deviation of 95th quantile of growth for species in early development stands |
            | early.upper.gr        | 97.5% of 95th quantile of growth for species in early development stands              |
            | early.lower.gr        | 2.5% of 95th quantile of growth for species in early development stands               |
