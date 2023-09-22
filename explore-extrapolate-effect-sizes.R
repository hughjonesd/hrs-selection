
# == Setup ====

library(ggplot2)
library(purrr)
library(dplyr)
library(stringr)
library(survey)
library(tidyr)

rand_orig <- readRDS("data/randhrs/rand-orig.rds")

# ragender: 1 male, 2 female
rand_orig <- rand_orig |>
  filter(
    (rabyear <= 1965 & as.numeric(ragender) == 1) |
    (rabyear <= 1970 & as.numeric(ragender) == 2)
  )

# african ancestry
pgen4a <- haven::read_dta("data/PGENSCORE4r3/pgenscore4a_r.dta")
# european ancestry
pgen4e <- haven::read_dta("data/PGENSCORE4r3/pgenscore4e_r.dta")

pgens <- list(white = pgen4e, black = pgen4a)

pcs <- c(paste0("PC1_5", LETTERS[1:5]), paste0("PC6_10", LETTERS[1:5]))

rands <- list()

eths <- c("black", "white")
for (eth in eths) {
  pgen <- pgens[[eth]]
  names(pgen) <- str_replace(names(pgen), "[EA]4_", "")
  names(pgen) <- str_replace(names(pgen), "^01", "GW01") # invalid names
  rands[[eth]] <- inner_join(rand_orig, pgen,
                     join_by(hhid, pn),
                     unmatched = "drop",
                     relationship = "one-to-one")
}

rand <- list_rbind(rands, names_to = "ethnicity")

trk_orig <- haven::read_dta("data/trk2020v3/trk2020tr_r.dta")
trk_weights <- trk_orig |> select(HHID, PN, matches("BIOWGTR"))

rand <- left_join(rand, trk_weights, by = join_by(hhid == HHID, pn == PN))
rand$weight <- ifelse(is.na(rand$MBIOWGTR), rand$NBIOWGTR, rand$MBIOWGTR)

all_pgs <- pgen4e |>
    rename_with(\(x) str_replace(x, "[EA]4_", "")) |>
    rename_with(\(x) str_replace(x, "^01", "GW01")) |>
    # weird that I have to operate on whole df to use select semantics
    select(-hhid, -pn, -version,
           -matches("^PC\\d"),
           ) |>
    names()

# scores and phenotypes:
# EDU2, EDU3, EA3: raedyrs -- educational attainment
# HEIGHT and HEIGHT2: r*height -- height
# BMI and BMI2: r*bmi -- bmi
# E4/A4_AD/AD2_IGAP13, GW/01AD2NA/WA_IGAP19: r*alzhee from wave 10, also
# related including r*memrye up to wave 9 ("memory-related disease")
# -- alzheimer's, nb each gen has diff scores
# GENCOG and GCOG2: r*cogtot -- general cognition/cognitive function
# EVRSMK_TAG10, SI: r*smokev -- smoking initiation ever/never. also r*smoken,
# smoke now
# LONGEVITY, radage_y -- outdated but apparently no replacement?
# CPD, CPD_GSCAN19, available for now and also (some years) "when you were
# smoking the most". But only in original file, not in RAND.
# See question concordance. -- cigarettes per day (average or max)
# ADHD_PGC10, ADHD_PGC17, ADHD
# MDD_PGC13, MDD2, r*cesd -- major depressive disorder
pgs_list <- list(
  # target phenotype = list of successive pgs names
  ph_education = c("EDU2_SSGAC16", "EDU3_SSGAC18", "EA3_W23_SSGAC18"),
  ph_height = c("HEIGHT_GIANT14", "HEIGHT2_GIANT18"),
  ph_bmi = c("BMI_GIANT15", "BMI2_GIANT18"),
  ph_alzheimers = c("AD_IGAP13", "AD2_IGAP13", "GWAD2NA_IGAP19", "GWAD2WA_IGAP19",
              "GW01AD2NA_IGAP19", "GW01AD2WA_IGAP19"),
  ph_cognition = c("GENCOG_CHARGE15", "GCOG2_CHARGE18"),
  ph_eversmoke = c("EVRSMK_TAG10", "SI_GSCAN19"),
  # ph_longevity = c("LONGEVITY"),
  ph_cigsperday = c("CPD_TAG10", "CPD_GSCAN19"),
  ph_adhd = c("ADHD_PGC10", "ADHD_PGC17"),
  ph_depression = c("MDD_PGC13", "MDD2_PGC18")
)

rand <- rand |>
  mutate(
    ph_education = as.numeric(raedyrs),
    ph_height = r10height, # r10 doesn't have many missing; not much gained from adding more
    ph_bmi = r10bmi,
    ph_alzheimers = pmax(as.numeric(r7memrye), as.numeric(r10alzhee), na.rm = TRUE),
    ph_cognition = pmax(r9cogtot, r10cogtot, r11cogtot, na.rm = TRUE),
    ph_eversmoke = as.numeric(r10smokev),
    # ph_cigsperday = NA_real_, # for now
    # ph_adhd = NA_real_, # nothing obvious
    ph_depression = r10cesd,
  ) |>
  mutate(.by = rabyear, # calculate this pooling ethnicity
    mean_raevbrn = mean(raevbrn, na.rm = TRUE),
    rlrs = raevbrn/mean_raevbrn
  )

rand_long <- rand
rand <- rand_long |>
  select(
    matches("^ph_"),
    all_of(all_pgs),
    all_of(pcs),
    hhid, pn,
    weight,
    ragender,
    raevbrn,
    rlrs,
    ethnicity,
    raehsamp,
    raestrat
  )

randsv <- svydesign(
    id = ~ raehsamp,
		strata = ~ raestrat,
		weights = ~ weight,
		nest = TRUE,
		data = rand |> drop_na(weight)
		)
randsv <- randsv |> subset(ethnicity == "white")


# == regressions ====


# from the man himself
# https://stackoverflow.com/a/73956738/10522567
rsquared <- function (dep_var, design, mod) {
  dep_var <- reformulate(dep_var)
  total_var <-svyvar(dep_var, design, na.rm = TRUE)
  resids <-
  resid_var <- summary(mod)$dispersion
  1 - resid_var/total_var
}


phenotypes <- intersect(names(pgs_list), names(randsv$variables))
rsquares_phen <- list()
rsquares_fert <- list()
for (phenotype in phenotypes) {
  pgs <- pgs_list[[phenotype]]

  fmls <- map(pgs, \(x) reformulate(c(x, pcs), response = phenotype))
  mods <- map(fmls, \(x) svyglm(x, design = randsv))
  rsqs <- map_dbl(mods, \(m) rsquared(phenotype, design = randsv, mod = m))

  fmls_fert <- map(pgs, \(x) reformulate(c(x, pcs), response = "raevbrn"))
  mods_fert <- map(fmls_fert, \(x) svyglm(x, design = randsv))
  rsqs_fert <- map_dbl(mods_fert, \(m) rsquared("raevbrn", design = randsv, mod = m))

  names(rsqs) <- names(rsqs_fert) <- pgs
  rsquares_phen[[phenotype]] <- rsqs
  rsquares_fert[[phenotype]] <- rsqs_fert
}

for (phenotype in phenotypes) {
  dfr <- data.frame(fert = rsquares_fert[[phenotype]], phen = rsquares_phen[[phenotype]])
  m1 <- lm(fert ~ phen, dfr)
  m2 <- lm(fert ~ 0 + phen, dfr)
  p1 <- predict(m1, newdata = data.frame(phen = 0.4))
  p2 <- predict(m2, newdata = data.frame(phen = 0.4))
  plot(rsquares_phen[[phenotype]], rsquares_fert[[phenotype]],
       main = phenotype,
       sub = paste(round(p1, 3), "/", round(p2, 3)),
       xlim = c(0, 0.5), ylim = c(0, 0.5),
       xlab = "R2 on phenotype", ylab = "R2 on fertility")
  points(c(0.4, 0.4), c(p1, p2), pch = 3)
  abline(m1, col = "red")
  abline(m2, col = "blue")
}

# == bootstraps on edu pgs ====

# from the PGI repository guide article:
# in the special case of a univariate regression,
# the bias is b_ghat = 1/rho b_g
# where rho = sqrt(h^2/R^2), h^2 is [chip-]heritability, R^2 is R-squared
# of regression on target phenotype;
# b_ghat is the observed regression using your pgs;
# b_g is the true effect of the "true pgs".

boot_edu <- function (w, data) {
  edu_pgs <- c("EDU2_SSGAC16", "EDU3_SSGAC18", "EA3_W23_SSGAC18")
  fs_fert <- map(edu_pgs, \(x) reformulate(c(x, pcs), response = "rlrs"))
  fs_edu <- map(edu_pgs, \(x) reformulate(c(x, pcs), response = "ph_education"))

  ms_fert <- map(fs_fert, \(x) lm(x, data = data, weights = w))
  ms_edu <- map(fs_edu, \(x) lm(x, data = data, weights = w))
  bhat_fert <- map2_dbl(ms_fert, edu_pgs, \(x, y) coef(x)[[y]])
  bhat_edu <- map2_dbl(ms_edu, edu_pgs, \(x, y) coef(x)[[y]])

  r2_edu <- map_dbl(ms_edu, \(x) {
    an <- car::Anova(x)
    an$`Sum Sq`[1]/sum(an$`Sum Sq`)  # incremental R2 of first term
  })

  snp_h2 <- 0.10 # From Abdel and Karen "Dissecting..."
  twin_h2 <- 0.43
  r2_edu <- c(r2_edu, snp_h2, twin_h2)
  # the correction factor
  rhos <- outer(r2_edu, r2_edu, \(x, y) sqrt(x/y)) # lower triangular part of this
                                                   # shows how much you multiply
                                                   # to get from COLNAME to ROWNAME
  rownames(rhos) <- colnames(rhos) <- c(edu_pgs, "Chip", "Twin")

  # multiply first column by bhat_edu[1] to give predictions from the
  # first PGS to the next ones; etc.
  # The diagonals of b_edu are exactly equal to bhat_edu
  # arow,column gives the predicted value of b by ROWNAME for COLNAME
  b_edu <- apply(rhos, 1, \(x) x * c(bhat_edu, NA, NA))
  b_fert <- apply(rhos, 1, \(x) x * c(bhat_fert, NA, NA))

  # don't include the fourth/fifth rows of NA predictions "from chip/twin h2"
  b_fert <- b_fert[1:3, ]

  res <- c(b_fert)
  res_names <- outer(rownames(b_fert), colnames(b_fert),
                      \(x, y) paste0(x, ".to.", y))
  names(res) <- c(res_names)
  res
}

rand_boot <- randsv |> as.svrepdesign(type = "bootstrap", replicates = 99)
results <- withReplicates(rand_boot, boot_edu)
results_df <- results |>
  as.data.frame() |>
  tibble::rownames_to_column(var = "rn") |>
  cbind(confint(results)) |>
  tibble::remove_rownames() |>
  setNames(c("rn", "mean", "SE", "conf.low", "conf.high")) |>
  mutate(
    source.pgs = str_remove(rn, ".to.*"),
    target.pgs = str_remove(rn, ".*.to."),
    target.pgs = forcats::fct_relevel(target.pgs, "Twin", "Chip", "EA3_W23_SSGAC18",
                                  "EDU3_SSGAC18", "EDU2_SSGAC16"),
    rn = NULL,
    .before = 1
  ) |>
  arrange(target.pgs, source.pgs)

results_df |>
  mutate(
    Target = target.pgs,
    Source = ifelse(source.pgs == target.pgs, "True", source.pgs)
  ) |>
  ggplot(aes(y = Target, color = Source)) +
    geom_pointrange(aes(x = mean, xmin = conf.low, xmax = conf.high,
                        linetype = Target == "Twin"),
                    position = position_dodge(.4)) +
    scale_color_manual(values = c(True = "red",
                                  EDU2_SSGAC16 = "steelblue",
                                  EDU3_SSGAC18 = "green4",
                                  EA3_W23_SSGAC18 = "darkblue")) +
    scale_linetype(guide = "none") +
    labs(x = "Corrected effect") +
    theme_minimal()



# the confidence intervals usually seem wider than if you use return.replicates = TRUE
# and calculate them yourself. Even though they appear to be just normal-theory c.i.s

# See this useful stackoverflow question and answer:
# https://stackoverflow.com/questions/75028201/how-does-the-r-survey-package-compute-bootstrap-confidence-intervals
# == old stuff ====

edu_pgs <- c("EDU2", "EDU3", "EA3")

edu_r2 <- edu_pgs |>
  map(\(x) reformulate(c(x, pcs), response = "raedyrs")) |>
  map(\(x) lm(x, data = rands$white, weights = weight)) |>
  map(broom::glance) |>
  setNames(edu_pgs) |>
  list_rbind(names_to = "pgs")

edu_r2_evb <- edu_pgs |>
  map(\(x) reformulate(c(x, pcs), response = "raevbrn")) |>
  map(\(x) lm(x, data = rands$white, weights = weight)) |>
  map(broom::glance) |>
  setNames(edu_pgs) |>
  list_rbind(names_to = "pgs")

edu_pgs_eff <- edu_pgs |>
  map(\(x)reformulate(c(x, pcs), response = "raevbrn")) |>
  map(\(x) lm(x, data = rands$white, weights = weight)) |>
  map(broom::tidy) |>
  list_rbind() |>
  filter(term %in% edu_pgs)

effs <- data.frame(edu_r2 = edu_r2$r.squared,
                   evb_r2 = edu_r2_evb$r.squared,
                   evb_eff = edu_pgs_eff$estimate)

# obviously the point about these extrapolations is they're ridiculous...
# also they don't even include the effects of sampling variation on the point
# estimates!

# also I should remind myself whether the errors-in-variables model expects
# that the R2 are linearly related or the effect size(s)? I think it's the R2...

ggplot(effs, aes(edu_r2, evb_r2)) +
  geom_smooth(method = "lm", fullrange = TRUE) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, color = "red", fill = "red",
              fullrange = TRUE) +
  geom_point(size = 3) +
  xlim(0, 0.45) + ylim(0, 0.025)

ggplot(effs, aes(edu_r2, evb_eff)) +
  geom_smooth(method = "lm", fullrange = TRUE) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, color = "red", fill = "red",
              fullrange = TRUE) +
  geom_point(size = 3) +
  xlim(0, 0.45)
