
rw <- rand |> subset(ethnicity == "white")
rw$variables$income_rs <- c(scale(rw$variables$income_resid)) # mean 0 sd 1

fml_inc <- map(pgs, \(x) reformulate(c("income_rs", pcs), response = x))
names(fml_inc) <- pgs
mod_inc <- map(fml_inc, \(x) svyglm(x, design = rw))
tidy_inc <- map(mod_inc, broom::tidy, conf.int = TRUE) |>
  list_rbind(names_to = "pgs") |>
  mutate(
    p.value.adj = stats::p.adjust(p.value, "fdr"),
    sig = p.value.adj < 0.05
  ) |>
  filter(term == "income_rs")



tidy_inc |>
  arrange(desc(estimate)) |>
  mutate(
    pgs = nice_names[pgs],
    pgs = fct_reorder(pgs, estimate),
    color = rep(rep(1:8, each = 4), 2)[1:nrow(tidy_inc)],
    color = factor(color)
  ) |>
  ggplot(aes(y = pgs, x = estimate)) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = color, shape = sig),
                  fill = "white") +
  geom_text(aes(label = pgs, x = conf.low, color = color), hjust = "right", size = 3) +
  scale_shape_manual(values = c("FALSE" = "circle filled", "TRUE" = "circle")) +
  scale_color_brewer(type = "qual", palette = 2) +
  theme_minimal() +
  coord_cartesian(xlim = c(-0.2, 0.15)) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank()) +
  labs(
    title = "Standardized coefficients of income on polygenic scores in HRS",
    caption = glue::glue("HRS white sample. Income is age-residualized and standardized.
    Filled circles are significant after multiple testing correction (false discovery rate < 0.05).
    Controls include 10 principal components of genetic data."),
    y = "",
    x = "Beta"
  )


rw_l <- rw$variables |>
  select(all_of(pgs), all_of(pcs), income_rs, income_resid, raedyrs) |>
  pivot_longer(cols = all_of(pgs), names_to = "pgs")

pgs_sig <- tidy_inc$pgs[tidy_inc$sig]

library(ggdark)
rw_l |>
  filter(pgs %in% pgs_sig) |>
  mutate(
    pgs = nice_names[pgs]
  ) |>
  ggplot(aes(income_rs, value)) +
  facet_wrap(vars(pgs)) +
  geom_point(alpha = 0.1, size = 0.1) +
  geom_smooth(aes(color = pgs), method = "loess") +
  coord_cartesian( ylim = c(-1, 1), xlim = c(-2, 2)) +
  dark_theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "Scatterplots of income against significant PGS (fdr < 0.05) with loess smoother",
    subtitle = "Plots show the centre of the distribution but smoother uses all data",
    caption = "Source: Health and Retirement Study (US)",
    x = "Income (residualized and standardized)",
    y = "PGS"
  )

