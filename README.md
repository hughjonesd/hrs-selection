# hrs-selection

Investigating natural selection in the US Health and Retirement Study. Now published in [Behavior Genetics](https://link.springer.com/article/10.1007/s10519-024-10189-8).

The paper is at [hrs-selection.pdf](https://github.com/hughjonesd/hrs-selection/blob/master/hrs-selection.pdf).

A quick informal intro is [here](https://open.substack.com/pub/wyclif/p/work-in-progress-on-natural-selection?r=2k7zr&utm_campaign=post&utm_medium=web).

# To reproduce

You'll need [quarto](https://quarto.org) and [R](https://r-project.org).

You'll need data from the RAND HRS, including genetic data, and the PGI repository: see the `setup` chunk in `hrs-selection.qmd` for details on where to put them. We can help if you are stuck.

Clone this repository with `git clone https://github.com/hughjonesd/hrs-selection`. 
Start R within the root directory. It will install [renv](https://rstudio.github.io/renv/articles/renv.html) 
and the required R libraries. Then from the shell, 
run `quarto render hrs-selection.qmd` to produce the PDF. 


