library(rmarkdown)
library(knitr)

args = commandArgs(trailingOnly=TRUE)

sapply(c("questions", "answers"), dir.create, showWarnings = FALSE)

has.answers <- FALSE
has.answer.figs <-'hide'
render(args[1], output_dir = "questions/",
			output_format = "all")

has.answers <- TRUE
has.answer.figs <-'asis'
render(args[1], output_dir = "answers/",
			output_format = "all")

