library(rmarkdown)
library(knitr)

args = commandArgs(trailingOnly=TRUE)

#sapply(c("questions", "answers"), dir.create, showWarnings = FALSE)


render(args[1], output_format = "all")
