#all: pdf word Rscript
all: pdf 

pdf: $(patsubst %.rmd, %.pdf, $(wildcard *.rmd))
#word: $(patsubst %.Rmd, %.docx, $(wildcard *.Rmd))
Rscript: $(patsubst %.Rmd, %.R, $(wildcard *.Rmd))

%.md: %.rmd
	Rscript -e "library(knitr); knit('$<')"	

%.pdf: %.md 
	pandoc -t beamer -V theme:UoB $< -o $@

#%.docx: %.md 
#	pandoc -H format.sty -V fontsize=12pt $< -o $@  

%.R: %.rmd
	Rscript -e "library(knitr); purl('$<')"	