pdf = $(patsubst %.rmd, %.pdf, $(wildcard *.rmd))
html = $(patsubst %.rmd, %.html, $(wildcard *.rmd))

all: $(pdf) 

%.pdf: %.rmd
	Rscript --vanilla build.r $< 

clean: 
	rm -rf $(pdf) $(html)