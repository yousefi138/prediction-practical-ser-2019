pdf = $(patsubst %.rmd, questions/%.pdf, $(wildcard *.rmd))

all: $(pdf)

questions/%.pdf: %.rmd
	Rscript --vanilla build.r $< 

clean: 
	rm -rf questions/ answers/