pdf = $(patsubst %.rmd, questions/%.pdf, $(wildcard *.rmd))
pdf.final = $(patsubst %.rmd, questions/%-final.pdf, $(wildcard *.rmd))

all: $(pdf) $(pdf.final)

questions/%.pdf: %.rmd
	Rscript --vanilla build.r $< 

## annoying way to trim out whitespace slides
questions/%-final.pdf: questions/%.pdf
	pdftk A=$< cat A1-18 A22-32 A39 output $@

clean: 
	rm -rf questions/ answers/