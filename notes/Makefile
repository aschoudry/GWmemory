USE_BIBTEX=yes

all: Paper.pdf

%.pdf: %.tex
	pdflatex $<
ifeq ($(USE_BIBTEX),yes)
	-bibtex $*
endif
	pdflatex $<
	pdflatex $<
	rm -f *.aux *.blg *.log

clean:
	rm -f *.log *.aux *.out *.toc *.snm *.nav *.bbl *.blg *~

veryclean:
	make clean
	rm -f *.dvi *.ps *.pdf

spellchk:
	ispell *.tex
