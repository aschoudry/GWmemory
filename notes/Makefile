#
# $Header: /home/ashok/gravitational_wave_memory_project/notes/Makefile,v 1.4 2019/01/01 10:57:22 ashok Exp $
#

DOCS = Paper.pdf

BYPROD = Paper.ps Paper.blg

TEXS = Paper.tex 

FIGS = 

default: ${DOCS}

#Paper.pdf: Paper.ps
#	ps2pdf Paper.ps Paper.pdf

#Paper.ps: Paper.dvi ${FIGS}
#	dvips -t a4 -P pdf -o Paper.ps Paper.dvi

Paper.pdf: Paper.tex ${FIGS} ${TEXS}
#	latex Paper && latex Paper && latex Paper && bibtex Paper
	pdflatex Paper  && bibtex Paper && bibtex Paper && pdflatex Paper && pdflatex Paper

Paper.ps: Paper.pdf ${FIGS}
	pdftops Paper.pdf

%.pdf: %.eps
	epstopdf --outfile=$@ ${@:.pdf=.eps}

clean:
	rm -f *.log *.out *.dvi *.end *~ ${DOCS} ${BYPROD}
