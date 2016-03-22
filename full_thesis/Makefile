SHELL = /bin/bash

main = full_thesis.pdf

.PHONY : main clean

main : $(main)

$(main): *.tex ../bibliography.bib ../thesis.sty
	pdflatex ${@:.pdf=} && bibtex ${@:.pdf=} && pdflatex ${@:.pdf=} && pdflatex ${@:.pdf=} && pdflatex ${@:.pdf=}

clean :
	rm -f $(main:.pdf=){.aux,.bbl,.blg,.log,.out,.pdf,.toc,Notes.bib} $(texfigs) $(texonlyfigs)
