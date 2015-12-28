SHELL = /bin/bash

main = introduction.pdf

.PHONY : main clean

main : $(main)

$(main): *.tex ../bibliography.bib ../style_file.sty
	pdflatex ${@:.pdf=} && bibtex ${@:.pdf=} && pdflatex ${@:.pdf=} && pdflatex ${@:.pdf=} && pdflatex ${@:.pdf=}

clean :
	rm -f $(main:.pdf=){.aux,.bbl,.blg,.log,.out,.pdf,Notes.bib} $(texfigs) $(texonlyfigs)
