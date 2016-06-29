SHELL = /bin/bash

main = full_thesis.pdf
final = ../GregoryAshton_Thesis_NoLinks.pdf

.PHONY : main clean final

main : $(main)
final : $(final)

$(main): *.tex ../bibliography.bib ../thesis.sty
	pdflatex ${@:.pdf=} && bibtex ${@:.pdf=} && pdflatex ${@:.pdf=} && pdflatex ${@:.pdf=} && pdflatex ${@:.pdf=}

$(final): *.tex ../bibliography.bib ../thesis.sty
	pdflatex ${main:.pdf=} && bibtex ${main:.pdf=} && pdflatex ${main:.pdf=} && pdflatex ${main:.pdf=} && pdflatex '\def\nocolouredlinks{1} \input{${main:.pdf=}}'
	mv ${main:.pdf=}.pdf ../GregoryAshton_Thesis_NoLinks.pdf

	pdflatex ${main:.pdf=} && bibtex ${main:.pdf=} && pdflatex ${main:.pdf=} && pdflatex ${main:.pdf=} && pdflatex ${main:.pdf=}
	mv ${main:.pdf=}.pdf ../GregoryAshton_Thesis_WithLinks.pdf

clean :
	rm -f $(main:.pdf=){.aux,.bbl,.blg,.log,.out,.pdf,.toc,Notes.bib} $(texfigs) $(texonlyfigs)
