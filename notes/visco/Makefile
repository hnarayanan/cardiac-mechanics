# Modified from the simple build system I used for managing my
# doctoral dissertation.
#
# Copyright (C) 2007--2011 Harish Narayanan

BASENAME = visco

SOURCES = Makefile ${BASENAME}.tex \
	  ${BASENAME}.bib ${BASENAME}.bst \
	  auxiliary/*.tex sections/*.tex \
          images/pdf/*.pdf

all: pristine ${BASENAME}.pdf
	make clean

clean:
	@echo "Cleaning cruft:"
	rm -f *~
	rm -f ${BASENAME}.4ct
	rm -f ${BASENAME}.4tc
	rm -f ${BASENAME}.aux
	rm -f ${BASENAME}.bbl
	rm -f ${BASENAME}.blg
	rm -f ${BASENAME}.brf
	rm -f ${BASENAME}.dvi
	rm -f ${BASENAME}.idv
	rm -f ${BASENAME}.lg
	rm -f ${BASENAME}.log
	rm -f ${BASENAME}.out
	rm -f ${BASENAME}.spl
	rm -f ${BASENAME}.tmp
	rm -f ${BASENAME}.xref
	(cd auxiliary; make clean)
	(cd sections; make clean)

pristine: clean
	@echo "Removing older output files:"
	rm -f *.css
	rm -f *.html
	rm -f *.pdf
	rm -f *.png

${BASENAME}.pdf: ${SOURCES}
	@echo "Creating the PDF file:"
	pdflatex -interaction=batchmode ${BASENAME} && \
	bibtex -terse ${BASENAME} && \
	pdflatex -interaction=batchmode ${BASENAME} && \
	pdflatex -interaction=batchmode ${BASENAME} && \
	pdflatex -interaction=batchmode ${BASENAME}
	egrep --color=auto -i '(warn|error|full)' ${BASENAME}.log