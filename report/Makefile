PANDOC_OPTIONS := --listings --include-in-header=header.tex --natbib --extract-media figs

.PHONY: tex2pdf clean
tex2pdf: report.tex
	latexmk -xelatex $< report.pdf

clean:
	latexmk -c report.tex
	rm -f *.pandoc.tex *.pandoc.pdf

report.pandoc.tex: report.md
	pandoc $(PANDOC_OPTIONS) $< -o $@

report.pandoc.pdf: report.md
	pandoc $(PANDOC_OPTIONS) --pdf-engine=xelatex $< -o $@
