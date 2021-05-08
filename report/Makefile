PANDOC_OPTIONS := --listings --include-in-header=header.tex


report.tex: report.md
	pandoc $(PANDOC_OPTIONS) $< -o $@

report.pdf: report.md
	pandoc $(PANDOC_OPTIONS) --pdf-engine=xelatex $< -o $@
