PANDOC_OPTIONS := --toc

report.tex: report.md
	pandoc $(PANDOC_OPTIONS) --standalone $< -o $@

report.pdf: report.md
	pandoc $(PANDOC_OPTIONS) --pdf-engine=xelatex $< -o $@
