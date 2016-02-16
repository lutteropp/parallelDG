BUILD=build
BASE=present
DEPS=$(wildcard *.tex) $(wildcard *.sty)

default: $(BASE).pdf

$(BUILD)/$(BASE).pdf: $(BASE).tex $(DEPS)
	mkdir -p $(BUILD)
	TEXINPUTS=style:$(TEXINPUTS) latexmk -pdf -halt-on-error -file-line-error -outdir='$(BUILD)' $<

$(BASE).pdf: $(BUILD)/$(BASE).pdf
	mv $< $@ # atomic update
	cp $@ $<

.PHONY: show clean

show: $(BASE).pdf
	gnome-open $< || open $<

clean:
	rm -rf build
