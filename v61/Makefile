TEXFILE = v61

HEADERS = $(CURDIR)/..

all: build/$(TEXFILE).pdf

# hier Python-Skripte:
build/plot.pdf: python_code/v61.py ../matplotlibrc ../header-matplotlib.tex | build
	# so that matplotlib can find the tex header when running
	# LaTeX in the tmp directory
	# and set the matplotlibrc
	TEXINPUTS=$$(pwd)/..: MATPLOTLIBRC=../matplotlibrc python python_code/v61.py

# hier weitere Abhängigkeiten für build/vXXX.pdf deklarieren:
build/$(TEXFILE).pdf: build/plot.pdf

build/$(TEXFILE).pdf: $(TEXFILE).tex FORCE | build
	TEXINPUTS=$(HEADERS): latexmk -r $(HEADERS)/latexmkrc $<

preview: $(TEXFILE).tex FORCE | build
	TEXINPUTS=$(HEADERS): latexmk -r $(HEADERS)/latexmkrc $<

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
