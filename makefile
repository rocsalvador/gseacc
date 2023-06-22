SRC_DIR := src
R_DIR := R
TARGET := gseacc

cc:
	g++ -O3 -Wall -c src/main.cc src/gsea.cc
	g++ -o gsea main.o gsea.o

build:
	rm -rf $(TARGET)
	R -e 'library("Rcpp");filenames <- c(Sys.glob("$(SRC_DIR)/*.cc"), Sys.glob("$(SRC_DIR)/*.hh"));Rcpp.package.skeleton("$(TARGET)", cpp_files = filenames, code_files = "$(R_DIR)/gseacc.R")'
	cp DESCRIPTION $(TARGET)
	cp man/* $(TARGET)/man
	cp  $(SRC_DIR)/Makevars $(TARGET)/src
	R CMD build $(TARGET)

install: build
	R CMD INSTALL $(TARGET)_1.0.tar.gz 

debug:
	g++ -g -c src/main.cc src/gsea.cc
	g++ -o gsea main.o gsea.o

clean:
	rm -rf $(TARGET)_1.0.tar.gz **.o gsea
