SRC_DIR := src
R_DIR := R
TARGET := gseacc

cc:
	g++ -O3 -Wall -c src/main.cc src/gsea.cc
	g++ -o gsea main.o gsea.o

install:
	R CMD build .
	R CMD INSTALL $(TARGET)_1.0.tar.gz

debug:
	g++ -g -c src/main.cc src/gsea.cc
	g++ -o gsea main.o gsea.o

clean:
	rm -rf $(TARGET)_1.0.tar.gz **.o gsea
