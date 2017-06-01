CC=gcc
CXX=g++
CFLACGS=-I.
PREFIX=/usr/local

all: fqextract fqparse

fqextract: fqextract.c
	$(CC) -o fqextract fqextract.c -I.

fqparse: fqparse.cpp
	$(CXX) -Wall -o fqparse fqparse.cpp


.PHONY: install
install: fqextract fqparse
	mkdir -p $(PREFIX)/bin
	cp $< $(PREFIX)/bin/fqextract
	cp $< $(PREFIX)/bin/fqparse

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/bin/fqextract
	rm -f $(PREFIX)/bin/fqparse

