CC=gcc
CFLACGS=-I.
PREFIX=/usr/local

fqextract: fqextract.c
        $(CC) -o fqextract fqextract.c -I.


.PHONY: install
install: fqextract
        mkdir -p $(PREFIX)/bin
        cp $< $(PREFIX)/bin/fqextract

.PHONY: uninstall
uninstall:
        rm -f $(PREFIX)/bin/fqextract

