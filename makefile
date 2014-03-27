CC=gcc
IDIR=../include
LIBS=-lm
CFLAGS=-DDEBUG -Wall -g $(LIBS) -I$(IDIR)
VPATH=src
BUILDDIR=build

all: tree_2d.o tree_hybrid.o main.o
	cd $(BUILDDIR) && $(CC) $(CFLAGS) tree_2d.o tree_hybrid.o main.o -o main

test: tree_2d.o tree_hybrid.o test.o
	cd $(BUILDDIR) && $(CC) $(CFLAGS) tree_2d.o tree_hybrid.o test.o -o test

.c.o: builddir
	mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $(BUILDDIR)/$@

clean:
	rm -r ./build/
