CC=gcc
CFLAGS=-O3 -Wall
LIBS=-lhts -lbam -lz -pthread -lbz2

TARGET=bamFilters

DIROBJ=obj
DIRSRC=src

all: $(DIROBJ)/getfa.o $(DIROBJ)/dust.o $(DIROBJ)/bamFilters.o
	$(CC) $(DIROBJ)/bamFilters.o $(DIROBJ)/dust.o $(DIROBJ)/getfa.o $(CFLAGS) $(SAMTOOLS_CFLAGS) -L$(SAMTOOLS_LIBDIR) $(LIBS) -o $(TARGET)

$(DIROBJ)/getfa.o: $(DIRSRC)/getfa.c
	$(CC)  $(CFLAGS) -c -o $(DIROBJ)/getfa.o $(DIRSRC)/getfa.c

$(DIROBJ)/dust.o: $(DIRSRC)/dust.c
	$(CC) $(CFLAGS) -c -o $(DIROBJ)/dust.o $(DIRSRC)/dust.c

$(DIROBJ)/bamFilters.o: $(DIRSRC)/bamFilters.c
	$(CC) $(DIRSRC)/bamFilters.c $(CFLAGS) $(SAMTOOLS_CFLAGS) -L$(SAMTOOLS_LIBDIR) $(LIBS) -c -o $(DIROBJ)/bamFilters.o

clean:
	rm -f $(DIROBJ)/getfa.o $(DIROBJ)/dust.o $(DIROBJ)/bamFilters.o
