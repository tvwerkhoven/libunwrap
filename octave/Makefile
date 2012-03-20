
CC=gcc

CC_FLAGS=-O3
#CC_FLAGS=-Wall #-O3
#DEBUGFLAGS=-ggdb -DDEBUG0

INCLUDE=../c/unwrap

all: unwrap_oct

unwrap_oct: $(OBJS)

	mkoctfile unwrap_oct.cc ../c/unwrap/libunwrap.c -I$(INCLUDE)

clean:
	rm -f *.o core *.oct
