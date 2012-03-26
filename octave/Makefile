
CC=gcc

CC_FLAGS=-O3
#CC_FLAGS=-Wall #-O3
#DEBUGFLAGS=-ggdb -DDEBUG0

INCLUDE=../c/unwrap

all: unwrap_flood  unwrap_qua


unwrap_flood: $(OBJS)

	mkoctfile unwrap_flood.cc ../c/unwrap/libunwrap.cc -I$(INCLUDE)

unwrap_qua: $(OBJS)

	mkoctfile unwrap_qua.cc ../c/unwrap/libunwrap.cc -I$(INCLUDE)


clean:
	rm -f *.o core *.oct
