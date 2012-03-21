
CC=gcc

#CC_FLAGS=-O3
CC_FLAGS=-Wall #-O3
DEBUGFLAGS=-ggdb -DHAVE_DEBUG
CC_FLAGS_BONUS=-std=c99

#INCLUDE=../c/unwrap

TESTOBJS=libunwrap.o test_unwrap.o
LIB=-lm

all: test_unwrap

test_unwrap: $(TESTOBJS)

	$(CC) $(TESTOBJS) $(LIB) -o test_unwrap



test_unwrap.o : test_unwrap.c  libunwrap.h
libunwrap.o : libunwrap.c  libunwrap.h


clean:
	rm -f *.o core *.oct
.c.o:
	$(CC) $(CC_FLAGS) $(DEBUGFLAGS) $(CC_FLAGS_BONUS) $(INCLUDE) -c $<
