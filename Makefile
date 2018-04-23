CC		= g++

CFLAGS		= -g -Wall -fopenmp

CLIBS		= -lm 

all: pmf

pmf: source/pmf.cpp sched.o
	$(CC) $(CFLAGS) -o pmf source/pmf.cpp sched.o

sched.o: source/sched.h source/sched.cpp
	$(CC) $(CFLAGS) -c source/sched.cpp

clean:		
	rm -rf *~ *.o pmf pmf.dSYM a.out core main  
