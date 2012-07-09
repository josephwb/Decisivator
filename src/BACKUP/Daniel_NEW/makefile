LIBS = -lm
CC=g++
CPPFLAGS = -g -O3 -m64
SRCS = GAoptimize.cpp
HDRS = GAoptimize.h
OBJS = GAoptimize.o

All: 
	make optimize

optimize: $(OBJS)
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJS) -o GAoptimize
	rm -f *.o *~


