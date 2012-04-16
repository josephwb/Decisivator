LIBS = -lm
CC=g++
CPPFLAGS = -g -O3 -m64
SRCS = optimize.cpp
HDRS = optimize.h
OBJS = optimize.o

All: 
	make optimize

optimize: $(OBJS)
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJS) -o optimize
	rm -f *.o *~


