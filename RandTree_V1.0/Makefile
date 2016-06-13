CC	= /usr/local/bin/gcc-5

#############################################################################
## NODETYPE specifies what type of node to use in the trees, use flag to change it

#CFLAGS = -g3 
CFLAGS = -g3 -pg

# INCDIRS = -I /usr/local/include
# DEPEND  = makedepend

# *** THE LINK FLAGS ***
#LINKFLAGS = -fast

LIBS = -lm

SRCS    = ptr_tree.c stack.c randtree.c 
OBJS	= ptr_tree.o stack.o randtree.o 
HDRS    = ptr_tree.h stack.h 

all: rt

clean:
	rm -f *.o a.out *.a rt

.c.o: $(SRCS) $(HDRS)
	$(CC) -c $(INCDIRS) $(CFLAGS) $(OPTFLAG) $(DEFS) $<

at_harness: at_harness.c
	$(CC) $(LINKFLAGS) -o $@ $(OBJS) $(LIBDIRS) $(LIBS)

rt: $(OBJS) $(LIBOBJS)
	 $(CC) $(LINKFLAGS) -pg -o $@ $(OBJS) $(LIBDIRS) $(LIBS)
