CC=             gcc
CFLAGS=         -O3 -fomit-frame-pointer -static
#CFLAGS=         -O3 -g -D_DEBUG
DFLAGS= 

OBJS=		arc.o attachPEinfo.o bubble.o check.o compactEdge.o \
		concatenateEdge.o connect.o contig.o cutTipPreGraph.o cutTip_graph.o \
		darray.o dfib.o dfibHeap.o fib.o fibHeap.o \
		hashFunction.o kmer.o lib.o loadGraph.o loadPath.o \
		loadPreGraph.o localAsm.o main.o map.o mem_manager.o \
		newhash.o node2edge.o orderContig.o output_contig.o output_pregraph.o \
		output_scaffold.o pregraph.o prlHashCtg.o prlHashReads.o prlRead2Ctg.o \
		prlRead2path.o prlReadFillGap.o read2scaf.o readInterval.o stack.o\
		readseq1by1.o ReadTrace.o scaffold.o searchPath.o seq.o splitReps.o #sortContig.o
PROG=           ../SOAPdenovo-Trans-31kmer
#CFLAGS += -DMER31
ifdef 127mer
CFLAGS += -DMER127
PROG = ../SOAPdenovo-Trans-127mer
else
ifdef 63mer
CFLAGS += -DMER63
PROG = ../SOAPdenovo-Trans-63mer
else
CFLAGS += -DMER31
PROG = ../SOAPdenovo-Trans-31mer
endif
endif
ifdef debug
CFLAGS += -DDEBUG
endif
INCLUDES=	-Iinc
SUBDIRS=        . 
LIBPATH=	
LIBS=      	-pthread -lm -lrt -lbam -lz -L./inc
EXTRA_FLAGS=

BIT_ERR = 0
ifeq (,$(findstring $(shell uname -m), x86_64 ppc64 ia64))
BIT_ERR = 1
endif


LINUX = 0
ifneq (,$(findstring Linux,$(shell uname)))
LINUX = 1
EXTRA_FLAGS += -Wl,--hash-style=both
endif

ifneq (,$(findstring $(shell uname -m), x86_64))
CFLAGS += -m64
endif

ifneq (,$(findstring $(shell uname -m), ia64))
CFLAGS += 
endif

ifneq (,$(findstring $(shell uname -m), ppc64))
CFLAGS += -mpowerpc64
endif

.SUFFIXES:.c .o

.c.o:
		@printf "Compiling $<...                             \r"; \
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< || echo "Error in command: $(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $<"

all:		clean SOAPdenovo

.PHONY:all clean install

envTest:
		@test $(BIT_ERR) != 1 || sh -c 'echo "Fatal: 64bit CPU and Operating System required!";false;'

SOAPdenovo:	envTest $(OBJS)
		@$(CC) $(CFLAGS) -o $(PROG) $(OBJS) $(LIBPATH) $(LIBS) $(ENTRAFLAGS)
		@printf "Linking...\r"
		@printf "$(PROG) compilation done.\n";

clean:
		@rm -fr gmon.out *.o a.out *.exe *.dSYM $(PROG) *~ *.a *.so.* *.so *.dylib
		@printf "$(PROG) cleaning done.\n";

install:
		@cp $(PROG) ../bin/
		@printf "$(PROG) installed at ../bin/$(PROG)\n"
