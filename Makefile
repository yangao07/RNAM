CC      = gcc
CFLAGS  = -O3 -c -Wall 
DFLAGS  = -g -Wall
LDFLAGS = -lz -lm -lpthread

SOURCE  = main.c rest_index.c rest_aln.c build_debwt.c bntseq.c kmer_hash.c utils.c kstring.c
OBJS    = $(SOURCE:.c=.o)

PROG    = rest
DEBUG   = gdb_rest
LIT     = lit_rest
DMARCRO = -D __DEBUG__
LMARCRO = -D __LIT__

all: $(SOURCE) $(PROG)

$(PROG): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 

$(DEBUG):
	$(CC) $(DFLAGS) $(SOURCE) $(DMARCRO) -o $@ $(LDFLAGS)
$(LIT):
	$(CC) $(DFLAGS) $(SOURCE) $(LMARCRO) $(DMARCRO) -o $@ $(LDFLAGS)

clean:
	rm -f *.o $(PROG) $(DEBUG)

