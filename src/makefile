EXEC   = c_potential.so

SRCS   = c_pairwise_dist.c c_potential.c multi_potential.c


OBJS   = $(SRCS:.c=.o)
INCL   = include.h #forcetree.h proto.h


CFLAGS =  -std=c99 -shared -fPIC -O3 -fopenmp  #  -Wall
LNKCMD =  ld -L/usr/lib -L/usr/local/lib  -shared

LIBS   =  -lm 

CC     =  gcc

$(EXEC): $(OBJS) 
	$(CC) -L/usr/lib -L/usr/local/lib  -shared $(OBJS) $(LIBS) $(CFLAGS)  -o $(EXEC)


$(OBJS): $(INCL) 



clean:
	rm $(OBJS) $(EXEC)
