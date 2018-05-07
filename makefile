EXEC   = c_potential.so

SRCS   = c_potential.c # main.c ngbtree3d.c selectb.c


OBJS   = $(SRCS:.c=.o)
INCL   = #forcetree.h proto.h


CFLAGS =  -std=c99 -shared -fPIC -O1 -g  #  -Wall
LNKCMD =  ld -L/usr/lib -L/usr/local/lib  -shared

LIBS   =  -lm 

CC     =  gcc

$(EXEC): $(OBJS) 
	$(CC) -L/usr/lib -L/usr/local/lib  -shared $(OBJS) $(LIBS) $(CFLAGS)  -o $(EXEC)


$(OBJS): $(INCL) 



clean:
	rm $(OBJS) $(EXEC)

