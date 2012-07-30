#CFLAGS = -O3 -DSAVE_EVERY_STEPS -g
CFLAGS = -O3
all: le2d le2d_cf le2d_soa le2d_sse

le2d: le_core.h le_core.c main.c
	gcc $(CFLAGS) -o le2d main.c le_core.c -lm

le2d_cf: le_core.h le_core.c main_cf.c
	gcc $(CFLAGS) -o le2d_cf main_cf.c le_core.c -lm

le2d_soa: le_core.h le_core.c main_soa.c
	gcc $(CFLAGS) -o le2d_soa main_soa.c le_core.c -lm

le2d_ispc: le_core.h le_ispc.ispc main_ispc.c le_core.c
	ispc -O3 le_ispc.ispc -o le_ispc.o -h le_ispc.h
	gcc $(CFLAGS) -o le2d_ispc main_ispc.c le_core.c le_ispc.o -lm

le2d_sse: le_core.h le_core.c le_sse.c main_sse.c
	gcc $(CFLAGS) -o le2d_sse le_sse.c le_core.c main_sse.c -lm

clean:
	rm -f le2d le2d_cf le2d_soa le2d_ispc le2d_sse *.o
