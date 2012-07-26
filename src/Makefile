CFLAGS = -O3 -DSAVE_EVERY_STEPS
#CFLAGS = -O3
all: le2d le2d_cf le2d_soa

le2d: le_core.h le_core.c main.c
	gcc $(CFLAGS) -o le2d main.c le_core.c -lm

le2d_cf: le_core.h le_core.c main_cf.c
	gcc $(CFLAGS) -o le2d_cf main_cf.c le_core.c -lm

le2d_soa: le_core.h le_core.c main_soa.c
	gcc $(CFLAGS) -o le2d_soa main_soa.c le_core.c -lm

le2d_ispc: le_core.h le_ispc.ispc
	ispc -O3 le_ispc.ispc -o le_ispc.o -h le_ispc.h
clean:
	rm -f le2d le2d_cf le2d_soa
