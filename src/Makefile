
OBJS1 = sort.o
OBJS2 = myrand.o sort.o slfast_2D.o leasq_sparse.o eigen.o dirplg.o switchsub.o stridip.o
OBJS4 = myrand.o sort.o slfast_4D.o leasq_sparse.o eigen.o dirplg.o switchsub.o stridip.o

CFLAGS= -O -g -m64

all:  libs programs

libs: slfast_2D.o slfast_4D.o leasq_sparse.o myrand.o dirplg.o eigen.o switchsub.o stridip.o sort.o

programs: boot_uncert bootmech_2D bootmech_4D bootmech_4D_test satsi_2D satsi_2D_tradeoff satsifast_2D satsi_4D satsi_4D_tradeoff satsifast_4D 

boot_uncert: ${OBJS1} boot_uncert.c
	gcc ${CFLAGS} boot_uncert.c -o boot_uncert ${OBJS1} -lm
bootmech_2D: ${OBJS2} bootmech_2D.c
	gcc ${CFLAGS} bootmech_2D.c -o bootmech_2D ${OBJS2} -lm
bootmech_4D: ${OBJS4} bootmech_4D.c
	gcc ${CFLAGS} bootmech_4D.c -o bootmech_4D ${OBJS4} -lm
bootmech_4D_test: ${OBJS4} bootmech_4D_test.c
	gcc ${CFLAGS} bootmech_4D_test.c -o bootmech_4D_test ${OBJS4} -lm
satsi_2D: ${OBJS2} satsi_2D.c
	gcc ${CFLAGS} satsi_2D.c -o satsi_2D ${OBJS2} -lm
satsi_2D_tradeoff: ${OBJS2} satsi_2D.c
	gcc ${CFLAGS} satsi_2D_tradeoff.c -o satsi_2D_tradeoff ${OBJS2} -lm
satsifast_2D: ${OBJS2} satsifast_2D.c
	gcc ${CFLAGS} satsifast_2D.c -o satsifast_2D ${OBJS2} -lm
satsi_4D: ${OBJS4} satsi_4D.c
	gcc ${CFLAGS} satsi_4D.c -o satsi_4D ${OBJS4} -lm
satsi_4D_tradeoff: ${OBJS4} satsi_4D.c
	gcc ${CFLAGS} satsi_4D_tradeoff.c -o satsi_4D_tradeoff ${OBJS4} -lm
satsifast_4D: ${OBJS4} satsifast_4D.c
	gcc ${CFLAGS} satsifast_4D.c -o satsifast_4D ${OBJS4} -lm

slfast_2D.o: slfast_2D.c 
	gcc -c ${CFLAGS} slfast_2D.c 
slfast_4D.o: slfast_4D.c 
	gcc -c ${CFLAGS} slfast_4D.c 
leasq_sparse.o: leasq_sparse.c 
	gcc -c ${CFLAGS} leasq_sparse.c 
myrand.o: myrand.c 
	gcc -c ${CFLAGS} myrand.c 
dirplg.o: dirplg.c 
	gcc -c ${CFLAGS} dirplg.c 
eigen.o: eigen.c 
	gcc -c ${CFLAGS} eigen.c 
switchsub.o: switchsub.c 
	gcc -c ${CFLAGS} switchsub.c 
stridip.o: stridip.c 
	gcc -c ${CFLAGS} stridip.c 
sort.o: sort.c 
	gcc -c ${CFLAGS} sort.c 

clean: 
	rm *.o
