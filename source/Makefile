CC=gcc
CFLAGS=-I.
DEPS=

multalin: afichseq.o aligne.o basproc1.o cluster.o coefdisk.o commande.o copright.o disk.o drivers.o fast.o gcgdisk.o lstfseq.o ma.o maglob.o mbgbdisk.o msfdisk.o msgerr.o muldisk.o parametr.o portab.o swapd.o util.o 
	gcc -lm -o mutalin *.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 
	
