echo "TODO: create a makefile for this"
find *.c | xargs -n 1 cc -c
cc -lm afichseq.o aligne.o basproc1.o cluster.o coefdisk.c commande.o copright.o disk.o drivers.o fast.o gcgdisk.o lstfseq.o ma.o maglob.o mbgbdisk.o msfdisk.o msgerr.o muldisk.o parametr.o portab.o swapd.o   util.o -o ma
echo "run: ./ma to test"
