DIRS=source

all: force_look
	cd source; $(MAKE) $(MFLAGS)

#clean:
#	cd source; $(MAKE) clean

force_look:
	true

install:
	cd source; $(MAKE) install

