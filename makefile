.PHONY: all cpu gpu clean

all: gpu

cpu:
	$(MAKE) -f makefile.cpu

gpu:
	$(MAKE) -f makefile.cuda

clean:
	$(MAKE) -f makefile.cpu clean
	$(MAKE) -f makefile.cuda clean
