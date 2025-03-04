include Makefile.in

default: all

all:
	(cd src1 ; make )
	(cd src2 ; make )
	(cd src3 ; make )
clean:
	(cd src1 ; make clean)
	(cd src2 ; make clean)
	(cd src3 ; make clean)
	-rm run/sol[123]-*
