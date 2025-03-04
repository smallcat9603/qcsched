# before compilation, do not forget "module load metis/4.0.3" for methis 
#INSTALL_DIR = /work/01/jh210022a/share/waitio/opt/
INSTALL_DIR = /work/share/

INCLUDE   = -I$(INSTALL_DIR)/waitio/include/ -I../waitio/
FFLAGS    =   -O3
CFLAGS    =   -O3
LDFLAGS    =  -O3  -L$(IMPIBASE)/lib/ -L$(IMPIBASE)/lib/release/ -L$(IMPIBASE)/../../compiler/2021.2.0/linux/compiler/lib/intel64_lin -L$(IMPIBASE)/../../compiler/2021.2.0/linux/lib /work/opt/local/x86_64/cores/intel/2021.2.0/compiler/2021.2.0/linux/compiler/lib/intel64_lin/libimf.so
LIB       = 
F90 = mpiifort
CC = mpiicc

IMPIBASE = /work/opt/local/x86_64/cores/intel/2021.2.0/mpi/2021.2.0
CC2     =  $(IMPIBASE)/bin/mpiicc
FC2     =  $(IMPIBASE)/bin/mpiifort
INCLUDE2 = -I$(IMPIBASE)/include/ -I$(INSTALL_DIR)/waitio/include/  -I../waitio/

F90OPTa= -align array64byte -O3 -axCORE-AVX512 -qopenmp $(INCLUDE2)
LOPT  =   -O3 -axCORE-AVX512 -qopenmp $(INCLUDE2)
CLOPT  =  -O3 -axCORE-AVX512 -qopenmp $(INCLUDE2)
FLIBS  =  /work/opt/local/x86_64/cores/intel/2021.2.0/compiler/2021.2.0/linux/compiler/lib/intel64_lin/libimf.so

CC3=/opt/FJSVxtclanga/tcsds-1.2.31/bin/mpifccpx
FC3=/opt/FJSVxtclanga/tcsds-1.2.31/bin/mpifrtpx
CFLAGS3    =   -Kfast -DTOFU -g -Kopenmp  -I$(INSTALL_DIR)/waitio/include/  -I../waitio/
FFLAGS3    =   -Kfast -DTOFU -g -Kopenmp -I$(INSTALL_DIR)/waitio/include/  -I../waitio/
LDFLAGS3    =  -Kfast -g -Kopenmp
FLIBS3  = 
LOPT3  =  -O3 -Kopenmp
INCLUDE3 =  -I$(INSTALL_DIR)/waitio/include/

WAITIOLIBS	= -L../waitio/ -L$(INSTALL_DIR)/waitio/lib -lwaitio_intel -lwaitio_mpi_intel
WAITIOLIBS2	= -L../waitio/ -L$(INSTALL_DIR)/waitio/lib -lwaitio_a64fx -lwaitio_mpi_a64fx

.f.o:
	$(F90) -c $(F90OPTa)  $(FFLAGS) $< 

.c.o:
	$(CC) -c $(CFLAGS) $< 

test-intel:	test.c
	$(CC2) -o $@ $<  $(LOPT) $(CLOPT) $(LIB) $(FLIBS) $(WAITIOLIBS) -L$(IMPIBASE)/lib/ $(LDFLAGS) $(INCLUDE2)

test-mpi-intel:	test-mpi.c
	$(CC2) -o $@ $<  $(LOPT) $(CLOPT) $(LIB) $(FLIBS) $(WAITIOLIBS) -L$(IMPIBASE)/lib/ $(LDFLAGS) $(INCLUDE2)

test-a64fx:	test.c
	$(CC3) -o $@ $< $(LOPT3) $(CFLAGS3) $(LIB3) $(FLIBS3) $(WAITIOLIBS2) -L$(IMPIBASE)/lib/ $(LDFLAGS3)

test-mpi-a64fx:	test-mpi.c
	$(CC3) -o $@ $< $(LOPT3) $(CFLAGS3) $(LIB3) $(FLIBS3) $(WAITIOLIBS2) -L$(IMPIBASE)/lib/ $(LDFLAGS3)

TARGET = test-intel test-a64fx test-mpi-intel test-mpi-a64fx

clean    :
	rm -f *.o *.mod *.log $(TARGET)  *~ *.lst *.done

all	: $(TARGET)
