include ../Makefile.in

TARGET= ../run/sol1
TARGET1  = ../run/sol1-intel
TARGET2  = ../run/sol1-a64fx

all	:  $(TARGET) $(TARGET1) $(TARGET2)


SRC     =         \
test1.o

SRC2     =         \
test1_waitio.f

OBJ     = ${SRC:.f=.o}
OBJ1     = ${SRC2:.f=.o}
OBJ2     = ${SRC2:.f=-a64fx.o}

.f.o:
	$(F90) -c -loglist $(F90OPT) $<

$(TARGET):$(OBJ)
	$(F90) -o $@ $(LOPT) $(OBJ) $(LIB)

$(TARGET1):$(OBJ1)
	$(F90) -o $@ $(LOPT) $(F90OPT) $(OBJ1) $(LIB) $(WAITIOLIBS)

$(TARGET2):$(OBJ2)
	$(FC2) -o $@ $(LOPT2) $(FFLAGS2) $(OBJ2) $(LIB2) $(WAITIOLIBS2)

clean    :
	rm -f *.o *.mod *.log $(TARGET) *~

test1_waitio-a64fx.o: test1_waitio.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<
# Makefile End
