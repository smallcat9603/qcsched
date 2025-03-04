include ../Makefile.in

TARGET  = ../run/sol3
TARGET1  = ../run/sol3-intel
TARGET2  = ../run/sol3-a64fx

all	:  $(TARGET) $(TARGET1) $(TARGET2)

SRC     =         \
hpcmw_util.o hpcmw_fem_cntl.o hpcmw_fem_util.o hpcmw_solver_cntl.o\
hpcmw_fem_mesh.o hpcmw_solver_matrix.o hpcmw_all.o\
solver_SR.o solver_CG.o \
solver11.o test1.o util.o hpcmw_init.o input_cntl.o input_grid.o\
mat_con0.o mat_con1.o mat_ass_main_1.o\
mat_ass_bc.o hpcmw_finalize.o rcm.o

SRC2     =         \
hpcmw_util_waitio.f hpcmw_fem_cntl.f hpcmw_fem_util.f hpcmw_solver_cntl.f\
hpcmw_fem_mesh.f hpcmw_solver_matrix.f hpcmw_all.f\
solver_SR.f solver_CG.f \
solver11.f test1_waitio.f util.f hpcmw_init_waitio.f input_cntl.f input_grid.f\
mat_con0.f mat_con1.f mat_ass_main_1.f\
mat_ass_bc.f hpcmw_finalize.f rcm.f


OBJ     = ${SRC:.f=.o}
OBJ1     = ${SRC2:.f=.o}
OBJ2     = ${SRC2:.f=-a64fx.o}

.f.o:
	$(F90) -c $(F90OPT) $<

$(TARGET):$(OBJ)
	$(F90) -o $@ $(LOPT) $(F90OPT) $(OBJ) $(LIB)

$(TARGET1):$(OBJ1)
	$(F90) -o $@ $(LOPT) $(F90OPT) $(OBJ1) $(LIB) $(WAITIOLIBS)

$(TARGET2):$(OBJ2)
	$(FC2) -o $@ $(LOPT2) $(FFLAGS2) $(OBJ2) $(LIB2) $(WAITIOLIBS2)

clean    :
	rm -f *.o *.mod *.log $(TARGET) *~ *.lst

hpcmw_util_waitio-a64fx.o: hpcmw_util_waitio.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

hpcmw_fem_cntl-a64fx.o: hpcmw_fem_cntl.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

hpcmw_fem_util-a64fx.o: hpcmw_fem_util.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

hpcmw_solver_cntl-a64fx.o : hpcmw_solver_cntl.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

hpcmw_fem_mesh-a64fx.o: hpcmw_fem_mesh.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

hpcmw_solver_matrix-a64fx.o: hpcmw_solver_matrix.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

hpcmw_all-a64fx.o: hpcmw_all.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

solver_SR-a64fx.o: solver_SR.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

solver_CG-a64fx.o: solver_CG.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

solver11-a64fx.o: solver11.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

test1_waitio-a64fx.o: test1_waitio.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

util2-a64fx.o: util2.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

hpcmw_init_waitio-a64fx.o: hpcmw_init_waitio.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

input_cntl-a64fx.o: input_cntl.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

input_grid-a64fx.o: input_grid.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

mat_con0-a64fx.o: mat_con0.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

mat_con1-a64fx.o: mat_con1.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

mat_ass_main_1-a64fx.o: mat_ass_main_1.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

mat_ass_bc-a64fx.o: mat_ass_bc.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

hpcmw_finalize-a64fx.o: hpcmw_finalize.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

output-a64fx.o: output.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

partitioner-a64fx.o: partitioner.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

initial_mesh-a64fx.o: initial_mesh.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

part-a64fx.o: part.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

calc_edgcut-a64fx.o: calc_edgcut.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

cre_local_data-a64fx.o: cre_local_data.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

interface_nodes-a64fx.o: interface_nodes.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

metis-a64fx.o: metis.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

neib_pe-a64fx.o: neib_pe.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

paraset-a64fx.o: paraset.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

proc_local-a64fx.o: proc_local.f  
	$(FC2) $(FFLAGS2) -c -o $@ $<

local_data-a64fx.o: local_data.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

util-a64fx.o: util.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

rcm-a64fx.o: rcm.f
	$(FC2) $(FFLAGS2) -c -o $@ $<

# Makefile End
