ZOD = ../$(ZOLTAN_OBJ_DIR)
ZOM = $(F90_MODULE_PREFIX)$(ZOD)

COMPILE = $(F90) -g
COMPILE_LIB = -lm

# To use purify with zfdrive (yes, it really works!), replace the COMPILE and 
# COMPILE_LIB macros with the two macros below and compile on paunchy.
# (FYI:  dbx told me which libraries are automatically linked using $(F90);
# these libraries must be explicitly specified when using cc.)
#COMPILE =  purify -best-effort -follow-child-processes=yes \
#                -cache-dir=/Net/local/proj/zoltan/tmp/purify \
#                -chain-length=12 /Net/local/bin/cc 
#COMPILE_LIB =   -lc -lm -lfui -lfai -lfai2 -lfsumai -lfprodai \
#	        -lfminlai -lfmaxlai -lfminvai -lfmaxvai \
#	        -lfsu -lsunmath -lmp -lm


zfdrive:	$(ZOD)/fdr_migrate.o \
		$(ZOD)/farg.o \
		$(ZOD)/mpistring.o \
		$(ZOD)/lampmpi.o \
		$(ZOD)/fdr_sort.o \
		$(ZOD)/fdr_main.o \
		$(ZOD)/fdr_loadbal.o \
		$(ZOD)/fdr_input.o \
		$(ZOD)/fdr_param_file.o \
		$(ZOD)/fdr_chaco_io.o \
		$(ZOD)/fdr_mm_io.o \
		$(ZOD)/fdr_const.o \
		$(ZOD)/mmio.o \
		$(ZOD)/libzoltan.a
		$(COMPILE) $(LDFLAGS) -o      zfdrive $(ZOD)/fdr_main.o \
			        $(ZOD)/fdr_loadbal.o $(ZOD)/fdr_input.o \
				$(ZOD)/fdr_sort.o $(ZOD)/fdr_param_file.o \
				$(ZOD)/fdr_mm_io.o $(ZOD)/mmio.o \
			        $(ZOD)/fdr_chaco_io.o $(ZOD)/fdr_const.o \
			        $(ZOD)/fdr_migrate.o $(ZOD)/farg.o \
			        $(ZOD)/mpistring.o $(ZOD)/lampmpi.o \
			        -L$(ZOD) -lzoltan $(LNK_LIBS) $(COMPILE_LIB)
		/bin/mv zfdrive $(ZOD)

$(ZOD)/farg.o: $(FARG).f
	$(F90) -o $(ZOD)/farg.o -c $(FARG).f

$(ZOD)/fdr_main.o: 	fdr_main.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_input.o $(ZOD)/fdr_chaco_io.o \
			$(ZOD)/fdr_mm_io.o $(ZOD)/mmio.o \
			$(ZOD)/fdr_loadbal.o $(ZOD)/fdr_sort.o
	$(F90) $(ZOM) -c fdr_main.f90
	/bin/mv *.o $(ZOD)

$(ZOD)/fdr_migrate.o: 	fdr_migrate.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_chaco_io.o \
			$(ZOD)/fdr_mm_io.o $(ZOD)/mmio.o 
	$(F90) $(ZOM) -c fdr_migrate.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_loadbal.o: 	fdr_loadbal.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_migrate.o $(ZOD)/fdr_input.o \
			$(ZOD)/fdr_param_file.o
	$(F90) $(ZOM) -c fdr_loadbal.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_sort.o: 	fdr_sort.f90 
	$(F90) $(ZOM) -c fdr_sort.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_input.o: 	fdr_input.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o 
	$(F90) $(ZOM) -c fdr_input.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_param_file.o: 	fdr_param_file.f90 $(ZOD)/mpi_h.o
	$(F90) $(ZOM) -c fdr_param_file.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_chaco_io.o: 	fdr_chaco_io.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_input.o 
	$(F90) $(ZOM) -c fdr_chaco_io.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_mm_io.o: 	fdr_mm_io.f90 $(ZOD)/mmio.o $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o $(ZOD)/fdr_sort.o \
			$(ZOD)/fdr_input.o 
	$(F90) $(ZOM) -c fdr_mm_io.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/mmio.o: 	mmio.f 
	$(F90) -o $(ZOD)/mmio.o -c mmio.f

$(ZOD)/fdr_const.o: 	fdr_const.f90 
	$(F90) $(ZOM) -c fdr_const.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/mpi_h.o: 	mpi_h.f
	$(F90) $(ZOLTAN_INC_PATH) -c mpi_h.f
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/lampmpi.o: 	lampmpi.f90
	$(F90) -c lampmpi.f90
	/bin/mv *.o $(ZOD)

$(ZOD)/mpistring.o: 	mpistring.c
	$(CC) -c mpistring.c
	/bin/mv *.o $(ZOD)

clean:
	@rm -f $(ZOD)/libzoltan.a $(ZOD)/fdr*.o $(ZOD)/farg.o \
	      $(ZOD)/mpi_h.o $(ZOD)/lampmpi.o $(ZOD)/mpistring.o\
	      $(ZOD)/*.mod debug V* *.vo *.dbg $(ZOD)/zfdrive
	@rm -Rf album
