ZOD = ../$(ZOLTAN_OBJ_DIR)
ZOM = $(F90_MODULE_PREFIX)$(ZOD)

zfdrive:	$(ZOD)/fdr_migrate.o \
		$(ZOD)/farg.o \
		$(ZOD)/mpistring.o \
		$(ZOD)/lampmpi.o \
		$(ZOD)/fdr_main.o \
		$(ZOD)/fdr_loadbal.o \
		$(ZOD)/fdr_input.o \
		$(ZOD)/fdr_chaco_io.o \
		$(ZOD)/fdr_const.o \
		$(ZOD)/libzoltan.a
	$(F90) -o 	zfdrive $(ZOD)/fdr_main.o \
				$(ZOD)/fdr_loadbal.o $(ZOD)/fdr_input.o \
        			$(ZOD)/fdr_chaco_io.o $(ZOD)/fdr_const.o \
				$(ZOD)/fdr_migrate.o $(ZOD)/farg.o \
				$(ZOD)/mpistring.o $(ZOD)/lampmpi.o \
				-L$(ZOD) -lzoltan $(LNK_LIBS) -lm
	/bin/mv zfdrive $(ZOD)

$(ZOD)/farg.o: $(FARG).f
	$(F90) -o $(ZOD)/farg.o -c $(FARG).f

$(ZOD)/fdr_main.o: 	fdr_main.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_input.o $(ZOD)/fdr_chaco_io.o \
			$(ZOD)/fdr_loadbal.o 
	$(F90) $(ZOM) -c fdr_main.f90
	/bin/mv *.o $(ZOD)

$(ZOD)/fdr_migrate.o: 	fdr_migrate.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_chaco_io.o
	$(F90) $(ZOM) -c fdr_migrate.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_loadbal.o: 	fdr_loadbal.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_migrate.o
	$(F90) $(ZOM) -c fdr_loadbal.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_input.o: 	fdr_input.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o 
	$(F90) $(ZOM) -c fdr_input.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_chaco_io.o: 	fdr_chaco_io.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_input.o 
	$(F90) $(ZOM) -c fdr_chaco_io.f90
	/bin/mv *.o *.mod $(ZOD)

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
	rm -f $(ZOD)/libzoltan.a $(ZOD)/fdr*.o $(ZOD)/farg.o \
	      $(ZOD)/mpi_h.o $(ZOD)/lampmpi.o $(ZOD)/mpistring.o\
	      $(ZOD)/*.mod debug V* *.vo *.dbg $(ZOD)/zfdrive
	rm -Rf album
