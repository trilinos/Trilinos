ZOD = $(ZOLTAN_OBJ_DIR)

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
	$(F90) $(FFLAGS) -o 	zfdrive $(ZOD)/fdr_main.o \
				$(ZOD)/fdr_loadbal.o $(ZOD)/fdr_input.o \
        			$(ZOD)/fdr_chaco_io.o $(ZOD)/fdr_const.o \
				$(ZOD)/fdr_migrate.o $(ZOD)/farg.o \
				$(ZOD)/mpistring.o $(ZOD)/lampmpi.o \
				-L$(ZOD) -lzoltan $(MPILIB) -lm
	/bin/mv zfdrive $(ZOD)

$(ZOD)/farg.o: $(FARG).f
	$(F90) $(FFLAGS) -o $(ZOD)/farg.o -c $(FARG).f

$(ZOD)/fdr_main.o: 	fdr_main.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_input.o $(ZOD)/fdr_chaco_io.o \
			$(ZOD)/fdr_loadbal.o $(ZOD)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c fdr_main.f90
	/bin/mv *.o $(ZOD)

$(ZOD)/fdr_migrate.o: 	fdr_migrate.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_chaco_io.o $(ZOD)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c fdr_migrate.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_loadbal.o: 	fdr_loadbal.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_migrate.o $(ZOD)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c fdr_loadbal.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_input.o: 	fdr_input.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c fdr_input.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_chaco_io.o: 	fdr_chaco_io.f90 $(ZOD)/mpi_h.o $(ZOD)/fdr_const.o \
			$(ZOD)/fdr_input.o $(ZOD)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c fdr_chaco_io.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fdr_const.o: 	fdr_const.f90 $(ZOD)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c fdr_const.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/mpi_h.o: 	mpi_h.f
	$(F90) $(FFLAGS) $(MPIINC) -c mpi_h.f
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/lampmpi.o: 	lampmpi.f90
	$(F90) $(FFLAGS) -c lampmpi.f90
	/bin/mv *.o $(ZOD)

$(ZOD)/mpistring.o: 	mpistring.c
	$(CC) -c mpistring.c
	/bin/mv *.o $(ZOD)

$(ZOD)/libzoltan.a:
	cd ..; $(MAKE) $(TARGET_ARCH)

clean:
	rm -f $(ZOD)/fdr*.o $(ZOD)/mpi_h.o $(ZOD)/lampmpi.o $(ZOD)/mpistring.o\
	      *.mod debug V* *.vo *.dbg $(ZOD)/zfdrive
	rm -Rf album
