zfdrive: dr_migrate.o farg.o mpistring.o lampmpi.o \
	dr_main.o dr_loadbal.o dr_input.o dr_chaco_io.o dr_const.o \
	$(ZOLTAN_OBJ_DIR)/libzoltan.a
	$(F90) $(FFLAGS) -o zfdrive dr_main.o dr_loadbal.o dr_input.o \
        dr_chaco_io.o dr_const.o dr_migrate.o farg.o mpistring.o lampmpi.o \
	-L$(ZOLTAN_OBJ_DIR) -lzoltan $(MPILIB) -lm

farg.o: $(FARG).f
	$(F90) $(FFLAGS) -o farg.o -c $(FARG).f

dr_main.o: dr_main.f90 mpi_h.o dr_const.o dr_input.o \
	dr_chaco_io.o dr_loadbal.o $(ZOLTAN_OBJ_DIR)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c dr_main.f90

dr_migrate.o: dr_migrate.f90 mpi_h.o dr_const.o \
	dr_chaco_io.o $(ZOLTAN_OBJ_DIR)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c dr_migrate.f90

dr_loadbal.o: dr_loadbal.f90 mpi_h.o dr_const.o \
	dr_migrate.o $(ZOLTAN_OBJ_DIR)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c dr_loadbal.f90

dr_input.o: dr_input.f90 mpi_h.o dr_const.o $(ZOLTAN_OBJ_DIR)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c dr_input.f90

dr_chaco_io.o: dr_chaco_io.f90 mpi_h.o dr_const.o dr_input.o \
	$(ZOLTAN_OBJ_DIR)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c dr_chaco_io.f90

dr_const.o: dr_const.f90 $(ZOLTAN_OBJ_DIR)/libzoltan.a
	$(F90) $(FFLAGS) $(ZOLTAN_MOD_DIR) -c dr_const.f90

mpi_h.o: mpi_h.f
	$(F90) $(FFLAGS) $(MPIINC) -c mpi_h.f

lampmpi.o: lampmpi.f90
	$(F90) $(FFLAGS) -c lampmpi.f90

mpistring.o: mpistring.c
	$(CC) -c mpistring.c

$(ZOLTAN_OBJ_DIR)/libzoltan.a:
	cd ..; $(MAKE) $(TARGET_ARCH)

clean:
	rm -f *.o *.mod debug V* *.vo *.dbg zfdrive
	rm -Rf album
