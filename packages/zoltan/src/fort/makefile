ZOD = ../$(ZOLTAN_OBJ_DIR)
ZOM = $(F90_MODULE_PREFIX)$(ZOD)

all: 			$(ZOD)/cwrap.o $(ZOD)/ztypes.o \
			$(ZOD)/lb_user_const.o $(ZOD)/fwrap.o

$(ZOD)/fwrap.o: 	$(ZOD)/fwrap.f90 \
			$(ZOD)/ztypes.o  \
			$(ZOD)/lb_user_const.o
	$(F90) -c $(ZOM) $(ZOD)/fwrap.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fwrap.f90: 	fwrap.fpp $(ZOD)/sppr $(ZOD)/sppr_header
	$(ZOD)/sppr < fwrap.fpp > $(ZOD)/fwrap.f90

$(ZOD)/lb_user_const.o:	lb_user_const.f90 \
			$(ZOD)/ztypes.o
	$(F90) -c  $(ZOM) lb_user_const.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/ztypes.o: 	$(ZOD)/ztypes.f90
	$(F90) -c  $(ZOM) $(ZOD)/ztypes.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/ztypes.f90:	ztypes.fpp $(ZOD)/sppr $(ZOD)/sppr_header
	$(ZOD)/sppr < ztypes.fpp > $(ZOD)/ztypes.f90

$(ZOD)/sppr_header: $(SPPR_HEAD)
	cp $(SPPR_HEAD) $(ZOD)/sppr_header

$(ZOD)/cwrap.o: cwrap.c
	$(CC) -c cwrap.c
	/bin/mv *.o $(ZOD)

$(ZOD)/sppr: sppr.f90
	$(LOCAL_F90) -o $(ZOD)/sppr sppr.f90
	/bin/mv *.mod $(ZOD)

clean:
	rm -f $(ZOD)/*.f90 $(ZOD)/*.mod debug V* $(ZOD)/*.vo \
	      $(ZOD)/cwrap.o \
	      $(ZOD)/fwrap.o \
	      $(ZOD)/ztypes.o \
	      $(ZOD)/lb_user_const.o \
              *.dbg $(ZOD)/zfdrive $(ZOD)/sppr $(ZOD)/sppr_header
	rm -Rf album
