all: cwrap.o ztypes.o lb_user_const.o fwrap.o

fwrap.o: fwrap.f90 ztypes.o lb_user_const.o
	$(F90) -c fwrap.f90

fwrap.f90: fwrap.fpp sppr sppr_header
	./sppr < fwrap.fpp > fwrap.f90

lb_user_const.o: lb_user_const.f90 ztypes.o
	$(F90) -c lb_user_const.f90

ztypes.o: ztypes.f90
	$(F90) -c ztypes.f90

ztypes.f90: ztypes.fpp sppr sppr_header
	./sppr < ztypes.fpp > ztypes.f90

sppr_header: $(SPPR_HEAD)
	cp $(SPPR_HEAD) sppr_header

cwrap.o: cwrap.c
	$(CC) -c cwrap.c

sppr: sppr.f90
	/usr/sbin/f90 -o sppr sppr.f90

clean:
	rm -f *.o fwrap.f90 ztypes.f90 *.mod debug V* *.vo *.dbg zfdriver sppr \
	      sppr_header fwrap.F90 ztypes.F90
	rm -Rf album
