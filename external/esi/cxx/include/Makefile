#
# makefile to configure ESI headers for use with MPI
#
CP=/bin/cp
LS=/bin/ls
RM=/bin/rm -f

include ../../make.vars

default:

clean:
	rm -f *.x *.o *.a *~ Templates.DB ti_files \
	ii_files SunWS_cache *.*~ .*~ *.log;
	cd esi; make clean;
