##############################################################################
# Zoltan Dynamic Load-Balancing Library for Parallel Applications            #
# Copyright (c) 2000, Sandia National Laboratories.                          #
# This software is distributed under the GNU Lesser General Public License.  #
# For more info, see the README file in the top-level Zoltan directory.      # 
##############################################################################
##############################################################################
# CVS File Information
#    $RCSfile$
#    $Author$
#    $Date$
#    $Revision$
##############################################################################

ZOD = ../$(ZOLTAN_OBJ_DIR)
ZOM = $(F90_MODULE_PREFIX)$(ZOD)

all: 			$(ZOD)/cwrap.o $(ZOD)/ztypes.o \
			$(ZOD)/lb_user_const.o $(ZOD)/fwrap.o

$(ZOD)/fwrap.o: 	$(ZOD)/fwrap.f90 \
			$(ZOD)/ztypes.o  \
			$(ZOD)/lb_user_const.o
	$(F90) -c $(ZOM) $(ZOD)/fwrap.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fwrap.f90: 	fwrap.fpp \
			set_numgeom.if		set_numgeom.fn \
			set_geom.if		set_geom.fn \
			set_numedges.if 	set_numedges.fn \
			set_edgelist.if		set_edgelist.fn \
			set_numobj.if 		set_numobj.fn \
			set_objlist.if		set_objlist.fn \
			set_firstobj.if		set_firstobj.fn \
			set_nextobj.if		set_nextobj.fn \
			set_numborderobj.if	set_numborderobj.fn \
			set_borderobjlist.if	set_borderobjlist.fn \
			set_firstborderobj.if	set_firstborderobj.fn \
			set_nextborderobj.if	set_nextborderobj.fn \
			set_premigrate.if	set_premigrate.fn \
			set_midmigrate.if	set_midmigrate.fn \
			set_postmigrate.if	set_postmigrate.fn \
			set_objsize.if		set_objsize.fn \
			set_packobj.if		set_packobj.fn \
			set_unpackobj.if	set_unpackobj.fn \
			set_numchild.if 	set_numchild.fn \
			set_childweight.if 	set_childweight.fn \
			set_childlist.if 	set_childlist.fn \
			set_numcoarseobj.if 	set_numcoarseobj.fn \
			set_coarseobjlist.if 	set_coarseobjlist.fn \
			set_firstcoarseobj.if	set_firstcoarseobj.fn \
			set_nextcoarseobj.if	set_nextcoarseobj.fn \
			$(ZOD)/sppr $(ZOD)/sppr_header
	/bin/cp $(ZOD)/sppr_header sppr_header
	$(ZOD)/sppr < fwrap.fpp > $(ZOD)/fwrap.f90

$(ZOD)/lb_user_const.o:	lb_user_const.f90 \
			$(ZOD)/ztypes.o
	$(F90) -c  $(ZOM) lb_user_const.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/ztypes.o: 	$(ZOD)/ztypes.f90
	$(F90) -c  $(ZOM) $(ZOD)/ztypes.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/ztypes.f90:	ztypes.fpp $(ZOD)/sppr $(ZOD)/sppr_header
	/bin/cp $(ZOD)/sppr_header sppr_header
	$(ZOD)/sppr < ztypes.fpp > $(ZOD)/ztypes.f90

$(ZOD)/sppr_header: $(SPPR_HEAD)
	cp $(SPPR_HEAD) $(ZOD)/sppr_header

$(ZOD)/cwrap.o: cwrap.c $(ZOD)/sppr_header
	/bin/cp $(ZOD)/sppr_header sppr_header
	$(CC) -c cwrap.c
	/bin/mv *.o $(ZOD)

$(ZOD)/sppr: sppr.f90
	$(LOCAL_F90) -o $(ZOD)/sppr sppr.f90

clean:
	rm -f $(ZOD)/*.f90 $(ZOD)/*.mod debug V* $(ZOD)/*.vo \
	      $(ZOD)/cwrap.o \
	      $(ZOD)/fwrap.o \
	      $(ZOD)/ztypes.o \
	      $(ZOD)/lb_user_const.o \
              sppr_header *.vo *.mod *.dbg \
	      $(ZOD)/zfdrive $(ZOD)/sppr $(ZOD)/sppr_header
	rm -Rf album
