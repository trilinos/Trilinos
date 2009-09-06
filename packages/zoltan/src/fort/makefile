##############################################################################
# Zoltan Library for Parallel Applications                                   #
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
			$(ZOD)/zoltan_user_data.o  $(ZOD)/lb_user_const.o \
			$(ZOD)/fwrap.o 

$(ZOD)/fwrap.o: 	$(ZOD)/fwrap.f90 \
			$(ZOD)/ztypes.o  \
			$(ZOD)/zoltan_user_data.o \
			$(ZOD)/lb_user_const.o 
	$(F90) -c $(ZOM) $(ZOD)/fwrap.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/fwrap.f90: 	fwrap.fpp lbif.fpp lbfn.fpp \
		set_numgeom.if			set_numgeom.fn \
		set_geommulti.if		set_geommulti.fn \
		set_geom.if			set_geom.fn \
		set_numedges.if 		set_numedges.fn \
		set_numedgesmulti.if 		set_numedgesmulti.fn \
		set_edgelist.if			set_edgelist.fn \
		set_edgelistmulti.if		set_edgelistmulti.fn \
		set_numobj.if 			set_numobj.fn \
		set_partition.if 		set_partition.fn \
		set_partitionmulti.if		set_partitionmulti.fn \
		set_objlist.if			set_objlist.fn \
		set_firstobj.if			set_firstobj.fn \
		set_nextobj.if			set_nextobj.fn \
		set_numborderobj.if		set_numborderobj.fn \
		set_borderobjlist.if		set_borderobjlist.fn \
		set_firstborderobj.if		set_firstborderobj.fn \
		set_nextborderobj.if		set_nextborderobj.fn \
		set_premigrate.if		set_premigrate.fn \
		set_midmigrate.if		set_midmigrate.fn \
		set_postmigrate.if		set_postmigrate.fn \
		set_objsize.if			set_objsize.fn \
		set_packobj.if			set_packobj.fn \
		set_unpackobj.if		set_unpackobj.fn \
		set_objsizemulti.if		set_objsizemulti.fn \
		set_packobjmulti.if		set_packobjmulti.fn \
		set_unpackobjmulti.if		set_unpackobjmulti.fn \
		set_numchild.if 		set_numchild.fn \
		set_childweight.if 		set_childweight.fn \
		set_childlist.if 		set_childlist.fn \
		set_numcoarseobj.if 		set_numcoarseobj.fn \
		set_coarseobjlist.if 		set_coarseobjlist.fn \
		set_firstcoarseobj.if		set_firstcoarseobj.fn \
		set_nextcoarseobj.if		set_nextcoarseobj.fn \
		set_numgeom.if.lbif		set_numgeom.fn.lbfn \
		set_geom.if.lbif		set_geom.fn.lbfn \
		set_numedges.if.lbif 		set_numedges.fn.lbfn \
		set_edgelist.if.lbif		set_edgelist.fn.lbfn \
		set_numobj.if.lbif 		set_numobj.fn.lbfn \
		set_objlist.if.lbif		set_objlist.fn.lbfn \
		set_firstobj.if.lbif		set_firstobj.fn.lbfn \
		set_nextobj.if.lbif		set_nextobj.fn.lbfn \
		set_numborderobj.if.lbif	set_numborderobj.fn.lbfn \
		set_borderobjlist.if.lbif	set_borderobjlist.fn.lbfn \
		set_firstborderobj.if.lbif	set_firstborderobj.fn.lbfn \
		set_nextborderobj.if.lbif	set_nextborderobj.fn.lbfn \
		set_premigrate.if.lbif		set_premigrate.fn.lbfn \
		set_midmigrate.if.lbif		set_midmigrate.fn.lbfn \
		set_postmigrate.if.lbif		set_postmigrate.fn.lbfn \
		set_objsize.if.lbif		set_objsize.fn.lbfn \
		set_packobj.if.lbif		set_packobj.fn.lbfn \
		set_unpackobj.if.lbif		set_unpackobj.fn.lbfn \
		set_numchild.if.lbif 		set_numchild.fn.lbfn \
		set_childweight.if.lbif 	set_childweight.fn.lbfn \
		set_childlist.if.lbif 		set_childlist.fn.lbfn \
		set_hgcs.if			set_hgcs.fn \
		set_hgsizecs.if			set_hgsizecs.fn \
		set_hgedgeweights.if		set_hgedgeweights.fn \
		set_hgsizeedgeweights.if	set_hgsizeedgeweights.fn \
                set_numfixedobj.if              set_numfixedobj.fn \
                set_fixedobjlist.if             set_fixedobjlist.fn \
		set_numcoarseobj.if.lbif 	set_numcoarseobj.fn.lbfn \
		set_coarseobjlist.if.lbif 	set_coarseobjlist.fn.lbfn \
		set_firstcoarseobj.if.lbif	set_firstcoarseobj.fn.lbfn \
		set_nextcoarseobj.if.lbif	set_nextcoarseobj.fn.lbfn \
		set_hiernumlevels.if		set_hiernumlevels.fn \
		set_hierpartition.if		set_hierpartition.fn \
		set_hiermethod.if		set_hiermethod.fn \
		$(ZOD)/sppr 
	$(ZOD)/sppr < fwrap.fpp > $(ZOD)/fwrap.f90

$(ZOD)/zoltan_user_data.o:	zoltan_user_data.f90 \
			$(ZOD)/ztypes.o
	$(F90) -c  $(ZOM) zoltan_user_data.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/lb_user_const.o:	lb_user_const.f90 \
			$(ZOD)/ztypes.o
	$(F90) -c  $(ZOM) lb_user_const.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/ztypes.o: 	$(ZOD)/ztypes.f90
	$(F90) -c  $(ZOM) $(ZOD)/ztypes.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/ztypes.f90:	ztypes.fpp $(ZOD)/sppr
	$(ZOD)/sppr < ztypes.fpp > $(ZOD)/ztypes.f90


$(ZOD)/cwrap.o: cwrap.c
	$(CC) -c cwrap.c
	/bin/mv *.o $(ZOD)

$(ZOD)/sppr: sppr.f90
	$(LOCAL_F90) -o $(ZOD)/sppr sppr.f90

clean:
	@rm -f $(ZOD)/*.f90 $(ZOD)/*.mod debug V* $(ZOD)/*.vo \
	      $(ZOD)/cwrap.o \
	      $(ZOD)/fwrap.o \
	      $(ZOD)/ztypes.o \
	      $(ZOD)/zoltan_user_data.o \
	      $(ZOD)/lb_user_const.o \
              sppr_header *.vo *.mod *.dbg \
	      $(ZOD)/zfdrive $(ZOD)/sppr $(ZOD)/sppr_header
	@rm -Rf album
