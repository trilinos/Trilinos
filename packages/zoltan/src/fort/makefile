# @HEADER
#
########################################################################
#
#  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
#                  Copyright 2012 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Karen Devine	kddevin@sandia.gov
#                    Erik Boman	        egboman@sandia.gov
#
########################################################################
#
# @HEADER

ZOD = ../$(ZOLTAN_OBJ_DIR)
ZOM = $(F90_MODULE_PREFIX)$(ZOD) -I../fort

all: 			$(ZOD)/cwrap.o $(ZOD)/ztypes.o \
			$(ZOD)/zoltan_user_data.o  \
			$(ZOD)/fwrap.o 

$(ZOD)/fwrap.o: 	$(ZOD)/ztypes.o  \
			$(ZOD)/zoltan_user_data.o \
	 	fwrap.f90 \
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
		set_hgcs.if			set_hgcs.fn \
		set_hgsizecs.if			set_hgsizecs.fn \
		set_hgedgeweights.if		set_hgedgeweights.fn \
		set_hgsizeedgeweights.if	set_hgsizeedgeweights.fn \
                set_numfixedobj.if              set_numfixedobj.fn \
                set_fixedobjlist.if             set_fixedobjlist.fn \
		set_hiernumlevels.if		set_hiernumlevels.fn \
		set_hierpartition.if		set_hierpartition.fn \
		set_hiermethod.if		set_hiermethod.fn 
	$(F90) -c $(ZOM) fwrap.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/zoltan_user_data.o:	zoltan_user_data.f90 \
			$(ZOD)/ztypes.o
	$(F90) -c  $(ZOM) zoltan_user_data.f90
	/bin/mv *.o *.mod $(ZOD)

$(ZOD)/ztypes.o: 	ztypes.f90
	$(F90) -c  $(ZOM) ztypes.f90
	/bin/mv *.o *.mod $(ZOD)


$(ZOD)/cwrap.o: cwrap.c
	$(CC) -c cwrap.c
	/bin/mv *.o $(ZOD)


clean:
	@rm -f $(ZOD)/*.f90 $(ZOD)/*.mod debug V* $(ZOD)/*.vo \
	      $(ZOD)/cwrap.o \
	      $(ZOD)/fwrap.o \
	      $(ZOD)/ztypes.o \
	      $(ZOD)/zoltan_user_data.o \
              *.vo *.mod *.dbg \
	      $(ZOD)/zfdrive 
	@rm -Rf album
