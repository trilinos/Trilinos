#@HEADER
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

EXEEXT = .exe

#
#  All of the AmeosShellScript and TestAmesosScript stuff is designed to 
#  copy some files form the source tree to the build tree
#
#  The idea behind the mechanism is to fool autotools into thinking that
#  AmesosShellScript is an executable and that the way that you create
#  that executable is the "cp" command.
#
#  The files are renamed as they are copied so that configure can work 
#  within the source directory (even though we recommend that users 
#  run configure from a separate build directory.
#

noinst_PROGRAMS = amesos_test \
                  TestBasic

TestBasic_SOURCES = $(srcdir)/In_TestBasic.csh

TestBasic$(EXEEXT): $(srcdir)/In_TestBasic.csh  $(srcdir)/In_AmesosSuperludist.csh  \
	$(srcdir)/In_AmesosDscpack.csh Matrices
	cp  $(srcdir)/In_TestBasic.csh TestBasic.csh
	cp  $(srcdir)/In_AmesosSuperludist.csh AmesosSuperludist.csh 
	cp  $(srcdir)/In_AmesosDscpack.csh AmesosDscpack.csh 
	cp  $(srcdir)/In_AmesosSuperlu.csh AmesosSuperlu.csh 
	cp  $(srcdir)/In_AmesosScalapack.csh AmesosScalapack.csh 
	cp  $(srcdir)/In_AmesosMumps.csh AmesosMumps.csh 

Matrices: $(srcdir)/662_bus.rsa
	cp  $(srcdir)/662_bus.rsa  662_bus_out.rsa                              
	cp  $(srcdir)/662_bus.rsa  662_bus_out.rsa                              
	cp  $(srcdir)/In_ImpcolA.rua  ImpcolA.rua                               
	cp  $(srcdir)/In_ImpcolA.rua  ImpcolA.rua                               
	cp  $(srcdir)/In_ImpcolB.rua  ImpcolB.rua                               
	cp  $(srcdir)/In_ImpcolB.rua  ImpcolB.rua                               
	cp  $(srcdir)/In_ImpcolC.rua  ImpcolC.rua                               
	cp  $(srcdir)/In_ImpcolC.rua  ImpcolC.rua                               
	cp  $(srcdir)/In_ImpcolD.rua  ImpcolD.rua                               
	cp  $(srcdir)/In_ImpcolD.rua  ImpcolD.rua                               
	cp  $(srcdir)/In_ImpcolE.rua  ImpcolE.rua                               
	cp  $(srcdir)/In_ImpcolE.rua  ImpcolE.rua                               
	cp  $(srcdir)/In_fidapm05.rua  fidapm05.rua                             
	cp  $(srcdir)/In_fidapm05.rua  fidapm05.rua                             
	cp  $(srcdir)/In_SuperLU.triU  SuperLU.triU                               
	cp  $(srcdir)/In_SuperLU.rua  SuperLU.rua                               
	cp  $(srcdir)/In_SuperLU.rua  SuperLU.rua                               
	cp  $(srcdir)/In_bcsstk01.mtx bcsstk01.mtx                              
	cp  $(srcdir)/In_bcsstk02.mtx bcsstk02.mtx                              
	cp  $(srcdir)/In_bcsstk04.mtx bcsstk04.mtx                              
	cp  $(srcdir)/In_nos1.mtx nos1.mtx                                      
	cp  $(srcdir)/In_nos1.mtx nos1.mtx                                      
	cp  $(srcdir)/In_nos5.mtx nos5.mtx                                      
	cp  $(srcdir)/In_nos5.mtx nos5.mtx                                      
	cp  $(srcdir)/In_FourByFour.mtx FourByFour.mtx                          
	cp  $(srcdir)/In_FourByFour.mtx FourByFour.mtx                          
	cp  $(srcdir)/In_Diagonal.mtx Diagonal.mtx                              
	cp  $(srcdir)/In_Diagonal.mtx Diagonal.mtx                              
	cp  $(srcdir)/In_KheadK.mtx KheadK.mtx                                  
	cp  $(srcdir)/In_KheadK.mtx KheadK.mtx                                  

SRC_CXX =  $(srcdir)/SparseDirectTimingVars.cpp \
           $(srcdir)/SparseDirectTimingVars.h \
           $(srcdir)/SparseSolverResult.cpp \
           SparseSolverResult.h 

amesos_test_SOURCES = $(srcdir)/Amesos_TestDriver.cpp \
	 $(srcdir)/Amesos_TestSolver.cpp $(srcdir)/Amesos_TestSolver.h \
	 $(srcdir)/Amesos_TestMultiSolver.cpp  $(srcdir)/Amesos_TestMrhsSolver.cpp \
	 $(SRC_CXX)

#Specify locations of include files within Trilinos, but outside amesos
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src \
	-I$(top_builddir)/../epetra/src \
	-I$(top_srcdir)/../teuchos/src -I$(top_builddir)/../teuchos/src \
	-I$(top_srcdir)/../epetraext/src/transform -I$(top_builddir)/../epetraext/src/transform \
	-I$(top_srcdir)/../epetraext/src -I$(top_builddir)/../epetraext/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../triutils/src 

amesos_test_DEPENDENCIES = \
	$(top_builddir)/../epetraext/src/libepetraext.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libamesos.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a 

amesos_test_CXXFLAGS = -I$(top_srcdir)/src

#
#  libamesos.a is loaded twice.  The second time allows umfpack to pick up AMD.
#
amesos_test_LDADD = $(top_builddir)/src/libamesos.a \
	$(top_builddir)/../epetraext/src/libepetraext.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a \
	$(top_builddir)/../triutils/src/libtriutils.a \
	$(AMESOS_LIBS)  $(top_builddir)/src/libamesos.a \
	$(LAPACK_LIBS) \
	$(BLAS_LIBS) $(FLIBS) $(FLIBS)

EXTRA_DIST = $(srcdir)/662_bus.rsa \
	$(srcdir)/In_ImpcolA.rua \
        $(srcdir)/In_ImpcolB.rua \
        $(srcdir)/In_ImpcolC.rua \
        $(srcdir)/In_ImpcolD.rua \
        $(srcdir)/In_ImpcolE.rua \
        $(srcdir)/In_fidapm05.rua \
        $(srcdir)/In_SuperLU.rua \
        $(srcdir)/In_bcsstk01.mtx \
        $(srcdir)/In_bcsstk02.mtx \
        $(srcdir)/In_bcsstk04.mtx \
        $(srcdir)/In_nos1.mtx \
        $(srcdir)/In_nos5.mtx \
        $(srcdir)/In_FourByFour.mtx \
        $(srcdir)/In_Diagonal.mtx \
        $(srcdir)/In_KheadK.mtx \
        $(srcdir)/In_KheadSorted.mtx
