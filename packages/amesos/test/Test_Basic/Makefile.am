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
include $(top_builddir)/Makefile.export.amesos

if USING_GNUMAKE
EXAMPLE_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(AMESOS_LIBS) $(AMESOS_EXTRA_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES))
else
EXAMPLE_LIBS = $(AMESOS_LIBS) $(AMESOS_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES)
endif

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

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
	$(srcdir)/In_AmesosDscpack.csh Matrices $(srcdir)/In_AmesosUmfpack.csh \
	$(srcdir)/In_AmesosKlu.csh $(srcdir)/In_AmesosSuperlu.csh \
	 $(srcdir)/In_AmesosScalapack.csh \
	 $(srcdir)/In_AmesosLapack.csh \
	  $(srcdir)/In_AmesosTaucs.csh \
	  $(srcdir)/In_AmesosPardiso.csh \
	  $(srcdir)/In_AmesosPastix.csh \
	  $(srcdir)/In_AmesosParaklete.csh \
	  $(srcdir)/In_AmesosMumps.csh \
	$(srcdir)/In_AmesosKluSerial.csh  $(srcdir)/In_AmesosSuperluSerial.csh \
	$(srcdir)/In_AmesosUmfpackSerial.csh 
	cp  $(srcdir)/In_TestBasic.csh TestBasic.csh
	cp  $(srcdir)/In_AmesosSuperludist.csh AmesosSuperludist.csh 
	cp  $(srcdir)/In_AmesosDscpack.csh AmesosDscpack.csh 
	cp  $(srcdir)/In_AmesosUmfpack.csh AmesosUmfpack.csh 
	cp  $(srcdir)/In_AmesosKlu.csh AmesosKlu.csh 
	cp  $(srcdir)/In_AmesosSuperlu.csh AmesosSuperlu.csh 
	cp  $(srcdir)/In_AmesosScalapack.csh AmesosScalapack.csh 
	cp  $(srcdir)/In_AmesosLapack.csh AmesosLapack.csh 
	cp  $(srcdir)/In_AmesosTaucs.csh AmesosTaucs.csh 
	cp  $(srcdir)/In_AmesosPardiso.csh AmesosPardiso.csh 
	cp  $(srcdir)/In_AmesosPastix.csh AmesosPastix.csh 
	cp  $(amesos_paraklete_test) AmesosParaklete.csh 
	cp  $(srcdir)/In_AmesosMumps.csh AmesosMumps.csh 
	cp  $(srcdir)/In_AmesosKluSerial.csh AmesosKluSerial.csh 
	cp  $(srcdir)/In_AmesosSuperluSerial.csh AmesosSuperluSerial.csh 
	cp  $(srcdir)/In_AmesosUmfpackSerial.csh AmesosUmfpackSerial.csh 
	cp  $(srcdir)/In_Suppressions Suppressions 

Matrices: $(srcdir)/662_bus.rsa
	cp  $(srcdir)/662_bus.rsa  662_bus_out.rsa                              
	cp  $(srcdir)/In_ImpcolA.rua  ImpcolA.rua                               
	cp  $(srcdir)/In_ImpcolB.rua  ImpcolB.rua                               
	cp  $(srcdir)/In_ImpcolC.rua  ImpcolC.rua                               
	cp  $(srcdir)/In_ImpcolD.rua  ImpcolD.rua                               
	cp  $(srcdir)/In_ImpcolE.rua  ImpcolE.rua                               
	cp  $(srcdir)/In_fidapm05.rua  fidapm05.rua                             
	cp  $(srcdir)/In_Dense.triU  Dense.triU                               
	cp  $(srcdir)/In_SuperLU.triU  SuperLU.triU                               
	cp  $(srcdir)/In_SuperLU.rua  SuperLU.rua                               
	cp  $(srcdir)/In_bcsstk01.mtx bcsstk01.mtx                              
	cp  $(srcdir)/In_bcsstk02.mtx bcsstk02.mtx                              
	cp  $(srcdir)/In_bcsstk04.mtx bcsstk04.mtx                              
	cp  $(srcdir)/In_nos1.mtx nos1.mtx                                      
	cp  $(srcdir)/In_nos5.mtx nos5.mtx                                      
	cp  $(srcdir)/In_FourByFour.mtx FourByFour.mtx                          
	cp  $(srcdir)/In_StructurallySingular.mtx StructurallySingular.mtx   
	cp  $(srcdir)/In_NumericallySingular.mtx NumericallySingular.mtx  
	cp  $(srcdir)/In_Diagonal.mtx Diagonal.mtx                              
	cp  $(srcdir)/In_MissingADiagonal.mtx MissingADiagonal.mtx             
	cp  $(srcdir)/In_KheadK.mtx KheadK.mtx                                  
	cp  $(srcdir)/In_KheadSorted.mtx  KheadSorted.mtx
	cp  $(srcdir)/In_Khead.triS Khead.triS  
	cp  $(srcdir)/In_count_coverage count_coverage  

SRC_CXX =  $(srcdir)/SparseDirectTimingVars.cpp \
           $(srcdir)/SparseDirectTimingVars.h \
           $(srcdir)/SparseSolverResult.cpp \
           SparseSolverResult.h 

if HAVE_AMESOS_TRIUTILS
amesos_test_SOURCES = $(srcdir)/Amesos_TestDriver.cpp \
	 $(srcdir)/Amesos_TestSolver.cpp $(srcdir)/Amesos_TestSolver.h \
	 $(srcdir)/Amesos_TestMultiSolver.cpp  $(srcdir)/Amesos_TestMrhsSolver.cpp \
	 $(SRC_CXX)
else
amesos_test_SOURCES = $(srcdir)/DoNothing.cpp
endif

if HAVE_AMESOS_EPETRAEXT
amesos_paraklete_test = $(srcdir)/In_AmesosParaklete.csh
else
amesos_paraklete_test = $(srcdir)/In_AmesosParaklete_Sans_EpetraExt.csh
endif


TestBasic_DEPENDENCIES = $(top_builddir)/src/libamesos.a
TestBasic_LDADD        = ${EXAMPLE_LIBS}

amesos_test_DEPENDENCIES = $(top_builddir)/src/libamesos.a
amesos_test_LDADD        = ${EXAMPLE_LIBS}

EXTRA_DIST = $(srcdir)/662_bus.rsa \
	$(srcdir)/In_ImpcolA.rua \
        $(srcdir)/In_ImpcolB.rua \
        $(srcdir)/In_ImpcolC.rua \
        $(srcdir)/In_ImpcolD.rua \
        $(srcdir)/In_ImpcolE.rua \
        $(srcdir)/In_fidapm05.rua \
        $(srcdir)/In_SuperLU.rua \
        $(srcdir)/In_Dense.triU \
        $(srcdir)/In_SuperLU.triU \
        $(srcdir)/In_bcsstk01.mtx \
        $(srcdir)/In_bcsstk02.mtx \
        $(srcdir)/In_bcsstk04.mtx \
        $(srcdir)/In_nos1.mtx \
        $(srcdir)/In_nos5.mtx \
        $(srcdir)/In_FourByFour.mtx \
        $(srcdir)/In_StructurallySingular.mtx \
        $(srcdir)/In_NumericallySingular.mtx \
        $(srcdir)/In_Diagonal.mtx \
        $(srcdir)/In_MissingADiagonal.mtx \
        $(srcdir)/In_KheadK.mtx \
        $(srcdir)/In_Khead.triS \
        $(srcdir)/In_KheadSorted.mtx \
	$(srcdir)/In_AmesosSuperludist.csh  \
	$(srcdir)/In_AmesosDscpack.csh $(srcdir)/In_AmesosUmfpack.csh \
	$(srcdir)/In_AmesosKlu.csh $(srcdir)/In_AmesosSuperlu.csh \
	$(srcdir)/In_AmesosScalapack.csh \
	$(srcdir)/In_AmesosLapack.csh \
	$(srcdir)/In_AmesosTaucs.csh \
	$(srcdir)/In_AmesosPardiso.csh \
	$(srcdir)/In_AmesosPastix.csh \
	$(srcdir)/In_AmesosParaklete.csh \
	$(srcdir)/In_AmesosParaklete_Sans_EpetraExt.csh \
	$(srcdir)/In_AmesosMumps.csh \
	$(srcdir)/In_AmesosKluSerial.csh  $(srcdir)/In_AmesosSuperluSerial.csh \
	$(srcdir)/In_AmesosUmfpackSerial.csh \
	$(srcdir)/In_Suppressions \
	$(srcdir)/Amesos_TestDriver.cpp \
	$(srcdir)/Amesos_TestSolver.cpp \
	$(srcdir)/Amesos_TestMultiSolver.cpp \
	$(srcdir)/Amesos_TestMrhsSolver.cpp \
	$(srcdir)/SparseDirectTimingVars.cpp \
	$(srcdir)/SparseDirectTimingVars.h \
	$(srcdir)/Amesos_TestSolver.h \
	$(srcdir)/SparseSolverResult.cpp \
	$(srcdir)/SparseSolverResult.h \
	$(srcdir)/In_count_coverage 

clean:  clean-am
	rm -f  662_bus_out.rsa \
	ImpcolA.rua \
        ImpcolB.rua \
        ImpcolC.rua \
        ImpcolD.rua \
        ImpcolE.rua \
        fidapm05.rua \
        SuperLU.rua \
        Dense.triU \
        SuperLU.triU \
        bcsstk01.mtx \
        bcsstk02.mtx \
        bcsstk04.mtx \
        nos1.mtx \
        nos5.mtx \
        FourByFour.mtx \
        StructurallySingular.mtx \
        NumericallySingular.mtx \
        Diagonal.mtx \
        MissingADiagonal.mtx \
        KheadK.mtx \
        Khead.triS \
        KheadSorted.mtx \
	AmesosSuperludist.csh  \
	TestBasic.csh  \
	AmesosDscpack.csh AmesosUmfpack.csh \
	AmesosKlu.csh AmesosSuperlu.csh \
	AmesosScalapack.csh \
	AmesosLapack.csh \
	AmesosTaucs.csh \
	AmesosPardiso.csh \
	AmesosPastix.csh \
	AmesosParaklete.csh \
	AmesosMumps.csh \
	AmesosKluSerial.csh  AmesosSuperluSerial.csh \
	AmesosUmfpackSerial.csh \
	count_coverage \
	Suppressions \
	SST.* \
	AME.* 

all-libs:
	cd $(top_builddir)/src ; $(MAKE) all-libs 

