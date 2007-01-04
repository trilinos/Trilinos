#@HEADER
# ************************************************************************
# 
#                 IFPACK
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

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package ifpack subdirectory example built successfully."
	@echo ""

include $(top_builddir)/Makefile.export

if HAVE_IFPACK_GALERI
include $(top_builddir)/../galeri/Makefile.export.galeri
endif

if HAVE_IFPACK_GALERI
EXAMPLE_LIBS     = $(IFPACK_EXTRA_LIBS)     $(IFPACK_LIBS)      $(GALERI_LIBS)
EXAMPLE_INCLUDES = $(IFPACK_EXTRA_INCLUDES) $(IFPACK_INCLUDES)  $(GALERI_INCLUDES)
else
EXAMPLE_LIBS     = $(IFPACK_EXTRA_LIBS)     $(IFPACK_LIBS)      
EXAMPLE_INCLUDES = $(IFPACK_EXTRA_INCLUDES) $(IFPACK_INCLUDES) 
endif 

EXEEXT = .exe

EXAMPLE_CORE =	Ifpack_ex_Reordering           \
		Ifpack_ex_Filtering           

if HAVE_IFPACK_AZTECOO
if HAVE_IFPACK_GALERI
EXAMPLE_EXT = Ifpack_ex_Factory              \
	      Ifpack_ex_ICT
if HAVE_IFPACK_AMESOS
EXAMPLE_EXT2 = Ifpack_ex_Amesos               \
	       Ifpack_ex_BlockRelaxation
endif
else
EXAMPLE_EXT =
EXAMPLE_EXT2 =
endif
endif

noinst_PROGRAMS = $(EXAMPLE_CORE) $(EXAMPLE_EXT) $(EXAMPLE_EXT2)

Ifpack_ex_Factory_SOURCES      = $(srcdir)/Ifpack_ex_Factory.cpp
Ifpack_ex_Factory_DEPENDENCIES = $(top_builddir)/src/libifpack.a 
Ifpack_ex_Factory_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Ifpack_ex_Factory_LDADD        = $(EXAMPLE_LIBS) 

Ifpack_ex_Amesos_SOURCES      = $(srcdir)/Ifpack_ex_Amesos.cpp
Ifpack_ex_Amesos_DEPENDENCIES = $(top_builddir)/src/libifpack.a 
Ifpack_ex_Amesos_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Ifpack_ex_Amesos_LDADD        = $(EXAMPLE_LIBS) 

Ifpack_ex_ICT_SOURCES      = $(srcdir)/Ifpack_ex_ICT.cpp
Ifpack_ex_ICT_DEPENDENCIES = $(top_builddir)/src/libifpack.a 
Ifpack_ex_ICT_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Ifpack_ex_ICT_LDADD        = $(EXAMPLE_LIBS)

Ifpack_ex_BlockRelaxation_SOURCES      = $(srcdir)/Ifpack_ex_BlockRelaxation.cpp
Ifpack_ex_BlockRelaxation_DEPENDENCIES = $(top_builddir)/src/libifpack.a 
Ifpack_ex_BlockRelaxation_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Ifpack_ex_BlockRelaxation_LDADD        = $(EXAMPLE_LIBS)

Ifpack_ex_Reordering_SOURCES      = $(srcdir)/Ifpack_ex_Reordering.cpp
Ifpack_ex_Reordering_DEPENDENCIES = $(top_builddir)/src/libifpack.a 
Ifpack_ex_Reordering_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Ifpack_ex_Reordering_LDADD        = $(EXAMPLE_LIBS) 

Ifpack_ex_Filtering_SOURCES      = $(srcdir)/Ifpack_ex_Filtering.cpp
Ifpack_ex_Filtering_DEPENDENCIES = $(top_builddir)/src/libifpack.a 
Ifpack_ex_Filtering_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Ifpack_ex_Filtering_LDADD        = $(EXAMPLE_LIBS)
