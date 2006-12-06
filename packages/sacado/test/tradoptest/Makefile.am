## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
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
# Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
# (etphipp@sandia.gov).
# 
# ************************************************************************
# @HEADER

EXEEXT = .exe

noinst_PROGRAMS = \
	tradoptesta1 \
	tradoptestb1 \
	tradoptestc1 \
	tradoptesta2 \
	tradoptestb2 \
	tradoptestc2 \
	tradoptesta3 \
	tradoptestb3 \
	tradoptestc3 \
	tradoptesta4 \
	tradoptestb4 \
	tradoptestc4 \
	tradoptesta5 \
	tradoptestb5 \
	tradoptestc5 \
	tradoptesta6 \
	tradoptestb6 \
	tradoptestc6

testall : $(noinst_PROGRAMS)
	./tradoptesta1$(EXEEXT)
	./tradoptesta2$(EXEEXT)
	./tradoptesta3$(EXEEXT)	
	./tradoptesta4$(EXEEXT)
	./tradoptesta5$(EXEEXT)
	./tradoptesta6$(EXEEXT)
	./tradoptestb1$(EXEEXT)
	./tradoptestb2$(EXEEXT)
	./tradoptestb3$(EXEEXT)	
	./tradoptestb4$(EXEEXT)
	./tradoptestb5$(EXEEXT)
	./tradoptestb6$(EXEEXT)
	./tradoptestc1$(EXEEXT)
	./tradoptestc2$(EXEEXT)
	./tradoptestc3$(EXEEXT)	
	./tradoptestc4$(EXEEXT)
	./tradoptestc5$(EXEEXT)
	./tradoptestc6$(EXEEXT)

tradoptesta1_SOURCES = tradoptesta.cpp 
tradoptesta1_CXXFLAGS = -USACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptesta1_DEPENDENCIES = $(top_builddir)/src/libsacado.a
tradoptesta1_LDADD = -L$(top_builddir)/src -lsacado $(LIBS) $(FLIBS)

tradoptesta2_SOURCES = tradoptesta.cpp 
tradoptesta2_CXXFLAGS = -USACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptesta2_DEPENDENCIES = $(tradoptesta1_DEPENDENCIES)
tradoptesta2_LDADD = $(tradoptesta1_LDADD)

tradoptesta3_SOURCES = tradoptesta.cpp 
tradoptesta3_CXXFLAGS = -USACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptesta3_DEPENDENCIES = $(tradoptesta1_DEPENDENCIES)
tradoptesta3_LDADD = $(tradoptesta1_LDADD)

tradoptesta4_SOURCES = tradoptesta.cpp 
tradoptesta4_CXXFLAGS = -DSACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptesta4_DEPENDENCIES = $(tradoptesta1_DEPENDENCIES)
tradoptesta4_LDADD = $(tradoptesta1_LDADD)

tradoptesta5_SOURCES = tradoptesta.cpp 
tradoptesta5_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptesta5_DEPENDENCIES = $(tradoptesta1_DEPENDENCIES)
tradoptesta5_LDADD = $(tradoptesta1_LDADD)

tradoptesta6_SOURCES = tradoptesta.cpp 
tradoptesta6_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptesta6_DEPENDENCIES = $(tradoptesta1_DEPENDENCIES)
tradoptesta6_LDADD = $(tradoptesta1_LDADD)


tradoptestb1_SOURCES = tradoptestb.cpp 
tradoptestb1_CXXFLAGS = -USACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestb1_DEPENDENCIES = $(top_builddir)/src/libsacado.a
tradoptestb1_LDADD = -L$(top_builddir)/src -lsacado $(LIBS) $(FLIBS)

tradoptestb2_SOURCES = tradoptestb.cpp 
tradoptestb2_CXXFLAGS = -USACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestb2_DEPENDENCIES = $(tradoptestb1_DEPENDENCIES)
tradoptestb2_LDADD = $(tradoptestb1_LDADD)

tradoptestb3_SOURCES = tradoptestb.cpp 
tradoptestb3_CXXFLAGS = -USACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestb3_DEPENDENCIES = $(tradoptestb1_DEPENDENCIES)
tradoptestb3_LDADD = $(tradoptestb1_LDADD)

tradoptestb4_SOURCES = tradoptestb.cpp 
tradoptestb4_CXXFLAGS = -DSACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestb4_DEPENDENCIES = $(tradoptestb1_DEPENDENCIES)
tradoptestb4_LDADD = $(tradoptestb1_LDADD)

tradoptestb5_SOURCES = tradoptestb.cpp 
tradoptestb5_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestb5_DEPENDENCIES = $(tradoptestb1_DEPENDENCIES)
tradoptestb5_LDADD = $(tradoptestb1_LDADD)

tradoptestb6_SOURCES = tradoptestb.cpp 
tradoptestb6_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestb6_DEPENDENCIES = $(tradoptestb1_DEPENDENCIES)
tradoptestb6_LDADD = $(tradoptestb1_LDADD)


tradoptestc1_SOURCES = tradoptestc.cpp 
tradoptestc1_CXXFLAGS = -USACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestc1_DEPENDENCIES = $(top_builddir)/src/libsacado.a
tradoptestc1_LDADD = -L$(top_builddir)/src -lsacado $(LIBS) $(FLIBS)

tradoptestc2_SOURCES = tradoptestc.cpp 
tradoptestc2_CXXFLAGS = -USACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestc2_DEPENDENCIES = $(tradoptestc1_DEPENDENCIES)
tradoptestc2_LDADD = $(tradoptestc1_LDADD)

tradoptestc3_SOURCES = tradoptestc.cpp 
tradoptestc3_CXXFLAGS = -USACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestc3_DEPENDENCIES = $(tradoptestc1_DEPENDENCIES)
tradoptestc3_LDADD = $(tradoptestc1_LDADD)

tradoptestc4_SOURCES = tradoptestc.cpp 
tradoptestc4_CXXFLAGS = -DSACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestc4_DEPENDENCIES = $(tradoptestc1_DEPENDENCIES)
tradoptestc4_LDADD = $(tradoptestc1_LDADD)

tradoptestc5_SOURCES = tradoptestc.cpp 
tradoptestc5_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestc5_DEPENDENCIES = $(tradoptestc1_DEPENDENCIES)
tradoptestc5_LDADD = $(tradoptestc1_LDADD)

tradoptestc6_SOURCES = tradoptestc.cpp 
tradoptestc6_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestc6_DEPENDENCIES = $(tradoptestc1_DEPENDENCIES)
tradoptestc6_LDADD = $(tradoptestc1_LDADD)

tradoptesta.cpp : testgen.as tradoptesta.c0
	$(SHELL) testgen.as $(srcdir)/tradoptesta.c0 > tradoptesta.cpp

tradoptestb.cpp : testgen.as tradoptestb.c0
	$(SHELL) testgen.as $(srcdir)/tradoptestb.c0 > tradoptestb.cpp

tradoptestc.cpp : testgen.as tradoptestc.c0
	$(SHELL) testgen.as $(srcdir)/tradoptestc.c0 > tradoptestc.cpp

CLEANFILES = tradoptesta.cpp tradoptestb.cpp tradoptestc.cpp
