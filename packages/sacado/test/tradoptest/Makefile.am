## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
#                 Copyright (2006) Sandia Corporation
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

# Broke the original tradoptest.c0 into 5 pieces, tradoptest[a-e],
# in hopes that small enough pieces would compile on beowulf.

EXEEXT = .exe

noinst_PROGRAMS = \
	tradoptesta1 \
	tradoptesta2 \
	tradoptesta3 \
	tradoptesta4 \
	tradoptesta5 \
	tradoptesta6 \
	tradoptestb1 \
	tradoptestb2 \
	tradoptestb3 \
	tradoptestb4 \
	tradoptestb5 \
	tradoptestb6 \
	tradoptestc1 \
	tradoptestc2 \
	tradoptestc3 \
	tradoptestc4 \
	tradoptestc5 \
	tradoptestc6 \
	tradoptestd1 \
	tradoptestd2 \
	tradoptestd3 \
	tradoptestd4 \
	tradoptestd5 \
	tradoptestd6 \
	tradopteste1 \
	tradopteste2 \
	tradopteste3 \
	tradopteste4 \
	tradopteste5 \
	tradopteste6 \
	tradoptestf1 \
	tradoptestf2 \
	tradoptestf3 \
	tradoptestf4 \
	tradoptestf5 \
	tradoptestf6

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
	./tradoptestd1$(EXEEXT)
	./tradoptestd2$(EXEEXT)
	./tradoptestd3$(EXEEXT)	
	./tradoptestd4$(EXEEXT)
	./tradoptestd5$(EXEEXT)
	./tradoptestd6$(EXEEXT)
	./tradopteste1$(EXEEXT)
	./tradopteste2$(EXEEXT)
	./tradopteste3$(EXEEXT)	
	./tradopteste4$(EXEEXT)
	./tradopteste5$(EXEEXT)
	./tradopteste6$(EXEEXT)
	./tradoptestf1$(EXEEXT)
	./tradoptestf2$(EXEEXT)
	./tradoptestf3$(EXEEXT)	
	./tradoptestf4$(EXEEXT)
	./tradoptestf5$(EXEEXT)
	./tradoptestf6$(EXEEXT)

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


tradoptestd1_SOURCES = tradoptestd.cpp 
tradoptestd1_CXXFLAGS = -USACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestd1_DEPENDENCIES = $(top_builddir)/src/libsacado.a
tradoptestd1_LDADD = -L$(top_builddir)/src -lsacado $(LIBS) $(FLIBS)

tradoptestd2_SOURCES = tradoptestd.cpp 
tradoptestd2_CXXFLAGS = -USACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestd2_DEPENDENCIES = $(tradoptestd1_DEPENDENCIES)
tradoptestd2_LDADD = $(tradoptestd1_LDADD)

tradoptestd3_SOURCES = tradoptestd.cpp 
tradoptestd3_CXXFLAGS = -USACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestd3_DEPENDENCIES = $(tradoptestd1_DEPENDENCIES)
tradoptestd3_LDADD = $(tradoptestd1_LDADD)

tradoptestd4_SOURCES = tradoptestd.cpp 
tradoptestd4_CXXFLAGS = -DSACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestd4_DEPENDENCIES = $(tradoptestd1_DEPENDENCIES)
tradoptestd4_LDADD = $(tradoptestd1_LDADD)

tradoptestd5_SOURCES = tradoptestd.cpp 
tradoptestd5_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestd5_DEPENDENCIES = $(tradoptestd1_DEPENDENCIES)
tradoptestd5_LDADD = $(tradoptestd1_LDADD)

tradoptestd6_SOURCES = tradoptestd.cpp 
tradoptestd6_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestd6_DEPENDENCIES = $(tradoptestd1_DEPENDENCIES)
tradoptestd6_LDADD = $(tradoptestd1_LDADD)


tradopteste1_SOURCES = tradopteste.cpp 
tradopteste1_CXXFLAGS = -USACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradopteste1_DEPENDENCIES = $(top_builddir)/src/libsacado.a
tradopteste1_LDADD = -L$(top_builddir)/src -lsacado $(LIBS) $(FLIBS)

tradopteste2_SOURCES = tradopteste.cpp 
tradopteste2_CXXFLAGS = -USACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradopteste2_DEPENDENCIES = $(tradopteste1_DEPENDENCIES)
tradopteste2_LDADD = $(tradopteste1_LDADD)

tradopteste3_SOURCES = tradopteste.cpp 
tradopteste3_CXXFLAGS = -USACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradopteste3_DEPENDENCIES = $(tradopteste1_DEPENDENCIES)
tradopteste3_LDADD = $(tradopteste1_LDADD)

tradopteste4_SOURCES = tradopteste.cpp 
tradopteste4_CXXFLAGS = -DSACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradopteste4_DEPENDENCIES = $(tradopteste1_DEPENDENCIES)
tradopteste4_LDADD = $(tradopteste1_LDADD)

tradopteste5_SOURCES = tradopteste.cpp 
tradopteste5_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradopteste5_DEPENDENCIES = $(tradopteste1_DEPENDENCIES)
tradopteste5_LDADD = $(tradopteste1_LDADD)

tradopteste6_SOURCES = tradopteste.cpp 
tradopteste6_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradopteste6_DEPENDENCIES = $(tradopteste1_DEPENDENCIES)
tradopteste6_LDADD = $(tradopteste1_LDADD)


tradoptestf1_SOURCES = tradoptestf.cpp 
tradoptestf1_CXXFLAGS = -USACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestf1_DEPENDENCIES = $(top_builddir)/src/libsacado.a
tradoptestf1_LDADD = -L$(top_builddir)/src -lsacado $(LIBS) $(FLIBS)

tradoptestf2_SOURCES = tradoptestf.cpp 
tradoptestf2_CXXFLAGS = -USACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestf2_DEPENDENCIES = $(tradoptestf1_DEPENDENCIES)
tradoptestf2_LDADD = $(tradoptestf1_LDADD)

tradoptestf3_SOURCES = tradoptestf.cpp 
tradoptestf3_CXXFLAGS = -USACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestf3_DEPENDENCIES = $(tradoptestf1_DEPENDENCIES)
tradoptestf3_LDADD = $(tradoptestf1_LDADD)

tradoptestf4_SOURCES = tradoptestf.cpp 
tradoptestf4_CXXFLAGS = -DSACADO_NAMESPACE \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestf4_DEPENDENCIES = $(tradoptestf1_DEPENDENCIES)
tradoptestf4_LDADD = $(tradoptestf1_LDADD)

tradoptestf5_SOURCES = tradoptestf.cpp 
tradoptestf5_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_AUTO_AD_Const \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestf5_DEPENDENCIES = $(tradoptestf1_DEPENDENCIES)
tradoptestf5_LDADD = $(tradoptestf1_LDADD)

tradoptestf6_SOURCES = tradoptestf.cpp 
tradoptestf6_CXXFLAGS = -DSACADO_NAMESPACE -DRAD_NO_EQ_ALIAS \
	-I$(top_srcdir)/src \
	-I$(top_builddir)/src
tradoptestf6_DEPENDENCIES = $(tradoptestf1_DEPENDENCIES)
tradoptestf6_LDADD = $(tradoptestf1_LDADD)

tradoptesta.cpp : testgen.as tradoptesta.c0
	$(SHELL) testgen.as $(srcdir)/tradoptesta.c0 > tradoptesta.cpp

tradoptestb.cpp : testgen.as tradoptestb.c0
	$(SHELL) testgen.as $(srcdir)/tradoptestb.c0 > tradoptestb.cpp

tradoptestc.cpp : testgen.as tradoptestc.c0
	$(SHELL) testgen.as $(srcdir)/tradoptestc.c0 > tradoptestc.cpp

tradoptestd.cpp : testgen.as tradoptestd.c0
	$(SHELL) testgen.as $(srcdir)/tradoptestd.c0 > tradoptestd.cpp

tradopteste.cpp : testgen.as tradopteste.c0
	$(SHELL) testgen.as $(srcdir)/tradopteste.c0 > tradopteste.cpp

tradoptestf.cpp : testgen.as tradoptestf.c0
	$(SHELL) testgen.as $(srcdir)/tradoptestf.c0 > tradoptestf.cpp

CLEANFILES = \
	tradoptesta.cpp \
	tradoptestb.cpp \
	tradoptestc.cpp \
	tradoptestd.cpp \
	tradopteste.cpp \
	tradoptestf.cpp
