## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                           Sacado Package
#                 Copyright (2006) Sandia Corporation
# 
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
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

if ENABLE_SACADO_TEUCHOS

# Include export Makefile for correct paths
include $(top_builddir)/Makefile.export.sacado

noinst_PROGRAMS = \
	fad_expr_depth \
	fad_expr \
	fad_lj_grad \
	taylor_expr \
	fad_fe_jac_fill \
	rad_fe_adj_fill

endif

fad_expr_depth_SOURCES = fad_expr_depth.cpp 
fad_expr_depth_CXXFLAGS = -O3  \
	$(SACADO_INCLUDES) \
	-I$(top_srcdir)/Fad \
	-I$(top_srcdir)/test/utils
fad_expr_depth_DEPENDENCIES = \
        $(top_builddir)/src/libsacado.a \
	$(top_builddir)/test/utils/libsacadoutils.a
fad_expr_depth_LDADD = \
	-L$(top_builddir)/test/utils -lsacadoutils \
	$(SACADO_LIBS) 

fad_expr_SOURCES = fad_expr.cpp
fad_expr_CXXFLAGS = $(fad_expr_depth_CXXFLAGS)
fad_expr_DEPENDENCIES = $(fad_expr_depth_DEPENDENCIES)
fad_expr_LDADD = $(fad_expr_depth_LDADD)

fad_lj_grad_SOURCES = fad_lj_grad.cpp
fad_lj_grad_CXXFLAGS = $(fad_expr_depth_CXXFLAGS)
fad_lj_grad_DEPENDENCIES = $(fad_expr_depth_DEPENDENCIES)
fad_lj_grad_LDADD = $(fad_expr_depth_LDADD)

if ENABLE_ADOLC
TAYLOR_CXXFLAGS = -I$(ADOLC_INCDIR)
TAYLOR_LDADD = -L$(ADOLC_LIBDIR) $(ADOLC_LIBS)
else
TAYLOR_CXXFLAGS = 
TAYLOR_LDADD = 
endif

taylor_expr_SOURCES = taylor_expr.cpp
taylor_expr_CXXFLAGS = $(TAYLOR_CXXFLAGS) $(fad_expr_depth_CXXFLAGS)
taylor_expr_DEPENDENCIES = $(fad_expr_depth_DEPENDENCIES)
taylor_expr_LDADD = $(TAYLOR_LDADD) $(fad_expr_depth_LDADD)

fad_fe_jac_fill_SOURCES = fad_fe_jac_fill.cpp
fad_fe_jac_fill_CXXFLAGS = $(TAYLOR_CXXFLAGS) $(fad_expr_depth_CXXFLAGS)
fad_fe_jac_fill_DEPENDENCIES = $(fad_expr_depth_DEPENDENCIES)
fad_fe_jac_fill_LDADD = $(TAYLOR_LDADD) $(fad_expr_depth_LDADD)

rad_fe_adj_fill_SOURCES = rad_fe_adj_fill.cpp
rad_fe_adj_fill_CXXFLAGS = $(TAYLOR_CXXFLAGS) $(fad_expr_depth_CXXFLAGS)
rad_fe_adj_fill_DEPENDENCIES = $(fad_expr_depth_DEPENDENCIES)
rad_fe_adj_fill_LDADD = $(TAYLOR_LDADD) $(fad_expr_depth_LDADD)