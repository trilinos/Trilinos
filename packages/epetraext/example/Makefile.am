#@HEADER
# ************************************************************************
# 
#     EpetraExt: Epetra Extended - Linear Algebra Services Package 
#                 Copyright (2001) Sandia Corporation
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
	@echo "Trilinos package epetraext subdirectory example built successfully."
	@echo ""

if BUILD_EXAMPLES

if BUILD_INOUT
if BUILD_EPETRAEXT_TRIUTILS
INOUT_SUBDIRS = inout
endif
endif

if BUILD_ZOLTAN
if BUILD_EPETRAEXT_TRIUTILS
ZOLTAN_SUBDIRS = zoltan
endif
endif

if BUILD_MATLAB

MATLAB_SUBDIRS = matlab

endif

if BUILD_COLORING

MAPCOLORING_SUBDIRS = MapColoring

endif

SUBDIRS = $(INOUT_SUBDIRS) $(ZOLTAN_SUBDIRS) $(MATLAB_SUBDIRS) \
	$(MAPCOLORING_SUBDIRS)

endif

EXTRA_DIST = model_evaluator/GLpApp/data/square/square.1.000 \
	model_evaluator/GLpApp/data/square/square.2.000 \
	model_evaluator/GLpApp/data/square/square.9x9.2.000 \
	model_evaluator/GLpApp/data/square/square.9x9.2.001 \
	model_evaluator/GLpApp/GenSQP_DataPool.hpp \
	model_evaluator/GLpApp/GenSQP_Vector.hpp \
	model_evaluator/GLpApp/GenSQP_YUEpetraVector.hpp \
	model_evaluator/GLpApp/GLpApp_AdvDiffReactOptModel.hpp \
	model_evaluator/GLpApp/GLpApp_AdvDiffReactOptModelCreator.hpp \
	model_evaluator/GLpApp/GLpApp_GLpYUEpetraDataPool.hpp \
	model_evaluator/GLpApp/rect2DMeshGenerator.hpp

