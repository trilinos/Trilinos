#@HEADER
# ************************************************************************
# 
#               Isorropia: Partitioning and Load Balancing Package
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
# Questions? Contact Alan Williams (william@sandia.gov) 
#                 or Erik Boman    (egboman@sandia.gov)
# 
# ************************************************************************
#@HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package isorropia subdirectory test built successfully."
	@echo ""

if BUILD_ZOLTAN
ZOLTAN_TEST_DIR = zoltan
ZOLTAN_TEST = rebalance_zoltan.exe
else
ZOLTAN_TEST_DIR =
ZOLTAN_TEST =
endif

#
# Note: utils must be built first, before the other subdirs
#
SUBDIRS = utils utest rebalance_1d_default $(ZOLTAN_TEST_DIR)

TESTS = \
utest/utest.exe \
rebalance_1d_default/rebalance_1d_default.exe \
$(ZOLTAN_TEST)

EXTRA_DIST = definition

