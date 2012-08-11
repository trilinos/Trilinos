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

# The following line helps the test harness recover from build errors.

all-local:

include $(top_builddir)/Makefile.export.zoltan

EXEEXT = .exe

noinst_PROGRAMS = simpleBLOCK simpleRCB simpleGRAPH simplePHG problemGRAPH

simpleBLOCK_SOURCES      = simpleBLOCK.c
simpleBLOCK_DEPENDENCIES = $(ZOLTAN_DEPS)
simpleBLOCK_CFLAGS     = $(ZOLTAN_INCLUDES)
simpleBLOCK_LDADD        = $(ZOLTAN_LIBS)

simpleRCB_SOURCES      = simpleRCB.c
simpleRCB_DEPENDENCIES = $(ZOLTAN_DEPS)
simpleRCB_CFLAGS     = $(ZOLTAN_INCLUDES)
simpleRCB_LDADD        = $(ZOLTAN_LIBS)

simpleGRAPH_SOURCES      = simpleGRAPH.c
simpleGRAPH_DEPENDENCIES = $(ZOLTAN_DEPS)
simpleGRAPH_CFLAGS     = $(ZOLTAN_LIBRARY_INCLUDES)
simpleGRAPH_LDADD        = $(ZOLTAN_LIBS)

simplePHG_SOURCES      = simplePHG.c
simplePHG_DEPENDENCIES = $(ZOLTAN_DEPS)
simplePHG_CFLAGS     = $(ZOLTAN_INCLUDES)
simplePHG_LDADD        = $(ZOLTAN_LIBS)

problemGRAPH_SOURCES      = problemGRAPH.c
problemGRAPH_DEPENDENCIES = $(ZOLTAN_DEPS)
problemGRAPH_CFLAGS     = $(ZOLTAN_LIBRARY_INCLUDES)
problemGRAPH_LDADD        = $(ZOLTAN_LIBS)

EXTRA_DIST = classicMakefile objects.txt hypergraph.txt graph.txt mesh.txt

FNAMES=$(EXTRA_DIST:%=$(srcdir)/%)

all-local:
	if test -e Makefile.am ; then : ; else cp $(FNAMES) $(builddir); fi
