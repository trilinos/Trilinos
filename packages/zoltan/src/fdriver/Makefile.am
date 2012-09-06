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

include $(top_builddir)/Makefile.export.zoltan

F77=$(FC)
AM_FFLAGS=$(FCFLAGS)

#if USING_GNUMAKE
#EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ZOLTAN_LIBS))
#EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ZOLTAN_INCLUDES))
#else
EXPORT_LIBS = $(ZOLTAN_LIBS)
EXPORT_INC_PATH = $(ZOLTAN_INCLUDES)
#endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

noinst_PROGRAMS = zfdrive

if NAG_F90_COMPILER
FARG_SRC = $(srcdir)/farg_nagf95.f
else
FARG_SRC = $(srcdir)/farg_typical.f
endif

zfdrive_SOURCES = \
	$(srcdir)/mpi_h.f \
	$(srcdir)/mpistring.c \
	$(FARG_SRC) \
	$(srcdir)/mmio.f \
	$(srcdir)/fdr_sort.f90 \
	$(srcdir)/fdr_const.f90 \
	$(srcdir)/fdr_input.f90 \
	$(srcdir)/fdr_chaco_io.f90 \
	$(srcdir)/fdr_param_file.f90 \
	$(srcdir)/fdr_mm_io.f90 \
	$(srcdir)/fdr_migrate.f90 \
	$(srcdir)/fdr_loadbal.f90 \
	$(srcdir)/fdr_main.f90 

zfdrive_DEPENDENCIES = \
			$(ZOLTAN_DEPS)

zfdrive_CPPFLAGS = $(EXPORT_INC_PATH) -I$(srcdir)/ -I$(top_srcdir)/src/ch/

zfdrive_LDADD = $(EXPORT_LIBS)

EXTRA_DIST = $(srcdir)/CMakeLists.txt \
	$(srcdir)/Makefile.am \
	$(srcdir)/Makefile.in \
	$(srcdir)/README.mpich \
	$(srcdir)/farg_nagf95.f \
	$(srcdir)/farg_typical.f \
	$(srcdir)/makefile \
	$(srcdir)/zoltan_user_data.f90.old

MOSTLYCLEANFILES = *.mod
