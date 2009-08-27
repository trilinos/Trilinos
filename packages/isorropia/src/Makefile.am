# @HEADER
# ************************************************************************
#
#               Isorropia: Partitioning Package
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
#                 or Erik Boman (egboman@sandia.gov)
#
# ************************************************************************
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Epetra. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

include $(top_builddir)/Makefile.export.isorropia

CORE = \
	$(srcdir)/Isorropia_Utils.cpp \
	$(srcdir)/Isorropia_Exception.cpp

CORE_H = \
	$(srcdir)/Isorropia_CostDescriber.hpp \
	$(srcdir)/Isorropia_Operator.hpp \
	$(srcdir)/Isorropia_Partitioner.hpp \
	$(srcdir)/Isorropia_Partitioner2D.hpp \
	$(srcdir)/Isorropia_Colorer.hpp \
	$(srcdir)/Isorropia_LevelScheduler.hpp \
	$(srcdir)/Isorropia_Orderer.hpp \
	$(srcdir)/Isorropia_Redistributor.hpp \
	$(srcdir)/Isorropia_Utils.hpp \
	$(srcdir)/Isorropia_Exception.hpp \
	$(srcdir)/Isorropia_ConfigDefs.hpp \
	$(srcdir)/Isorropia_Version.hpp

EPETRA = \
	$(srcdir)/epetra/Isorropia_EpetraCostDescriber.cpp \
	$(srcdir)/epetra/Isorropia_EpetraOperator.cpp \
	$(srcdir)/epetra/Isorropia_EpetraPartitioner.cpp \
	$(srcdir)/epetra/Isorropia_EpetraColorer.cpp \
	$(srcdir)/epetra/Isorropia_EpetraLevelScheduler.cpp \
	$(srcdir)/epetra/Isorropia_EpetraOrderer.cpp \
	$(srcdir)/epetra/Isorropia_Epetra.cpp \
	$(srcdir)/epetra/Isorropia_EpetraRedistributor.cpp \
	$(srcdir)/epetra/Isorropia_EpetraLibrary.cpp \
	$(srcdir)/epetra/Isorropia_EpetraProber.cpp

EPETRA_H = \
	$(srcdir)/epetra/Isorropia_EpetraCostDescriber.hpp \
	$(srcdir)/epetra/Isorropia_EpetraPartitioner.hpp \
	$(srcdir)/epetra/Isorropia_EpetraOperator.hpp \
	$(srcdir)/epetra/Isorropia_EpetraColorer.hpp \
	$(srcdir)/epetra/Isorropia_EpetraLevelScheduler.hpp \
	$(srcdir)/epetra/Isorropia_EpetraOrderer.hpp \
	$(srcdir)/epetra/Isorropia_Epetra.hpp \
	$(srcdir)/epetra/Isorropia_EpetraRedistributor.hpp \
	$(srcdir)/epetra/Isorropia_EpetraLibrary.hpp \
	$(srcdir)/epetra/Isorropia_EpetraProber.hpp

EPETRA_ZOLTAN = \
	$(srcdir)/epetra/zoltanLib/QueryObject.cpp \
	$(srcdir)/epetra/Isorropia_EpetraZoltanLib.cpp
#	$(srcdir)/epetra/zoltanLib/Isorropia_Zoltan_Repartition.cpp

EPETRA_ZOLTAN_H = \
	$(srcdir)/epetra/zoltanLib/QueryObject.hpp \
	$(srcdir)/epetra/Isorropia_EpetraZoltanLib.hpp
#	$(srcdir)/epetra/zoltanLib/Isorropia_Zoltan_Repartition.hpp

DEFAULT_INCLUDES = $(ISORROPIA_INCLUDES)

#Isorropia now depends on Epetra
#if BUILD_EPETRA
EPETRA_SOURCES = $(EPETRA)
EPETRA_HEADERFILES = $(EPETRA_H)
EPETRA_INC=-I$(srcdir)/epetra

#Isorropia now depends on Zoltan
#if BUILD_ZOLTAN
EPETRA_ZOLTAN_SOURCES = $(EPETRA_ZOLTAN)
EPETRA_ZOLTAN_HEADERFILES = $(EPETRA_ZOLTAN_H)
ZOLTAN_INC=-I$(srcdir)/epetra/zoltanLib
#endif

#endif

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ISORROPIA_INCLUDES))
else
EXPORT_INCLUDES = $(ISORROPIA_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES) $(EPETRA_INC) $(ZOLTAN_INC)

# ------------------------------------------------------------------------
# Isorropia library specifications
# ------------------------------------------------------------------------

ISORROPIA_LIB = libisorropia.a


libisorropia_a_SOURCES = \
	$(CORE) $(EPETRA_SOURCES) $(EPETRA_ZOLTAN_SOURCES)


#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libisorropia_a_AR = $(ALTERNATE_AR)
else

libisorropia_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libisorropia_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

lib_LIBRARIES = $(ISORROPIA_LIB)

include_HEADERS = $(CORE_H) $(EPETRA_HEADERFILES) $(EPETRA_ZOLTAN_HEADERFILES)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Isorropia_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Isorropia_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in
