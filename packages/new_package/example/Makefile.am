#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = hello_test

hello_test_SOURCES = $(srcdir)/cxx_main.cpp

#Specify locations of include files within Trilinos, but outside New_Package
# REQUIRED CHANGE - Add packages that new_package is dependent on
# OPTIONAL CHANGE - Remove epetra if new_package does not need it
#  
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

# REQUIRED CHANGE - Change the nam eof the library
hello_test_DEPENDENCIES = \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/src/libnew_package.a 

hello_test_CXXFLAGS = -I$(top_srcdir)/src

# REQUIRED CHANGE - Change the name of the library
hello_test_LDADD = -L$(top_builddir)/src -lnew_package \
	-L$(top_builddir)/../epetra/src -lepetra \
	$(BLAS_LIBS) $(FLIBS)

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 
