# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------
ABSTRACT_H = \
	$(srcdir)/NOX_Abstract_Vector.H  \
	$(srcdir)/NOX_Abstract_Group.H

# Common is **not** an installed header because it is created by configure.
COMMON = \
	$(srcdir)/common/NOX_Config.h

LINESEARCH_H = \
	$(srcdir)/NOX_Linesearch_Generic.H \
	$(srcdir)/NOX_Linesearch_Fullstep.H \
	$(srcdir)/NOX_Linesearch_Halving.H \
	$(srcdir)/NOX_Linesearch_Manager.H \
	$(srcdir)/NOX_Linesearch_MoreThunte.H \
	$(srcdir)/NOX_Linesearch_Polynomial.H 

LINESEARCH = \
	$(srcdir)/NOX_Linesearch_Fullstep.C \
	$(srcdir)/NOX_Linesearch_Halving.C \
	$(srcdir)/NOX_Linesearch_Manager.C \
	$(srcdir)/NOX_Linesearch_MoreThunte.C \
	$(srcdir)/NOX_Linesearch_Polynomial.C 

PARAMETER_H = \
	$(srcdir)/NOX_Parameter_Entry.H \
	$(srcdir)/NOX_Parameter_List.H 

PARAMETER = \
	$(srcdir)/NOX_Parameter_Entry.C \
	$(srcdir)/NOX_Parameter_List.C

SOLVER_H = \
	$(srcdir)/NOX_Solver_Generic.H \
	$(srcdir)/NOX_Solver_Newton.H 

SOLVER = \
	$(srcdir)/NOX_Solver_Newton.C 

STATUS_H = \
	$(srcdir)/NOX_Status_Test.H \
	$(srcdir)/NOX_Status_All_Tests.H \
	$(srcdir)/NOX_Status_AbsResid.H \
	$(srcdir)/NOX_Status_Combo.H \
	$(srcdir)/NOX_Status_MaxResid.H \
	$(srcdir)/NOX_Status_RelResid.H 

STATUS = \
	$(srcdir)/NOX_Status_AbsResid.C \
	$(srcdir)/NOX_Status_Combo.C \
	$(srcdir)/NOX_Status_MaxResid.C \
	$(srcdir)/NOX_Status_RelResid.C 

UTILS_H = \
	$(srcdir)/NOX_Utils.H 

UTILS = \
	$(srcdir)/NOX_Utils.C


EPETRA_H = \
	$(srcdir)/NOX_Epetra_Vector.H \
	$(srcdir)/NOX_Epetra_Group.H \
	$(srcdir)/NOX_Epetra_SharedJacobian.H \
	$(srcdir)/NOX_Epetra_Interface.H 

EPETRA = \
	$(srcdir)/NOX_Epetra_Vector.C \
	$(srcdir)/NOX_Epetra_Group.C \
	$(srcdir)/NOX_Epetra_SharedJacobian.C 

# ------------------------------------------------------------------------
# Nox library specifications
# ------------------------------------------------------------------------

NOX_LIB = libnox.a

NOX_H = \
	$(ABSTRACT_H) \
	$(LINESEARCH_H) \
	$(PARAMETER_H) \
	$(SOLVER_H) \
	$(STATUS_H) \
	$(UTILS_H) 


libnox_a_SOURCES = \
	$(COMMON) \
	$(LINESEARCH) \
	$(PARAMETER) \
	$(SOLVER) \
	$(STATUS) \
	$(UTILS) 


# ------------------------------------------------------------------------
# Noxepetra library specifications
# ------------------------------------------------------------------------

if BUILD_NOXEPETRA

NOXEPETRA_LIB = libnoxepetra.a

NOXEPETRA_H = $(EPETRA_H)

else

NOXEPETRA =

endif

libnoxepetra_a_SOURCES = $(EPETRA)

libnoxepetra_a_CXXFLAGS = \
	-I$(TRILINOS_HOME)/packages/epetra/src \
	-I$(TRILINOS_HOME)/packages/aztecoo/src


# ------------------------------------------------------------------------
# List of all libraries to install in package directory
# ------------------------------------------------------------------------

lib_LIBRARIES = $(NOX_LIB) $(NOXEPETRA_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in package directory
# ------------------------------------------------------------------------

include_HEADERS = $(NOX_H) $(NOXEPETRA_H)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






