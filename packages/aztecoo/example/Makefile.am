#@HEADER
#@HEADER

#SUBDIRS = azoo_iterate_hb aztec_app aztec_hb AztecOO AztecOO_UserOpUserMat \
#	IfpackAztecOO simple_epetra_aztecoo
# note that simple_epetra_aztecoo, any AztecOO*, and azoo should not be built when "core Aztec" is 
# being built

#One example is dependent on Triutils.  Rather than making AztecOO dependent on
# Triutils, the example is built only if Triutils is built.
if BUILD_TRIUTILS
USER_OP_USER_MAT_SUBDIR=AztecOO_UserOpUserMat
endif

SUBDIRS = $(USER_OP_USER_MAT_SUBDIR) simple_epetra_aztecoo
