#@HEADER
#@HEADER
#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose tests are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# test1/ rather than this file.

if HAVE_TEUCHOS_EXPERIMENTAL

SUBDIRS = scripts BLAS DenseMatrix LAPACK ParameterList Containers Time Random RefCountPtr XML

else

SUBDIRS = scripts BLAS DenseMatrix LAPACK ParameterList

endif
