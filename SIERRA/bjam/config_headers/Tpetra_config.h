#ifndef TPETRA_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define TPETRA_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define TPETRA_DEPRECATED
#  endif
#endif


