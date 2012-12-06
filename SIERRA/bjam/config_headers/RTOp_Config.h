#define HAVE_MPI

/* #undef HAVE_RTOP_DEBUG */

/* #undef HAVE_RTOP_EXPLICIT_INSTANTIATION */

#ifndef RTOP_DEPRECATED
#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define RTOP_DEPRECATED  __attribute__((__deprecated__))
#  else
#    define RTOP_DEPRECATED
#  endif
#endif

