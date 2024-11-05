/* -*- C -*- */

#if defined(FC_FN_LOWER)
#  if defined(FC_FN_UNDER)
#    define DGEMV dgemv_
#  elif defined(FC_FN_NO_UNDER)
#    define DGEMV dgemv
#  else
#    error "Unknown setting for underscore convention."
#  endif
#elif defined(FC_FN_UPPER)
#  if defined(FC_FN_UNDER)
#    define DGEMV DGEMV_
#  elif defined(FC_FN_NO_UNDER)
#    define DGEMV DGEMV
#  else
#    error "Unknown setting for underscore convention."
#  endif
#else
#  error "Unknown setting for function case."
#endif

extern void DGEMV(void);

int main (void)
{
    DGEMV();
}
