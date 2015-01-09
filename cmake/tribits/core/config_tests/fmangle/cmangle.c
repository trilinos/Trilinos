/* -*- C -*- */

#if defined(FC_FN_LOWER)
#  if defined(FC_FN_UNDER)
#    define FN1 f_
#    define FN2 f_f_
#  elif defined(FC_FN_SECOND_UNDER)
#    define FN1 f_
#    define FN2 f_f__
#  elif defined(FC_FN_NO_UNDER)
#    define FN1 f
#    define FN2 f_f
#  else
#    error "Unknown setting for underscore convention."
#  endif
#elif defined(FC_FN_UPPER)
#  if defined(FC_FN_UNDER)
#    define FN1 F_
#    define FN2 F_F_
#  elif defined(FC_FN_SECOND_UNDER)
#    define FN1 F_
#    define FN2 F_F__
#  elif defined(FC_FN_NO_UNDER)
#    define FN1 F
#    define FN2 F_F
#  else
#    error "Unknown setting for underscore convention."
#  endif
#else
#  error "Unknown setting for function case."
#endif

extern void FN1 (void);
extern void FN2 (void);

int main (void)
{
    FN1 ();
    FN2 ();
}
