#ifndef _Epetra_ESI_CHK_ERR_h_
#define _Epetra_ESI_CHK_ERR_h_

/*
Define a macro for checking error returns when a function is called,
and produce something intelligible if the error return is non-zero.
*/

#ifdef CHK_ERR
#undef CHK_ERR
#endif

#ifdef EPETRA_ESI_ABORT_ON_ERROR
#define EPETRA_ESI_ERR_BEHAVIOR(a) { if (a != 0) abort(); else return(a); }
#else
#define EPETRA_ESI_ERR_BEHAVIOR(a) return(a);
#endif

#define CHK_ERR(a) { int pesi_err; if ((pesi_err = a) != 0) { \
                      cerr << "Epetra_ESI ERROR " << pesi_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; \
                      EPETRA_ESI_ERR_BEHAVIOR(pesi_err) \
                   } }
#endif

