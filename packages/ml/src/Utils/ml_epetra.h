#ifndef ML_EPETRA_H
#define ML_EPETRA_H

// prints out an error message if variable is not zero,
// and return this value. This macro always returns.
#define ML_RETURN(ml_err) \
  { if (ml_err != 0) { \
    cerr << "ML::ERROR:: " << ml_err << ", " \
      << __FILE__ << ", line " << __LINE__ << endl; } \
      return(ml_err);  } 

// prints out an error message if variable is not zero,
// and return this value.
#define ML_CHK_ERR(ml_err) \
  { if (ml_err != 0) { \
    cerr << "ML::ERROR:: " << ml_err << ", " \
      << __FILE__ << ", line " << __LINE__ << endl; \
      return(ml_err);  } }

// prints out an error message if variable is not zero
// and returns.
#define ML_CHK_ERRV(ml_err) \
  { if (ml_err != 0) { \
    cerr << "ML::ERROR:: " << ml_err << ", " \
      << __FILE__ << ", line " << __LINE__ << endl; \
    return; } }

#define ML_EXIT(ml_err) \
  { if (ml_err != 0) { \
    cerr << "ML::FATAL ERROR:: " << ml_err << ", " \
      << __FILE__ << ", line " << __LINE__ << endl; } \
    exit(ml_err); }

#endif
