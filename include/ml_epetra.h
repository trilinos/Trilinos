/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_EPETRA_H
#define ML_EPETRA_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

// prints out an error message if variable is not zero,
// and return this value. This macro always returns.
#define ML_RETURN(ml_err) \
  { if (ml_err != 0) { \
    std::cerr << "ML::ERROR:: " << ml_err << ", " \
      << __FILE__ << ", line " << __LINE__ << std::endl; } \
      return(ml_err);  }

// prints out an error message if variable is not zero,
// and return this value.
#define ML_CHK_ERR(ml_err) \
  { if (ml_err != 0) { \
    std::cerr << "ML::ERROR:: " << ml_err << ", " \
      << __FILE__ << ", line " << __LINE__ << std::endl; \
      return(ml_err);  } }

// prints out an error message if variable is not zero
// and returns.
#define ML_CHK_ERRV(ml_err) \
  { if (ml_err != 0) { \
    std::cerr << "ML::ERROR:: " << ml_err << ", " \
      << __FILE__ << ", line " << __LINE__ << std::endl; \
    return; } }

#define ML_EXIT(ml_err) \
  { if (ml_err != 0) { \
    std::cerr << "ML::FATAL ERROR:: " << ml_err << ", " \
      << __FILE__ << ", line " << __LINE__ << std::endl; } \
    exit(ml_err); }

#endif
