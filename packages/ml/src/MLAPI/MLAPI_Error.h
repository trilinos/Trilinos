#ifndef MLAPI_ERROR_H
#define MLAPI_ERROR_H

#ifndef ML_THROW
#ifdef HAVE_ML_CFUNC
// some old compilers do not have __func__
#define ML_THROW(str,val) { \
  std::cerr << "ERROR: In function/method " << __func__ << "()" << endl; \
  std::cerr << "ERROR: File " << __FILE__ << ", line " << __LINE__ << endl; \
  std::cerr << "ERROR: " << str << endl; \
  throw(val); \
  }
#else
#define ML_THROW(str,val) { \
  std::cerr << "ERROR: File " << __FILE__ << ", line " << __LINE__ << endl; \
  std::cerr << "ERROR: " << str << endl; \
  throw(val); \
  }
#endif // HAVE_ML_CFUNC
#endif // ndef ML_THROW
#endif // MLAPI_ERROR_H
