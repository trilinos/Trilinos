#ifndef MLAPI_ERROR_H
#define MLAPI_ERROR_H

#include <string>

namespace MLAPI {

typedef struct StackEntry {
  int line;
  std::string FileName;
  std::string FuncName;
} Entry;

#define StackPush() \
  StackPush_(__PRETTY_FUNCTION__, __FILE__, __LINE__)

void StackPush_(std::string FuncName, std::string FileName, int line);

void StackPop();

void StackPrint();

} // namespace MLAPI

#ifndef ML_THROW
#ifdef HAVE_ML_CFUNC
// some old compilers do not have __func__
#define ML_THROW(str,val) { \
  std::cerr << "ERROR: In " << __PRETTY_FUNCTION__ << "()" << endl; \
  std::cerr << "ERROR: File " << __FILE__ << ", line " << __LINE__ << endl; \
  std::cerr << "ERROR: " << str << endl; \
  StackPrint(); \
  throw(val); \
  }
#else
#define ML_THROW(str,val) { \
  std::cerr << "ERROR: File " << __FILE__ << ", line " << __LINE__ << endl; \
  std::cerr << "ERROR: " << str << endl; \
  StackPrint(); \
  throw(val); \
  }
#endif // HAVE_ML_CFUNC
#endif // ndef ML_THROW

#endif // MLAPI_ERROR_H
