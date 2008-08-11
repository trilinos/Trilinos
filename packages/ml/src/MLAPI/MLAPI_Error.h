/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#ifndef MLAPI_ERROR_H
#define MLAPI_ERROR_H

#include <string>
#include <iostream>

namespace MLAPI {

typedef struct StackEntry {
  int line;
  std::string FileName;
  std::string FuncName;
} Entry;

#ifdef MLAPI_CHECK
#ifdef HAVE_ML_CFUNC
#define StackPush() \
  StackPush_(__PRETTY_FUNCTION__, __FILE__, __LINE__)
#else
#define StackPush() \
  StackPush_("function not available", __FILE__, __LINE__)
#endif
#else
#define StackPush()
#endif

#ifdef MLAPI_CHECK
void StackPush_(std::string FuncName, std::string FileName, int line);

void StackPop();

void StackPrint();
#else
inline void StackPop() {}
inline void StackPrint() {std::cout << "Compile with -DMLAPI_CHECK to get the function stack" << std::endl;}
#endif


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
