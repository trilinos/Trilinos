#pragma once
#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <stdio.h>
#include <string.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <memory>

#include <cmath>


/// \file util.hpp
/// \brief Utility functions and constant integer class like an enum class.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// This provides utility functions for implementing mini-app for incomplete
/// sparse matrix factorization with task-data parallelism e.g., parameter
/// classes, error handling, ostream << overloading.
///
/// Note: The reference of the "static const int" members in the enum-like 
/// classes should not be used as function arguments but their values only. 


using namespace std;

namespace Example {
  
#undef CHKERR
#define CHKERR(ierr)                                                    \
  if (ierr != 0) { cout << "Error in " << __FILE__ << ", " << __LINE__ << endl; /* return ierr; */ }


#define MSG_NOT_YET_IMPLEMENTED ">> Not yet implemented"
#define MSG_INVALID_INPUT(what) ">> Invaid input argument: " #what
#define ERROR(msg)                              \
  { std::runtime_error(msg); }

  /// \class Partition
  /// \brief Matrix partition parameters.
  class Partition { 
  public:
    static const int Top         = 101;
    static const int Bottom      = 102;

    static const int Left        = 201;
    static const int Right       = 202;

    static const int TopLeft     = 401;
    static const int TopRight    = 402;
    static const int BottomLeft  = 403;
    static const int BottomRight = 404;
  };


  /// \class Uplo
  /// \brief Matrix upper/lower parameters.
  class Uplo {
  public:
    static const int Upper = 501;
    static const int Lower = 502;
    
    template<int uplo, typename OrdinalType> 
    bool 
    is(const OrdinalType i, const OrdinalType j) {
      return true;
    }

    // template<typename OrdinalType> 
    // bool 
    // is<Uplo::Lower,OrdinalType>(const OrdinalType i, const OrdinalType j) {
    //   return (i<=j);      
    // }

  };



  /// \class Side
  /// \brief Matrix left/right parameters.
  class Side {
  public:
    static const int Left  = 601;
    static const int Right = 602;
  };

  /// \class Diag
  /// \brief Matrix unit/non-unit diag parameters.
  class Diag {
  public:
    static const int Unit    = 701;
    static const int NonUnit = 702;
  };

  /// \brief Interface for overloaded stream operators.
  template<typename T> 
  inline 
  ostream& operator<<(ostream &os, const auto_ptr<T> &p) {
    return p->showMe(os);
  }

  /// \class Disp
  /// \brief Interface for the stream operator.
  class Disp {
    friend ostream& operator<<(ostream &os, const Disp &disp);
  public:
    Disp() { }
    virtual ostream& showMe(ostream &os) const {
      return os;
    }
  };

  /// \brief Implementation of the overloaded stream operator.
  inline 
  ostream& operator<<(ostream &os, const Disp &disp) {
    return disp.showMe(os);
  }  

}

#endif
