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

using namespace std;

namespace Example {
  
#undef CHKERR
#define CHKERR(ierr)                                                    \
  if (ierr != 0) { cout << "Error in " << __FILE__ << ", " << __LINE__ << endl; /* return ierr; */ }


#define MSG_NOT_YET_IMPLEMENTED ">> Not yet implemented"
#define MSG_INVALID_INPUT(what) ">> Invaid input argument: " #what
#define ERROR(msg)                              \
  { std::runtime_error(msg); }

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

  class Uplo {
  public:
    static const int Upper = 501;
    static const int Lower = 502;
  };

  class Side {
  public:
    static const int Left  = 601;
    static const int Right = 602;
  };

  class Diag {
  public:
    static const int Unit    = 701;
    static const int NonUnit = 702;
  };

  template<typename T> 
  inline 
  ostream& operator<<(ostream &os, const auto_ptr<T> &p) {
    return p->showMe(os);
  }

  class Disp {
    friend ostream& operator<<(ostream &os, const Disp &disp);
  public:
    Disp() { }
    virtual ostream& showMe(ostream &os) const {
      return os;
    }
  };
  
  inline 
  ostream& operator<<(ostream &os, const Disp &disp) {
    return disp.showMe(os);
  }  

  static map<string,string> g_graphviz_color = {
    { "ichol/scalar", "indianred2"},
    { "ichol/trsm",   "orange2"   },
    { "ichol/gemm",   "lightblue2"} };

}

#endif
