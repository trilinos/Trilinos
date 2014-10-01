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
#include <algorithm>

#include <cmath>

using namespace std;

namespace Example {
  
#undef CHKERR
#define CHKERR(ierr)                                                    \
  if (ierr != 0) { cout << "Error in " << __FILE__ << ", " << __LINE__ << endl; /* return ierr; */ }
}

#endif
