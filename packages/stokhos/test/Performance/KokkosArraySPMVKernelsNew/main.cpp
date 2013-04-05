// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include <string>
#include <iostream>
#include <cstdlib>

template <typename scalar>
int mainHost(bool test_flat, bool test_orig, bool test_block, bool symmetric,
	     bool mkl);

template <typename scalar>
int mainCuda(bool test_flat, bool test_orig, bool test_block, bool symmetric,
	     int device_id);

int main(int argc, char *argv[])
{
  // Defaults
  bool test_host = true;
  bool test_cuda = true;
  bool test_block = true;
  bool test_flat = true;
  bool test_orig = true;
  bool symmetric = true;
  bool single = false;
  bool mkl = false;
  int device = 0;

  // Parse command line arguments
  bool print_usage = false;
  int i=1;
  while (i<argc) {
    std::string s(argv[i]);
    if (s == "cuda")
      test_cuda = true;
    else if (s == "no-cuda")
      test_cuda = false;
    else if (s == "host")
      test_host = true;
    else if (s == "no-host")
      test_host = false;
    else if (s == "block")
      test_block = true;
    else if (s == "no-block")
      test_block = false;
    else if (s == "flat")
      test_flat = true;
    else if (s == "no-flat")
      test_flat = false;
    else if (s == "orig")
      test_orig = true;
    else if (s == "no-orig")
      test_orig = false;
    else if (s == "symmetric")
      symmetric = true;
    else if (s == "no-symmetric")
      symmetric = false;
    else if (s == "mkl")
      mkl = true;
    else if (s == "no-mkl")
      mkl = false;
    else if (s == "single")
      single = true;
    else if (s == "double")
      single = false;
    else if (s == "device") {
      ++i;
      device = std::atoi(argv[i]);
    }
    else if (s == "-h" || s == "--help")
      print_usage = true;
    else {
      std::cout << "Invalid argument:  " << s << std::endl;
      print_usage = true;
    }
    ++i;
  }
  if (print_usage) {
    std::cout << "Usage:" << std::endl
	      << "\t" << argv[0] 
	      << " [no-][cuda|host|block|flat|orig|symmetric] [single|double] [device device_id]" 
	      << std::endl << "Defaults are all enabled." << std::endl;
    return -1;
  }

  if (test_host) {
    if (single)
      mainHost<float>(test_flat, test_orig, test_block, symmetric, mkl);
    else
      mainHost<double>(test_flat, test_orig, test_block, symmetric, mkl);
  }
  if (test_cuda) {
    if (single)
      mainCuda<float>(test_flat, test_orig, test_block, symmetric, device);
    else
      mainCuda<double>(test_flat, test_orig, test_block, symmetric, device);
  }

  return 0 ;
}

