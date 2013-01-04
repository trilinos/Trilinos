/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <string>
#include <iostream>
#include <cstdlib>

namespace unit_test {
void product_tensor_legendre();
}

template <typename scalar>
int mainHost(bool test_flat, bool test_orig, bool test_block, bool check);

template <typename scalar>
int mainCuda(bool test_flat, bool test_orig, bool test_block, bool check, 
	     int device_id);

int main(int argc, char *argv[])
{
  // Defaults
  bool test_host = true;
  bool test_cuda = true;
  bool test_block = true;
  bool test_flat = true;
  bool test_orig = true;
  bool check = true;
  bool single = false;
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
    else if (s == "check")
      check = true;
    else if (s == "no-check")
      check = false;
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
	      << " [no-][cuda|host|block|flat|orig|check] [single|double] [device device_id]" 
	      << std::endl << "Defaults are all enabled." << std::endl;
    return -1;
  }

  unit_test::product_tensor_legendre();

#if 1

  if (test_host) {
    if (single)
      mainHost<float>(test_flat, test_orig, test_block, check);
    else
      mainHost<double>(test_flat, test_orig, test_block, check);
  }
  if (test_cuda) {
    if (single)
      mainCuda<float>(test_flat, test_orig, test_block, check, device);
    else
      mainCuda<double>(test_flat, test_orig, test_block, check, device);
  }

#endif

  return 0 ;
}

