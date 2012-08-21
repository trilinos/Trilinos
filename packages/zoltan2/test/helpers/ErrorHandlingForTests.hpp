// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef ERRORHANDLINGFORTESTS_HPP
#define ERRORHANDLINGFORTESTS_HPP

#include <Zoltan2_config.h>
#include <iostream>
#include <string>
#include <exception>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::reduceAll;

#ifdef HAVE_ZOLTAN2_MPI

#define TEST_FAIL_AND_THROW(comm, ok, s){ \
int gval, lval=( (ok) ? 0 : 1);       \
reduceAll<int,int>(comm, Teuchos::REDUCE_SUM, 1, &lval, &gval);\
if (gval){ \
  throw std::runtime_error(std::string(s)); \
} \
}

#define TEST_FAIL_AND_EXIT(comm, ok, s, code){ \
int gval, lval=( (ok) ? 0 : 1);       \
reduceAll<int,int>(comm, Teuchos::REDUCE_SUM, 1, &lval, &gval);\
if (gval){ \
  if ((comm).getRank() == 0){\
    std::cerr << "Error: " << s << std::endl;\
    std::cout << "FAIL" << std::endl;\
  } \
  exit(code);\
} \
}

#define TEST_FAIL_AND_RETURN(comm, ok, s){ \
int gval, lval=( (ok) ? 0 : 1);       \
reduceAll<int,int>(comm, Teuchos::REDUCE_SUM, 1, &lval, &gval);\
if (gval){ \
  if ((comm).getRank() == 0){\
    std::cerr << "Error: " << s << std::endl;\
    std::cout << "FAIL" << std::endl;\
  } \
  return; \
} \
}

#define TEST_FAIL_AND_RETURN_VALUE(comm, ok, s, rc){ \
int gval, lval=( (ok) ? 0 : 1);       \
reduceAll<int,int>(comm, Teuchos::REDUCE_SUM, 1, &lval, &gval);\
if (gval){ \
  if ((comm).getRank() == 0){\
    std::cerr << "Error: " << s << std::endl;\
    std::cout << "FAIL" << std::endl;\
  } \
  return (rc); \
} \
}

#else

#define TEST_FAIL_AND_THROW(comm, ok, s) \
if ((!ok)){ \
  throw std::runtime_error(std::string(s)); \
} 

#define TEST_FAIL_AND_EXIT(comm, ok, s, code) \
if (!(ok)){ \
  std::cerr << "Error: " << s << std::endl;\
  std::cout << "FAIL" << std::endl;\
  exit(code);\
} 

#define TEST_FAIL_AND_RETURN(comm, ok, s) \
if (!(ok)){ \
  std::cerr << "Error: " << s << std::endl;\
  std::cout << "FAIL" << std::endl;\
  return;\
} 

#define TEST_FAIL_AND_RETURN_VALUE(comm, ok, s, rc) \
if (!(ok)){ \
  std::cerr << "Error: " << s << std::endl;\
  std::cout << "FAIL" << std::endl;\
  return (rc);\
} 
#endif

int globalFail(const RCP<const Comm<int> > &comm, int fail)
{
  int gfail=0;
  reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &fail, &gfail);
  return gfail;
}

void printFailureCode(const RCP<const Comm<int> > &comm, int fail)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank)
      std::cout << rank << ": " << fail << std::endl;
    comm->barrier();
  }
  comm->barrier();
  if (rank==0) std::cout << "FAIL" << std::endl;
  exit(1);
}

#endif
