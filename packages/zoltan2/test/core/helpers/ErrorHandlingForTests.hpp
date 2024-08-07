// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    std::cout << "Error: " << s << std::endl;\
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
    std::cout << "Error: " << s << std::endl;\
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
    std::cout << "Error: " << s << std::endl;\
    std::cout << "FAIL" << std::endl;\
  } \
  return (rc); \
} \
}

int globalFail(const Comm<int> &comm, int fail)
{
  int gfail=0;
  reduceAll<int,int>(comm, Teuchos::REDUCE_SUM, 1, &fail, &gfail);
  return gfail;
}

void printFailureCode(const Comm<int> &comm, int fail)
{
  int rank = comm.getRank();
  int nprocs = comm.getSize();
  comm.barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank)
      std::cout << rank << ": " << fail << std::endl;
    comm.barrier();
  }
  comm.barrier();
  if (rank==0) std::cout << "FAIL" << std::endl;
  exit(1);
}

#endif
