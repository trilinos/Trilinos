/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_mpi.h>
#include <fei_iostream.hpp>
#include <fei_ReverseMapper.hpp>
#include <fei_VectorSpace.hpp>

#include <vector>
#include <cmath>

namespace {

TEUCHOS_UNIT_TEST(ReverseMapper, test1)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs > 1) return;

  fei::VectorSpace vspace(comm);

  int fieldIDs[2] = {1, 2};
  int fieldSizes[2] = {1, 2};

  vspace.defineFields(2, &fieldIDs[0], &fieldSizes[0]);

  int idTypes[2] = {1, 2};
  vspace.defineIDTypes(2, &idTypes[0]);

  std::vector<int> ids(10);
  for(size_t i=0; i<ids.size(); ++i) ids[i] = i;

  vspace.addDOFs(fieldIDs[0], idTypes[0], ids.size(), &ids[0]);
  vspace.addDOFs(fieldIDs[1], idTypes[0], ids.size(), &ids[0]);
  vspace.addDOFs(fieldIDs[0], idTypes[1], ids.size(), &ids[0]);

  vspace.initComplete();

  fei::ReverseMapper rm(vspace);

  //fei::VectorSpace produces a set of equation-numbers grouped by
  //field, then by id, then by id-type. In other words, all components of
  //a field are contiguous in the equation-space, then all equations for
  //fields at an id (such as a node) are contiguous, etc.
  //The VectorSpace initialized above was given both fields for ids with
  //type idTypes[0], and one field for ids with type idTypes[1].
  //There should be 50 equations total, numbered 0 .. 49. ids[0] with
  //idTypes[0] should have 3 equations. The 5th equation, eqn==4, should
  //be the first scalar component of the 2nd field on ids[1].

  fei::EqnRecord er1 = rm.getEqnRecord(4);

  TEUCHOS_TEST_EQUALITY(er1.IDType, idTypes[0], out, success);
  TEUCHOS_TEST_EQUALITY(er1.ID, ids[1], out, success);

  TEUCHOS_TEST_EQUALITY(er1.fieldID, fieldIDs[1], out, success);

  TEUCHOS_TEST_EQUALITY(er1.offset, 0, out, success);
}

}//namespace <anonymous>

