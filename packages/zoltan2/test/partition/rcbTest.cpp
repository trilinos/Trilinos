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

/*! \file rcb.cpp
    \brief An example of partitioning coordinates with RCB.
    \todo add more cases to this test.
*/

#include <Zoltan2_config.h>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> myTypes_t;
// Need to use zgno_t as zzgid_t in basic types, since zzgid_t=zgno_t in MultiVector


/*! \test rcbTest.cpp
    An example of the use of the RCB algorithm to partition coordinate data.
    \todo error handling
    \todo write some examples that don't use teuchos
    \todo check the solution, visualize it somehow
*/

void testFromDataFile(
  const RCP<const Teuchos::Comm<int> > & comm,
  int nParts,
  string &filename,
  bool doRemap
)
{
  int me = comm->getRank();
  if (me == 0)
    std::cout << "Parallel partitioning of " << filename << ".mtx: "
         << nParts << " parts." << std::endl;

  std::string fname(filename);
  UserInputForTests uinput(testDataFilePath, fname, comm, true);

  RCP<tMVector_t> coords = uinput.getUICoordinates();
  if (me == 0)
    std::cout << "Multivector length = " << coords->getGlobalLength()
         << " Num vectors = " << coords->getNumVectors() << std::endl;

  RCP<const tMVector_t> coordsConst = rcp_const_cast<const tMVector_t>(coords);

  typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;
  inputAdapter_t ia(coordsConst);
  if (me == 0)
    std::cout << "Adapter constructed" << std::endl;

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("num_global_parts", nParts);
  params.set("algorithm", "rcb");
  params.set("imbalance_tolerance", 1.1);
  if (doRemap) params.set("remap_parts", true); // bool parameter

#ifdef HAVE_ZOLTAN2_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params,
    MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif
  if (me == 0)
    std::cout << "Problem constructed" << std::endl;


  problem.solve();
  if (me == 0)
    std::cout << "Problem solved" << std::endl;
}

void serialTest(int numParts, bool doRemap)
{
  int numCoords = 1000;
  numParts *= 8;

  std::cout << "Serial partitioning: " << numParts << " parts." << std::endl;

  zgno_t *ids = new zgno_t [numCoords];
  if (!ids)
    throw std::bad_alloc();
  for (int i=0; i < numCoords; i++)
    ids[i] = i;
  ArrayRCP<zgno_t> globalIds(ids, 0, numCoords, true);

  Array<ArrayRCP<zscalar_t> > randomCoords(3);
  UserInputForTests::getUIRandomData(555, numCoords, 0, 10,
    randomCoords.view(0,3));

  typedef Zoltan2::BasicVectorAdapter<myTypes_t> inputAdapter_t;

  inputAdapter_t ia(numCoords, ids,
    randomCoords[0].getRawPtr(), randomCoords[1].getRawPtr(),
     randomCoords[2].getRawPtr(), 1,1,1);

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("num_global_parts", numParts);
  params.set("algorithm", "rcb");
  params.set("imbalance_tolerance", 1.1);
  if (doRemap) params.set("remap_parts", true); // bool parameter

#ifdef HAVE_ZOLTAN2_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> serialProblem(
    &ia, &params, MPI_COMM_SELF);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> serialProblem(&ia, &params);
#endif

  serialProblem.solve();
}

void meshCoordinatesTest(const RCP<const Teuchos::Comm<int> > & comm)
{
  int xdim = 40;
  int ydim = 60;
  int zdim = 20;
  UserInputForTests uinput(xdim, ydim, zdim, string("Laplace3D"), comm, true, true);

  RCP<tMVector_t> coords = uinput.getUICoordinates();

  size_t localCount = coords->getLocalLength();

  zscalar_t *x=NULL, *y=NULL, *z=NULL;
  x = coords->getDataNonConst(0).getRawPtr();
  y = coords->getDataNonConst(1).getRawPtr();
  z = coords->getDataNonConst(2).getRawPtr();

  const zgno_t *globalIds = coords->getMap()->getNodeElementList().getRawPtr();
  typedef Zoltan2::BasicVectorAdapter<tMVector_t> inputAdapter_t;

  inputAdapter_t ia(localCount, globalIds, x, y, z, 1, 1, 1);

  Teuchos::ParameterList params("test params");
  params.set("rectilinear", true); // bool parameter

#ifdef HAVE_ZOLTAN2_MPI
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params, MPI_COMM_WORLD);
#else
  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params);
#endif

  problem.solve();
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > tcomm = Tpetra::getDefaultComm();

  int rank = tcomm->getRank();
  int nParts = tcomm->getSize();
  bool doRemap = false;
  string filename = "USAir97";

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("file", &filename, "Name of the Matrix Market file to read");
  cmdp.setOption("nparts", &nParts, "Number of parts.");
  cmdp.setOption("remap", "no-remap", &doRemap, "Remap part numbers.");
  cmdp.parse(narg, arg);

  meshCoordinatesTest(tcomm);

  testFromDataFile(tcomm, nParts, filename, doRemap);

  if (rank == 0)
    serialTest(nParts, doRemap);

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}
