/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_XMLParameterListReader.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_HybridPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

template <class Node, class Scalar, class Ordinal>
void some_method(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node)
{
  typedef Tpetra::Map<Ordinal,Ordinal,Node>           Map;
  typedef Tpetra::Vector<Scalar,Ordinal,Ordinal,Node> Vector;
  using Teuchos::RCP;
  const Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  // create a Map
  const size_t numLocal   = 1000;
  RCP<const Map> map = Tpetra::createContigMapWithNode<Ordinal,Ordinal>(INVALID,numLocal,comm,node);
  RCP<Vector>    vec = Tpetra::createVector<Scalar>(map);
  // randomize the vector
  vec->randomize();
  // get the norm
  Scalar nrm2 = vec->norm2();
  if (comm->getRank() == 0) std::cout << "Norm of rand vector: " << nrm2 << std::endl;
  // zero the vector
  vec->putScalar(ZERO);
  // take norms; they should be zero
  nrm2 = vec->norm2();
  if (comm->getRank() == 0) std::cout << "Norm of zero vector: " << nrm2 << std::endl;
}

template <class Node>
class runTest {
  public:
  static void run(Teuchos::ParameterList &myMachPL, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) 
  {
    std::cout << "Running test with Node==" << Teuchos::typeName(*node) << " on rank " << comm->getRank() << "/" << comm->getSize() << std::endl;
    some_method<Node,double,int>(comm,node);
  }
};

int main(int argc, char **argv) {
  Teuchos::GlobalMPISession mpisess(&argc,&argv,&std::cout);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));

  // 
  // read machine file and initialize platform
  // 
  std::string machine_list(
    "<ParameterList>                                                               "
    "  <ParameterList name='%1=0'>                                                 "
    "    <Parameter name='NodeType' type='string' value='Kokkos::SerialNode'/>     "
    "  </ParameterList>                                                            "
    "  <ParameterList name='=-1'>                                                  "
    "    <Parameter name='NodeType' type='string' value='Kokkos::OpenMPNode'/>     "
    "    <Parameter name='Verbose' type='int' value='0'/>                          "
    "    <Parameter name='Num Threads' type='int' value='-1'/>                     "
    "  </ParameterList>                                                            "
    "  <ParameterList name='=-2'>                                                  "
    "    <Parameter name='NodeType' type='string' value='Kokkos::TBBNode'/>        "
    "    <Parameter name='Verbose' type='int' value='0'/>                          "
    "    <Parameter name='Num Threads' type='int' value='-1'/>                     "
    "  </ParameterList>                                                            "
    "  <ParameterList name='=-3'>                                                  "
    "    <Parameter name='NodeType' type='string' value='Kokkos::TPINode'/>        "
    "    <Parameter name='Verbose' type='int' value='0'/>                          "
    "    <Parameter name='Num Threads' type='int' value='0'/>                      "
    "  </ParameterList>                                                            "
    "  <ParameterList name='=-4'>                                                  "
    "    <Parameter name='NodeType' type='string' value='Kokkos::ThrustGPUNode'/>  "
    "    <Parameter name='Verbose' type='int' value='0'/>                          "
    "    <Parameter name='Device Number' type='int' value='0'/>                    "
    "  </ParameterList>                                                            "
    "</ParameterList>                                                              "
  );
  Teuchos::ParameterList machPL;
  Teuchos::updateParametersFromXmlString(machine_list, inOutArg(machPL));
  Tpetra::HybridPlatform platform(comm,machPL);
  platform.runUserCode<runTest>();

  if (comm->getRank() == 0) std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
}
