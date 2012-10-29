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
#include <Tpetra_RTI.hpp>

#include <functional>

// Function object that takes two arguments and returns the second
template <class S>
class AssignSecond {
  public:
  AssignSecond() {}
  S operator()(const S &s1, const S &s2) {return s2;}
};

// A few examples of RTI capability
template <class Node>
void simple_rti_examples(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node)
{
  using Teuchos::RCP;
  using Tpetra::RTI::reduce;
  using Tpetra::RTI::ZeroOp;
  using Tpetra::RTI::reductionGlob;
  using Tpetra::RTI::unary_transform;
  using Tpetra::RTI::binary_transform;
  const Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  // Create a simple map with 1000 local entries per node
  const size_t numLocal   = 1000*comm->getSize();
  auto map = Tpetra::createContigMapWithNode<int,int>(INVALID,numLocal,comm,node);
  auto   dx = Tpetra::createVector<double>(map),
         dy = Tpetra::createVector<double>(map);

  // Set dx to random
  dx->randomize();

  // Assign dy = dx, multiple ways!
  // Via functor
  binary_transform( *dy, *dx, AssignSecond<double>() );
  // Via C++11 lambda expression
  binary_transform( *dy, *dx, [](double, double xx){return xx;} );
  // Via convenient macro
  TPETRA_BINARY_TRANSFORM( dy, dx,   dx ); 

  // Perform multi-precision inner product...
  // floating point inner product with double precision accumulator
  float fresult; double dresult;
  auto fx = Tpetra::createVector<float>(map),
       fy = Tpetra::createVector<float>(map);
  TPETRA_BINARY_TRANSFORM( fx, dx, (float)dx );
  TPETRA_BINARY_TRANSFORM( fy, dy, (float)dy );
  // ... using a composite adaptor and standard functors
  fresult = reduce(*fx, *fy, reductionGlob<ZeroOp<double>>(std::multiplies<float>(), std::plus<double>()) );
  // ... using a convenience macro to generate all of that
  fresult = TPETRA_REDUCE2( fx, fy,   fx*fy, ZeroOp<double>, std::plus<double>() );

  // compare against double precision approach
  dresult = TPETRA_REDUCE2( dx, dy,   dx*dy, ZeroOp<double>, std::plus<double>() );
  if (comm->getRank() == 0) {
    std::cout << " Float product/double accumulator: " << fresult << std::endl;
    std::cout << "Double product/double accumulator: " << dresult << std::endl;
  }
}

// User-specified run object; think of this as main(), but templated on the node type
template <class Node>
class runTest {
  public:
    // HybridNode requires a method with this signature
    static void run(Teuchos::ParameterList &myMachPL, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) 
    {
      std::cout << "Running test with Node==" << Teuchos::typeName(*node) << " on rank " << comm->getRank() << "/" << comm->getSize() << std::endl;
      simple_rti_examples<Node>(comm,node);
    }
};

// Actual main(): configure the HybridPlatform
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

  // Load the machine ParameterList from the string above; this will typically 
  // be loaded from an XML file
  Teuchos::ParameterList machPL;
  Teuchos::updateParametersFromXmlString(machine_list, inOutArg(machPL));
  Tpetra::HybridPlatform platform(comm,machPL);
  platform.runUserCode<runTest>();

  if (comm->getRank() == 0) std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
}
