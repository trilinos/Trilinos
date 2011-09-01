// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER


#include "RTOpPack_version.hpp"
#include "RTOpPack_SPMD_apply_op.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Tuple.hpp"


template<class Scalar>
bool testRTOp(
  const Teuchos::Comm<Teuchos_Ordinal> &comm,
  const Teuchos::RCP<Teuchos::FancyOStream> &out 
  )
{

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::Ptr;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::tuple;
  using Teuchos::broadcast;

  typedef Teuchos_Ordinal Ordinal;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  bool success = true;

  *out << "\n" << RTOpPack::version() << "\n";

  *out
    << "\n***"
    << "\n*** Testing RTOp SPMD support using scalar type " << ST::name()
    << "\n***\n";

  const int
    procRank = rank(comm);

  const Ordinal localDim = 10; // ToDo: Make a commandline argument!
  Ordinal localOffset = 0;
  scan(
    comm, Teuchos::REDUCE_SUM,
    as<Ordinal>(1),&localDim,&localOffset
    );
  localOffset -= localDim;

  *out
    << "\nlocalDim = " << localDim
    << "\nlocalOffset = " << localOffset
    << "\n";

  Teuchos::Array<Scalar>
    _x(localDim);

  std::fill_n(&_x[0], localDim, Scalar(1.0));

  RTOpPack::SubVectorView<Scalar>
    x(localOffset, localDim, Teuchos::arcpFromArray(_x), 1);

  *out << "\nComputing the sum of x ...\n";

  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget>
    sumTarget = sumOp.reduct_obj_create();
  RTOpPack::SPMD_apply_op<Scalar>(&comm, sumOp, 1, &x, 0, 0, &*sumTarget);
  Scalar sum = sumOp(*sumTarget);
  *out << "\nsum(x) = " << sum << "\n";
  // ToDo: Test sum(x) == ???

  *out << "\nBroadcasting the sum of x to all processes ...\n";

  RCP<RTOpPack::ReductTarget>
    sumTarget2 = sumOp.reduct_obj_create();
  if(procRank==0)
    sumOp.reduce_reduct_objs(*sumTarget, sumTarget2());
 
  RTOpPack::ReductTargetSerializer<Scalar>
    sumTargetSerializer(Teuchos::rcp(&sumOp,false));

  broadcast<Ordinal, RTOpPack::ReductTarget>(
    comm, sumTargetSerializer, 0,
    tuple<Ptr<RTOpPack::ReductTarget> >(sumTarget2.ptr())()
    );

  *out << "\nbroadcast value = " << sumOp(*sumTarget2) << "\n";

  return success;

}

int main(int argc, char* argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::CommandLineProcessor;

  typedef Teuchos_Ordinal Ordinal;
  typedef Teuchos::OrdinalTraits<Ordinal> OT;

  bool success = true, result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  try {

		CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    
		CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    RCP<FancyOStream>
      out = VerboseObjectBase::getDefaultOStream();

    //*out << std::endl << Teuchos::Teuchos_Version() << std::endl << std::endl;

    RCP<const Teuchos::Comm<Ordinal> >
      comm = Teuchos::DefaultComm<Ordinal>::getComm();
    
    result = testRTOp<double>(*comm, out);
    if(!result) success = false;
    
    if(success)
      *out << "\nEnd Result: TEST PASSED\n";
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
    
  return ( success ? 0 : 1 );
  
}
