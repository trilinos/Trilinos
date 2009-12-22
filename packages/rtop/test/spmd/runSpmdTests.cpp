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
#include "Teuchos_arrayArg.hpp"


template<class Scalar>
bool testRTOp(
  const Teuchos::Comm<Teuchos_Ordinal> &comm,
  const Teuchos::RCP<Teuchos::FancyOStream> &out 
  )
{

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::arrayArg;
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

  broadcast<Ordinal>(
    comm, sumTargetSerializer, 0, 1,
    arrayArg<RTOpPack::ReductTarget*>(&*sumTarget2)()
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
