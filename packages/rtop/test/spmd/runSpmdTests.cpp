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
  const Teuchos::Comm<Teuchos_Index>                   &comm
  ,const Teuchos::RefCountPtr<Teuchos::FancyOStream>   &out 
  )
{

  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::arrayArg;

  typedef Teuchos_Index Ordinal;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  bool success = true, result;

  *out
    << "\n***"
    << "\n*** Testing RTOp using scalar type " << ST::name()
    << "\n***\n";

  const int
    procRank = rank(comm);

  const Ordinal localDim = 10; // ToDo: Make a commandline argument!
  Ordinal localOffset = 0;
  scan(
    comm,Teuchos::REDUCE_SUM
    ,Ordinal(1),&localDim,&localOffset
    );
  localOffset -= localDim;

  *out
    << "\nlocalDim = " << localDim
    << "\nlocalOffset = " << localOffset
    << "\n";

  Teuchos::Array<Scalar>
    _x(localDim);

  std::fill_n(&_x[0],localDim,Scalar(1.0));

  RTOpPack::SubVectorView<Scalar>
    x(localOffset,localDim,&_x[0],1);

  *out << "\nComputing the sum of x ...\n";

  RTOpPack::ROpSum<Scalar> sumOp;
  RefCountPtr<RTOpPack::ReductTarget>
    sumTarget = sumOp.reduct_obj_create();
  RTOpPack::SPMD_apply_op<Scalar>(
    &comm,sumOp,1,&x,0,0,&*sumTarget
    );
  Scalar sum = sumOp(*sumTarget);
  *out << "\nsum(x) = " << sum << "\n";
  // ToDo: Test sum(x) == ???

  *out << "\nBroadcasting the sum of x to all processes ...\n";

  RefCountPtr<RTOpPack::ReductTarget>
    sumTarget2 = sumOp.reduct_obj_create();
  if(procRank==0)
    sumOp.reduce_reduct_objs(*sumTarget,&*sumTarget2);
 
  RTOpPack::ReductTargetSerializer<Scalar>
    sumTargetSerializer(Teuchos::rcp(&sumOp,false));

  broadcast(
    comm,sumTargetSerializer,0,1
    ,arrayArg<RTOpPack::ReductTarget*>(&*sumTarget2)()
    );

  *out << "\nbroadcast value = " << sumOp(*sumTarget2) << "\n";

  return success;

}

int main(int argc, char* argv[])
{

  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::CommandLineProcessor;

  typedef Teuchos_Index Ordinal;
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

    RefCountPtr<FancyOStream>
      out = VerboseObjectBase::getDefaultOStream();

    //*out << std::endl << Teuchos::Teuchos_Version() << std::endl << std::endl;

    RefCountPtr<const Teuchos::Comm<Ordinal> >
      comm = Teuchos::DefaultComm<Ordinal>::getComm();
    
    result = testRTOp<double>(*comm,out);
    if(!result) success = false;
    
    if(success)
      *out << "\nEnd Result: TEST PASSED\n";
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
    
  return ( success ? 0 : 1 );
  
}
