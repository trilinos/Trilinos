#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_Version.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <Tpetra_MatrixIO.hpp>
#include <Tpetra_RTI.hpp>
#include <TpetraExt_TypeStack.hpp>

#include <iostream>

#ifdef HAVE_TPETRA_QD
#include <qd/qd_real.h>
#endif

#include "MultiPrecCG.hpp"

/** \file MultiPrecExample_qd_dd_double.cpp
    \brief An example of a multi-precision algorithm, using a flexible preconditioned CG with recursive precision.
 */

int main(int argc, char *argv[])
{
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using std::pair;
  using std::make_pair;
  using std::plus;
  using std::cout;
  using std::endl;
  using TpetraExamples::make_pair_op;
  using Tpetra::RTI::reductionGlob;
  using Tpetra::RTI::ZeroOp;
  using Tpetra::RTI::binary_pre_transform_reduce;
  using Tpetra::RTI::binary_transform;

#ifdef HAVE_TPETRA_QD
  TPETRAEXT_TYPESTACK3(MPStack, qd_real, dd_real, double )
#else
  TPETRAEXT_TYPESTACK2(MPStack, double, float )
#endif

  // 
  // Get the default communicator and node
  //
  auto &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  auto comm = platform.getComm();
  auto node = platform.getNode();
  const int myImageID = comm->getRank();

  // Static types
  typedef typename MPStack::type   S;
  typedef int                     LO;
  typedef int                     GO;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
  typedef Tpetra::CrsMatrix<S,LO,GO,Node> CrsMatrix;
  typedef Tpetra::Vector<S,LO,GO,Node>       Vector;

  //
  // Get example parameters from command-line processor
  //  
  bool verbose = (myImageID==0);
  std::string matfile;
  std::string xmlfile;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("matrix-file",&matfile,"Filename for matrix");
  cmdp.setOption("param-file", &xmlfile,"XML file for solver parameters");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  Teuchos::RCP<Teuchos::FancyOStream> out; 
  if (verbose) out = Teuchos::getFancyOStream(Teuchos::rcp(&std::cout,false));
  else         out = Teuchos::getFancyOStream(Teuchos::rcp(new Teuchos::oblackholestream()));

  // 
  // Say hello, print some communicator info
  //
  *out << "\n" << Tpetra::version() << endl
       << "Comm info: " << *comm
       << "Node type: " << Teuchos::typeName(*node) << endl
       << "Outermost scalar: " << Teuchos::TypeNameTraits<S>::name() << endl
       << endl;

  // read the matrix
  RCP<CrsMatrix> A;
  Tpetra::Utils::readHBMatrix(matfile,comm,node,A);

  // get the solver parameters
  Teuchos::ParameterList stackPL;
  // default solver stack parameters
  std::string xmlString(
    " <ParameterList>                                                       \n"
    "   <Parameter name='tolerance' value='1e-60' type='double'/>           \n"
    "   <Parameter name='numIters' value='25' type='int'/>                  \n"
    "   <Parameter name='verbose' value='2' type='int'/>                    \n"
    "   <ParameterList name='child'>                                        \n"
    "     <Parameter name='tolerance' value='1e-28' type='double'/>         \n"
    "     <Parameter name='numIters' value='25' type='int'/>                \n"
    "     <Parameter name='verbose' value='2' type='int'/>                  \n"
    "     <ParameterList name='child'>                                      \n"
    "       <Parameter name='tolerance' value='1e-14' type='double'/>       \n"
    "       <Parameter name='numIters' value='1000' type='int'/>            \n"
    "       <Parameter name='verbose' value='0' type='int'/>                \n"
    "       <Parameter name='Extract Diagonal' value='true' type='bool'/>   \n"
    "     </ParameterList>                                                  \n"
    "   </ParameterList>                                                    \n"
    " </ParameterList>                                                      \n"
  );
  Teuchos::updateParametersFromXmlString(xmlString,&stackPL);
  if (xmlfile != "") Teuchos::updateParametersFromXmlFile(xmlfile,&stackPL);

  // init the solver stack
  TpetraExamples::RFPCGInit<S,LO,GO,Node> init(A);
  RCP<ParameterList> db = Tpetra::Ext::initStackDB<MPStack>(stackPL,init);
  
  // choose a solution, compute a right-hand-side
  auto x = Tpetra::createVector<S>(A->getRowMap()),
       b = Tpetra::createVector<S>(A->getRowMap());
  x->randomize();
  A->apply(*x,*b);
  {
    auto bx = db->get<RCP<Vector>>("bx");
    binary_transform( *bx, *b, [](S, S bi) {return bi;}); // bx = b
  }

  // call the solve
  TpetraExamples::recursiveFPCG<MPStack,LO,GO,Node>(out,*db);

  //
  // Print timings
  //
  Teuchos::TimeMonitor::summarize( *out );

  //
  // Test the result
  //
  auto xhat = db->get<RCP<Vector>>("bx"),
       bhat = Tpetra::createVector<S>(A->getRowMap());
  A->apply(*xhat,*bhat);
  // compute bhat-b, while simultaneously computing |bhat-b|^2 and |b|^2
  auto nrms = binary_pre_transform_reduce(*bhat, *b, 
                                    reductionGlob<ZeroOp<pair<S,S>>>( 
                                    [](S bhati, S bi){ return bi-bhati;}, // bhati = bi-bhat
                                    [](S bhati, S bi){ return make_pair(bhati*bhati, bi*bi); },
                                    make_pair_op<S,S>(plus<S>())) );
  const S enrm = Teuchos::ScalarTraits<S>::squareroot(nrms.first),
          bnrm = Teuchos::ScalarTraits<S>::squareroot(nrms.second);
  // check that residual is as requested
  *out << "|b - A*x|/|b|: " << enrm / bnrm << endl;
  const double tolerance = db->get<double>("tolerance");
  *out << "tolerance: " << tolerance << endl;
  if (enrm / bnrm < tolerance) {
    *out << "End Result: TEST PASSED" << endl;
  }
  return 0;
}

/** \example MultiPrecExample_qd_dd_double.cpp 
    Demonstrate using Tpetra::RTI and a multi-precision flexible preconditioned CG, Tpetra::TypeStack and related utilities.
  */
