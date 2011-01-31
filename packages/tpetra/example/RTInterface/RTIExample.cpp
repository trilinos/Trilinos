#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_Version.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_RTI.hpp"

#include <iostream>
#include <iomanip>
#include <functional>

/** \file RTIExample.cpp
    \brief An example of the Tpetra Vector/MultiVector Reduction-Transformation Interface.
 */

namespace TpetraExamples {

  /** \class mprec_mult 
    \brief mprec_mult is a function object that takes two arguments of a specified precision and multiplies them using a different precision.
   */
  template <class Arg1, class Arg2, class MultPrec>
  class mprec_mult : public std::binary_function<Arg1,Arg2,MultPrec> {
  public:
    mprec_mult() {}
    MultPrec operator()(const Arg1& x, const Arg2& y) const {return (MultPrec)(x) * (MultPrec)(y);}
  };

}

int main(int argc, char *argv[])
{
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  using Teuchos::TimeMonitor;
  using Tpetra::RTI::ZeroOp;
  using Tpetra::RTI::reductionGlob;
  using TpetraExamples::mprec_mult;

  // 
  // Get the default communicator and node
  //
  auto &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  auto comm = platform.getComm();
  auto node = platform.getNode();
  const int myImageID = comm->getRank();

  //
  // Get example parameters from command-line processor
  //  
  bool verbose = (myImageID==0);
  int numGlobal_user = 100*comm->getSize();
  int numTimeTrials = 3;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("global-size",&numGlobal_user,"Global test size.");
  cmdp.setOption("num-time-trials",&numTimeTrials,"Number of trials in timing loops.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // Say hello, print some communicator info
  //
  if (verbose) {
    std::cout << "\n" << Tpetra::version() << std::endl;
    std::cout << "Comm info: " << *comm;
    std::cout << "Node type: " << Teuchos::typeName(*node) << std::endl;
    std::cout << std::endl;
  }


  // 
  // Create a simple map with 5 local entries per node
  // 
  Tpetra::global_size_t numGlobalRows = numGlobal_user;
  auto map = Tpetra::createUniformContigMapWithNode<int,int>(numGlobalRows, comm, node);
  const size_t numLocalRows = map->getNodeNumElements();
  auto x = Tpetra::createVector<float>(map),
       y = Tpetra::createVector<float>(map);
  auto z = Tpetra::createVector<double>(map),
       w = Tpetra::createVector<double>(map);

  //
  // Initialization and simple reduction
  //

  // sets x[i] = 1.0f
  Tpetra::RTI::unary_transform( *x, [](float xi){return 1.0f;} );
  // sets y[i] = x[i]
  Tpetra::RTI::binary_transform( *y, *x, [](float /*yi*/, float xi) {return xi;} );
  // sets y[i] = plus(y[i], x[i]) = y[i] + x[i]
  Tpetra::RTI::binary_transform( *y, *x, std::plus<float>() );
  // compute dot(x,y) in single precision, sum all x[i]*y[i]
  float fresult = Tpetra::RTI::reduce( *x, *y, 
                                       reductionGlob<
                                         ZeroOp<float>>( 
                                         std::multiplies<float>(), 
                                         std::plus<float>()
                                       ) );
  if (verbose) {
    std::cout << std::left << "dot( ones, twos ) result == " << fresult << ", expected value is " << numGlobalRows*2.0f << std::endl;
  }

  //
  // Single precision testing
  //

  // set x = [1, 1e-4, ..., 1e-4]
  {
    Teuchos::ArrayRCP<float> xdata = x->get1dViewNonConst();
    for (size_t i=1; i < numLocalRows; ++i) {
      xdata[i] = 1.0f / 4096.0f;
    }
  }
  // set y[i] = x while computing square of norm2(x) (using y): multiply x[i] * y[i] in float, sum them in float, 0.0f is the additive identity
  fresult = Tpetra::RTI::binary_pre_transform_reduce(*y, *x, reductionGlob<ZeroOp<float>>( [](float /*yi*/, float xi){return xi;}, std::multiplies<float>(), std::plus<float>()) );

  // compute pure float inner product alone, timing it for comparison
  auto timePF = TimeMonitor::getNewTimer("Pure float dot()");
  {
    TimeMonitor lcltimer(*timePF);    
    for (int l=0; l != numTimeTrials; ++l) 
       fresult = Tpetra::RTI::reduce(*x, *y, reductionGlob<ZeroOp<float>>(std::multiplies<float>(), std::plus<float>()) );
  }
  if (verbose) {
    std::cout << std::left << std::endl << std::setw(25) << "pure float result" << " == " << std::setprecision(12) << std::scientific << fresult << std::endl;
  }

  //
  // Mixed precision testing
  // 

  // compute dot(x,y) with double accumulator: multiply x[i] * y[i] in float, sum them in double, 0.0 is the additive identity
  auto timeAD = TimeMonitor::getNewTimer("Double acc. dot()");
  {
    TimeMonitor lcltimer(*timeAD);
    for (int l=0; l != numTimeTrials; ++l) 
       fresult = Tpetra::RTI::reduce(*x, *y, reductionGlob<ZeroOp<double>>(std::multiplies<float>(), std::plus<double>()) );
  }
  if (verbose) {
    std::cout << std::left << std::setw(25) << "double acc. result" << " == " << std::setprecision(12) << std::scientific << fresult << std::endl;
  }

  // compute dot(x,y) in full double precision: multiply x[i] * y[i] in double, sum them in double, 0.0 is the additive identity
  double dresult = 0.0;
  auto timeMAD = TimeMonitor::getNewTimer("Double mult./acc. dot()");
  {
    TimeMonitor lcltimer(*timeMAD);
    for (int l=0; l != numTimeTrials; ++l) 
       dresult = Tpetra::RTI::reduce(*x, *y, reductionGlob< ZeroOp<double>>(
                                                            mprec_mult<float,float,double>(), 
                                                            std::plus<double>()) );
  }
  if (verbose) {
    std::cout << std::left << std::setw(25) << "double mult/add result" << " == " << std::setprecision(12) << std::scientific << dresult << std::endl;
  }

  // compute dot(x,y) in full double precision: multiply x[i] * y[i] in double, sum them in double, 0.0 is the additive identity
  auto timePD = TimeMonitor::getNewTimer("Pure double dot()");
  // set [w,z] = x
  Tpetra::RTI::binary_transform( *w, *x, [](double /*wi*/, float xi) -> double {return xi;} );
  Tpetra::RTI::binary_transform( *z, *x, [](double /*zi*/, float xi) -> double {return xi;} );
  {
    TimeMonitor lcltimer(*timePD);
    for (int l=0; l != numTimeTrials; ++l) 
       dresult = Tpetra::RTI::reduce(*w, *z, reductionGlob< ZeroOp<double>>(
                                                            std::multiplies<double>(), 
                                                            std::plus<double>()) );
  }
  if (verbose) {
    std::cout << std::left << std::setw(25) << "pure double result" << " == " << std::setprecision(12) << std::scientific << dresult << std::endl
              << std::endl;
  }

  //
  // Compute z = alpha*x, where z is double precision and x is single precision
  //
  Tpetra::RTI::binary_transform( *z, *x, [](double /*zi*/, float xi) -> double {return 2.0*xi;} );

  //
  // Print timings
  //
  if (verbose) {
    TimeMonitor::summarize( std::cout );
  }

  if (verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}

/** \example RTIExample.cpp
    Demonstrate using Tpetra::RTI methods for transforming/reducing Tpetra::Vector objects.
  */

