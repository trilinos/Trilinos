#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_Version.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_RTI.hpp"
#include <functional>
#include <Teuchos_Tuple.hpp>

/** \file RTIExample.cpp
    \brief An example of the Tpetra Vector/MultiVector Reduction-Transformation Interface.
 */


namespace TpetraExamples {

  /** \class project2nd 
    \brief project2nd is a function object that takes two arguments and returns its first argument; the second argument is unused. It is essentially a generalization of identity to the case of a Binary Function.
   */
  template <class Arg1, class Arg2=Arg1>
  class project2nd : public std::binary_function<Arg1,Arg2,Arg2> {
  public:
    project2nd() {}
    Arg1 operator()(const Arg1& /* x */, const Arg2& y) const {return y;}
  };

  /** \class mprec_mult 
    \brief mprec_mult is a function object that takes two arguments of a specified precision and multiplies them using a different precision.
   */
  template <class Arg1, class Arg2, class MultPrec>
  class mprec_mult : public std::binary_function<Arg1,Arg2,MultPrec> {
  public:
    mprec_mult() {}
    MultPrec operator()(const Arg1& x, const Arg2& y) const {return (MultPrec)(x) * (MultPrec)(y);}
  };

  /** \class scale2nd 
    \brief Ignore the first, scale the second by alpha
   */
  template <class Arg1, class Arg2=Arg1, class ResPrec=Arg1>
  class scale2nd : public std::binary_function<Arg1,Arg2,ResPrec> {
  protected:
    ResPrec alpha_;
  public:
    scale2nd(const ResPrec &alpha) : alpha_(alpha) {}
    ResPrec operator()(const Arg1& x, const Arg2& y) const {return alpha_ * (ResPrec)(y);}
  };

  /** \class Tconstant 
    \brief Tconstant is a unary function object that always returns the same value.
   */
  template <class T, class Arg=T>
  class Tconstant : public std::unary_function<Arg,T> {
  protected:
    T c_;    
  public:
    Tconstant(const T &c) : c_(c) {}
    T operator()(const Arg& /* x */) const {return c_;}
  };

  /** \class TId
      \brief TId returns the absolute value of the input
  */
  template <class Arg>
  class TId : public std::unary_function<Arg,Arg> {
  public:
    TId() {}
    Arg operator()(const Arg& x) const {return x;}
  };

}

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  typedef Tpetra::DefaultPlatform::DefaultPlatformType           Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
  typedef Tpetra::Map<int,int,Node>                              Map;
  typedef Tpetra::Vector<float,int,int,Node>                     Vector;
  typedef Tpetra::Vector<double,int,int,Node>                    DVector;

  using Tpetra::RTI::ZeroOp;
  using Tpetra::RTI::reductionGlob;
  using TpetraExamples::Tconstant;
  using TpetraExamples::project2nd;
  using TpetraExamples::project2nd;
  using TpetraExamples::mprec_mult;
  using TpetraExamples::scale2nd;

  // 
  // Get the default communicator and node
  //
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();
  const int myRank = comm->getRank();

  //
  // Get example parameters from command-line processor
  //  
  bool verbose = (myRank==0);
  int numGlobal_user = 5*comm->getSize();
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("global-size",&numGlobal_user,"Global test size.");
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
  Teuchos::RCP<const Map> map = Tpetra::createUniformContigMapWithNode<int,int,Node>(numGlobalRows, comm, node);
  Teuchos::RCP<Vector>  x = Tpetra::createVector<float>(map),
                        y = Tpetra::createVector<float>(map);
  Teuchos::RCP<DVector> z = Tpetra::createVector<double>(map);

  //
  // Solitary precision testing
  //

  // sets x[i] = Tconstant<float>(1.0f)(x[i]) = 1.0f
  Tpetra::RTI::unary_transform(*x, Tconstant<float>(1.0f) );
  // sets y[i] = project2nd(y[i], x[i]) = x[i]
  Tpetra::RTI::binary_transform(*y, *x, project2nd<float,float>() );
  // sets y[i] = plus<float>(y[i],x[i]) = y[i] + x[i]
  Tpetra::RTI::binary_transform(*y, *x, std::plus<float>() );
  // compute dot(x,y) in single precision
  // returns plus<float>_i  multiplies<float>( x[i], y[i] ) = \sum_i x[i]*y[i]
  float fresult = Tpetra::RTI::reduce(*x, *y, reductionGlob<ZeroOp<float> >(std::multiplies<float>(), std::plus<float>()) );
  if (verbose) {
    std::cout << "result == " << fresult << ", expected value is " << numGlobalRows*2.0f << std::endl;
  }
  // set y[i] = x while computing square of norm2(x) (using y)
  fresult = Tpetra::RTI::binary_pre_transform_reduce(*y, *x, reductionGlob<ZeroOp<float> >(project2nd<float>(), std::multiplies<float>(), std::plus<float>()) );
  if (verbose) {
    std::cout << "result == " << fresult << ", expected value is " << numGlobalRows*1.0f << std::endl;
  }

  //
  // Mixed precision testing
  // 

  // compute dot(x,y) with double accumulator
  double dresult = Tpetra::RTI::reduce(*x, *y, reductionGlob<ZeroOp<double> >(std::multiplies<float>(), std::plus<double>()) );
  if (verbose) {
    std::cout << "result == " << dresult << ", expected value is " << numGlobalRows*1.0 << std::endl;
  }
  // compute dot(x,y) in full double precision
  dresult = Tpetra::RTI::reduce(*x, *y, reductionGlob<ZeroOp<double> >(mprec_mult<float,float,double>(), std::plus<double>()) );
  if (verbose) {
    std::cout << "result == " << dresult << ", expected value is " << numGlobalRows*1.0 << std::endl;
  }
  // compute z = alpha*x, where z is double precision and x is single precision
  // sets z[i] = scale2nd(z[i], x[i]) = alpha*x[i]
  Tpetra::RTI::binary_transform(*z, *x, scale2nd<float,double>(2.0) );

  if (verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}
