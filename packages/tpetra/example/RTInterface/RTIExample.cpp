#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_RTI.hpp"
#include <functional>

/** \file RTIExample.cpp
    \brief An example of the Tpetra Vector/MultiVector Reduction-Transformation Interface.
 */


namespace TpetraExamples {

  /** \class project2nd 
    \brief project2nd is a function object that takes two arguments and returns its first argument; the second argument is unused. It is essentially a generalization of identity to the case of a Binary Function.
   */
  template <class Arg1, class Arg2>
  class project2nd : public std::binary_function<Arg1,Arg2,Arg1> {
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

  /** \class Tconstant 
    \brief Tconstant is a unary function object that always returns the same value.
   */
  template <class Arg, class Result>
  class Tconstant : public std::unary_function<Arg,Result> {
  protected:
    Result c_;    
  public:
    Tconstant(const Result &c) : c_(c) {}
    Result operator()(const Arg& /* x */) const {return c_;}
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

  // 
  // Get the default communicator and node
  //
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();

  // 
  // Create a simple map with 5 local entries per node
  // 
  const size_t numLocal = 5;
  const Tpetra::global_size_t numGlobalRows = numLocal*comm->getSize();
  Teuchos::RCP<const Map> map = Tpetra::createUniformContigMapWithNode<int,int,Node>(numGlobalRows, comm, node);
  Teuchos::RCP<Vector> x = Tpetra::createVector<float>(map),
                       y = Tpetra::createVector<float>(map),
                       z = Tpetra::createVector<double>(map);

  //
  // Solitary precision testing
  //

  // sets x[i] = Tconstant<float>(1.0f)(x[i]) = 1.0f
  Tpetra::RTI::unary_transform(*x, TpetraExamples::Tconstant<float>(1.0f) );
  // sets y[i] = project2nd(y[i], x[i]) = x[i]
  Tpetra::RTI::binary_transform(*y, (const Vector &)*x, TpetraExamples::project2nd<float,float>() );
  // sets y[i] = plus<float>(y[i],x[i]) = y[i] + x[i]
  Tpetra::RTI::binary_transform(*y, (const Vector &)*x, std::plus<float>() );
  // compute dot(x,y) in single precision
  // returns plus<float>_i  multiplies<float>( x[i], y[i] ) = \sum_i x[i]*y[i]
  Tpetra::RTI::reduce(tuple<const Vector&>(*x,*y), std::multiplies<float>(), std::plus<float>() );
  // set y[i] = x while computing norm2(x)
  Tpetra::RTI::transform_reduce(*y, (const Vector &)*x, TpetraExamples::identity<float>(), std::mutliplies<float(), std::plus<float>() );

  //
  // Mixed precision testing
  // 

  // compute dot(x,y) with full double accumulator
  Tpetra::RTI::reduce(tuple<const Vector&>(*x,*y), std::multiplies<float>(), std::plus<double>() );
  // compute dot(x,y) in double precision
  Tpetra::RTI::reduce(tuple<const Vector&>(*x,*y), TpetraExamples::mprec_mult<float,float,double>(), std::plus<double>() );
  // compute z = alpha*x, where z is double precision and x is single precision
  // sets z[i] = binder1st<multiplies,alpha>( x[i] ) = multiplies(alpha,x[i]) = alpha*x[i]
  const double alpha = 2.0;
  Tpetra::RTI::unary_transform(*z, (const Vector &)*x, std::bind1st( std::multiplies<double>(), alpha ) );

  std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
}
