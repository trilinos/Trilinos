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

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Tpetra_Version.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_RTI.hpp>

#include <iostream>
#include <iomanip>
#include <functional>
#include <utility>

/** \file RTIExample.cpp
    \brief An example of the Tpetra Vector/MultiVector Reduction-Transformation Interface.
 */

namespace TpetraExamples {

  /// \class Pair
  /// \tparam T1 Type of the first element in the pair.
  /// \tparam T2 Type of the second element in the pair.
  ///
  /// This class exists because std::pair does not have a volatile
  /// overload of operator=.  Kokkos::Compat::KokkosDeviceWrapperNode
  /// needs that.  Otherwise, this class looks almost exactly like
  /// std::pair.  Note the use of makePair instead of std::make_pair
  /// (see declaration and definition of makePair below).
  template<class T1, class T2>
  struct Pair {
    T1 first;
    T2 second;

    //! Default constructor (invokes default constructors of T1 and T2).
    Pair () {}

    //! The usual constructor (pass in a T1 and a T2).
    Pair (const T1& theFirst, const T2& theSecond) :
      first (theFirst), second (theSecond)
    {}

    //! Copy constructor (invokes copy constructors of T1 and T2).
    Pair (const Pair& p) : first (p.first), second (p.second) {}

    //! Copy constructor that takes a volatile input.
    Pair (const volatile Pair& p) : first (p.first), second (p.second) {}

    /// \brief Nonvolatile assignment operator.
    ///
    /// Invokes nonvolatile assignment operators of T1 and T2.
    Pair& operator= (const Pair& p) {
      first = p.first;
      second = p.second;
      return *this;
    }

    /// \brief Volatile assignment operator that takes a volatile input.
    ///
    /// This operator, and the other volatile overloads in this class,
    /// is why this class exists, and why we can't use std::pair
    /// (because it lacks those volatile overloads).  Invokes volatile
    /// assignment operators of T1 and T2.  Built-in types like int
    /// and double have those.  User-defined T1 and T2 might not, in
    /// which case the compiler will give an error if you try to
    /// compile this class.
    volatile Pair& operator= (const volatile Pair& p) volatile {
      first = p.first;
      second = p.second;
      return *this;
    }

    //! Volatile assignment operator that takes a nonvolatile input.
    volatile Pair& operator= (const Pair& p) volatile {
      first = p.first;
      second = p.second;
      return *this;
    }
  };

  //! Equivalent of std::make_pair for the Pair class above.
  template<class T1, class T2>
  Pair<T1, T2> makePair (const T1& first, const T2& second) {
    return Pair<T1, T2> (first, second);
  }

  /** \class mprec_mult
   * \brief Function object multiplying two different types in a third type's precision.
   *
   * mprec_mult's operator() converts its two arguments, of types Arg1
   * resp. Arg2, into MultPrec objects.  It then uses MultPrec's
   * operator* to multiply them, and returns the result as a MultPrec
   * object.
   *
   * \tparam Arg1 Type of the operator's first argument.
   * \tparam Arg2 Type of the operator's second argument.
   * \tparam MultPrec Type in which the multiplication is performed,
   *   and the return type.
   */
  template <class Arg1, class Arg2, class MultPrec>
  class mprec_mult : public std::binary_function<Arg1,Arg2,MultPrec> {
  public:
    mprec_mult() {}
    inline MultPrec
    operator() (const Arg1& x, const Arg2& y) const {
      return (MultPrec)(x) * (MultPrec)(y);
    }
  };

  /** \class pair_op
   * \brief Reduction function object that combines pairs of T1, T2.
   *
   * A pair_op takes a binary operator "op" of type Op, whose
   * arguments may be of types T1 or T2.  It uses op to combine the
   * elements of two pairs of (T1,T2) pairwise: the first element of
   * the first pair with the first element of the second pair, and the
   * second element of the first pair with the second element of the
   * second pair.
   *
   * \tparam T1 Type of the first element in the input pairs.
   * \tparam T2 Type of the second element in the input pairs.
   * \tparam Op A binary operator for which operator()(T1,T1) and
   *   operator()(T2,T2) are syntactically correct.
   */
  template <class T1, class T2, class Op>
  class pair_op :
    public std::binary_function<Pair<T1,T2>, Pair<T1,T2>, Pair<T1,T2> > {
  private:
    Op op_;
  public:
    pair_op (Op op) : op_ (op) {}

    Pair<T1,T2> operator () (const Pair<T1,T2>& a, const Pair<T1,T2>& b) const {
      return makePair (op_ (a.first, b.first), op_ (a.second, b.second));
    }
  };

  /// \fn make_pair_op
  /// \brief Nonmember constructor for \c pair_op.
  template <class T1, class T2, class Op>
  pair_op<T1,T2,Op> make_pair_op (Op op) {
    return pair_op<T1,T2,Op> (op);
  }

} // namespace TpetraExamples


// Partial specialization of Teuchos::SerializationTraits for the Pair
// class declared above.  This teaches Teuchos::Comm how to send and
// receive Pair instances.
//
// FIXME (mfh 22 Feb 2014) This only works if T1 and T2 are directly
// serializable.
namespace Teuchos {
  template<typename Ordinal, typename T1, typename T2>
  class SerializationTraits<Ordinal, ::TpetraExamples::Pair<T1, T2> >
    : public DirectSerializationTraits<Ordinal, ::TpetraExamples::Pair<T1, T2> >
  {};
} // namespace Teuchos


namespace Tpetra {
  namespace RTI {
    //
    // Specialization of ZeroOp for a pair of two different types.
    //
    // Note: "<::" is not a legal character combination for beginning
    // a template argument list, even with C++11.  This is why we
    // leave a space between the less-than sign and the double colon.
    // We use the initial double colon to tell the compiler that we
    // specifically want the TpetraExamples namespace declared above,
    // and not some other TpetraExamples namespace (if one should
    // exist) nested in the Tpetra namespace.
    template <class T1, class T2>
    class ZeroOp< ::TpetraExamples::Pair<T1,T2> > {
    public:
      static inline ::TpetraExamples::Pair<T1,T2> identity () {
        return ::TpetraExamples::makePair (Teuchos::ScalarTraits<T1>::zero (),
                                           Teuchos::ScalarTraits<T2>::zero ());
      }
    };
  }
}

int main(int argc, char *argv[])
{
  //
  // Set up MPI, if applicable.
  //
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  using Teuchos::TimeMonitor;
  using Tpetra::RTI::ZeroOp;
  using Tpetra::RTI::reductionGlob;
  using TpetraExamples::mprec_mult;
  using TpetraExamples::makePair;
  using TpetraExamples::Pair;

  //
  // Get the default communicator
  //
  auto comm = Tpetra::getDefaultComm();
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
    std::cout << std::endl;
  }


  //
  // Create a simple map with 5 local entries per process
  //
  Tpetra::global_size_t numGlobalRows = numGlobal_user;
  auto map = Tpetra::createUniformContigMapWithNode<int,int>(numGlobalRows, comm);
  const size_t numLocalRows = map->getNodeNumElements();
  auto x = Tpetra::createVector<float>(map),
       y = Tpetra::createVector<float>(map);
  auto z = Tpetra::createVector<double>(map),
       w = Tpetra::createVector<double>(map),
       v = Tpetra::createVector<double>(map);

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
    std::cout << std::left << "dot( ones, twos ) result == " << fresult
              << ", expected value is " << numGlobalRows*2.0f << std::endl;
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
  // Two simultaneous dot products: z'*w and z'*v at the same time, saves a pass through z
  //
  Tpetra::RTI::unary_transform( *z, [](double xi){return 1.0f;} );
  Tpetra::RTI::unary_transform( *w, [](double xi){return 2.0f;} );
  Tpetra::RTI::unary_transform( *v, [](double xi){return 3.0f;} );
  auto timeTwoDotInd = TimeMonitor::getNewTimer("dot(z,w) and dot(z,v) independently");
  {
    TimeMonitor lcltimer(*timeTwoDotInd);
    for (int l=0; l != numTimeTrials; ++l)
    {
       dresult = Tpetra::RTI::reduce(*z, *w, reductionGlob< ZeroOp<double>>(
                                                            std::multiplies<double>(),
                                                            std::plus<double>()) );
       dresult = Tpetra::RTI::reduce(*z, *v, reductionGlob< ZeroOp<double>>(
                                                            std::multiplies<double>(),
                                                            std::plus<double>()) );
    }
  }
  auto timeTwoDotSim = TimeMonitor::getNewTimer("dot(z,w) and dot(z,v) simultaneously");
  {
    TimeMonitor lcltimer(*timeTwoDotSim);
    Pair<double, double> tdres;
    for (int l=0; l != numTimeTrials; ++l) {
      tdres = Tpetra::RTI::reduce (*z, *w, *v,
                                   reductionGlob<ZeroOp<Pair<double,double> > > ([] (double zi, double wi, double vi) {
                                       return makePair (zi*wi, zi*vi);
                                     },
                                     TpetraExamples::make_pair_op<double,double> (std::plus<double> ()))
                                  );
    }
  }

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

