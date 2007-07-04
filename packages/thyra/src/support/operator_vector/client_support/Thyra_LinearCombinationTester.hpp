// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_LINEARCOMBINATIONTESTER_HPP
#define THYRA_LINEARCOMBINATIONTESTER_HPP

#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_LinearCombinationImpl.hpp"
#include "Thyra_VectorImpl.hpp"
#include "Thyra_TestSpecifier.hpp"
#include "Thyra_TesterBase.hpp"
#include "Teuchos_ScalarTraits.hpp"

#define TESTER(form1, form2)\
  {\
    Vector<Scalar> _val1 = form1;\
    TEST_FOR_EXCEPT(_val1.constPtr().get()==0);\
    Vector<Scalar> _val2 = form2;\
    TEST_FOR_EXCEPT(_val2.constPtr().get()==0);\
    ScalarMag err = norm2(_val1-_val2);\
    if (!checkTest(spec_, err, "[" #form1 "] == [" #form2 "]")) pass = false;\
  }

namespace Thyra
{
  using Teuchos::RCP;
  using Teuchos::ScalarTraits;

  /** */
  template <class Scalar>
  class LinearCombinationTester : public TesterBase<Scalar>
  {
  public:

    /** */
    LinearCombinationTester(const RCP<const Comm<int> >& comm,
                            const VectorSpace<Scalar>& vecSpace, 
                            Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const TestSpecifier<Scalar>& spec);

    /** */
    bool runAllTests() const ;


  private:

    /** */
    bool nonModifyingOpTests() const ;

    /** */
    bool selfModifyingOpTests() const ;

    TestSpecifier<Scalar> spec_;

  };

  template <class Scalar> 
  inline LinearCombinationTester<Scalar>
  ::LinearCombinationTester(const RCP<const Comm<int> >& comm,
                            const VectorSpace<Scalar>& vecSpace,
                            Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const TestSpecifier<Scalar>& spec)
    : TesterBase<Scalar>(comm, vecSpace, vecSpace.dim(), out), spec_(spec)
  {;}

  template <class Scalar> 
  inline bool LinearCombinationTester<Scalar>
  ::runAllTests() const
  {
    bool pass = true;

    pass = this->nonModifyingOpTests() && pass;
    pass = this->selfModifyingOpTests() && pass;

    return pass;
  }

  template <class Scalar> 
  inline bool LinearCombinationTester<Scalar>
  ::nonModifyingOpTests() const
  {
    bool pass = true;
    typedef typename Teuchos::ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType ScalarMag;

    LinearOperator<Scalar> A = this->randomDenseOp();
    LinearOperator<Scalar> B = this->randomDenseOp();
    LinearOperator<Scalar> C = this->randomDenseOp();

    Vector<Scalar> x = A.domain().createMember();
    Vector<Scalar> y = A.domain().createMember();
    Vector<Scalar> z = A.domain().createMember();

    randomizeVec(x);
    randomizeVec(y);
    randomizeVec(z);

    this->out() << "starting linear combination tests" << std::endl;

    TESTER(x*Scalar(2.0), Scalar(2.0)*x);

    TESTER(Scalar(2.0)*(x + y), Scalar(2.0)*x + Scalar(2.0)*y);

    TESTER(Scalar(2.0)*(x - y), Scalar(2.0)*x - Scalar(2.0)*y);

    TESTER((x + y)*Scalar(2.0), Scalar(2.0)*x + Scalar(2.0)*y);

    TESTER((x - y)*Scalar(2.0), Scalar(2.0)*x - Scalar(2.0)*y);

    TESTER(Scalar(2.0)*(x - y), -Scalar(2.0)*(y - x));

    TESTER(Scalar(0.25)*(Scalar(2.0)*(x + y) - Scalar(2.0)*(x - y)), y);

    TESTER((Scalar(2.0)*A)*x, Scalar(2.0)*(A*x));

    TESTER(Scalar(2.0)*(A*x), (A*x)*Scalar(2.0));

    TESTER(A*(B*x), (A*B)*x);

    TESTER(Scalar(2.0)*A*(B*x), A*(B*(Scalar(2.0)*x)));

    TESTER(Scalar(3.0)*(Scalar(2.0)*A)*x, Scalar(6.0)*(A*x));

    TESTER(A*x + y, y + A*x);

    TESTER(z + (A*x + B*y), (B*y + A*x) + z);

    TESTER(z - (A*x + B*y), Scalar(-1.0)*((B*y + A*x) - z));

    TESTER(C*z + (A*x + B*y), (B*y + A*x) + C*z);

    TESTER(C*z - (A*x + B*y), Scalar(-1.0)*((B*y + A*x) - C*z));

    TESTER(Scalar(2.0)*z + (A*x + B*y), (B*y + A*x) + Scalar(2.0)*z);

    TESTER(Scalar(2.0)*z - (A*x + B*y), Scalar(-1.0)*((B*y + A*x) - Scalar(2.0)*z));

    TESTER(A*x - y, Scalar(-1.0)*(y - A*x));

    TESTER(A*x + B*y, B*y + A*x);

    TESTER(A*x - B*y - A*x + B*y +z, z);

    TESTER(Scalar(2.0)*(A*x + y), Scalar(2.0)*A*x + Scalar(2.0)*y);

    TESTER(Scalar(2.0)*(A*x + B*y), A*x + B*y + A*x + B*y);

    TESTER(Scalar(2.0)*(y + A*x), Scalar(2.0)*y + Scalar(2.0)*(A*x));

    TESTER(x + Scalar(2.0)*A*y, x + Scalar(2.0)*(A*y));

    TESTER(Scalar(2.0)*A*y + x, Scalar(2.0)*(A*y) + x);

    TESTER(Scalar(2.0)*(A*x + B*y), Scalar(2.0)*A*x + Scalar(2.0)*B*y);

    TESTER(Scalar(2.0)*(A*x - B*y), Scalar(2.0)*A*x - Scalar(2.0)*B*y);

    TESTER(Scalar(2.0)*(A*x + B*y + z), Scalar(2.0)*A*x + Scalar(2.0)*B*y + Scalar(2.0)*z);

    TESTER(Scalar(2.0)*(A*x + Scalar(3.0)*B*y), Scalar(2.0)*A*x + Scalar(6.0)*B*y);

    TESTER(Scalar(2.0)*(A*x + Scalar(3.0)*(z + B*y)), Scalar(2.0)*A*x + Scalar(6.0)*B*y + Scalar(6.0)*z);

    TESTER(Scalar(2.0)*(z + A*x + B*y + z), Scalar(2.0)*A*x + Scalar(2.0)*B*y + Scalar(4.0)*z);

    TESTER(Scalar(2.0)*(Scalar(3.0)*(z + A*x) + B*y), Scalar(6.0)*z + Scalar(6.0)*A*x + Scalar(2.0)*B*y);

    TESTER(Scalar(2.0)*(Scalar(3.0)*(z + A*x) + Scalar(4.0)*(B*y + z)), Scalar(6.0)*z + Scalar(6.0)*A*x + Scalar(8.0)*B*y + Scalar(8.0)*z);
    
    TESTER((A*x + B*y) + (A*y + B*x), (A + B)*x + (A+B)*y);
    TESTER((A*x + B*y) - (A*y + B*x), A*x - A*y + B*y - B*x);

    TESTER((A*x + B*y) + Scalar(2.0)*(A*y + B*x), A*(x + Scalar(2.0)*y) + B*(Scalar(2.0)*x + y));
    TESTER((A*x + B*y) - Scalar(2.0)*(A*y + B*x), A*(x - Scalar(2.0)*y) + B*(y - Scalar(2.0)*x));

    return pass;
  }

  template <class Scalar> 
  inline bool LinearCombinationTester<Scalar>
  ::selfModifyingOpTests() const
  {
    bool pass = true;
    typedef typename Teuchos::ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType ScalarMag;

    LinearOperator<Scalar> A = this->randomDenseOp();
    LinearOperator<Scalar> B = this->randomDenseOp();
    LinearOperator<Scalar> C = this->randomDenseOp();

    Vector<Scalar> x = A.domain().createMember();
    Vector<Scalar> y = A.domain().createMember();
    Vector<Scalar> z = A.domain().createMember();

    randomizeVec(x);
    randomizeVec(y);
    randomizeVec(z);

    Vector<Scalar> a = copy(x);
    Vector<Scalar> b = copy(y);
    Vector<Scalar> c = copy(z);
    
    TestSpecifier<Scalar>
      specLooser(spec_.doTest(),1e+1*spec_.warningTol(),1e+2*spec_.errorTol());

    this->out() << "starting linear combination tests" << std::endl;

    x = Scalar(2.0)*A*x;
    ScalarMag err = norm2(x - Scalar(2.0)*A*a);
    if (!checkTest(spec_, err, "x=Scalar(2.0)*A*x")) pass = false;

    a = copy(x);
    x = x + y;
    err = norm2(x - (a + y));
    if (!checkTest(spec_, err, "x=x+y")) pass = false;

    a = copy(x);
    x = y + x;
    err = norm2(x - (y + a));
    if (!checkTest(spec_, err, "x=y+x")) pass = false;

    a = copy(x);
    x = A*x + B*y;
    err = norm2(x - (A*a + B*y));
    if (!checkTest(spec_, err, "x=A*x+B*y")) pass = false;

    a = copy(x);
    x = B*y + A*x;
    err = norm2(x - (B*y + A*a));
    if (!checkTest(spec_, err, "x=B*y+A*x")) pass = false;

    a = copy(x);
    x = A*x + (B*y + C*x);
    err = norm2(x - (A*a + (B*y + C*a)));
    if (!checkTest(specLooser, err, "x=A*x + (B*y + C*x)")) pass = false;

    a = copy(x);
    x = (A*x + B*y) + C*x;
    err = norm2(x - ((A*a + B*y) + C*a));
    if (!checkTest(specLooser, err, "x=(A*x + B*y) + C*x")) pass = false;

    /* test assignment of OpTimesLC into empty and non-empty vectors */
    Vector<Scalar> u;
    u = Scalar(2.0)*A*B*x;
    err = norm2(u - Scalar(2.0)*A*B*x);

    u = Scalar(2.0)*A*x;
    err = norm2(u - Scalar(2.0)*A*B*x);

    /* test assignment of LC2 into empty and non-empty vectors */
    Vector<Scalar> v;
    v = Scalar(2.0)*x + Scalar(3.0)*y;
    err = norm2(v - (Scalar(2.0)*x + Scalar(3.0)*y));

    v = Scalar(2.0)*x + Scalar(3.0)*y;
    err = norm2(v - (Scalar(2.0)*x + Scalar(3.0)*y));

    return pass;
  }
 
}

// We had better not leave a global macro named TESTSER just lying around
// awating disaster!
#undef TESTER

#endif
