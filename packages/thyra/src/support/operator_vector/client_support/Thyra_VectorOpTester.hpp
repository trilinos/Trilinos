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

#ifndef THYRA_VECTOROPTESTER_HPP
#define THYRA_VECTOROPTESTER_HPP

#include "Thyra_VectorImpl.hpp"
#include "Thyra_VectorSpaceImpl.hpp"
#include "Thyra_TesterBase.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_LinearCombinationImpl.hpp"

namespace Thyra
{
using namespace Teuchos;
using std::ostream;

/** 
 * Run comparisons between element-wise calculations and Vector member
 * functions.
 */
template <class Scalar>
class VectorOpTester : public TesterBase<Scalar>
{
public:
  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** */
  VectorOpTester(const RCP<const Comm<Ordinal> >& comm,
    const VectorSpace<Scalar>& space,
    Teuchos::RCP<Teuchos::FancyOStream>& out,
    const TestSpecifier<Scalar>& spec);

  /** */
  bool runAllTests() const ;

  /** */
  bool sumTest() const ;

  /** */
  bool setElementUsingBracketTest() const ;

  /** */
  bool dotStarTest() const ;

  /** */
  bool scalarMultTest() const ;

  /** */
  bool overloadedUpdateTest() const ;

private:
  TestSpecifier<Scalar> spec_;
};

template <class Scalar> 
inline VectorOpTester<Scalar>
::VectorOpTester(const RCP<const Comm<Ordinal> >& comm,
  const VectorSpace<Scalar>& space,
  Teuchos::RCP<Teuchos::FancyOStream>& out,
  const TestSpecifier<Scalar>& spec)
  : TesterBase<Scalar>(comm, space, 1, out),
    spec_(spec)
{;}

template <class Scalar> 
inline bool VectorOpTester<Scalar>
::runAllTests() const
{
  bool pass = true;

  pass = setElementUsingBracketTest() && pass;
  pass = sumTest() && pass;
  pass = dotStarTest() && pass;
  pass = scalarMultTest() && pass;
  pass = overloadedUpdateTest() && pass;

  return pass;
}

template <class Scalar> 
inline bool VectorOpTester<Scalar>
::sumTest() const 
{

  using std::endl;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  ScalarMag err = ST::magnitude(ST::zero());
    
  this->out() << "running vector addition test..." << endl;

  Vector<Scalar> a = this->space().createMember();
  Vector<Scalar> b = this->space().createMember();
  Vector<Scalar> x = this->space().createMember();
  Vector<Scalar> y = this->space().createMember();
  randomizeVec(a);
  randomizeVec(b);
    
  /* do the operation with member functions */
  x = a + b ;
    
  /* do the operation elementwise. For off-processor elements, 
   * this is a no-op */
  for (int i=0; i<this->space().dim(); i++)
  {
    Scalar a_i = a[i];
    Scalar b_i = b[i];
    y[i] = a_i + b_i;
  }
  err = normInf(x-y);
  return this->checkTest(spec_, err, "vector addition");
}

 


template <class Scalar> 
inline bool VectorOpTester<Scalar>
::setElementUsingBracketTest() const 
{

  using std::endl;

  typedef Teuchos::ScalarTraits<Ordinal> IST;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  /* this test requires a comparison operator to make sense */
  if (ST::isComparable && spec_.doTest())
  {
    ScalarMag err = ST::magnitude(ST::zero());

    Vector<Scalar> a = this->space().createMember();

    /* We will load a vector with zeros, then set a randomly 
     * selected element to 1.0 and 
     * another to -1.0. Run minloc and maxloc to test that the minimum and
     * maximum values are at the expected locations. */

    zeroOut(a);

    int N = this->space().dim();
    int nTrials = 50;

    /* pick a seed in a way that is deterministic across processors */
    unsigned int seed = 107*N % 101;
    ST::seedrandom( seed );

    for (int t=0; t<nTrials; t++)
    {
      zeroOut(a);
      Ordinal minIndex = IST::random() % N;
      Ordinal maxIndex = IST::random() % N;
      /* skip cases where we've generated identical indices */
      if (minIndex == maxIndex) continue;
      a[minIndex] = -ST::one();
      a[maxIndex] = ST::one();
      Ordinal findMin = -1;
      Ordinal findMax = -1;
      Scalar minVal = Thyra::minloc(a, findMin);
      Scalar maxVal = Thyra::maxloc(a, findMax);
      err += ST::magnitude(-ST::one() - minVal);
      err += ST::magnitude(ST::one() - maxVal);
      err += ST::magnitude(findMin - minIndex);
      err += ST::magnitude(findMax - maxIndex);
    }

    return this->checkTest(spec_, err, "bracket operator");
  }
  else
  {
    this->out() << "skipping bracket operator test..." << endl;
    return true;
  }
}

/* specialize the set element test for complex types to a no-op. 
 * This is because minloc and maxloc do not exist for complex vectors 
 */
#ifdef HAVE_THYRA_COMPLEX


template <> 
inline bool VectorOpTester<std::complex<double> >
::setElementUsingBracketTest() const 
{
  using std::endl;
  this->out() << "skipping vector bracket operator test..." << endl;
  return true;
}


#if defined(HAVE_THYRA_FLOAT)
template <> 
inline bool VectorOpTester<std::complex<float> >
::setElementUsingBracketTest() const 
{
  using std::endl;
  this->out() << "skipping vector bracket operator test..." << endl;
  return true;
}
#endif


#endif // HAVE_THYRA_COMPLEX
  

template <class Scalar> 
inline bool VectorOpTester<Scalar>
::dotStarTest() const 
{

  using std::endl;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  ScalarMag err = ST::magnitude(ST::zero());

  this->out() << "running vector dotStar test..." << endl;
    
  Vector<Scalar> a = this->space().createMember();
  Vector<Scalar> b = this->space().createMember();
  Vector<Scalar> x = this->space().createMember();
  Vector<Scalar> y = this->space().createMember();
  Vector<Scalar> z = this->space().createMember();
  randomizeVec(a);
  randomizeVec(b);
    
  x = dotStar(a,b);
  z = dotSlash(x, b);
    
  err = normInf(a-z);
    
  this->out() << "|dotStar error|=" << err << endl;
  return this->checkTest(spec_, err, "dot-star");
}



template <class Scalar> 
inline bool VectorOpTester<Scalar>
::scalarMultTest() const 
{

  using std::endl;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  ScalarMag err = ST::magnitude(ST::zero());

  this->out() << "running vector scalarMult test..." << endl;
    
  Vector<Scalar> a = this->space().createMember();
  Vector<Scalar> x = this->space().createMember();
  Vector<Scalar> y = this->space().createMember();
  randomizeVec(a);
    
    
  /* do the operation with member functions */
  Scalar alpha = 3.14;
  x = alpha * a;
    
  /* do the operation elementwise */
    
  for (int i=0; i<this->space().dim(); i++)
  {
    Scalar a_i = a[i];
    y[i]= alpha *a_i;
  }
    
  err = normInf(x-y);
    
  this->out() << "|scalarMult error|=" << err << endl;
  return this->checkTest(spec_, err, "scalar multiplication");

}
 
template <class Scalar> 
inline bool VectorOpTester<Scalar>
::overloadedUpdateTest() const 
{

  using std::endl;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  ScalarMag err = ST::magnitude(ST::zero());

  this->out() << "running vector overloadedUpdate test..." << endl;
    
  Vector<Scalar> a = this->space().createMember();
  Vector<Scalar> b = this->space().createMember();
  Vector<Scalar> x = this->space().createMember();
  Vector<Scalar> y = this->space().createMember();
  randomizeVec(a);
  randomizeVec(b);
    
  /* do the operation with member functions */
  Scalar alpha = 3.14;
  Scalar beta = 1.4;
  x = alpha*a + beta*b;
    
  /* do the operation elementwise */
  for (int i=0; i<this->space().dim(); i++)
  {
    Scalar a_i = a[i];
    Scalar b_i = b[i];
    y[i] = alpha*a_i + beta*b_i;
  }
    
  err = normInf(x-y);
    
  this->out() << "|overloadedUpdate error|=" << err << endl;
  return this->checkTest(spec_, err, "overloaded update");

}


} // namespace Thyra

#endif
