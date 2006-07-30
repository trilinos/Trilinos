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
//#include "Thyra_SUNDIALS_Ops.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_MPIComm.hpp"
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
    VectorOpTester(const VectorSpace<Scalar>& space,
                   Teuchos::RefCountPtr<Teuchos::FancyOStream>& out,
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

    /** */
    bool reciprocalTest() const ;

    /** */
    bool minQuotientTest() const ;

    /** */
    bool addScalarTest() const ;

    /** */
    bool compareToScalarTest() const ;

    /** */
    bool constraintMaskTest() const ;

    /** */
    bool indexTest() const ;

  private:
    TestSpecifier<Scalar> spec_;
  };

  template <class Scalar> 
  inline VectorOpTester<Scalar>
  ::VectorOpTester(const VectorSpace<Scalar>& space,
                   Teuchos::RefCountPtr<Teuchos::FancyOStream>& out,
                   const TestSpecifier<Scalar>& spec)
    : TesterBase<Scalar>(space, 1, out),
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
    pass = reciprocalTest() && pass;
    pass = minQuotientTest() && pass;
    pass = constraintMaskTest() && pass;
    pass = compareToScalarTest() && pass;
    pass = indexTest() && pass;

    return pass;
  }

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::sumTest() const 
  {

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
    typedef Teuchos::ScalarTraits<Index> IST;
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
            Index minIndex = IST::random() % N;
            Index maxIndex = IST::random() % N;
            /* skip cases where we've generated identical indices */
            if (minIndex == maxIndex) continue;
            a[minIndex] = -ST::one();
            a[maxIndex] = ST::one();
            Index findMin = -1;
            Index findMax = -1;
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
#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
  template <> 
  inline bool VectorOpTester<std::complex<double> >
  ::setElementUsingBracketTest() const 
  {
    this->out() << "skipping vector bracket operator test..." << endl;
    return true;
  }

  template <> 
  inline bool VectorOpTester<std::complex<float> >
  ::setElementUsingBracketTest() const 
  {
    this->out() << "skipping vector bracket operator test..." << endl;
    return true;
  }
#endif


  

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::dotStarTest() const 
  {
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

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::reciprocalTest() const 
  {
#ifdef TRILINOS_6
    typedef Teuchos::ScalarTraits<Scalar> ST;
    ScalarMag err = ST::magnitude(ST::zero());

    if (spec_.doTest())
      {
        this->out() << "running vector reciprocal test..." << endl;

        Vector<Scalar> a = this->space().createMember();
        randomizeVec(a);

        Vector<Scalar> y = this->space().createMember();

        /* load the operation elementwise */
        int low = lowestLocallyOwnedIndex(this->space());
        int high = low + numLocalElements(this->space());

        int denomsAreOK = true;
        for (int i=low; i<high; i++)
          {
            Scalar a_i = a[i];
            if (a_i != Teuchos::ScalarTraits<Scalar>::zero()) 
              {
                y[i] = Teuchos::ScalarTraits<Scalar>::one()/a_i;
              }
            else
              {
                denomsAreOK=false;
                y[i] = a_i;
              }
          }
        
        Vector<Scalar> x = this->space().createMember();
        int tDenomsAreOK = Thyra::VInvTest(*(a.constPtr()), x.ptr().get());
        err = (x - y).norm2();

#ifdef HAVE_MPI
        int localDenomsAreOK = denomsAreOK;
        MPI_Allreduce( (void*) &localDenomsAreOK, (void*) &denomsAreOK, 
                       1, MPI_INT, MPI_LAND, comm_.getComm());
#endif

        this->out() << "|reciprocal error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            this->out() << "vector reciprocal test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (tDenomsAreOK != denomsAreOK)
          {
            this->out() << "vector reciprocal test FAILED: trilinosDenomsOK="
                 << tDenomsAreOK << ", denomsOK=" << denomsAreOK << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            this->out() << "WARNING: vector reciprocal test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        this->out() << "skipping vector reciprocal test..." << endl;
      }
    this->out() << "vector reciprocal test PASSED: error=" << err << ", tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::minQuotientTest() const 
  {
#ifdef TRILINOS_6
    typedef Teuchos::ScalarTraits<Scalar> ST;
    ScalarMag err = ST::magnitude(ST::zero());

    if (spec_.doTest())
      {
        this->out() << "running vector minQuotient test..." << endl;

        Vector<Scalar> a = this->space().createMember();
        Vector<Scalar> b = this->space().createMember();

        randomizeVec(a);
        randomizeVec(b);

        /* perform the operation elementwise */
        int low = lowestLocallyOwnedIndex(this->space());
        int high = low + numLocalElements(this->space());

        Scalar minQLocal = Teuchos::ScalarTraits<Scalar>::rmax();
        for (int i=low; i<high; i++)
          {
            Scalar a_i = a[i];
            Scalar b_i = b[i];
            if (b_i != Teuchos::ScalarTraits<Scalar>::zero())
              {
                Scalar q = a_i/b_i;
                if (q < minQLocal) minQLocal = q;
              }
          }

        Scalar minQ = minQLocal;
        comm_.allReduce((void*) &minQLocal, (void*) &minQ, 1, Teuchos::MPIComm::DOUBLE,
                        Teuchos::MPIComm::MIN);
	

        Scalar tMinQ = Thyra::VMinQuotient(*(a.constPtr()), *(b.constPtr()));
        this->out() << "trilinos minQ = " << tMinQ << endl;
        this->out() << "elemwise minQ = " << minQ << endl;
        err = ST::magnitude(minQ - tMinQ);
        
        this->out() << "min quotient error=" << err << endl;
        if (err > spec_.errorTol())
          {
            this->out() << "min quotient test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            this->out() << "WARNING: min quotient test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        this->out() << "skipping min quotient test..." << endl;
      }
    this->out() << "min quotient test PASSED: error=" << err << ", tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }



  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::constraintMaskTest() const 
  {
#ifdef TRILINOS_6
    typedef Teuchos::ScalarTraits<Scalar> ST;
    ScalarMag err = ST::magnitude(ST::zero());

    if (spec_.doTest())
      {
        this->out() << "running vector constraintMask test..." << endl;

        Vector<Scalar> a = this->space().createMember();
        Vector<Scalar> c = this->space().createMember();
        randomizeVec(a);

        Vector<Scalar> y = this->space().createMember();
        Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

        /* load the operation elementwise */
        int low = lowestLocallyOwnedIndex(this->space());
        int high = low + numLocalElements(this->space());

        int allFeasible = true;
        for (int i=low; i<high; i++)
          {
            int feasible = true;
            Scalar a_i = a[i];
            switch(i%4)
              {
              case 0:
                c[i] = -2.0;
                feasible = a_i < zero;
                break;
              case 1:
                c[i] = -1.0;
                feasible = a_i <= zero;
                break;
              case 2:
                c[i] = 1.0;
                feasible = a_i > zero;
                break;
              case 3:
                c[i] = 2.0;
                feasible = a_i >= zero;
                break;
              default:
                TEST_FOR_EXCEPTION(true, logic_error, "impossible!");
              }
            y[i] = (Scalar) !feasible;
            allFeasible = allFeasible && feasible;
          }
	
        Vector<Scalar> m = this->space().createMember();
        int tAllFeasible = Thyra::VConstrMask(*(a.constPtr()), *(c.constPtr()), m.ptr().get());
        err = (m - y).norm2();

#ifdef HAVE_MPI
        int localAllFeas = allFeasible;
        this->out() << "local all feas=" << localAllFeas << endl;
        MPI_Allreduce( (void*) &localAllFeas, (void*) &allFeasible, 
                       1, MPI_INT, MPI_LAND, comm_.getComm());
        this->out() << "globalal all feas=" << allFeasible << endl;
#endif

        this->out() << "|constraintMask error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            this->out() << "vector constraintMask test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (allFeasible != tAllFeasible)
          {
            this->out() << "vector constraintMask test FAILED: trilinosFeas="
                 << tAllFeasible << ", feas=" << allFeasible << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            this->out() << "WARNING: vector constraintMask test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        this->out() << "skipping vector constraintMask test..." << endl;
      }
    this->out() << "vector constraintMask test PASSED: error=" << err << ", tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }
  

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::compareToScalarTest() const 
  {
#ifdef TRILINOS_6
    typedef Teuchos::ScalarTraits<Scalar> ST;
    ScalarMag err = ST::magnitude(ST::zero());

    if (spec_.doTest())
      {
        this->out() << "running vector compare-to-scalar test..." << endl;

        Vector<Scalar> a = this->space().createMember();
        Vector<Scalar> x = this->space().createMember();
        Vector<Scalar> y = this->space().createMember();
        randomizeVec(a);

        /* do the operation with member functions */
        Scalar s = 0.5;
        Thyra::VCompare(s, *(a.constPtr()), x.ptr().get());

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(this->space());
        int high = low + numLocalElements(this->space());

        for (int i=low; i<high; i++)
          {
            Scalar a_i = a[i];
            y[i] = fabs(a_i) >= s ;
          }
	
        err = normInf(x-y);

        this->out() << "|compare-to-scalar error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            this->out() << "vector compare-to-scalar test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            this->out() << "WARNING: vector compare-to-scalar test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        this->out() << "skipping vector compare-to-scalar test..." << endl;
      }
    this->out() << "vector compare-to-scalar test PASSED: error=" << err << ", tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }


  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::indexTest() const 
  {
#ifdef TRILINOS_6
    typedef Teuchos::ScalarTraits<Scalar> ST;
    ScalarMag err = ST::magnitude(ST::zero());

    if (spec_.doTest())
      {
        this->out() << "running vector index test..." << endl;
        Vector<Scalar> a = this->space().createMember();
        Vector<Scalar> x = this->space().createMember();
        Vector<Scalar> y = this->space().createMember();
        randomizeVec(a);
        /* do the operation with member functions */
        Scalar s = 0.5;
        Thyra::VCompare(s, *(a.constPtr()), x.ptr().get());

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(this->space());
        int high = low + numLocalElements(this->space());

        for (int i=low; i<high; i++)
          {
            Scalar a_i = a[i];
            y[i] =  fabs(a_i) >= s;
          }
	
        err = normInf(x-y);

        this->out() << "|index error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            this->out() << "vector index test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            this->out() << "WARNING: vector index test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        this->out() << "skipping vector index test..." << endl;
      }
    this->out() << "vector index test PASSED: error=" << err << ", tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }
  


}
#endif
