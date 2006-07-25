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
#include "Thyra_TestSpecifier.hpp"
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
  class VectorOpTester
  {
  public:
    /** \brief Local typedef for promoted scalar magnitude */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** */
    VectorOpTester(const VectorSpace<Scalar>& space,
                   const TestSpecifier<Scalar>& spec,
                   const Teuchos::MPIComm& comm = Teuchos::MPIComm::world());

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

    /** */
    void randomizeVec(Vector<Scalar>& x) const ;

    TestSpecifier<Scalar> spec_;

    VectorSpace<Scalar> space_;

    Teuchos::MPIComm comm_;

    Teuchos::RefCountPtr<Teuchos::FancyOStream> out_;

  };

  template <class Scalar> 
  inline VectorOpTester<Scalar>
  ::VectorOpTester(const VectorSpace<Scalar>& space,
                   const TestSpecifier<Scalar>& spec,
                   const Teuchos::MPIComm& comm)
    : spec_(spec), space_(space), comm_(comm), 
      out_(Teuchos::VerboseObjectBase::getDefaultOStream())
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
  inline void VectorOpTester<Scalar>
  ::randomizeVec(Vector<Scalar>& x) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Thyra::randomize(Scalar(-ST::one()),Scalar(+ST::one()),x.ptr().get());
  }

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::sumTest() const 
  {

    typedef Teuchos::ScalarTraits<Scalar> ST;
    
    if (spec_.doTest())
      {
        *out_ << "running vector addition test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);

        /* do the operation with member functions */
        x = a + b ;
        
        /* do the operation elementwise. For off-processor elements, this is a no-op */
        for (int i=0; i<space_.dim(); i++)
          {
            *out_ << "i=" << i << " of N=" << space_.dim() << endl;
            Scalar a_i = a[i];
            Scalar b_i = b[i];
            y[i] = a_i + b_i;
          }

        ScalarMag err = normInf(x-y);

        *out_ << "|sum error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector sum test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector sum test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        *out_ << "skipping vector addition test..." << endl;
      }
    *out_ << "vector addition test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

 


  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::setElementUsingBracketTest() const 
  {
    
    typedef Teuchos::ScalarTraits<int> ST;

    if (spec_.doTest())
      {
        Vector<Scalar> a = space_.createMember();

        /* We will load a vector with zeros, then set a randomly selected element to 1.0 and 
         * another to -1.0. Run minloc and maxloc to test that the minimum and
         * maximum values are at the expected locations. */

        zeroOut(a);

        int N = space_.dim();
        int nTrials = 50;

        /* pick a seed in a way that is deterministic across processors */
        unsigned int seed = 107*N % 101;
        ST::seedrandom( seed );

        ScalarMag err = 0.0;

        for (int t=0; t<nTrials; t++)
          {
            zeroOut(a);
            Index minIndex = ST::random() % N;
            Index maxIndex = ST::random() % N;
            /* skip cases where we've generated identical indices */
            if (minIndex == maxIndex) continue;
            a[minIndex] = -1.0;
            a[maxIndex] = 1.0;
            Index findMin = -1;
            Index findMax = -1;
            Scalar minVal = Thyra::minloc(a, findMin);
            Scalar maxVal = Thyra::maxloc(a, findMax);
            err += ST::magnitude(-1.0 - minVal);
            err += ST::magnitude(1.0 - maxVal);
            err += ST::magnitude(findMin - minIndex);
            err += ST::magnitude(findMax - maxIndex);
          }
        
        *out_ << "|bracket error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector bracket test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector setElementUsingBracket test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }

      }
    else
      {
        *out_ << "skipping vector setElementUsingBracket test..." << endl;
      }
    *out_ << "vector setElementUsingBracket test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }


  

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::dotStarTest() const 
  {
    if (spec_.doTest())
      {
        *out_ << "running vector dotStar test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        Vector<Scalar> z = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);

        x = dotStar(a,b);
        z = dotSlash(x, b);
        
        ScalarMag err = normInf(a-z);

        *out_ << "|dotStar error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector dotStar test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector dotStar test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        *out_ << "skipping vector dotStar test..." << endl;
      }
    *out_ << "vector dotStar test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::scalarMultTest() const 
  {
    if (spec_.doTest())
      {
        *out_ << "running vector scalarMult test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);


        /* do the operation with member functions */
        x = 3.14*a;

        /* do the operation elementwise */

        for (int i=0; i<space_.dim(); i++)
          {
            Scalar a_i = a[i];
            y[i]= 3.14*a_i;
          }
	
        ScalarMag err = normInf(x-y);

        *out_ << "|scalarMult error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector scalarMult test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector scalarMult test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        *out_ << "skipping vector scalarMult test..." << endl;
      }
    *out_ << "vector scalarMult test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }
 
  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::overloadedUpdateTest() const 
  {
    if (spec_.doTest())
      {
        *out_ << "running vector overloadedUpdate test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = 3.14*a + 1.4*b;

        /* do the operation elementwise */
        for (int i=0; i<space_.dim(); i++)
          {
            Scalar a_i = a[i];
            Scalar b_i = b[i];
            y[i] = 3.14*a_i + 1.4*b_i;
          }
	
        ScalarMag err = normInf(x-y);

        *out_ << "|overloadedUpdate error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector overloadedUpdate test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector overloadedUpdate test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        *out_ << "skipping vector overloadedUpdate test..." << endl;
      }
    *out_ << "vector overloadedUpdate test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::reciprocalTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        *out_ << "running vector reciprocal test..." << endl;

        Vector<Scalar> a = space_.createMember();
        randomizeVec(a);

        Vector<Scalar> y = space_.createMember();

        /* load the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        int denomsAreOK = true;
        for (int i=low; i<high; i++)
          {
            Scalar a_i = a[i];
            if (a_i != Teuchos::ScalarTraits<Scalar>::zero()) 
              {
                y[i] = 1.0/a_i;
              }
            else
              {
                denomsAreOK=false;
                y[i] = a_i;
              }
          }
        
        Vector<Scalar> x = space_.createMember();
        int tDenomsAreOK = Thyra::VInvTest(*(a.ptr()), x.ptr().get());
        ScalarMag err = (x - y).norm2();

#ifdef HAVE_MPI
        int localDenomsAreOK = denomsAreOK;
        MPI_Allreduce( (void*) &localDenomsAreOK, (void*) &denomsAreOK, 
                       1, MPI_INT, MPI_LAND, comm_.getComm());
#endif

        *out_ << "|reciprocal error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector reciprocal test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (tDenomsAreOK != denomsAreOK)
          {
            *out_ << "vector reciprocal test FAILED: trilinosDenomsOK="
                 << tDenomsAreOK << ", denomsOK=" << denomsAreOK << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector reciprocal test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        *out_ << "skipping vector reciprocal test..." << endl;
      }
    *out_ << "vector reciprocal test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::minQuotientTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        *out_ << "running vector minQuotient test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();

        randomizeVec(a);
        randomizeVec(b);

        /* perform the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

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
	

        Scalar tMinQ = Thyra::VMinQuotient(*(a.ptr()), *(b.ptr()));
        *out_ << "trilinos minQ = " << tMinQ << endl;
        *out_ << "elemwise minQ = " << minQ << endl;
        ScalarMag err = ST::magnitude(minQ - tMinQ);
        
        *out_ << "min quotient error=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "min quotient test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: min quotient test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        *out_ << "skipping min quotient test..." << endl;
      }
    *out_ << "min quotient test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }



  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::constraintMaskTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        *out_ << "running vector constraintMask test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> c = space_.createMember();
        randomizeVec(a);

        Vector<Scalar> y = space_.createMember();
        Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

        /* load the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

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
	
        Vector<Scalar> m = space_.createMember();
        int tAllFeasible = Thyra::VConstrMask(*(a.ptr()), *(c.ptr()), m.ptr().get());
        ScalarMag err = (m - y).norm2();

#ifdef HAVE_MPI
        int localAllFeas = allFeasible;
        *out_ << "local all feas=" << localAllFeas << endl;
        MPI_Allreduce( (void*) &localAllFeas, (void*) &allFeasible, 
                       1, MPI_INT, MPI_LAND, comm_.getComm());
        *out_ << "globalal all feas=" << allFeasible << endl;
#endif

        *out_ << "|constraintMask error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector constraintMask test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (allFeasible != tAllFeasible)
          {
            *out_ << "vector constraintMask test FAILED: trilinosFeas="
                 << tAllFeasible << ", feas=" << allFeasible << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector constraintMask test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        *out_ << "skipping vector constraintMask test..." << endl;
      }
    *out_ << "vector constraintMask test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }
  

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::compareToScalarTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        *out_ << "running vector compare-to-scalar test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);

        /* do the operation with member functions */
        Scalar s = 0.5;
        Thyra::VCompare(s, *(a.ptr()), x.ptr().get());

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            Scalar a_i = a[i];
            y[i] = fabs(a_i) >= s ;
          }
	
        ScalarMag err = normInf(x-y);

        *out_ << "|compare-to-scalar error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector compare-to-scalar test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector compare-to-scalar test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        *out_ << "skipping vector compare-to-scalar test..." << endl;
      }
    *out_ << "vector compare-to-scalar test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }


  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::indexTest() const 
  {
#ifdef TRILINOS_6
    if (spec_.doTest())
      {
        *out_ << "running vector index test..." << endl;
        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        /* do the operation with member functions */
        Scalar s = 0.5;
        Thyra::VCompare(s, *(a.ptr()), x.ptr().get());

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            Scalar a_i = a[i];
            y[i] =  fabs(a_i) >= s;
          }
	
        ScalarMag err = normInf(x-y);

        *out_ << "|index error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            *out_ << "vector index test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            *out_ << "WARNING: vector index test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        *out_ << "skipping vector index test..." << endl;
      }
    *out_ << "vector index test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }
  


}
#endif
