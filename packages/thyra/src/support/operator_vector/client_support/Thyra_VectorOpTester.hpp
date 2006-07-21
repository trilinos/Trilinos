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
  using Thyra::TestSpecifier;

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
    bool setElementTest() const ;

    /** */
    bool setElementUsingBracketTest() const ;

    /** */
    bool dotStarTest() const ;

    /** */
    bool dotSlashTest() const ;

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

  };

  template <class Scalar> 
  inline VectorOpTester<Scalar>
  ::VectorOpTester(const VectorSpace<Scalar>& space,
                   const TestSpecifier<Scalar>& spec,
                   const Teuchos::MPIComm& comm)
    : spec_(spec), space_(space), comm_(comm)
  {;}

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::runAllTests() const
  {
    bool pass = true;

    pass = sumTest() && pass;
    pass = setElementTest() && pass;
#ifdef BRACKET_TEST
    pass = setElementUsingBracketTest() && pass;
#endif
    pass = dotStarTest() && pass;
    pass = dotSlashTest() && pass;
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
    if (spec_.doTest())
      {
        cout << "running vector addition test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);

        /* do the operation with member functions */
        x = a + b ;

        /* do the operation with member functions */

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            double a_i = ((ConstVector<Scalar>) a)[i];
            double b_i = ((ConstVector<Scalar>) b)[i];
            y[i] = a_i + b_i;
          }

        double err = normInf(x-y);

        cout << "|sum error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector sum test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector sum test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector addition test..." << endl;
      }
    cout << "vector addition test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::setElementTest() const 
  {
    if (spec_.doTest())
      {
        cout << "running setElement test..." << endl;

        Vector<Scalar> a = space_.createMember();
	
        /* we will load a vector with a_i = i, and then do
         * the sum of all elements. If done correctly, the sum will equal 
         * N*(N+1)*(2N+1)/6.
         */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            a[i]=i;
          }
        Vector<double> b = copy(a);

        b = dotStar(a, b);

        double sum = 0.0;
        for (int i=low; i<high; i++)
          {
            sum += i * a[i];
          }


#ifdef HAVE_MPI
        Scalar localSum = sum;
        MPI_Allreduce( (void*) &localSum, (void*) &sum, 
                       1, MPI_DOUBLE, MPI_SUM, comm_.getComm());
#endif
        double thyraSum = Thyra::sum(*(b.rawPtr()));

        double err = ::fabs(sum - thyraSum);

        cout << "|setElement error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector setElement test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector setElement test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector setElement test..." << endl;
      }
    cout << "vector setElement test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }


  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::setElementUsingBracketTest() const 
  {
    if (spec_.doTest())
      {
        cout << "running setElementUsingBracket test..." << endl;

        VectorSpace<Scalar> prodSp = 
          productSpace(tuple(space_, space_));
        Vector<Scalar> prod = prodSp.createMember();
        Vector<Scalar> a = prod.getBlock(0);
        Vector<Scalar> ab = prod.getBlock(1);
	
        /* we will load a vector with a_i = i, and then do
         * the sum of all elements. If done correctly, the sum will equal 
         * N*(N+1)*(2N+1)/6.
         */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            prod[i] = i;
            prod[i + space_.dim()] = i + space_.dim();
          }

        cout << "a = " << endl << a << endl;
        cout << "ab = " << endl << ab << endl;
        cout << "prod = " << endl << prod.getBlock(0) << prod.getBlock(1) << endl;
        Vector<double> b = copy(a);
        Vector<double> prodB = copy(prod);
        b = dotStar(a,b);
        prodB = dotStar(prodB, prod);

        double sum = 0.0;
        double sumP = 0.0;
        for (int i=low; i<high; i++)
          {
            cout << i << " " << prod[i] << " " << i*prod[i]
                 << endl;
            cout << i << " " << a[i] << " " << i*a[i]
                 << endl;
            cout << i << " " << prod[i] << " " << i*prod[i] << " "
                 << prod[i + space_.dim()] << endl;
            //sum += i * a[i];
            sum += i * a[i];
            sumP += i * prod[i] + (i) * prod[i + space_.dim()];
          }

#ifdef HAVE_MPI
        Scalar localSum = sum;
        MPI_Allreduce( (void*) &localSum, (void*) &sum, 
                       1, MPI_DOUBLE, MPI_SUM, comm_.getComm());
        Scalar localSumP = sumP;
        MPI_Allreduce( (void*) &localSumP, (void*) &sumP, 
                       1, MPI_DOUBLE, MPI_SUM, comm_.getComm());
#endif
	
        double thyraSum = Thyra::sum(*(b.ptr()));
        cout << "elemwise sum = " << sum << endl;
        cout << "thyra sum = " << thyraSum << endl;
        double thyraSumP = Thyra::sum(*(prodB.ptr()));
        cout << "elemwise sum = " << sumP << endl;
        cout << "thyra sum = " << thyraSumP << endl;

        double err = ::fabs(sum - thyraSum);
        double errP = ::fabs(sumP - thyraSumP);

        cout << "|setElement error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector setElement test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector setElementUsingBracket test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }

        cout << "|setElement errorP|=" << errP << endl;
        if (errP > spec_.errorTol())
          {
            cout << "product vector setElement test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (errP > spec_.warningTol())
          {
            cout << "WARNING: product vector setElementUsingBracket test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector setElementUsingBracket test..." << endl;
      }
    cout << "vector setElementUsingBracket test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }


  

  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::dotStarTest() const 
  {
    if (spec_.doTest())
      {
        cout << "running vector dotStar test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = dotStar(a,b);

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            double a_i = a[i];
            double b_i = b[i];
            y[i] = a_i * b_i;
          }
	
        double err = normInf(x-y);

        cout << "|dotStar error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector dotStar test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector dotStar test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector dotStar test..." << endl;
      }
    cout << "vector dotStar test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }


  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::dotSlashTest() const 
  {
    if (spec_.doTest())
      {
        cout << "running vector dotSlash test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = dotSlash(a,b);

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            double a_i = a[i];
            double b_i = b[i];
            y[i]=  a_i / b_i;
          }
	
        double err = normInf(x-y);

        cout << "|dotSlash error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector dotSlash test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector dotSlash test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector dotSlash test..." << endl;
      }
    cout << "vector dotSlash test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }

  
  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::scalarMultTest() const 
  {
    if (spec_.doTest())
      {
        cout << "running vector scalarMult test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);


        /* do the operation with member functions */
        x = 3.14*a;

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            double a_i = a[i];
            y[i]= 3.14*a_i;
          }
	
        double err = normInf(x-y);

        cout << "|scalarMult error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector scalarMult test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector scalarMult test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector scalarMult test..." << endl;
      }
    cout << "vector scalarMult test PASSED: tol = " 
         << spec_.errorTol() << endl;
    return true;
  }
 
  template <class Scalar> 
  inline bool VectorOpTester<Scalar>
  ::overloadedUpdateTest() const 
  {
    if (spec_.doTest())
      {
        cout << "running vector overloadedUpdate test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = 3.14*a + 1.4*b;

        /* do the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        for (int i=low; i<high; i++)
          {
            double a_i = a[i];
            double b_i = b[i];
            y[i] = 3.14*a_i + 1.4*b_i;
          }
	
        double err = normInf(x-y);

        cout << "|overloadedUpdate error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector overloadedUpdate test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector overloadedUpdate test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector overloadedUpdate test..." << endl;
      }
    cout << "vector overloadedUpdate test PASSED: tol = " 
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
        cout << "running vector reciprocal test..." << endl;

        Vector<Scalar> a = space_.createMember();
        randomizeVec(a);

        Vector<Scalar> y = space_.createMember();

        /* load the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        int denomsAreOK = true;
        for (int i=low; i<high; i++)
          {
            double a_i = a[i];
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
        double err = (x - y).norm2();

#ifdef HAVE_MPI
        int localDenomsAreOK = denomsAreOK;
        MPI_Allreduce( (void*) &localDenomsAreOK, (void*) &denomsAreOK, 
                       1, MPI_INT, MPI_LAND, comm_.getComm());
#endif

        cout << "|reciprocal error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector reciprocal test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (tDenomsAreOK != denomsAreOK)
          {
            cout << "vector reciprocal test FAILED: trilinosDenomsOK="
                 << tDenomsAreOK << ", denomsOK=" << denomsAreOK << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector reciprocal test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        cout << "skipping vector reciprocal test..." << endl;
      }
    cout << "vector reciprocal test PASSED: tol = " 
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
        cout << "running vector minQuotient test..." << endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();

        randomizeVec(a);
        randomizeVec(b);

        /* perform the operation elementwise */
        int low = lowestLocallyOwnedIndex(space_);
        int high = low + numLocalElements(space_);

        double minQLocal = Teuchos::ScalarTraits<Scalar>::rmax();
        for (int i=low; i<high; i++)
          {
            double a_i = a[i];
            double b_i = b[i];
            if (b_i != Teuchos::ScalarTraits<Scalar>::zero())
              {
                double q = a_i/b_i;
                if (q < minQLocal) minQLocal = q;
              }
          }

        double minQ = minQLocal;
        comm_.allReduce((void*) &minQLocal, (void*) &minQ, 1, Teuchos::MPIComm::DOUBLE,
                        Teuchos::MPIComm::MIN);
	

        double tMinQ = Thyra::VMinQuotient(*(a.ptr()), *(b.ptr()));
        cout << "trilinos minQ = " << tMinQ << endl;
        cout << "elemwise minQ = " << minQ << endl;
        double err = fabs(minQ - tMinQ);
        
        cout << "min quotient error=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "min quotient test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: min quotient test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        cout << "skipping min quotient test..." << endl;
      }
    cout << "min quotient test PASSED: tol = " 
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
        cout << "running vector constraintMask test..." << endl;

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
            double a_i = a[i];
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
        double err = (m - y).norm2();

#ifdef HAVE_MPI
        int localAllFeas = allFeasible;
        cout << "local all feas=" << localAllFeas << endl;
        MPI_Allreduce( (void*) &localAllFeas, (void*) &allFeasible, 
                       1, MPI_INT, MPI_LAND, comm_.getComm());
        cout << "globalal all feas=" << allFeasible << endl;
#endif

        cout << "|constraintMask error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector constraintMask test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (allFeasible != tAllFeasible)
          {
            cout << "vector constraintMask test FAILED: trilinosFeas="
                 << tAllFeasible << ", feas=" << allFeasible << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector constraintMask test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
      }
    else
      {
        cout << "skipping vector constraintMask test..." << endl;
      }
    cout << "vector constraintMask test PASSED: tol = " 
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
        cout << "running vector compare-to-scalar test..." << endl;

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
            double a_i = a[i];
            y[i] = fabs(a_i) >= s ;
          }
	
        double err = normInf(x-y);

        cout << "|compare-to-scalar error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector compare-to-scalar test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector compare-to-scalar test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector compare-to-scalar test..." << endl;
      }
    cout << "vector compare-to-scalar test PASSED: tol = " 
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
        cout << "running vector index test..." << endl;
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
            double a_i = a[i];
            y[i] =  fabs(a_i) >= s;
          }
	
        double err = normInf(x-y);

        cout << "|index error|=" << err << endl;
        if (err > spec_.errorTol())
          {
            cout << "vector index test FAILED: tol = " 
                 << spec_.errorTol() << endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            cout << "WARNING: vector index test could not beat tol = " 
                 << spec_.warningTol() << endl;
          }
	
      }
    else
      {
        cout << "skipping vector index test..." << endl;
      }
    cout << "vector index test PASSED: tol = " 
         << spec_.errorTol() << endl;
#endif
    return true;
  }
  


}
#endif
