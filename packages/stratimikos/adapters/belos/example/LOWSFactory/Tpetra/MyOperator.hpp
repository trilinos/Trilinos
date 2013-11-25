/*
// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef MY_OPERATOR_HPP
#define MY_OPERATOR_HPP

#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Teuchos_BLAS.hpp"

//! Simple example of a user's defined Tpetra::Operator class.
/*! 
 * This is a simple, single processor example of user's defined
 * Tpetra::Operator-derived class. The class is templated
 * on OrdinalType and ScalarType; possible choices are, 
 * for example, "float", "double", or "std::complex<double>".
 *
 * \warning This class is almost certainly broken;
 *   Tpetra::Operator has a different interface now.
 */
template <class OrdinalType, class ScalarType>
class MyOperator : public Tpetra::Operator<OrdinalType,ScalarType> 
{

public:

  //! Constructor
  MyOperator(const Tpetra::VectorSpace<OrdinalType,ScalarType>& vs,
	     const int nrows, const int *colptr,
	     const int nnz, const int *rowin, const ScalarType *vals)
    : _vs(vs), _nr(nrows), _nnz(nnz), _cptr(nrows+1), _rind(nnz), _vals(nnz)
  {
    std::copy<const int*,IntIter>(colptr,colptr+nrows+1,_cptr.begin());
    std::copy<const int*,IntIter>(rowin,rowin+nnz,_rind.begin());
    std::copy<const ScalarType*,STIter>(vals,vals+nnz,_vals.begin());
  }
  
  //! Deconstructor
  ~MyOperator()
  {
  }

  /** \name Functions Overridden from Tpetra::Operator. */
  //@{
  
  //! Returns the VectorSpace associated with the domain of this linear operator.
  Tpetra::VectorSpace<OrdinalType,ScalarType> const& getDomainDist() const { return _vs; };
  
  //! Returns the VectorSpace associated with the range of this linear operator.
  Tpetra::VectorSpace<OrdinalType,ScalarType> const& getRangeDist() const { return _vs; };
  
  //! Computes the matrix-std::vector multiplication y = Ax.
  void apply(Tpetra::Vector<OrdinalType,ScalarType> const& x, 
	     Tpetra::Vector<OrdinalType, ScalarType> & y, 
	     bool transpose=false) const 
  {
    // Get the indexes of the rows on this processor
    const int numMyElements = _vs.getNumMyEntries();
    const std::vector<int> &myGlobalElements = _vs.elementSpace().getMyGlobalElements();
    
    // Initialize output std::vector to zero.
    y.setAllToScalar( Teuchos::ScalarTraits<ScalarType>::zero() );

    assert (x.getNumGlobalEntries() == y.getNumGlobalEntries());
    assert (x.getNumGlobalEntries() == y.getNumGlobalEntries());
    
    // Apply operator
    int IA1, IA2, ri;
    int i,j;
    for (int myRow = 0; myRow < numMyElements; ++myRow ) {

      // For each row this processor owns, compute the value of A*x[myRow]
      const int rowIndex = myGlobalElements[myRow];
      for (j=0; j<_nr; j++) {		
	IA1 = _cptr[j]-1;
	IA2 = _cptr[j+1]-1;
	for (i=IA1; i<IA2; i++) {
	  ri = _rind[i]-1;
	  if (ri == rowIndex) 
	    y[rowIndex] += _vals[i]*x[j];
	} // end for (i= ...)
      } // end for (j= ...)
    } // end for (myRow= ...)
    
  };
  
  //@}
  
private:

  typedef typename std::vector<ScalarType>::iterator STIter;
  typedef std::vector<int>::iterator        IntIter;

  //! Tpetra std::vector space 
  Tpetra::VectorSpace<OrdinalType,ScalarType> _vs;

  //! Number of rows and columns
  int _nr, _nnz;
  //! Column pointers 
  std::vector<int> _cptr;
  //! Row indices
  std::vector<int> _rind;
  //! Values
  std::vector<ScalarType> _vals;
};

#endif //MY_OPERATOR_HPP
