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
 * for example, "float", "double", or "complex<double>".
 *
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
  
  //! Computes the matrix-vector multiplication y = Ax.
  void apply(Tpetra::Vector<OrdinalType,ScalarType> const& x, 
	     Tpetra::Vector<OrdinalType, ScalarType> & y, 
	     bool transpose=false) const 
  {
    // Get the indexes of the rows on this processor
    const int numMyElements = _vs.getNumMyEntries();
    const std::vector<int> &myGlobalElements = _vs.elementSpace().getMyGlobalElements();
    
    // Initialize output vector to zero.
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

  //! Tpetra vector space 
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
