// Kris
// 06.18.03 -- Removed comments/documentation; file too hard to edit otherwise. Will replace later.
//          -- Begin conversion from <ScalarType> template to <OrdinalType, ScalarType>
// 06.23.03 -- Finished conversion from <ScalarType> to <OrdinalType, ScalarType>
//          -- Tpetra_DenseMatrix.cpp is now obsolete
//          -- Added new constructor to allow construction of a submatrix
//          -- Altered copyMat to enable its use in new constructor
//          -- Commented out broken print() function
//          -- Fixed oneNorm() (uninitialized return variable was causing erroneous results)
// 06.24.03 -- Minor formatting changes
// 07.01.03 -- Added TempPrint() function to temporarily take the place of print() and operator<< while I figure out how to fix them
// 07.02.03 -- Added operator== and operator!= to make testing programs easier to write/read. Implementation of == isn't the most
//             efficient/robust, but it works. Will consider optimizing later.
//          -- Warning! Constructor DenseMatrix(DataAccess, const DenseMatrix<OrdinalType, ScalarType> &, int, int, int, int) (the
//             "submatrix grabber" constructor) does not work correctly when used with CV == View (always grabs submatrix from top
//             left corner).
// 07.07.03 -- Constructor bug detailed above (07.02) is now corrected (hopefully).
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_DENSEMATRIX_HPP_
#define _TEUCHOS_DENSEMATRIX_HPP_

#include "Teuchos_CompObject.cpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_DataAccess.hpp"

namespace Teuchos
{
  template<typename OrdinalType, typename ScalarType>
  class DenseMatrix : public CompObject, public Object, public BLAS<OrdinalType, ScalarType>
  {
  public:
    DenseMatrix();
    DenseMatrix(DataAccess CV, ScalarType* A, int stride, int numRows, int numCols);
    DenseMatrix(const DenseMatrix<OrdinalType, ScalarType> &Source);
    DenseMatrix(DataAccess CV, const DenseMatrix<OrdinalType, ScalarType> &Source, int rows, int cols, int startRow=0, int startCol=0);
    virtual ~DenseMatrix();
    int shape(int numRows, int numCols);
    int reshape(int numRows, int numCols);
    typename ScalarTraits<ScalarType>::magnitudeType DenseMatrix<OrdinalType, ScalarType>::oneNorm();
    int  multiply (char TransA, char TransB, ScalarType alpha, const DenseMatrix<OrdinalType, ScalarType> &A, const DenseMatrix<OrdinalType, ScalarType> &B, ScalarType beta);
    bool operator== (const DenseMatrix<OrdinalType, ScalarType> &Operand);
    bool operator!= (const DenseMatrix<OrdinalType, ScalarType> &Operand);
    ScalarType &operator () (int RowIndex, int ColIndex);
    const ScalarType &operator () (int RowIndex, int ColIndex) const;
    ScalarType *operator [] (int ColIndex);
    const ScalarType *operator [] (int ColIndex) const;
    int numRows() const { return(numRows_); };
    int numCols() const { return(numCols_); };
    ScalarType* values() const { return(values_); };
    int stride() const { return(stride_); };
    // virtual void print(ostream& os) const;

    void TempPrint();

  protected:
    void copyMat(ScalarType* inputMatrix, int strideInput, int numRows, int numCols, ScalarType* outputMatrix, int strideOutput, int startRow, int startCol);
    void deleteArrays();
    int numRows_;
    int numCols_;
    int stride_;
    bool valuesCopied_;
    ScalarType* values_;
  }; // class Teuchos_DenseMatrix

template<typename OrdinalType, typename ScalarType>
DenseMatrix<OrdinalType, ScalarType>::DenseMatrix() : CompObject(), numRows_(0), numCols_(0), stride_(0), valuesCopied_(false), values_(0) { }

template<typename OrdinalType, typename ScalarType>
DenseMatrix<OrdinalType, ScalarType>::DenseMatrix(DataAccess CV, ScalarType* A, int stride, int numRows, int numCols) : CompObject(), numRows_(numRows), numCols_(numCols), stride_(stride), valuesCopied_(false), values_(A)
{
  if(CV == Copy)
  {
    stride_ = numRows_;
    values_ = new ScalarType[stride_*numCols_];
    copyMat(A, stride, numRows_, numCols_, values_, stride_, 0, 0);
    valuesCopied_ = true;
  }
}

template<typename OrdinalType, typename ScalarType>
DenseMatrix<OrdinalType, ScalarType>::DenseMatrix(const DenseMatrix<OrdinalType, ScalarType> &Source) : CompObject(), numRows_(Source.numRows_), numCols_(Source.numCols_), stride_(Source.stride_), valuesCopied_(true), values_(Source.values_)
{
  stride_ = numRows_;
  values_ = new ScalarType[stride_*numCols_];
  copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_, 0, 0);
}

template<typename OrdinalType, typename ScalarType>
DenseMatrix<OrdinalType, ScalarType>::DenseMatrix(DataAccess CV, const DenseMatrix<OrdinalType, ScalarType> &Source, int rows, int cols, int startRow, int startCol) : CompObject(), numRows_(rows), numCols_(cols), stride_(Source.stride_), valuesCopied_(false), values_(Source.values_)
{
  if(CV == Copy)
    {
      stride_ = rows;
      values_ = new ScalarType[stride_ * cols];
      copyMat(Source.values_, Source.stride_, rows, cols, values_, stride_, startRow, startCol);
      valuesCopied_ = true;
    }
  else // CV == View
    {
      values_ = values_ + (stride_ * startRow) + startCol;
    }
}

template<typename OrdinalType, typename ScalarType>
int DenseMatrix<OrdinalType, ScalarType>::reshape(int numRows, int numCols)
{
  // Allocate space for new matrix
  ScalarType* values_tmp = new ScalarType[numRows * numCols];
  ScalarType zero = ScalarTraits<ScalarType>::zero();
  for(int k = 0; k < numRows * numCols; k++)
    {
      values_tmp[k] = zero;
    }
  int numRows_tmp = TEUCHOS_MIN(numRows_, numRows);
  int numCols_tmp = TEUCHOS_MIN(numCols_, numCols);
  if(values_ != 0)
    {
      copyMat(values_, stride_, numRows_tmp, numCols_tmp, values_tmp, numRows, 0, 0); // Copy principal submatrix of A to new A
    }
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows;
  numCols_ = numCols;
  stride_ = numRows_;
  values_ = values_tmp; // Set pointer to new A
  valuesCopied_ = true;
  return(0);
}

template<typename OrdinalType, typename ScalarType>
int DenseMatrix<OrdinalType, ScalarType>::shape(int numRows, int numCols)
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows;
  numCols_ = numCols;
  stride_ = numRows_;
  ScalarType zero = ScalarTraits<ScalarType>::zero();
  values_ = new ScalarType[stride_*numCols_];
  for(int k = 0; k < stride_ * numCols_; k++)
    {
      values_[k] = zero; // Zero out values
    }
  valuesCopied_ = true;
  return(0);
}

template<typename OrdinalType, typename ScalarType>
DenseMatrix<OrdinalType, ScalarType>::~DenseMatrix()
{
  deleteArrays();
}

template<typename OrdinalType, typename ScalarType>
void DenseMatrix<OrdinalType, ScalarType>::deleteArrays(void)
{
  if (valuesCopied_)
  {
    delete [] values_;
    values_ = 0;
    valuesCopied_ = false;
  }
}

template<typename OrdinalType, typename ScalarType>
void DenseMatrix<OrdinalType, ScalarType>::copyMat(ScalarType* inputMatrix, int strideInput, int numRows, int numCols, ScalarType* outputMatrix, int strideOutput, int startRow, int startCol)
{
  int i, j;
  ScalarType* ptr1;
  ScalarType* ptr2;
  for(j = 0; j < numCols; j++)
  {
    ptr1 = outputMatrix + (j * strideOutput);
    ptr2 = inputMatrix + (j + startCol) * strideInput + startRow;
    for(i = 0; i < numRows; i++)
      {
	*ptr1++ = *ptr2++;
      }
  }
}

template<typename OrdinalType, typename ScalarType>
typename ScalarTraits<ScalarType>::magnitudeType DenseMatrix<OrdinalType, ScalarType>::oneNorm()
{
  int i, j;
  typename ScalarTraits<ScalarType>::magnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  typename ScalarTraits<ScalarType>::magnitudeType absSum = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  ScalarType* ptr;
  for(j = 0; j < numCols_; j++)
    {
      ScalarType sum = 0;
      ptr = values_ + j * stride_;
      for(i = 0; i < numRows_; i++)
	{
	  sum += ScalarTraits<ScalarType>::magnitude(*ptr++);
	}
      absSum = ScalarTraits<ScalarType>::magnitude(sum);
      if(absSum > anorm)
	{
	  anorm = absSum;
	}
    }
  updateFlops(numCols_ * numCols_);
  return(anorm);
}

template<typename OrdinalType, typename ScalarType>
bool DenseMatrix<OrdinalType, ScalarType>::operator== (const DenseMatrix<OrdinalType, ScalarType> &Operand)
{
  bool result = 1;
  if((numRows_ != Operand.numRows_) || (numCols_ != Operand.numCols_))
    {
      result = 0;
    }
  else
    {
      int i, j;
      for(i = 0; i < numRows_; i++)
	{
	  for(j = 0; j < numCols_; j++)
	    {
	      if((*this)(i, j) != Operand(i, j))
		{
		  return 0;
		}
	    }
	}
    }
  return result;
}

template<typename OrdinalType, typename ScalarType>
bool DenseMatrix<OrdinalType, ScalarType>::operator!= (const DenseMatrix<OrdinalType, ScalarType> &Operand)
{
  return !((*this) == Operand);
}

template<typename OrdinalType, typename ScalarType>
ScalarType& DenseMatrix<OrdinalType, ScalarType>::operator () (int RowIndex, int ColIndex)
{
  if (RowIndex >= numRows_)
  {
    cout << "Row index = " << RowIndex << " Out of Range 0 - " << numRows_-1 << endl;
    TEUCHOS_CHK_REF(*values_); // Return reference to values_[0]
  }
  if (ColIndex >= numCols_)
  {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << numCols_-1 << endl;
    TEUCHOS_CHK_REF(*values_); // Return reference to values_[0]
  }
  return(values_[ColIndex * stride_ + RowIndex]);
}

template<typename OrdinalType, typename ScalarType>
const ScalarType& DenseMatrix<OrdinalType, ScalarType>::operator () (int RowIndex, int ColIndex) const
{
  if (RowIndex >= numRows_)
  {
    cout << "Row index = " << RowIndex << " Out of Range 0 - " << numRows_ - 1 << endl;
    TEUCHOS_CHK_REF(values_[0]); // Return reference to values_[0]
  }
  if (ColIndex >= numCols_)
  {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << numCols_ - 1 << endl;
    TEUCHOS_CHK_REF(values_[0]); // Return reference to values_[0]
  }
  return(values_[ColIndex * stride_ + RowIndex]);
}

template<typename OrdinalType, typename ScalarType>
const ScalarType* DenseMatrix<OrdinalType, ScalarType>::operator [] (int ColIndex) const
{
  if (ColIndex >= numCols_)
  {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << numCols_ - 1 << endl;
    TEUCHOS_CHK_PTR(0); // Return zero pointer
  }
  return(values_ + ColIndex * stride_);
}

template<typename OrdinalType, typename ScalarType>
ScalarType* DenseMatrix<OrdinalType, ScalarType>::operator [] (int ColIndex)
{
  if (ColIndex >= numCols_)
  {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << numCols_ - 1 << endl;
    TEUCHOS_CHK_PTR(0); // Return zero pointer
  }
  return(values_ + ColIndex * stride_);
}

template<typename OrdinalType, typename ScalarType>
int  DenseMatrix<OrdinalType, ScalarType>::multiply(char TransA, char TransB, ScalarType alpha, const DenseMatrix<OrdinalType, ScalarType> &A, const DenseMatrix<OrdinalType, ScalarType> &B, ScalarType beta)
{
  // Check for compatible dimensions
  int A_nrows = (TransA=='T') ? A.numCols() : A.numRows();
  int A_ncols = (TransA=='T') ? A.numRows() : A.numCols();
  int B_nrows = (TransB=='T') ? B.numCols() : B.numRows();
  int B_ncols = (TransB=='T') ? B.numRows() : B.numCols();
  if ((numRows_ != A_nrows) || (A_ncols != B_nrows) || (numCols_ != B_ncols))
    {
      TEUCHOS_CHK_ERR(-1); // Return error
    }
  // Call GEMM function
  GEMM(TransA, TransB, numRows_, numCols_, A_ncols, alpha, A.values(), A.stride(), B.values(), B.stride(), beta, values_, stride_);
  double nflops = 2 * numRows_;
  nflops *= numCols_;
  nflops *= A_ncols;
  updateFlops(nflops);
  return(0);
}

template<typename OrdinalType, typename ScalarType>
void DenseMatrix<OrdinalType, ScalarType>::TempPrint()
{
  int i, j;
  for(i = 0; i < numRows_; i++)
    {
      for(j = 0; j < numCols_; j++)
	{
	  cout << (*this)(i, j) << " ";
	}
      cout << endl;
    }
  cout << endl;
}

// template<typename OrdinalType, typename ScalarType>
// inline ostream& operator<<(ostream& os, const Teuchos::DenseMatrix<OrdinalType, ScalarType>& obj)
// {
//   // os << obj.Label();
//   //  obj.print(os);
//   // return os;
// }

} // namespace Teuchos

// #include "Teuchos_DenseMatrix.cpp"

#endif /* _TEUCHOS_DENSEMATRIX_HPP_ */
