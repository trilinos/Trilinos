/*Paul
27-May-2002 General cleanup. Changed method names to fit namingConvention (already done). Moved error-reporting macros to Tpetra_Object.h
06-August-2002 Changed to images (nothing changed).
*/

namespace Tpetra
{
//=============================================================================
template<class scalarType>
DenseMatrix<scalarType>::DenseMatrix(void)
  : CompObject(),
    numRows_(0),
    numCols_(0),
    stride_(0),
    valuesCopied_(false),
    values_(0)
{
}

//=============================================================================
template<class scalarType>
DenseMatrix<scalarType>::DenseMatrix(Tpetra_DataAccess CV, scalarType *A, int stride, 
					     int numRows, int numCols)
  : CompObject(),
    numRows_(numRows),
    numCols_(numCols),
    stride_(stride),
    valuesCopied_(false),    
    values_(A)

{
  if (CV==Copy)
  {
    stride_ = numRows_;
    values_ = new scalarType[stride_*numCols_];
    CopyMat(A, stride, numRows_, numCols_, values_, stride_);
    valuesCopied_ = true;
  }

}
//=============================================================================
template<class scalarType>
DenseMatrix<scalarType>::DenseMatrix(const DenseMatrix<scalarType>& Source)
  : CompObject(),  
    numRows_(Source.numRows_),
    numCols_(Source.numCols_),
    stride_(Source.stride_),
    valuesCopied_(true),
    values_(Source.values_)

{

  int i;
  stride_ = numRows_;
  values_ = new scalarType[stride_*numCols_];
  copyMat(Source.values_, Source.stride_, numRows_, numCols_, values_, stride_);
}
//=============================================================================
template<class scalarType>
int DenseMatrix<scalarType>::reshape(int numRows, int numCols)
{

  // Allocate space for new matrix
  scalarType * values_tmp = new scalarType[numRows*numCols];
  for (int k = 0; k < numRows*numCols; k++) values_tmp[k] = 0.0; // Zero out values
  int numRows_tmp = TPETRA_MIN(numRows_, numRows);
  int numCols_tmp = TPETRA_MIN(numCols_, numCols);
  if (values_ != 0) CopyMat(values_, stride_, numRows_tmp, numCols_tmp, values_tmp, numRows); // Copy principal submatrix of A to new A
  
  DeleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows;
  numCols_ = numCols;
  stride_ = numRows_;
  values_ = values_tmp; // Set pointer to new A
  valuesCopied_ = true;

  return(0);
}
//=============================================================================
template<class scalarType>
int DenseMatrix<scalarType>::shape(int numRows, int numCols)
{
  deleteArrays(); // Get rid of anything that might be already allocated
  numRows_ = numRows;
  numCols_ = numCols;
  stride_ = numRows_;
  scalarType zero = 0; // NOTE:  This should be handled by the Traits mechanism
  values_ = new scalarType[stride_*numCols_];
  for (int k = 0; k < stride_*numCols_; k++) values_[k] = zero; // Zero out values
  valuesCopied_ = true;

  return(0);
}
//=============================================================================
template<class scalarType>
DenseMatrix<scalarType>::~DenseMatrix()
{
  deleteArrays();
}
//=============================================================================
template<class scalarType>
void DenseMatrix<scalarType>::deleteArrays(void)
{
  if (valuesCopied_)
  {
    delete [] values_;
    values_ = 0;
    valuesCopied_ = false;
  }
}
//=============================================================================
template<class scalarType>
void DenseMatrix<scalarType>::copyMat(scalarType * inputMatrix, int strideInput, int numRows, int numCols, 
					scalarType * outputMatrix, int strideOutput)
{

  int i, j;
  scalarType * ptr1;
  scalarType * ptr2;
  for (j=0; j<numCols; j++)
  {
    ptr1 = outputMatrix + j*strideOutput;
    ptr2 = inputMatrix + j*strideInput;
    for (i=0; i<numRows; i++)
      *ptr1++ = *ptr2++;
  }
  return;
}
//=============================================================================
template<class scalarType>
ScalarTraits<scalarType>::magnitudeType DenseMatrix<scalarType>::oneNorm(void)
{
  int i, j;

  ScalarTraits<scalarType>::magnitudeType anorm = ScalarTraits<scalarType>::magnitude(ScalarTraits<scalarType>::zero());
  ScalarTraits<scalarType>::magnitudeType absSum;
  scalarType * ptr;
  for (j=0; j<numCols_; j++)
  {
    scalarType sum=0;
    ptr = values_ + j*stride_;
    for (i=0; i<numRows_; i++)
      sum += ScalarTraits<scalarType>::magnitude(*ptr++);
    absSum = ScalarTraits<scalarType>::magnitude(sum);
    if (absSum>anorm)
      anorm = absSum;
  }
  updateFlops(numCols_*numCols_);
  return(anorm);
}
//=========================================================================
template<class scalarType>
scalarType& DenseMatrix<scalarType>::operator () (int RowIndex, int ColIndex)
{
  if (RowIndex>=numRows_)
  {
    cout << "Row index = " << RowIndex << " Out of Range 0 - " << numRows_-1 << endl;
    TPETRA_CHK_REF(*values_); // Return reference to values_[0]
  }
  if (ColIndex>=numCols_)
  {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << numCols_-1 << endl;
    TPETRA_CHK_REF(*values_); // Return reference to values_[0]
  }
  return(values_[ColIndex*stride_ + RowIndex]);
}

//=========================================================================
template<class scalarType>
const scalarType& DenseMatrix<scalarType>::operator () (int RowIndex, int ColIndex) const
{
  if (RowIndex>=numRows_)
  {
    cout << "Row index = " << RowIndex << " Out of Range 0 - " << numRows_-1 << endl;
    TPETRA_CHK_REF(values_[0]); // Return reference to values_[0]
  }
  if (ColIndex>=numCols_)
  {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << numCols_-1 << endl;
    TPETRA_CHK_REF(values_[0]); // Return reference to values_[0]
  }
  return(values_[ColIndex*stride_ + RowIndex]);
}
//=========================================================================
template<class scalarType>
const scalarType* DenseMatrix<scalarType>::operator [] (int ColIndex) const
{
  if (ColIndex>=numCols_)
  {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << numCols_-1 << endl;
    TPETRA_CHK_PTR(0); // Return zero pointer
  }
  return(values_ + ColIndex*stride_);
}
//=========================================================================
template<class scalarType>
scalarType* DenseMatrix<scalarType>::operator [] (int ColIndex)
{
  if (ColIndex>=numCols_)
  {
    cout << "Column index = " << ColIndex << " Out of Range 0 - " << numCols_-1 << endl;
    TPETRA_CHK_PTR(0); // Return zero pointer
  }
  return(values_+ ColIndex*stride_);
}
//=========================================================================
template<class scalarType>
int  DenseMatrix<scalarType>::multiply (char TransA, char TransB, scalarType ScalarAB, 
				      const DenseMatrix<scalarType>& A, 
				      const DenseMatrix<scalarType>& B,
				      scalarType Scalar )
{
  // Check for compatible dimensions
  
  int A_nrows = (TransA=='T') ? A.numCols() : A.numRows();
  int A_ncols = (TransA=='T') ? A.numRows() : A.numCols();
  int B_nrows = (TransB=='T') ? B.numCols() : B.numRows();
  int B_ncols = (TransB=='T') ? B.numRows() : B.numCols();
  
  if (numRows_  != A_nrows     ||
      A_ncols   != B_nrows     ||
      numCols_  != B_ncols  ) TPETRA_CHK_ERR(-1); // Return error

    
  // Call GEMM function
  GEMM(TransA, TransB, numRows_, numCols_, A_ncols, ScalarAB, A.values(), A.stride(), 
       B.values(), B.stride(), Scalar, values_, stride_);
  double nflops = 2*numRows_;
  nflops *= numCols_;
  nflops *= A_ncols;
  updateFlops(nflops);

  return(0);
}

template<class scalarType>
void DenseMatrix<scalarType>::print(ostream& os) const
{
  for (int i=0; i<numRows_; i++)
  {
    for (int j=0; j<numCols_; j++)
    {
      os << (*this)(i,j) << "\t";
    }
    os << endl;
  }
  return;
}

}  // namespace Tpetra
