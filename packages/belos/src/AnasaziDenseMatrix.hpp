// AnasaziDenseMatrix.hpp: interface for the AnasaziDenseMatrix class.
#ifndef ANASAZI_DENSE_MATRIX_HPP
#define ANASAZI_DENSE_MATRIX_HPP

#include "BelosConfigDefs.hpp"

namespace Anasazi {

/*!	\class Anasazi::DenseMatrix

	\brief Anasazi's templated class for constructing Fortran-style dense matrices that
	are used by the eigensolver.

	\author Rich Lehoucq, Heidi Thornquist
*/

template <class TYPE>
class DenseMatrix  {
public:
	//@{ \name Constructors/Destructor.

	//! %Anasazi::DenseMatrix default constructor, creates nondimensional matrix.
	DenseMatrix();

	/*! \brief %Anasazi::DenseMatrix constructor, creates a dense matrix of 
		dimension \c rows by \c cols.  The elements of this matrix are initialized
		to zero.
	*/
	DenseMatrix(int rows, int cols);

	/*! \brief %Anasazi::DenseMatrix constructor, creates a dense matrix of 
		dimension \c rows by \c cols with \c ld leading dimension.  The leading
		dimension \c ld is expected to be larger than, or equal to, \c rows.  The
		elements of this matrix are initialized to zero.
	*/
	DenseMatrix(int rows, int ld, int cols);

	/*! \brief Creates a new %Anasazi::DenseMatrix that is an exact replica of \c A.
	*/
	DenseMatrix(const DenseMatrix<TYPE> & A);

	/*! \brief Creates a new %Anasazi::DenseMatrix of dimension \c rows by \c cols and
		copies the leading \c rows by \c cols subblock of \c A into it.
	*/
	DenseMatrix(const DenseMatrix<TYPE> & A, int rows, int cols);

	/*! \brief Creates a new %Anasazi::DenseMatrix of dimension \c rows by \c cols and
		copies a selected submatrix from \c A into it.  The first entry of this 
		submatrix is determined by \c start_rows and \c start_cols.
	*/
	DenseMatrix(const DenseMatrix<TYPE> & A, int start_rows, 
		int start_cols, int rows, int cols);

	//! %Anasazi::DenseMatrix destructor.
	~DenseMatrix();
	//@}

	//@{ \name Dimension update methods.		

	//! Set the rows of \c *this matrix to \c rows, which is expected to be greater than zero.
	void setrows(const int rows);

	//! Set the columns of \c *this matrix to \c cols, which is expected to be greater than zero.
	void setcols(const int cols);
	//@}

	//@{ \name Dimension information methods.	

	//! Returns the number of rows in \c *this matrix.
	int getrows() const;

	//! Returns the number of columns in \c *this matrix.
	int getcols() const;

	//! Returns the leading dimension of \c *this matrix.
	int getld() const;
	//@}

	//@{ \name Element access methods.

	//! Return the array containing all the elements of \c *this matrix.
	TYPE* getarray() const;
	//@}

	//@{ \name Norm method.

	//! Compute the Frobenius norm of \c *this matrix.	
	TYPE getfronorm() const;
	//@}

	//@{ \name Initialization methods.

	/*! \brief Replace the values of \c *this with the values in \c array.
	*/
	void setvalues(TYPE* array, const int ldx);

	/*! \brief Replace each element of \c *this matrix with \c value.  
		The default is zero.
	*/
	void init( const TYPE value = 0.0 );
	//@}

	//@{ \name Print method.

	//! Print \c *this dense matrix
	void print();
	//@}
private:
	const int _rows, _ld, _cols;
	const bool _created;
	TYPE *_array,*_endp1;
};

template<class TYPE> 
DenseMatrix<TYPE>::DenseMatrix(): _rows(0), _ld(0), _cols(0), 
						_created(false), _array(0), _endp1(0) { 
//	cout << "ctor1:Anasazi::DenseMatrix " << this << endl;
}

template<class TYPE> 
DenseMatrix<TYPE>::DenseMatrix(int rows, int cols): 
							_rows(rows), _ld(rows), _cols(cols), 
						_created(true), _array(0), _endp1(0) { 
//	cout << "ctor2:Anasazi::DenseMatrix " << this << endl;
	if (_ld*_cols > 0) {
		_array = new TYPE[_ld*_cols];
		_endp1 = _array + _ld*_cols;
	}
	assert(_array!=NULL);
	assert(_ld*_rows>0);
	init();  // Initialize array values to zero
}

template<class TYPE> 
DenseMatrix<TYPE>::DenseMatrix(int rows, int ld, int cols): 
							_rows(rows), _ld(ld), _cols(cols), 
						_created(true), _array(0), _endp1(0) { 
//	cout << "ctor3:Anasazi::DenseMatrix " << this << endl;
	if (_ld*_cols > 0 ) {
		_array = new TYPE[_ld*_cols];
		_endp1 = _array + _ld*_cols;
	}
	assert(_array!=NULL);
	assert(_ld*_rows>0);
	assert(_rows <= _ld);
	init();  // Initialize array values to zero
}

//
// Copy constructor
//
template<class TYPE>
DenseMatrix<TYPE>::DenseMatrix(const DenseMatrix& A):
							_rows(A.getrows()), _ld(A.getrows()), _cols(A.getcols()),
							_created(true), _array(0), _endp1(0) {
//	cout << "copy_ctor:Anasazi::DenseMatrix " << this << endl;
	if (_ld*_cols > 0 ) {
		_array = new TYPE[_ld*_cols]; assert(_array!=NULL);
		_endp1 = _array + _ld*_cols;
		if (_array) {
			int i,j;
			TYPE* ptr=A.getarray();
			for (j=0;j<_cols;j++) {
				for (i=0;i<_ld;i++) {
					_array[i+j*_ld] = ptr[i+j*A.getld()];
				}
			}
		}
	assert(_ld*_rows>0);
	assert(_rows <= _ld);
	}
}

// create a view
template<class TYPE>
DenseMatrix<TYPE>::DenseMatrix(const DenseMatrix& A, int rows, int cols):
							_rows(rows), _ld(A.getld()), _cols(cols),
							_created(false), _array(0), _endp1(0) {
//	cout << "view_ctor1:Anasazi::DenseMatrix " << this << endl;
	if (_rows*_cols > 0 && 
		_rows <= A.getrows() && 
		_cols <= A.getcols()) {
		_array = A.getarray(); 
		_endp1 = _array+(_cols-1)*_ld + _rows; // check
	}
	assert(_array!=NULL);
	assert(_rows <= A.getrows());
	assert(_cols <= A.getcols());
	assert(_cols*_rows>0);
	assert(_rows <= _ld);
}

// create a partial view
template<class TYPE>
DenseMatrix<TYPE>::DenseMatrix(const DenseMatrix& A, int start_row, 
											 int start_col, int rows, int cols):
							_rows(rows), _ld(A.getld()), _cols(cols),
							_created(false), _array(0), _endp1(0) {
//	cout << "view_ctor1:Anasazi::DenseMatrix " << this << endl;
	// start_col and start_row are zero based, i.e., rows and columns are indexed
	// from 0.
	if (_rows*_cols > 0 && 
		start_row*start_col >= 0 &&
		_rows <= A.getrows()-start_row && 
		_cols <= A.getcols()-start_col) {
		_array = A.getarray() + start_col*_ld + start_row; 
		_endp1 = _array+(_cols-1)*_ld + _rows; // check
	}
	assert(_array!=NULL);
	assert(start_row*start_col >= 0);
	assert(_rows <= A.getrows()-start_row);
	assert(_cols <= A.getcols()-start_col);
	assert(_cols*_rows>0);
	assert(_rows <= _ld);
}

template<class TYPE>
DenseMatrix<TYPE>::~DenseMatrix() {
	if (_created) { 
		delete [] _array;
	}
//	cout << "dtor:Anasazi::DenseMatrix " << this << endl;
}

template<class TYPE>
void DenseMatrix<TYPE>::setrows(const int rows) {
	if (rows >= 0) {
		_rows = rows;
	}
	else {
		_rows = 0;
	}
	assert(rows>=0);
}

template<class TYPE>
void DenseMatrix<TYPE>::setcols(const int cols) {
	if (cols >= 0) {
		_cols = cols;
	}
	else {
		_cols = 0;
	}
	assert(cols>=0);
}

template<class TYPE>
int DenseMatrix<TYPE>::getrows() const {
	return _rows;
}

template<class TYPE>
int DenseMatrix<TYPE>::getld() const {
	return _ld;
}

template<class TYPE>
int DenseMatrix<TYPE>::getcols() const {
	return _cols;
}

template<class TYPE>
TYPE DenseMatrix<TYPE>::getfronorm() const {
	int i,j;
	TYPE norm=0.0;
	for (j=0; j<_cols; j++) {
		for (i=0; i<_rows; i++) {
			norm += pow(_array[i+j*_ld],2);
		}
	}
	norm = sqrt(norm);
	return norm;
}

template<class TYPE>
TYPE* DenseMatrix<TYPE>::getarray() const {
	return _array;
}

template<class TYPE>
void DenseMatrix<TYPE>::setvalues(TYPE* array, const int ld) {
	int i,j;
	for (j=0; j<_cols; j++ ) {
		for (i=0; i<_rows; i++ ) {
			_array[i+j*_ld] = array[i+j*ld];
		}
	}
}

template<class TYPE>
void DenseMatrix<TYPE>::init( const TYPE value ) {
	int i,j;
	for (j=0; j<_cols; j++ ) {
		for (i=0; i<_rows; i++ ) {
			_array[i+j*_ld] = value;
		}
	}
}

template<class TYPE>
void DenseMatrix<TYPE>::print() {
	int i,j;
//	cout.setf(ios::scientific); //use scientific format
//	cout.precision(4);
	for (i=0; i<_rows; i++ ) {
		for (j=0; j<_cols; j++ ) {
			cout <<_array[i+j*_ld]<<'\t';
		}
		cout << endl;
	}
}

} //end namespace Anasazi
#endif
// End of file AnsaziDenseMatrix.hpp

