// AnasaziDenseMatrix.hpp: interface for the AnasaziDenseMatrix class.
#ifndef ANASAZI_DENSE_MATRIX_HPP
#define ANASAZI_DENSE_MATRIX_HPP
#include <iostream>

template <class TYPE>
class AnasaziDenseMatrix  {
public:
	AnasaziDenseMatrix();
	AnasaziDenseMatrix(int rows, int cols);
	AnasaziDenseMatrix(int rows, int ld, int cols);
	// copy constructor
	AnasaziDenseMatrix(const AnasaziDenseMatrix<TYPE> &);
	// create a view
	AnasaziDenseMatrix(const AnasaziDenseMatrix<TYPE> &, int, int);
	// create a view from a submatrix
	AnasaziDenseMatrix(const AnasaziDenseMatrix<TYPE> &, int, int, int, int);
	~AnasaziDenseMatrix();
	void setrows(const int);
	void setcols(const int);
	int getrows() const;
	int getcols() const;
	int getld() const;
	TYPE getfronorm() const;
	TYPE* getarray() const;
	void setvalues(TYPE* array, const int ldx);
	void DisplayMat();
private:
	const int _rows, _ld, _cols;
	const bool _created;
	TYPE *_array,*_endp1;
};

template<class TYPE> 
AnasaziDenseMatrix<TYPE>::AnasaziDenseMatrix(): _rows(0), _ld(0), _cols(0), 
						_created(false), _array(0), _endp1(0) { 
//	std::cout << "ctor1:AnasaziDenseMatrix " << this << std::endl;
}

template<class TYPE> 
AnasaziDenseMatrix<TYPE>::AnasaziDenseMatrix(int rows, int cols): 
							_rows(rows), _ld(rows), _cols(cols), 
						_created(true), _array(0), _endp1(0) { 
//	std::cout << "ctor2:AnasaziDenseMatrix " << this << std::endl;
	if (_ld*_cols > 0) {
		_array = new TYPE[_ld*_cols];
		_endp1 = _array + _ld*_cols;
	}
	assert(_array);
	assert(_ld*_rows>0);
}

template<class TYPE> 
AnasaziDenseMatrix<TYPE>::AnasaziDenseMatrix(int rows, int ld, int cols): 
							_rows(rows), _ld(ld), _cols(cols), 
						_created(true), _array(0), _endp1(0) { 
//	std::cout << "ctor3:AnasaziDenseMatrix " << this << std::endl;
	if (_ld*_cols > 0 ) {
		_array = new TYPE[_ld*_cols];
		_endp1 = _array + _ld*_cols;
	}
	assert(_array);
	assert(_ld*_rows>0);
	assert(_rows <= _ld);
}

//
// Copy constructor
//
template<class TYPE>
AnasaziDenseMatrix<TYPE>::AnasaziDenseMatrix(const AnasaziDenseMatrix& A):
							_rows(A.getrows()), _ld(A.getrows()), _cols(A.getcols()),
							_created(true), _array(0), _endp1(0) {
//	std::cout << "copy_ctor:AnasaziDenseMatrix " << this << std::endl;
	if (_ld*_cols > 0 ) {
		_array = new TYPE[_ld*_cols]; assert(_array);
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
AnasaziDenseMatrix<TYPE>::AnasaziDenseMatrix(const AnasaziDenseMatrix& A, int rows, int cols):
							_rows(rows), _ld(A.getld()), _cols(cols),
							_created(false), _array(0), _endp1(0) {
//	std::cout << "view_ctor1:AnasaziDenseMatrix " << this << std::endl;
	if (_rows*_cols > 0 && 
		_rows <= A.getrows() && 
		_cols <= A.getcols()) {
		_array = A.getarray(); 
		_endp1 = _array+(_cols-1)*_ld + _rows; // check
	}
	assert(_array);
	assert(start_row*start_col >= 0);
	assert(_rows <= A.getrows());
	assert(_cols <= A.getcols());
	assert(_cols*_rows>0);
	assert(_rows <= _ld);
}

// create a partial view
template<class TYPE>
AnasaziDenseMatrix<TYPE>::AnasaziDenseMatrix(const AnasaziDenseMatrix& A, int start_row, 
											 int start_col, int rows, int cols):
							_rows(rows), _ld(A.getld()), _cols(cols),
							_created(false), _array(0), _endp1(0) {
//	std::cout << "view_ctor1:AnasaziDenseMatrix " << this << std::endl;
	// start_col and start_row are zero based, i.e., rows and columns are indexed
	// from 0.
	if (_rows*_cols > 0 && 
		start_row*start_col >= 0 &&
		_rows <= A.getrows()-start_row && 
		_cols <= A.getcols()-start_col) {
		_array = A.getarray() + start_col*_ld + start_row; 
		_endp1 = _array+(_cols-1)*_ld + _rows; // check
	}
	assert(_array);
	assert(start_row*start_col >= 0);
	assert(_rows <= A.getrows()-start_row);
	assert(_cols <= A.getcols()-start_col);
	assert(_cols*_rows>0);
	assert(_rows <= _ld);
}

template<class TYPE>
AnasaziDenseMatrix<TYPE>::~AnasaziDenseMatrix() {
	if (_created) { 
		delete [] _array;
	}
//	std::cout << "dtor:AnasaziDenseMatrix " << this << std::endl;
}

template<class TYPE>
void AnasaziDenseMatrix<TYPE>::setrows(const int rows) {
	if (rows >= 0) {
		_rows = rows;
	}
	else {
		_rows = 0;
	}
	assert(rows>=0);
}

template<class TYPE>
void AnasaziDenseMatrix<TYPE>::setcols(const int cols) {
	if (cols >= 0) {
		_cols = cols;
	}
	else {
		_cols = 0;
	}
	assert(cols>=0);
}

template<class TYPE>
int AnasaziDenseMatrix<TYPE>::getrows() const {
	return _rows;
}

template<class TYPE>
int AnasaziDenseMatrix<TYPE>::getld() const {
	return _ld;
}

template<class TYPE>
int AnasaziDenseMatrix<TYPE>::getcols() const {
	return _cols;
}

template<class TYPE>
TYPE AnasaziDenseMatrix<TYPE>::getfronorm() const {
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
TYPE* AnasaziDenseMatrix<TYPE>::getarray() const {
	return _array;
}

template<class TYPE>
void AnasaziDenseMatrix<TYPE>::setvalues(TYPE* array, const int ld) {
	int i,j;
	for (j=0; j<_cols; j++ ) {
		for (i=0; i<_rows; i++ ) {
			_array[i+j*_ld] = array[i+j*ld];
		}
	}
}

template<class TYPE>
void AnasaziDenseMatrix<TYPE>::DisplayMat() {
	//
	int i, j, numcols, numrows;
	TYPE* ary = getarray();
	numcols = getcols();
	numrows = getrows();

	for (i=0; i<getrows(); i++) {
		for (j=0; j<getcols(); j++){
			cout << ary[i+j*numrows] << " ";
		}
		cout << endl;
	}
}
//
#endif
