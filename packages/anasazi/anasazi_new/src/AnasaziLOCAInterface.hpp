// File AnasaziLOCA.hpp: interface for the AnasaziLOCA class.
//
#ifndef ANASAZI_LOCA_HPP
#define ANASAZI_LOCA_HPP

#include "AnasaziMatrix.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"
#include "NOX_Abstract_Vector.H"
#include "LOCA_Abstract_Group.H"
#include "NOX_Parameter_List.H"

namespace Anasazi {

enum DataAccess {Copy, View};

//
//--------template class AnasaziLOCAMat-----------------------
//
template <class TYPE> 
class LOCAMat : public Matrix<TYPE> {
public:
	LOCAMat( NOX::Parameter::List&, LOCA::Continuation::AnasaziGroup& );
	~LOCAMat();
	ReturnType ApplyMatrix ( const MultiVec<TYPE>& x, 
				       MultiVec<TYPE>& y ) const;
private:
	NOX::Parameter::List& locaParams;
	LOCA::Continuation::AnasaziGroup& locaGroup;
};
//------------------------------------------------------------
//
//--------template class AnasaziLOCAVec-------------------------------------
//
template <class TYPE>
class LOCAVec : public MultiVec<TYPE> {
public:
	friend class LOCAMat<TYPE>;
// constructors
	LOCAVec(const NOX::Abstract::Vector& N_vec, int NumVecs );
	LOCAVec(const vector< NOX::Abstract::Vector *> N_vecPtrs, 
			DataAccess type = Copy );
	LOCAVec(const LOCAVec<TYPE>& source, 
			DataAccess type = Copy );
	LOCAVec(DataAccess type, const LOCAVec<TYPE>& source, 
			int index[], int NumVecs); 

	~LOCAVec();
	//
	//  member functions inherited from MultiVec
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  not copied; instead a new MultiVec is created containing
	//  a non-zero amount of columns.
	//
	virtual MultiVec<TYPE> * Clone ( const int );
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied and a new stand-alone MultiVector is created.
	//  (deep copy).
	//
	virtual MultiVec<TYPE> * CloneCopy ();
	//
	//  Selective deep copy (or copy) constructor.
	//
	virtual MultiVec<TYPE> * CloneCopy ( int [], int );
	//
	//  the following is a virtual view constructor returning
	//  a pointer to the pure virtual class. vector values are 
	//  shared and hence no memory is allocated for the columns.
	//
	virtual MultiVec<TYPE> * CloneView ( int [], int );
	//
	virtual int GetNumberVecs () const;
	virtual int GetVecLength () const;
	//
	//  set a block of this multivec with the multivecs specified by
	//  the index.
	//
	virtual void SetBlock ( MultiVec<TYPE>& A, int index[], 
		int NumVecs ); 
	//
	// *this <- alpha * A * B + beta * (*this)
	//
	virtual void MvTimesMatAddMv ( TYPE alpha, MultiVec<TYPE>& A, 
		Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta );
	//
	// *this <- alpha * A + beta * B
	//
	virtual void MvAddMv ( TYPE alpha , MultiVec<TYPE>& A, TYPE beta,
		MultiVec<TYPE>& B);
	//
	// B <- alpha * A^T * (*this)
	//
	virtual void MvTransMv ( TYPE alpha, MultiVec<TYPE>& A, 
		Teuchos::SerialDenseMatrix<int,TYPE>& B );
	//
	// alpha[i] = norm of i-th column of (*this)
	//
	virtual void MvNorm ( TYPE* normvec );
	//
	// random vectors in i-th column of (*this)
	//
	virtual	void MvRandom();
        //
        // initializes each element of (*this) with alpha
        //
        virtual void MvInit ( TYPE alpha );
	//
	// print (*this)
	//
	virtual void MvPrint();
	//
	// Return a pointer to specific NOX::Abstract::Vector for LOCA computation.
	// If index is not a valid index for this multivec, then nothing is done.
	// NOTE:  This method is not included in the AnasaziMultiVec virtual base class.
	//
	virtual void GetNOXVector( NOX::Abstract::Vector& Vec, int index );
	// 
private:
// Data container
 	vector< NOX::Abstract::Vector* > mvPtrs;
	DataAccess CV;
};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

template<class TYPE>
LOCAVec<TYPE>::LOCAVec( const NOX::Abstract::Vector& N_vec, int NumVecs ) :
					mvPtrs(NumVecs), CV(Copy)
{
	for (int i=0; i<NumVecs; i++) {
		mvPtrs[i] = N_vec.clone(NOX::ShapeCopy);
		mvPtrs[i]->init(0.0);
//		cout<<"nox_vec_init "<<i<<"\t"<<
//			typeid(*(mvPtrs[i])).name()<<endl;
	}
}

template<class TYPE>
LOCAVec<TYPE>::LOCAVec( const vector< NOX::Abstract::Vector *> N_vecPtrs,
					DataAccess type ) : 
					mvPtrs(N_vecPtrs.size()), CV(type)
{
	int i;
	if (type == Copy) {
		for ( i=0; i<mvPtrs.size(); i++) {
			mvPtrs[i] = N_vecPtrs[i]->clone(NOX::DeepCopy);
		}
	} 
	else {
		for ( i=0; i<mvPtrs.size(); i++) {
			mvPtrs[i] = N_vecPtrs[i];
		}
	}
}

template<class TYPE>
LOCAVec<TYPE>::LOCAVec( const LOCAVec<TYPE>& source, 
			DataAccess type ) : mvPtrs(source.mvPtrs.size()),
			CV(type)
{
	int i;

	if (type == Copy) {
		for (i=0; i<mvPtrs.size(); i++) {
			mvPtrs[i] = source.mvPtrs[i]->clone(NOX::DeepCopy);
		}
	}
	else {
		for (i=0; i<mvPtrs.size(); i++) {
			mvPtrs[i] = source.mvPtrs[i];
		}
	}
}

template<class TYPE>
LOCAVec<TYPE>::LOCAVec( DataAccess type, const LOCAVec<TYPE>& source, 
			       	int index[], int NumVecs ): mvPtrs(NumVecs), CV(type)
{
	int i;

	if (type == Copy) {
		for ( i=0; i<NumVecs; i++ ) {
			mvPtrs[i] = source.mvPtrs[ index[i] ]->clone(NOX::DeepCopy);
//			cout<<"ALV_copy_init "<<i<<"\t"<<
//				typeid(*(mvPtrs[i])).name()<<endl;
		}
	} 
	else {
		for ( i=0; i<NumVecs; i++ ) {
			mvPtrs[i] = source.mvPtrs[ index[i] ];
//			cout<<"ALV_view_init "<<i<<"\t"<<
//				typeid(*(mvPtrs[i])).name()<<endl;
		}
	}
}

template<class TYPE>
LOCAVec<TYPE>::~LOCAVec()
{
	if (CV == Copy) {
		for (unsigned int i=0; i<mvPtrs.size(); i++) {
			delete mvPtrs[i];
		}
	}
}
//
//  member functions inherited from MultiVec
//
//
//  Simulating a virtual copy constructor. If we could rely on the co-variance
//  of virtual functions, we could return a pointer to LOCAVec<TYPE>
//  (the derived type) instead of a pointer to the pure virtual base class.
//
template<class TYPE>
MultiVec<TYPE>* LOCAVec<TYPE>::Clone ( const int NumVecs ) {
	LOCAVec *ptr_alv = new LOCAVec(*(mvPtrs[0]),NumVecs);
	return ptr_alv; // safe upcast.
}
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied.
	//
template<class TYPE>
MultiVec<TYPE>* LOCAVec<TYPE>::CloneCopy() {
	LOCAVec *ptr_alv = new LOCAVec(*this);
	return ptr_alv; // safe upcast
}

template<class TYPE>
MultiVec<TYPE>* LOCAVec<TYPE>::CloneCopy ( int index[], int NumVecs ) {
	
	LOCAVec *ptr_alv = new LOCAVec( Copy, *this, index, NumVecs );
	return ptr_alv; // safe upcast.
}

template<class TYPE>
MultiVec<TYPE>* LOCAVec<TYPE>::CloneView ( int index[], int NumVecs ) {
	
	LOCAVec *ptr_alv = new LOCAVec( View, *this, index, NumVecs );
	return ptr_alv; // safe upcast.
}

template<class TYPE>
int LOCAVec<TYPE>::GetNumberVecs () const {
	return mvPtrs.size();
}

template<class TYPE>
int LOCAVec<TYPE>::GetVecLength () const {
	return mvPtrs[0]->length();
}

template<class TYPE>
void LOCAVec<TYPE>::SetBlock( MultiVec<TYPE>& A, int index[], int NumVecs ) {

	int i, ind;
	LOCAVec *A_vec = dynamic_cast<LOCAVec *>(&A); assert(A_vec!=NULL);
	int MyNumVecs = mvPtrs.size();
	for (i=0; i<NumVecs; i++) {
		ind = index[i];
		if (ind < MyNumVecs) {
			delete mvPtrs[ind];
			mvPtrs[ind] = A_vec->mvPtrs[i]->clone(NOX::DeepCopy);				
		}
	}
}
//
// *this <- alpha * A * B + beta * (*this)
//
template<class TYPE>
void LOCAVec<TYPE>::MvTimesMatAddMv ( TYPE alpha, MultiVec<TYPE>& A, 
				      Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta ) 
{
	int i,j;
	LOCAVec *A_vec = dynamic_cast<LOCAVec *>(&A); assert(A_vec!=NULL);
	int m = B.numRows();
	int n = B.numCols();
	int ldb = B.stride();
	TYPE *Bvals = B.values();  	
	LOCAVec *temp_vec = new LOCAVec(*(mvPtrs[0]),n);
	temp_vec->MvInit(0.0);
	TYPE one = 1.0;
//
//	*this <- alpha * A * B + beta *(*this)
//
	for (j=0; j<n; j++) {
		for (i=0; i<m; i++) {
			temp_vec->mvPtrs[j]->update(Bvals[j*ldb+i], *(A_vec->mvPtrs[i]),one);
		}				
		mvPtrs[j]->update(alpha,*(temp_vec->mvPtrs[j]),beta);
	}
	delete temp_vec;
}
//
// *this <- alpha * A + beta * B
//
template<class TYPE>
void LOCAVec<TYPE>::MvAddMv ( TYPE alpha , MultiVec<TYPE>& A, 
						   TYPE beta, MultiVec<TYPE>& B) {
	const TYPE zero = 0.0;
	LOCAVec *A_vec = dynamic_cast<LOCAVec *>(&A); assert(A_vec!=NULL);
	LOCAVec *B_vec = dynamic_cast<LOCAVec *>(&B); assert(B_vec!=NULL);

	for (int i=0; i<mvPtrs.size(); i++) {
		mvPtrs[i]->update(alpha, *(A_vec->mvPtrs[i]), beta, *(B_vec->mvPtrs[i]), zero);
	}		
}
//
// dense B <- alpha * A^T * (*this)
//
template<class TYPE>
void LOCAVec<TYPE>::MvTransMv ( TYPE alpha, MultiVec<TYPE>& A,
				Teuchos::SerialDenseMatrix<int,TYPE>& B) {
	int i,j;
	LOCAVec *A_vec = dynamic_cast<LOCAVec *>(&A); assert(A_vec!=NULL);
	int m = B.numRows();
	int n = B.numCols();
	int ldb = B.stride();
	TYPE *Bvals = B.values();  	

	for (j=0; j<n; j++) {
		for (i=0; i<m; i++) {
			Bvals[j*ldb + i] = alpha * mvPtrs[j]->dot(*(A_vec->mvPtrs[i]));
		}
	}
}
//
// array[i] = norm of i-th column of (*this)
//
template<class TYPE>
void LOCAVec<TYPE>::MvNorm ( TYPE * normvec ) 
{
	if (normvec) {
		for (int i=0; i<mvPtrs.size(); i++) {
//			cout<<i<<"\t"<<typeid(*(mvPtrs[i])).name()<<endl;
			normvec[i] = mvPtrs[i]->norm();
		}
	}
}
//
// random vectors in i-th column of (*this)
//
template<class TYPE>
void LOCAVec<TYPE>::MvRandom () 
{
	for (unsigned int i=0; i<mvPtrs.size(); i++) {
		mvPtrs[i]->random();
	}	
}
//
// initializes each element of (*this) with alpha
//
template<class TYPE>
void LOCAVec<TYPE>::MvInit ( TYPE alpha ) 
{
	for (int i=0; i<mvPtrs.size(); i++) {
		mvPtrs[i]->init( alpha );
	}	
}
//
//  print multivectors
//
template<class TYPE>
void LOCAVec<TYPE>::MvPrint() {
//	cout << *this << endl;
}
//
//  return individual NOX Vector
//
template<class TYPE>
void LOCAVec<TYPE>::GetNOXVector( NOX::Abstract::Vector& Vec, int index ) 
{
	if (index < mvPtrs.size()) { Vec = *(mvPtrs[index]); }
}

///////////////////////////////////////////////////////////////
//
// implementation of the AnasaziLOCAMat class.
//
////////////////////////////////////////////////////////////////////
//
// Anasazi::Matrix constructors
//
template <class TYPE>
LOCAMat<TYPE>::LOCAMat(NOX::Parameter::List& params, 
					LOCA::Continuation::AnasaziGroup& group) :
					locaParams(params), locaGroup(group) {
//	cout << "ctor:Anasazi::LOCAMat " << this << endl;
	}

template <class TYPE>
LOCAMat<TYPE>::~LOCAMat() {
//	cout << "dtor:Anasazi::LOCAMat " << this << endl;
	}
//
// Matrix matrix multiply
//
template <class TYPE>
ReturnType LOCAMat<TYPE>::ApplyMatrix ( const MultiVec<TYPE>& x, 
					MultiVec<TYPE>& y ) const {
	
	NOX::Abstract::Group::ReturnType res;
	MultiVec<TYPE> &temp_x = const_cast<MultiVec<TYPE> &>(x);
	LOCAVec<TYPE> *x_vec = dynamic_cast<LOCAVec<TYPE> *>(&temp_x); assert(x_vec!=NULL);
	LOCAVec<TYPE> *y_vec = dynamic_cast<LOCAVec<TYPE> *>(&y); assert(y_vec!=NULL);

	int NumVecs = x_vec->GetNumberVecs();
	for (int i=0; i<NumVecs; i++) {
		res = locaGroup.applyAnasaziOperator(locaParams, *(x_vec->mvPtrs[i]), 
						*(y_vec->mvPtrs[i])); 
//		res = locaGroup.applyJacobian(*(x_vec->mvPtrs[i]), *(y_vec->mvPtrs[i])); 
	}
	if (res == NOX::Abstract::Group::Ok) {
	    return Ok;
	} else {
	    return Failed;
	}
}

} // end Anasazi namespace
#endif 
 // end of file ANASAZI_LOCA_HPP
