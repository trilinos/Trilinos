
#ifndef _TEUCHOS_SERIALDENSEVECTOR_HPP_
#define _TEUCHOS_SERIALDENSEVECTOR_HPP_

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Object.hpp" 
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Teuchos {

  template<typename OrdinalType, typename ScalarType>
  class SerialDenseVector : public SerialDenseMatrix<OrdinalType,ScalarType> {
    
  public:
    
    SerialDenseVector();
    SerialDenseVector(int length);
    SerialDenseVector(DataAccess CV, ScalarType* values, int length);
    SerialDenseVector(const SerialDenseVector<OrdinalType,ScalarType>& Source);
    int size(int length) {return(SerialDenseMatrix<OrdinalType, ScalarType>::shape(length, 1));};
    int resize(int length) {return(SerialDenseMatrix<OrdinalType,ScalarType>::reshape(length, 1));};
    virtual ~SerialDenseVector ();
    bool operator == (const SerialDenseVector<OrdinalType, ScalarType> &Operand);
    bool operator != (const SerialDenseVector<OrdinalType, ScalarType> &Operand);
    SerialDenseVector<OrdinalType,ScalarType>& operator = (const SerialDenseVector<OrdinalType,ScalarType>& Source);
    ScalarType& operator () (int index);
    const ScalarType& operator () (int index) const;
    ScalarType& operator [] (int index);
    const ScalarType& operator [] (int index) const;
    int length() const {return(numRows_);};
    virtual void print(ostream& os) const;
};

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::SerialDenseVector() : SerialDenseMatrix<OrdinalType,ScalarType>() {}

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::SerialDenseVector( int length ) : SerialDenseMatrix<OrdinalType,ScalarType>( length, 1 ) {}

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::SerialDenseVector(DataAccess CV, ScalarType* values, int length) : 
    SerialDenseMatrix<OrdinalType,ScalarType>( CV, values, length, length, 1 ) {}

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::SerialDenseVector(const SerialDenseVector<OrdinalType, ScalarType> &Source) :
    SerialDenseMatrix<OrdinalType,ScalarType>( Source ) {}

  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>::~SerialDenseVector() {}
  
  template<typename OrdinalType, typename ScalarType>
  SerialDenseVector<OrdinalType, ScalarType>& SerialDenseVector<OrdinalType,ScalarType>::operator = (const SerialDenseVector<OrdinalType, ScalarType>& Source) 
  {
    SerialDenseMatrix<OrdinalType,ScalarType>::operator=(Source); 
    return(*this);
  }

  template<typename OrdinalType, typename ScalarType>
  bool SerialDenseVector<OrdinalType, ScalarType>::operator == (const SerialDenseVector<OrdinalType, ScalarType> &Operand) 
  {
    bool result = 1;
    if(numRows_ != Operand.numRows_)
      {
	result = 0;
      }
    else
      {
	int i;
	for(i = 0; i < numRows_; i++) {
	  if((*this)(i) != Operand(i))
	    {
	      return 0;
	    }
	}
      }
    return result;
  }

  template<typename OrdinalType, typename ScalarType>
  bool SerialDenseVector<OrdinalType, ScalarType>::operator != (const SerialDenseVector<OrdinalType, ScalarType> &Operand)
  {
    return !((*this)==Operand);
  }

  template<typename OrdinalType, typename ScalarType>
  void SerialDenseVector<OrdinalType, ScalarType>::print(ostream& os) const
  {
    os << endl;
    if(valuesCopied_)
      os << "Values_copied : yes" << endl;
    else
      os << "Values_copied : no" << endl;
      os << "Length : " << numRows_ << endl;
    if(numRows_ == 0) {
      os << "(vector is empty, no values to display)" << endl;
    } else {
      for(int i = 0; i < numRows_; i++) {
	  os << (*this)(i) << " ";
      }
      os << endl;
    }
  }

  //----------------------------------------------------------------------------------------------------
  //   Accessor methods 
  //----------------------------------------------------------------------------------------------------

  template<typename OrdinalType, typename ScalarType>
  inline ScalarType& SerialDenseVector<OrdinalType, ScalarType>::operator () (int index)
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    if (index >= numRows_)
      {
	cout << "Row index = " << index << " Out of Range 0 - " << numRows_-1 << endl;
	TEUCHOS_CHK_REF(*values_); // Return reference to values_[0]
      }
#endif
    return(values_[index]);
  }
  
  template<typename OrdinalType, typename ScalarType>
  inline const ScalarType& SerialDenseVector<OrdinalType, ScalarType>::operator () (int index) const
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    if (index >= numRows_)
      {
	cout << "Row index = " << index << " Out of Range 0 - " << numRows_ - 1 << endl;
	TEUCHOS_CHK_REF(values_[0]); // Return reference to values_[0]
      }
#endif
    return(values_[index]);
  }
  
  template<typename OrdinalType, typename ScalarType>
  inline const ScalarType& SerialDenseVector<OrdinalType, ScalarType>::operator [] (int index) const
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    if (index >= numRows_)
      {
	cout << "Row index = " << index << " Out of Range 0 - " << numRows_ - 1 << endl;
	TEUCHOS_CHK_PTR(0); // Return zero pointer
      }
#endif
    return(values_[index]);
  }
  
  template<typename OrdinalType, typename ScalarType>
  inline ScalarType& SerialDenseVector<OrdinalType, ScalarType>::operator [] (int index)
  {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    if (index >= numRows_)
      {
    cout << "Column index = " << index << " Out of Range 0 - " << numRows_ - 1 << endl;
    TEUCHOS_CHK_PTR(0); // Return zero pointer
      }
#endif
    return(values_[index]);
  }

} // namespace Teuchos

#endif /* _TEUCHOS_SERIALDENSEVECTOR_HPP_ */
