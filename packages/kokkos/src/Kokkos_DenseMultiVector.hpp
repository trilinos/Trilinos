//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2003) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef KOKKOS_DENSEMULTIVECTOR_H
#define KOKKOS_DENSEMULTIVECTOR_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_MultiVector.hpp"

namespace Kokkos {

//! Kokkos::DenseMultiVector: Kokkos default implementation of the abstract Kokkos::MultiVector.

/*! The Kokkos::DenseMultiVector provides a default implementation of Kokkos::MultiVector interface.

  At this time, the primary functions provided by Kokkos::DenseMultiVector is wrapping of a
  set of dense vector and providing access to multivector entries.  Two basic
  categories of data structures are supported:
  <ol>
  <li> MultiVector is described by an array of pointers:  In this situation, the 
       ith entry of this array of pointers
       is the starting address of a contiguous array of values for the ith vector 
       in the multivector.  The storage mode
       will be assumed if the getIsStrided() method returns false.
  <li> MultiVector is a regular strided two-dimensional array of values:  In this situation, 
       the increment between 
       elements in the row and column dimensions is specified by the getRowInc() 
       and getColInc() methods. This is
       a very general mechanism for describing strided access.  Typical situations include:
       <ul>
       <li> getRowInc() = getNumCols(), getColInc() = 1 - column entries are contiguous.  
       <li> getRowInc() = 1, getColInc() = getNumRows() - row entries are contiguous.  
       </ul>
       However, this mechanism also allows extraction of array subsections, real or
       imaginary parts from 
       complex-valued arrays.

       This storage mode will be assumed if getIsStrided() returns true.  The base 
       address for the 2D array
       will be obtain by call getValues() with the argument equal to 0.

  </ol>

*/    

  template<typename OrdinalType, typename ScalarType>
  class DenseMultiVector: public virtual MultiVector<OrdinalType,ScalarType> {
  public:

    //@{ \name Constructors/Destructor.

    //! Default constructor
    DenseMultiVector(void):
      dataInitialized_(false),
      numRows_(0),
      numCols_(0),
      rowInc_(0),
      colInc_(0),
      values_(0),
      allValues_(0),
      isStrided_(false) {};
  
    //! Copy constructor.
    DenseMultiVector(const DenseMultiVector& source):
      dataInitialized_(source.dataInitialized_),
      numRows_(source.numRows_),
      numCols_(source.numCols_),
      rowInc_(source.rowInc_),
      colInc_(source.colInc_),
      values_(source.values_),
      allValues_(source.allValues_),
      isStrided_(source.isStrided_) {

    if (isStrided_ && numCols_>0) {
      ScalarType ** tmpValues = new ScalarType*[numCols_];
      for (OrdinalType i=0;i<numCols_; i++) tmpValues[i] = values_[i];
      values_ = tmpValues;
    }
};

    //! DenseMultiVector Destructor
    virtual ~DenseMultiVector(){

    if (isStrided_ && numCols_>0) delete [] values_;
    
    };
    //@}

    //@{ \name Initialization methods
	
    //! Initialize using an array of pointers
    /*!
      This interface supports multivectors that are stored as an array of pointers.
      \param numRows (In)  Number of rows in multivector (length of each vector).
      \param numCols (In)  Number of columns in multivector (number of vectors).
      \param values (In)  Array of pointers of length getNumCols().  values[i] is a
      vector of length numRows.
      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(OrdinalType numRows, OrdinalType numCols, ScalarType ** values) {
      numRows_ = numRows;
      numCols_ = numCols;
      rowInc_ = 0;
      colInc_ = 1;
      values_ = values;
      allValues_ = 0;
      isStrided_ = false;
      dataInitialized_ = true;
      return(0);
      };
	
    //! Initialize using a two-dimensional array
    /*!
      This interface supports multivectors that are stored as 2D arrays, or subsections of one.
      \param numRows (In)  Number of rows in multivector (length of each vector).
      \param numCols (In)  Number of columns in multivector (number of vectors).
      \param values (In)  Pointer to the first entry in the multivector.  Subsequent column 
                          entries are spaced a distance of getColInc().  Subsequent row entries
			  are spaced by getRowInc() increments.
      \param rowInc (In) The increment between two elements in a row of the multivector.  
                         Typically this value should be set to numRows.
      \param colInc (In) The increment between two elements in a column of the multivector.  
                         Typically this value should be set to 1, which is the default value.

      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(OrdinalType numRows, OrdinalType numCols, ScalarType * values,
			 OrdinalType rowInc, OrdinalType colInc = 1) {
      numRows_ = numRows;
      numCols_ = numCols;
      rowInc_ = rowInc;
      colInc_ = colInc;
      values_ = 0;
      allValues_ = values;
      isStrided_ = true;
      dataInitialized_ = true;

    if (numCols_>0) {
      values_ = new ScalarType*[numCols_];
      for (OrdinalType i=0;i<numCols_; i++) values_[i] = allValues_+i*rowInc;
    }
      return(0);
      };
	
    //@}

    //@{ \name Multivector entry access methods.

    //! Returns an array of pointers such that the ith pointer points to an array of values in the ith column of the multivector.
    /*! Extract an array of pointers such that the ith pointer points to the values in the ith column of the multivector.  Note that
        the values are not copied by this method.  Memory allocation is 
	handled by the multivector object itself.  
    */
    virtual ScalarType ** getValues() const {return(values_);};

    //! Returns a pointer to an array of values for the ith column of the multivector.
    /*! Extract a pointer to the values in the ith column of the multivector.  Note that
        the values are not copied by this method.  Memory allocation is 
	handled by the multivector object itself.  Also, if the getIsStrided() method returns
	true, then the getColInc() should be used to access values in the ith column
	of the multivector, especially if getColInc() != 1.

	\param i (In) The column that should be returned.
    */
    virtual ScalarType * getValues(OrdinalType i) const {
      if (!dataInitialized_ || // No data to return
	  i<0 || // Out of range
	  i>=numRows_ // Out of range
	  ) return(0);
      return(values_[i]);
    };
	
	
    //@}

    //@{ \name DenseMultiVector Attribute access methods.
	
    //! Number of rows
    virtual OrdinalType getNumRows() const {return(numRows_);};
	
    //! Number of columns
    virtual OrdinalType getNumCols() const{return(numCols_);};
	
    //! Indicates whether or not array is strided
    virtual bool getIsStrided() const {return(isStrided_);};
	
    //! Increment between entries in a row of the multivector, normally = numRows().
    virtual OrdinalType getRowInc() const {return(rowInc_);};
	
    //! Increment between entries in a column of the multivector, normally = 1.
    virtual OrdinalType getColInc() const {return(colInc_);};
	
    //@}

  protected:
    bool dataInitialized_;
    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType rowInc_;
    OrdinalType colInc_;

    ScalarType ** values_;
    ScalarType * allValues_;

    bool isStrided_;
  };

} // namespace Kokkos
#endif /* KOKKOS_DENSEMULTIVECTOR_H */
