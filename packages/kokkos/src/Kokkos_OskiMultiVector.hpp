//@HEADER
// ************************************************************************
// 
//                Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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

// Author: Jim Willenbring jmwille@sandia.gov 08-17-2005

#ifndef KOKKOS_OSKIMULTIVECTOR_HPP
#define KOKKOS_OSKIMULTIVECTOR_HPP

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_MultiVector.hpp"
extern "C" {
  #include <oski/oski.h>
}
namespace Kokkos {

//! Kokkos::OskiMultiVector: Implementation of the abstract Kokkos::MultiVector  that wraps an oski_vecview_t object.

/*! 
  At this time, the primary function provided by Kokkos::OskiMultiVector is wrapping 
  an oski_vecview_t object and providing access to its entries.

*/    
// Add additional info above about striding, etc.
  
  template<typename OrdinalType, typename ScalarType>
  class OskiMultiVector: public virtual MultiVector<OrdinalType,ScalarType> {
  public:

    //! @name Constructors/Destructor

    //@{

    //! Default constructor
    OskiMultiVector(void):
      dataInitialized_(false),
      values_(0),
      allValues_(0),
      isStrided_(false) {};
  
    //! Copy constructor.
    OskiMultiVector(const OskiMultiVector& source):
      dataInitialized_(source.dataInitialized_),
      values_(source.values_),
      allValues_(source.allValues_),
      isStrided_(source.isStrided_) { 
      if (dataInitialized_) {
	x_view_ = oski_CopyVecView(source.x_view_);
	// Look at this piece -  is it necessary?
        OrdinalType numCols = x_view_->num_cols;
        if (isStrided_ &&  numCols>0) {
          ScalarType ** tmpValues = new ScalarType*[numCols];
          for (OrdinalType i=0;i<numCols; i++) tmpValues[i] = values_[i];
          values_ = tmpValues;
        }
      }
};

    //! OskiMultiVector Destructor
    virtual ~OskiMultiVector(){

    // **This is causing huge problems, see if "newing" arrays is the problem
    // Fixed? 
    if (dataInitialized_) oski_DestroyVecView(x_view_);
    if (isStrided_ && x_view_->num_cols>0) delete [] values_;
    
    };
    //@}

    //! @name Initialization methods

    //@{

    //! Initialize using an array of pointers
    /*!
      This interface supports multivectors that are stored as an array of pointers.
      \param numRows (In)  Number of rows in multivector (length of each vector).
      \param numCols (In)  Number of columns in multivector (number of vectors).
      \param values (In)  Array of pointers of length getNumCols().  values[i] is a
      vector of length numRows.  Vector must be stored with a constant stride
      between rows.  (May change in the future.)
      \return Integer error code, set to 0 if successful.
    */
    int initializeValues(OrdinalType numRows, OrdinalType numCols, ScalarType ** values) {
      values_ = values;
      allValues_ = 0;
      dataInitialized_ = true;
      //isStrided_ = false; don't need this unless we can use such a MultiVec
      for (int i = 0; i < numCols-1; i++) if (values_[i+1] - values_[i] != numRows) return(-1); // doesn't allow nonstrided vectors at this point
      isStrided_ = true;
      if (numCols==1) x_view_=oski_CreateVecView(*values,numRows,STRIDE_UNIT);
      else x_view_=oski_CreateMultiVecView(*values,numRows,numCols,LAYOUT_COLMAJ,numRows);

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
      if (numCols==1) x_view_=oski_CreateVecView(values,numRows,colInc);
      else x_view_=oski_CreateMultiVecView(values,numRows,numCols,LAYOUT_COLMAJ,numRows);
      values_ = 0;// **Figure out if I can/ need to accept this array of values
      allValues_ = values;
      isStrided_ = true;
      dataInitialized_ = true;

    if (numCols>0) {
      values_ = new ScalarType*[numCols];
      for (OrdinalType i=0;i<numCols; i++) values_[i] = allValues_+i*rowInc;
    }
      return(0);
      };
	
    //@}

    //! @name Multivector entry access methods

    //@{

    //! Returns an array of pointers such that the ith pointer points to an array of values in the ith column of the multivector.
    /*! Extract an array of pointers such that the ith pointer points to the values in the ith column of the multivector.  Note that
        the values are not copied by this method.  Memory allocation is 
	handled by the multivector object itself.  
    */
//    virtual ScalarType ** getValues() const {return(values_);};
    virtual ScalarType ** getValues() const {
	if (dataInitialized_) return(&(x_view_->val));
	else return(values_);
    };

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
	  i>=x_view_->num_rows // Out of range
	  ) return(0);
      return(values_[i]);
    };
	
	
    //@}

    //! @name DenseMultiVector Attribute access methods

    //@{
	
    //! Number of rows
    virtual OrdinalType getNumRows() const {
      if (dataInitialized_) return(x_view_->num_rows);
      else return(0);
    };
	
    //! Number of columns
    virtual OrdinalType getNumCols() const{
      if (dataInitialized_) return(x_view_->num_cols);
      else return(0);
    };
	
    //! Indicates whether or not array is strided
    virtual bool getIsStrided() const {
      if (dataInitialized_) return(isStrided_);
      else return(false);
    };
	
    //! Increment between entries in a row of the multivector, normally = numRows().
    virtual OrdinalType getRowInc() const {
      if (!dataInitialized_) return(0);
      OrdinalType leadDim = x_view_->stride;
      OrdinalType rowInc = 0;
      if (x_view_->orient == LAYOUT_COLMAJ) rowInc = leadDim;
      else rowInc = leadDim/x_view_->num_cols;
      return(rowInc);    
    };
	
    //! Increment between entries in a column of the multivector, normally = 1.
    virtual OrdinalType getColInc() const {
      if (!dataInitialized_) return(0);
      OrdinalType leadDim = x_view_->stride;
      OrdinalType colInc = 0;
      if (x_view_->orient == LAYOUT_ROWMAJ) colInc = leadDim;
      else colInc = leadDim/x_view_->num_rows;
      return(colInc);
    };

    //! Underlying OSKI MultiVector
    virtual oski_vecview_t getX_view() const{return(x_view_);};
	
    //@}

  protected:
    oski_vecview_t x_view_;

    bool dataInitialized_;

    ScalarType ** values_;
    ScalarType * allValues_;

    bool isStrided_;
  };

} // namespace Kokkos
#endif /* KOKKOS_OSKIMULTIVECTOR_HPP */

