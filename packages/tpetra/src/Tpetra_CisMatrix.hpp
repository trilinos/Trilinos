// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef _TPETRA_CISMATRIX_HPP_
#define _TPETRA_CISMATRIX_HPP_

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_CompObject.hpp>
#include <Kokkos_HbMatrix.hpp>
#include <Kokkos_DenseVector.hpp>
#include <Kokkos_BaseSparseMultiply.hpp>
#include "Tpetra_Object.hpp"
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_VectorSpace.hpp"

namespace Tpetra {

// forward declaration of CisMatrixData, needed to prevent circular inclusions
// actual #include statement at the end of this file
template<typename OrdinalType, typename ScalarType> class CisMatrixData;

//! Tpetra::CisMatrix: A class for constructing and using sparse compressed index matrices.
/*! CisMatrix enables the piecewise construction and use of sparse matrices
    where matrix entries are intended for either row or column access.

    At this time, the primary function provided by CisMatrix is matrix times vector and multiplication.
    It is also possible to extract a Kokkos::HBMatrix so that a Tpetra::CisMatrix can be used in Kokkos kernels.

<b>Constructing CisMatrix objects</b>

Constructing CisMatrix objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create a CisMatrix instance, using one of the constructors.
  <li> Enter values using the submitEntries methods.
  <li> Complete construction by calling fillComplete.
</ol>

<b>Primary and Secondary Distributions</b>
CisMatrix stores data using Compressed Index Space Storage. This is a generalization of Compressed Row Storage,
but allows the matrix to be either row-oriented or column oriented. Accordingly, CisMatrix refers to data using
a generalized vocabulary. The two main terms are the primary distribution and the secondary distribution:

In a row-oriented matrix, the primary distribution refers to rows, and the secondary distribution refers to columns.
In a column-oriented matrix, the primary distribution refers to columns, and the secondary distribution refers to rows.

These distributions are specified by Tpetra::VectorSpace objects. If both the primary and secondary distributions are
specified at construction, information about the secondary distribution is available from then on. But if only the primary
distribution is specified at construction, then CisMatrix will analyze the structure of the matrix and generate a 
VectorSpace that matches the distribution found in the matrix. Note that this is done when fillComplete is called, and so
information about the secondary distribution is not available prior to that.

<b> FillComplete </b>
A Tpetra::CisMatrix exists in one of two states. Prior to calling fillComplete, data is stored in a format optimized for
easy and efficient insertions/deletions/modifications. It is not a very efficient form for doing matvec operations though.
After calling fillComplete, the data is transformed into a format optimized for matvec operations. It is not very efficient
at modifying data though. So after fillComplete has been called, the matrix should be viewed as const, and cannot be modified.
Trying to do so will result in an exception being thrown.

<b> Counting Floating Point Operations </b>

CisMatrix inherits from Teuchos::CompObject, and keeps track of the number of floating point operations
associated with the \e this object. In conjunction with a timing object (such as Teuchos::Time), accurate
parallel performance information can be obtained. Further information can be found in the Teuchos documentation.

CisMatrix error codes (positive for non-fatal, negative for fatal):
<ol>
<li> +1  That global ID is not owned by this image.
<li> +2  Unsupported combine mode specified.
<li> +3  Requested distribution is not currently defined.
<li> +4  Cannot call that method until after fillComplete.
<li> +5  Cannot call that method after fillComplete.
<li> -99 Internal CisMatrix error. Contact developer.
</ol>

*/

template<typename OrdinalType, typename ScalarType>
class CisMatrix : public Object, public Teuchos::CompObject {

public:
  
	//@{ \name Constructor/Destructor Methods
  
  //! Constructor specifying the primary distribution only.
	CisMatrix(VectorSpace<OrdinalType, ScalarType> const& primaryDist, bool rowOriented = true)
		: Object("Tpetra::CisMatrix")
		, CisMatrixData_()
  {
    CisMatrixData_ = Teuchos::rcp(new CisMatrixData<OrdinalType, ScalarType>(primaryDist, rowOriented));
  }
  
  //! Constructor specifying both primary and secondary distributions.
	CisMatrix(VectorSpace<OrdinalType, ScalarType> const& primaryDist, 
							 VectorSpace<OrdinalType, ScalarType> const& secondaryDist, 
							 bool rowOriented = true)
		: Object("Tpetra::CisMatrix")
		, CisMatrixData_()
  {
		CisMatrixData_ = Teuchos::rcp(new CisMatrixData<OrdinalType, ScalarType>(primaryDist, secondaryDist, rowOriented));
  }
	
	//! copy constructor.
	CisMatrix(CisMatrix<OrdinalType, ScalarType> const& Source)
		: Object(Source.label())
		, CisMatrixData_(Source.CisMatrixData_)
	{}

	//! destructor.
	~CisMatrix() {}
  
	//@}
  
  //@{ \name Post-Construction Modification Routines
  
  //! Set all matrix entries equal to scalarThis
  void setAllToScalar(ScalarType scalarThis) {
    if(isFillCompleted())
      throw reportError("Cannot change matrix values after fillComplete has been called.", 5);
    typedef std::map<OrdinalType, ScalarType> OrdScalMap;
    typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;
    MapOfMaps& outermap = CisMatrixData_->indicesAndValues_;
    for(typename MapOfMaps::iterator i = outermap.begin(); i != outermap.end(); i++) {
      OrdScalMap& innermap = (*i).second;
      for(typename OrdScalMap::iterator j = innermap.begin(); j != innermap.end(); j++)
        (*j).second = scalarThis;
    }
  }
  
  //! Scale the current values of a matrix, \e this = scalarThis*\e this.
  void scale(ScalarType scalarThis) {
    if(isFillCompleted())
      throw reportError("Cannot change matrix values after fillComplete has been called.", 5);
    typedef std::map<OrdinalType, ScalarType> OrdScalMap;
    typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;
    MapOfMaps& outermap = CisMatrixData_->indicesAndValues_;
    for(typename MapOfMaps::iterator i = outermap.begin(); i != outermap.end(); i++) {
      OrdScalMap& innermap = (*i).second;
      for(typename OrdScalMap::iterator j = innermap.begin(); j != innermap.end(); j++)
        (*j).second *= scalarThis;
    }

    // update flops counter: nnz
    updateFlops(data().numMyNonzeros_);
  }

  //! Submit multiple entries, using global IDs.
  /*! All index values must be in the global space. Behavoir is defined by the CombineMode passed in. */
  void submitEntries(CombineMode CM, OrdinalType myRowOrColumn, OrdinalType numEntries, ScalarType const* values, OrdinalType const* indices) {
		if(isFillCompleted())
      throw reportError("Cannot change matrix values after fillComplete has been called.", 5);
    OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
    for(OrdinalType i = zero; i < numEntries; i++)
      submitEntry(CM, myRowOrColumn, *values++, *indices++);
  }

  //! Submit a single entry, using global IDs.
  /*! All index values must be in the global space. Behavoir is defined by the CombineMode passed in. */
  void submitEntry(CombineMode CM, OrdinalType myRowOrColumn, ScalarType value, OrdinalType index) {
    if(isFillCompleted())
      throw reportError("Cannot change matrix values after fillComplete has been called.", 5);
    // first check for proper index values
    if(!getPrimaryDist().isMyGlobalIndex(myRowOrColumn))
      throw reportError("Global primary index " + toString(myRowOrColumn) + " is not owned by this image.", 1);
    if(CisMatrixData_->haveSecondary_)
      if(!getSecondaryDist().isMyGlobalIndex(index))
        throw reportError("Global secondary index " + toString(myRowOrColumn) + "is not owned by this image.", 1);

    // create a map for that row/column if it doesn't exist
    if(CisMatrixData_->indicesAndValues_.find(myRowOrColumn) == CisMatrixData_->indicesAndValues_.end()) {
      std::map<OrdinalType, ScalarType> temp;
      CisMatrixData_->indicesAndValues_[myRowOrColumn] = temp;
    }
    
    // then submit the actual value
    std::map<OrdinalType, ScalarType>& innerMap = CisMatrixData_->indicesAndValues_[myRowOrColumn];
    if(CM == Add) {
      // innerMap[index] += value; // will this work??
      if(innerMap.find(index) != innerMap.end())
        innerMap[index] = innerMap[index] + value;
      else
        innerMap[index] = value;

      // update flops counter: 1
      updateFlops(Teuchos::OrdinalTraits<OrdinalType>::one());
    }
    else if(CM == Replace) {
      innerMap[index] = value;
    }
    else if(CM == Insert) {
      innerMap[index] = value; // change this to a call to insert
    }
    else
      throw reportError("Unknown Combine Mode.", 2);
      
    data().numMyNonzeros_++;
  }
  
  //! Signals that data entry is complete. Matrix data is converted into a more optimized form.
  /*! The domain distribution and range distribution will be set equal to the primary distribution. */
  void fillComplete() {
    fillComplete(getPrimaryDist(), getPrimaryDist());
  }
  
  //! Signals that data entry is complete. Matrix data is converted into a more optimized form.
  /*! The VectorSpaces passed in will be used for the domain and range distributions. */
  void fillComplete(VectorSpace<OrdinalType, ScalarType> const& domainSpace, VectorSpace<OrdinalType, ScalarType> const& rangeSpace) {
    if(isFillCompleted())
      throw reportError("Already fillCompleted.", -99);

		OrdinalType const ordinalZero = Teuchos::OrdinalTraits<OrdinalType>::zero();
		OrdinalType const ordinalOne = Teuchos::OrdinalTraits<OrdinalType>::one();

    // set domain and range distributions
    data().domain_ = domainSpace;
    data().range_ = rangeSpace;
    data().haveDomain_ = true;
    data().haveRange_ = true;
    
    // fill pntr_, indx_ and values_ arrays from map values
    OrdinalType const numPrimaryIndices = getPrimaryDist().getNumMyEntries();
    OrdinalType const numNonzeros = getNumMyNonzeros();
    data().pntr_.reserve(numPrimaryIndices + ordinalOne);
    data().indx_.reserve(numNonzeros);   // we need to allocate space in the vectors
    data().values_.reserve(numNonzeros); // since they were created with a size of zero

    OrdinalType pntr_value = ordinalZero; // value to put into pntr_[pntr_loc]
    
    typedef std::map<OrdinalType, ScalarType> OrdScalMap;
    typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;
    MapOfMaps& outermap = CisMatrixData_->indicesAndValues_;
    for(typename MapOfMaps::iterator i = outermap.begin(); i != outermap.end(); i++) {
      OrdScalMap& innermap = (*i).second;
      data().pntr_.push_back(pntr_value);
      for(typename OrdScalMap::iterator j = innermap.begin(); j != innermap.end(); j++) {
        data().indx_.push_back((*j).first);
        data().values_.push_back((*j).second);
        pntr_value++;
      }
    }

    data().pntr_.push_back(pntr_value); // pntr_ has an extra element on the end
    
    // create secondary distribution if we need to
    if(!data().haveSecondary_) {
			/* !!!!!!!!!!!!! Temporary fix:  leave secondary ES equal to primary ES !!!!!!!!!!!!!!!!!!!!!!!
			// create elementspace first
			OrdinalType numGlobalElements = ordinalZero - ordinalOne; // set to -1
			OrdinalType numMyElements = data().indx_.size();
			OrdinalType* elementList = &data().indx_.front(); // address of first element in indx_
			OrdinalType indexBase = getPrimaryDist().getIndexBase();
			Platform<OrdinalType, OrdinalType> const& platformO = getPrimaryDist().elementSpace().platform();
			ElementSpace<OrdinalType> elementspace(numGlobalElements, numMyElements, elementList, indexBase, platformO);
			//cout << "[CIS] nGE=" << numGlobalElements << " nME=" << numMyElements << " iB=" << indexBase;
			// then create vectorspace using it
			Platform<OrdinalType, ScalarType> const& platformS = platform();
			VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformS);*/
      
			data().haveSecondary_ = true;
			if(isRowOriented())
				data().haveCol_ = true;
			else
				data().haveRow_ = true;
    }

    // initialize Kokkos::HbMatrix (Classical form)
		int errorcode = data().HbMatrix_.initializeStructure(getRowDist().getNumMyEntries(), 
																												 getColumnDist().getNumMyEntries(), 
																												 isRowOriented(), 
																												 &data().pntr_.front(), 
																												 &data().indx_.front());
		if(errorcode) 
			throw reportError("HbMatrix_.initializeStructure returned non-zero. code = " + toString(errorcode) + ".", -99);

    errorcode = data().HbMatrix_.initializeValues(&data().values_.front());
		if(errorcode) 
			throw reportError("HbMatrix_.initializeValues returned non-zero. code = " + toString(errorcode) + ".", -99);

    data().fillCompleted_ = true;
  }
  
  //@}
  
  //@{ \name Computational Methods
  
  //! Computes the matrix-vector multiplication y = Ax
  void apply(Vector<OrdinalType, ScalarType> const& x, Vector<OrdinalType, ScalarType>& y) const {
    if(!isFillCompleted())
      throw reportError("Cannot apply until after fillComplete.", 4);
    
    // setup kokkos x vector
    Kokkos::DenseVector<OrdinalType, ScalarType> kx;
    ScalarType* scalarArray = const_cast<ScalarType*>(x.scalarPointer()); // x is logically const but not bitwise const
    int errorcode = kx.initializeValues(x.getNumMyEntries(), scalarArray);// we will not be modifying it, 
		                                                                      // but initializeValues is not a const function.
		if(errorcode) 
			throw reportError("kx.initializeValues returned non-zero. code = " + toString(errorcode) + ".", -99);
		

    // setup kokkos y vector
    Kokkos::DenseVector<OrdinalType, ScalarType> ky;
    errorcode = ky.initializeValues(y.getNumMyEntries(), y.scalarPointer());
		if(errorcode) 
			throw reportError("ky.initializeValues returned non-zero. code = " + toString(errorcode) + ".", -99);

    // setup kokkos sparsemultiply object
    Kokkos::BaseSparseMultiply<OrdinalType, ScalarType> axy;
    errorcode = axy.initializeStructure(data().HbMatrix_);
		if(errorcode) 
			throw reportError("axy.initializeStructure returned non-zero. code = " + toString(errorcode) + ".", -99);
    errorcode = axy.initializeValues(data().HbMatrix_);
		if(errorcode) 
			throw reportError("ky.initializeValues returned non-zero. code = " + toString(errorcode) + ".", -99);

		/*cout << "=BEFORE APPLY================================================" << endl;
		cout << "x:" << endl << x << endl;
		cout << "y:" << endl << y << endl;
		cout << "=============================================================" << endl;*/

    // do Kokkos apply operation
    errorcode = axy.apply(kx, ky);
		if(errorcode) 
			throw reportError("axy.apply returned non-zero. code = " + toString(errorcode) + ".", -99);

		/*cout << "=AFTER APPLY=================================================" << endl;
		cout << "x:" << endl << x << endl;
		cout << "y:" << endl << y << endl;
		cout << "=============================================================" << endl;*/

    // update flops counter: 2 * nnz
    updateFlops(data().numMyNonzeros_ + data().numMyNonzeros_);
  }
  
  //! Returns the global one norm of the matrix
  ScalarType normOne() const {
    if(rowOriented())
      if(!data().haveCol_)
        throw reportError("Cannot compute one-norm until column distribution is specified.", 3);
      else
        return(secondaryNorm());
    else
      return(primaryNorm());
  }

  //! Returns the global infinity norm of the matrix
  ScalarType normInf() const {
    if(rowOriented())
      return(primaryNorm());
    else
      if(!data().haveRow_)
        throw reportError("Cannot compute infinity-norm until row distribution is specified.", 3);
      else
        return(secondaryNorm());
  }
  
  //@}
  
  //@{ \name Attribute Access Methods (Most of these can only be called after fillComplete has been called.)
  
	//! Returns the number of nonzero entries in the global matrix.
	OrdinalType getNumGlobalNonzeros() const {
    return(globalSum(getNumMyNonzeros()));
  }
  
  //! Returns the number of nonzero entries in the calling image's portion of the matrix.
	OrdinalType getNumMyNonzeros() const {
    return(data().numMyNonzeros_);
  }
	
	//! Returns the number of global matrix rows.
	OrdinalType getNumGlobalRows() const {
    if(!data().haveRow_)
      throw reportError("Row distribution not specified.", 3);
    return(getRowDist().getNumGlobalEntries());
  }
	
	//! Returns the number of global matrix columns.
	OrdinalType getNumGlobalCols() const {
    if(!data().haveCol_)
      throw reportError("Column distribution not specified.", 3);
    return(getColDist().getNumGlobalEntries());
  }
  
  //! Returns the number of matrix rows owned by the calling image.
	OrdinalType getNumMyRows() const {
    if(!data().haveRow_)
      throw reportError("Row distribution not specified.", 3);
    return(getRowDist().getNumMyEntries());
  }
	
	//! Returns the number of matrix columns owned by the calling image.
	OrdinalType getNumMyCols() const {
    if(!data().haveCol_)
      throw reportError("Column distribution not specified.", 3);
    return(getColDist().getNumMyEntries());
  }
	
	//! Returns the number of global nonzero diagonal entries, based on global row/column index comparisons.
	OrdinalType getNumGlobalDiagonals() const {
    return(globalSum(getNumMyDiagonals()));
  }
	
	//! Returns the number of local nonzero diagonal entries, based on global row/column index comparisons.
	OrdinalType getNumMyDiagonals() const {
    OrdinalType numDiagonals = Teuchos::OrdinalTraits<OrdinalType>::zero();
    typedef std::map<OrdinalType, ScalarType> OrdScalMap;
    typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;
    MapOfMaps& outermap = CisMatrixData_->indicesAndValues_;
    for(typename MapOfMaps::iterator i = outermap.begin(); i != outermap.end(); i++) {
      OrdScalMap& innermap = (*i).second;
      for(typename OrdScalMap::iterator j = innermap.begin(); j != innermap.end(); j++)
        if((*i).first == (*j).first)
          numDiagonals++;
    }
    return(numDiagonals);
  }
  
  //! Returns the current number of nonzero entries in specified global index on this image.
	OrdinalType getNumEntries(OrdinalType index) const {
    std::map<OrdinalType, ScalarType>& innermap = data().indicesAndValues_[index];
    return(innermap.size());
  }
	
	//! Returns the maximum number of nonzero entries across all rows/columns on all images.
	OrdinalType getGlobalMaxNumEntries() const {
    OrdinalType localMax = getMyMaxNumEntries();
    OrdinalType globalMax = Teuchos::OrdinalTraits<OrdinalType>::zero();
    ordinalComm().maxAll(&localMax, &globalMax, Teuchos::OrdinalTraits<OrdinalType>::one());
    return(globalMax);
  }
  
	//! Returns the maximum number of nonzero entries across all rows/columns on this image.
	OrdinalType getMyMaxNumEntries() const {
    OrdinalType currentMax = Teuchos::OrdinalTraits<OrdinalType>::zero();
    typedef std::map<OrdinalType, ScalarType> OrdScalMap;
    typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;
    MapOfMaps& outermap = CisMatrixData_->indicesAndValues_;
    for(typename MapOfMaps::iterator i = outermap.begin(); i != outermap.end(); i++) {
      OrdinalType newMax = (*i).second.size();
      if(newMax > currentMax)
        currentMax = newMax;
    }
    return(currentMax);
  }
	
	//! Returns the index base for global indices for this matrix.
	OrdinalType getIndexBase() const {
    return(getPrimaryDist().getIndexBase());
  }
  
  //! Returns false if this matrix shares data with another CisMatrix instance (or instances).
  //bool isSoleOwner() const;

  //! Returns true if this matrix is row-oriented, and false if this matrix is column-oriented.
  bool isRowOriented() const {
    return(CisMatrixData_->rowOriented_);
  }

	//! Returns true if the state of this matrix is post-fillComplete, and false if it is pre-fillComplete.
	bool isFillCompleted() const {
		return(data().fillCompleted_);
	}
  
  //! Returns the VectorSpace that describes the primary distribution in this matrix.
  /*! In a row-oriented matrix, this will be the row VectorSpace. In a column-oriented matrix, this will be the column VectorSpace. */
  VectorSpace<OrdinalType, ScalarType> const& getPrimaryDist() const {
    return(CisMatrixData_->primary_);
	}
  
  //! Returns the VectorSpace that describes the secondary distribution in this matrix.
  /* In a row-oriented matrix, this will be the column VectorSpace. In a column-oriented matrix, this will be the row VectorSpace. */
  VectorSpace<OrdinalType, ScalarType> const& getSecondaryDist() const {
    if(!CisMatrixData_->haveSecondary_)
      throw reportError("Secondary distribution is not currently defined.", 3);
    return(CisMatrixData_->secondary_);
  }

  //! Returns the VectorSpace that describes the row distribution in this matrix.
  VectorSpace<OrdinalType, ScalarType> const& getRowDist() const {
    if(!CisMatrixData_->haveRow_)
      throw reportError("Row distribution is not currently defined.", 3);
    if(isRowOriented())
      return(CisMatrixData_->primary_);
    else
      return(CisMatrixData_->secondary_);
  }

  //! Returns the VectorSpace that describes the column distribution in this matrix.
  VectorSpace<OrdinalType, ScalarType> const& getColumnDist() const {
    if(!CisMatrixData_->haveCol_)
      throw reportError("Column distribution is not currently defined.", 3);
    if(isRowOriented())
      return(CisMatrixData_->secondary_);
    else
      return(CisMatrixData_->primary_);
  }
  
  //! Returns the VectorSpace associated with the domain of this matrix.
  VectorSpace<OrdinalType, ScalarType> const& getDomainDist() const {
    if(!CisMatrixData_->haveDomain_)
      throw reportError("Domain distribution is not currently defined.", 3);
    return(CisMatrixData_->domain_);
  }
  
  //! Returns the VectorSpace associated with the range of this matrix.
  VectorSpace<OrdinalType, ScalarType> const& getRangeDist() const {
    if(!CisMatrixData_->haveRange_)
      throw reportError("Range distribution is not currently defined.", 3);
    return(CisMatrixData_->range_);
  }
  
  //! Returns the Platform object used by this matrix
  Platform<OrdinalType, ScalarType> const& platform() const {
		return(*data().platform_);
	}
  
  //@}
  
  //@{ \name I/O Methods
  
  // Print method, used by the overloaded << operator
  void print(ostream& os) const {
    os << "Orientation: " << (data().rowOriented_ ? "Row" : "Column") << "-oriented" << endl;
		os << "State: " << (isFillCompleted() ? "Post-fillComplete" : "Pre-fillComplete") << endl;
    os << "Secondary distribution defined? " << (data().haveSecondary_ ? "yes" : "no") << endl;
    os << "Domain distribution defined? " << (data().haveDomain_ ? "yes" : "no") << endl;
    os << "Range distribution defined? " << (data().haveRange_ ? "yes" : "no") << endl;

    os << "Contents:" << endl;
    typedef std::map<OrdinalType, ScalarType> OrdScalMap;
    typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;
    MapOfMaps& outermap = CisMatrixData_->indicesAndValues_;
    cout << setw(20) << "Primary Index" << setw(20) << "Secondary Index" << setw(10) << "Value" << endl;
    for(typename MapOfMaps::iterator i = outermap.begin(); i != outermap.end(); i++) {
      OrdScalMap& innermap = (*i).second;
      for(typename OrdScalMap::iterator j = innermap.begin(); j != innermap.end(); j++)
        cout << setw(15) << (*i).first << setw(18) << (*j).first << setw(15) << (*j).second << endl;
    }

		if(isFillCompleted()) {
			os << "---HB data---" << endl;
			os << "pntr_: ";
			for(typename std::vector<OrdinalType>::const_iterator i = data().pntr_.begin(); i != data().pntr_.end(); i++)
				os << *i << " ";
			os << endl << "indx_: ";
			for(typename std::vector<OrdinalType>::const_iterator i = data().indx_.begin(); i != data().indx_.end(); i++)
				os << *i << " ";
			os << endl << "values_: ";
			for(typename std::vector<ScalarType>::const_iterator i = data().values_.begin(); i != data().values_.end(); i++)
				os << *i << " ";
			os << endl;
		}

		os << "---Defined Distributions---" << endl;
		os << "Primary:" << endl << getPrimaryDist();
		if(data().haveSecondary_) os << "Secondary:" << endl << getSecondaryDist();
		if(data().haveDomain_) os << "Domain:" << endl << getDomainDist();
		if(data().haveRange_) os << "Range:" << endl << getRangeDist();
	}
  
  //@}
  
  //@{ \name Miscellaneous Methods
  
  //! Inlined bracket operator, non-const version.
  //ScalarType* operator[] (OrdinalType index);
  
  //! Inlined bracket operator, const version.
  //ScalarType const* operator[] (OrdinalType index) const;
  
	//! Assignment operator
	CisMatrix<OrdinalType, ScalarType>& operator = (CisMatrix<OrdinalType, ScalarType> const& Source) {
    CisMatrixData_ = Source.CisMatrixData_;
    return(*this);
  }
  
  //@}
	
private:

	Teuchos::RefCountPtr< CisMatrixData<OrdinalType, ScalarType> > CisMatrixData_;
  
  // convenience functions for returning inner data class, both const and nonconst versions.
  CisMatrixData<OrdinalType, ScalarType>& data() {return(*CisMatrixData_);}
  CisMatrixData<OrdinalType, ScalarType> const& data() const {return(*CisMatrixData_);}

  // convenience functions for comm instances
  Comm<ScalarType, OrdinalType> const& comm() const {return(*data().comm_);}
  Comm<OrdinalType, OrdinalType> const& ordinalComm() const {return(*data().ordinalComm_);}

  // convenience function for doing Comm::sumAll on one OT variable
  OrdinalType globalSum(OrdinalType localNum) const {
    OrdinalType globalNum = Teuchos::OrdinalTraits<OrdinalType>::zero();
    ordinalComm().sumAll(&localNum, &globalNum, Teuchos::OrdinalTraits<OrdinalType>::one());
    return(globalNum);
  }

  // internal functions for doing norms
  // primaryNorm computes the infinity-norm in a row-oriented matrix, and the one-norm in a column-oriented matrix.
  // secondaryNorm computes the one-norm in a row-oriented matrixm and the infinity-norm in a column-oriented matrix.
  // otherwise, both of these functions would be duplicated inside both normOne() and normInf()
  ScalarType primaryNorm() const {
    typedef std::map<OrdinalType, ScalarType> OrdScalMap;
    typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;
    MapOfMaps& outermap = CisMatrixData_->indicesAndValues_;
    ScalarType maxSum = Teuchos::ScalarTraits<ScalarType>::zero();

    // iterate through all rows/columns indices
    for(typename MapOfMaps::iterator i = outermap.begin(); i != outermap.end(); i++) {
      OrdScalMap& innermap = (*i).second;
      ScalarType currentSum = Teuchos::ScalarTraits<ScalarType>::zero();
      // for each row/column, sum up all the values
      for(typename OrdScalMap::iterator j = innermap.begin(); j != innermap.end(); j++)
        currentSum += (*j).second;
      // if we have a new maxsum, set it
      if(currentSum > maxSum)
        maxSum = currentSum;
    }
    // now do Comm call to get global maxSum
    ScalarType globalMax = Teuchos::ScalarTraits<ScalarType>::zero();
    comm().maxAll(&maxSum, &globalMax, Teuchos::OrdinalTraits<OrdinalType>::one());
    
    // update flops counter: length-1 for each row/column
    updateFlops(getNumMyNonzeros() - getPrimaryDist().getNumMyEntries());

    return(globalMax);
  }

  ScalarType secondaryNorm() const {
    // create temporary vector to store partial sums in
    Vector<OrdinalType, ScalarType> sums(getSecondaryDist());
    typedef std::map<OrdinalType, ScalarType> OrdScalMap;
    typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;
    MapOfMaps& outermap = CisMatrixData_->indicesAndValues_;
    // iterate through all rows/columns indices
    for(typename MapOfMaps::iterator i = outermap.begin(); i != outermap.end(); i++) {
      OrdScalMap& innermap = (*i).second;
      // for each element, add value to its column/row's partial sum
      for(typename OrdScalMap::iterator j = innermap.begin(); j != innermap.end(); j++)
        sums[(*j).first] += (*j).second;
    }
    // get local max
    ScalarType localMax = sums.maxValue();

    // now do Comm call to get global max
    ScalarType globalMax = Teuchos::ScalarTraits<ScalarType>::zero();
    comm().maxAll(&localMax, &globalMax, Teuchos::OrdinalTraits<OrdinalType>::one());
    
    // update flops counter: length-1 for each column/row
    updateFlops(getNumMyNonzeros() - getSecondaryDist().getNumMyEntries());

    return(globalMax);
  }

}; // CisMatrix class

} // Tpetra namespace

#include "Tpetra_CisMatrixData.hpp"

#endif // _TPETRA_CISMATRIX_HPP_
