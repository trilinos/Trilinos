// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_FieldContainer.hpp
\brief  Header file for utility class to provide lexicographical containers.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_FIELDCONTAINER_HPP
#define INTREPID_FIELDCONTAINER_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
  /** \class Intrepid::FieldContainer
  \brief Implementation of a templated lexicographical container for a multi-indexed scalar quantity. 
  FieldContainer object stores a multi-indexed scalar value using the lexicographical index ordering: the
  rightmost index changes first and the leftmost index changes last. FieldContainer can be viewed
  as a dynamic multidimensional array whose values can be accessed in two ways: by their
  multi-index or by their enumeration, using an overloaded [] operator. The enumeration of a value
  gives the sequential order of the multi-indexed value in the container. The number of indices, i.e.,
  the rank of the container is unlimited. For containers with ranks up to 5 many of the methods are 
  optimized for faster execution. An overloaded () operator is also provided for such low-rank containers
  to allow element access by multi-index without having to create an auxiliary array for the multi-index.
  */
  template<class Scalar>
  class FieldContainer {
  protected:
    
    /** \brief Array to store the multi-indexed quantity 
    */
    Teuchos::Array<Scalar> data_;
    
    /**\brief Array to store dimensions (dimensions) for the multi-indices. Admissible range (dimension)
     for  the k-th index is <var>0 <= i_k < dimensions_[k]</var>. Size of this array defines the rank of 
     the multi-indexed quantity, i.e., the number of its indices.  
    */
    Teuchos::Array<int> dimensions_;
    
    /** \brief Store the first 5 dimensions explicitely for faster access operations
    */
    int dim0_;
    int dim1_;
    int dim2_;
    int dim3_;
    int dim4_;
    
  public:
    
    /** \brief Default destructor.
    */
    ~FieldContainer() {};
    
    
    /** \brief Default constructor.
    */
    FieldContainer() {
      data_.resize(0);
      dimensions_.resize(0);
    };
    
    
    /** \brief Copy constructor.
      */
    FieldContainer(const FieldContainer& right);
 
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                              Constructors of FieldContainer class                          //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    /** \brief Creates an empty rank-1 FieldContainer with the specified dimension.
      
      \param dim0    [in]      - dimension for the only index 
      */
    FieldContainer(const int dim0);
    
    
    /** \brief Creates an empty rank-2 FieldContainer with the specified dimensions.
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      */
    FieldContainer(const int dim0,
                   const int dim1);

    
    /** \brief Creates an empty rank-3 FieldContainer with the specified dimensions.
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
      */
    FieldContainer(const int dim0,
                   const int dim1,
                   const int dim2);

    
    /** \brief Creates an empty rank-4 FieldContainer with the specified dimensions.
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
      \param dim3    [in]      - dimension for the 4th index 
      */
    FieldContainer(const int dim0,
                   const int dim1,
                   const int dim2,
                   const int dim3);
    
    
    /** \brief Creates an empty rank-5 FieldContainer with the specified dimensions.
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
      \param dim3    [in]      - dimension for the 4th index 
      \param dim4    [in]      - dimension for the 5th index 
      */
    FieldContainer(const int dim0,
                   const int dim1,
                   const int dim2,
                   const int dim3,
                   const int dim4);
    
    
    /** \brief Creates an empty FieldContainer of arbitrary rank, using dimensions specified in an
      array. The size of the input array implicitely defines the rank of the container and its
      capacity is defined by the specified dimensions. 
      
      \param dimensions[in]           - array with container dimensions
      */
    FieldContainer(const Teuchos::Array<int>& dimensions);
    
    
    /** \brief Creates a FieldContainer of arbitrary rank, using dimensions specified in an
      array, and fills it with data from another array. The size of the input array implicitely 
      defines the rank of the container and its capacity is defined by the specified dimensions. 
      
      \param dimensions[in]           - array with container dimensions
      \param data[in]                 - array with container values
      */
    FieldContainer(const Teuchos::Array<int>&    dimensions,
                   const Teuchos::Array<Scalar>& data);
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                            Access methods of FieldContainer class                          //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    /** \brief Return rank of the FieldContainer = number of indices used to tag the multi-indexed value
      */
    int getRank() const;
    
    
    /** \brief Returns size of the FieldContainer defined as the product of its dimensions.
      */
    int getSize() const;
    
    
    /** \brief Returns array with the dimensions of the container
      */
    void getAllDimensions(Teuchos::Array<int>& dimensions) const;
    
    
    /** \brief Returns the specified dimension
      
      \param whichDim     [in]      - order of the dimension we want to get
      */
    int getDimension(const int whichDim) const;
    
    
    /** \brief Returns enumeration of a value (its order relative to the container), based on its 
      multi-index, for rank-2 containers. 
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      */
    int getEnumeration(const int i0,
                       const int i1) const;
    
    
    /** \brief Returns enumeration of a value (its order relative to the container), based on its 
      multi-index, for rank-3 containers. 
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      */
    int getEnumeration(const int i0,
                       const int i1,
                       const int i2) const;
 
    
    /** \brief Returns enumeration of a value (its order relative to the container), based on its 
      multi-index, for rank-4 containers. 
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      \param i3         [in]        - 4th index
      */
    int getEnumeration(const int i0,
                       const int i1,
                       const int i2,
                       const int i3) const;
    
    
    /** \brief Returns enumeration of a value (its order relative to the container), based on its 
      multi-index, for rank-5 containers. 
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      \param i3         [in]        - 4th index
      \param i4         [in]        - 5th index
      */
    int getEnumeration(const int i0,
                       const int i1,
                       const int i2,
                       const int i3,
                       const int i4) const;
    
        
    /** \brief Returns enumeration of a value (its order relative to the container), based on its 
      multi-index, for containers of arbitrary rank. 
      
      \param multiIndex   [in]      - array representing a multi-index
    */
    int getEnumeration(const Teuchos::Array<int>& multiIndex) const;
    
    
    /** \brief Returns the multi-index corresponding to the specified enumeration of a value
      
      \param multiIndex   [out]       - array containg multi-index of the specified enumeration
      \param valueEnum    [in]        - enumeration of the value (its order relative to the container)
    */
    void getMultiIndex(Teuchos::Array<int>& multiIndex,
                       const int            valueEnum) const;
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                          Methods to shape (resize) a field container                       //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    
    /** \brief Resets FieldContainer to trivial container (one with rank = 0 and size = 0)
    */
    void empty();
    
    
    /** \brief Resizes FieldContainer to a rank-1 container with the specified dimension. 
      
      \param dim0    [in]      - dimension for the 1st index 
    */
    void resize(const int dim0);
    
    
    /** \brief Resizes FieldContainer to a rank-2 container with specified dimensions. 

      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
    */
    void resize(const int dim0,
                const int dim1);
    
    
    /** \brief Resizes FieldContainer to a rank-3 container with specified dimensions. 
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
    */
    void resize(const int dim0,
                const int dim1,
                const int dim2);

    
    /** \brief Resizes FieldContainer to a rank-4 container with specified dimensions. 
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
      \param dim3    [in]      - dimension for the 4th index 
    */
    void resize(const int dim0,
                const int dim1,
                const int dim2,
                const int dim3);
    
    
    /** \brief Resizes FieldContainer to a rank-5 container with specified dimensions. 
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
      \param dim3    [in]      - dimension for the 4th index 
      \param dim4    [in]      - dimension for the 5th index 
    */
    void resize(const int dim0,
                const int dim1,
                const int dim2,
                const int dim3,
                const int dim4);
    
    
    /** \brief Resizes FieldContainer to arbitrary rank container with dimensions specified in the
      input array. The size of this array implicitely defined the rank of the FieldContainer.
      
      \param newDimensions[in]          - new upper values for index ranges
    */
    void resize(const Teuchos::Array<int>& newDimensions);
    
    
    /** \brief Resizes FieldContainer to have the same rank and dimensions as another FieldContainer
      
      \param anotherContainer[in]          - a FieldContainer
      */
    void resize(const FieldContainer<Scalar>& anotherContainer);
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                     Methods to read and write values to FieldContainer                     //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    /** \brief Fills a FieldContainer with zeroes.
      */
    void storeZero(); 
    
    
    /** \brief Retrieve value by its multi-index. To retrieve it by enumeration use the overloaded [].
      
      \param multiIndex [in]            - array containing multi-index of the desired value
      */
    Scalar getValue(const Teuchos::Array<int>& multiIndex) const;
    
    
    /** \brief Assign value by its multi-index. 
      
      \param dataValue [in]             - value to be assigned
      \param multiIndex [in]            - multi-index of the value
      */
    void setValue(const Scalar               dataValue,
                  const Teuchos::Array<int>& multiIndex);
    
    
    /** \brief Assign value by its enumeration (order relative to the FieldContainer)
      
      \param dataValue [in]             - value to be assigned
      \param order     [in]             - enumeration of the value 
      */
    void setValue(const Scalar  dataValue,
                  const int     order);
    
    
    /**\brief Assign values from Teuchos::Array without changing rank and index ranges of FieldContainer. 
      Size of the input array must match the size of the container
      
      \param dataArray[in]               - new values
      */
    void setValues(const Teuchos::Array<Scalar>& dataArray);


    /** \brief Exposes data of FieldContainer, data can be modified.
    */
    Teuchos::Array<Scalar> & getData() {
      return data_;
    }    


    /** \brief Exposes data of FieldContainer, data cannot be modified.
    */
    const Teuchos::Array<Scalar> & getData() const {
      return data_;
    }    

    
    /** \brief Overloaded () operators for rank-2 containers.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
    */
    const Scalar& operator () (const int i0,
                               const int i1) const;
    
    Scalar&       operator () (const int i0,
                               const int i1);
    
    
    
    /** \brief Overloaded () operator for rank-3 containers.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
    */
    const Scalar& operator () (const int i0,
                               const int i1,
                               const int i2) const;
    
    Scalar&       operator () (const int i0,
                               const int i1,
                               const int i2);
    
    
    /** \brief Overloaded () operator for rank-4 containers.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      \param i3         [in]        - 4th index
    */
    const Scalar& operator () (const int i0,
                               const int i1,
                               const int i2,
                               const int i3) const;
    
    Scalar&       operator () (const int i0,
                               const int i1,
                               const int i2,
                               const int i3);
    
    
    /** \brief Overloaded () operator for rank-5 containers.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      \param i3         [in]        - 4th index
      \param i4         [in]        - 5th index
    */
    const Scalar& operator () (const int i0,
                               const int i1,
                               const int i2,
                               const int i3,
                               const int i4) const;
    
    Scalar&       operator () (const int i0,
                               const int i1,
                               const int i2,
                               const int i3,
                               const int i4);
    
        
    /** \brief   Overloaded [] operator. Returns value based on its enumeration
    */
    const Scalar & operator [] (const int address) const;
    
    
    /** \brief Assignment operator <var>*this = right</var>.
    */
    FieldContainer& operator  = (const FieldContainer& right);
    
  }; // end class FieldContainer
  
  
  //--------------------------------------------------------------------------------------------//
  //                                                                                            //
  //                        Function declarations related to FieldContainer                     //
  //                                                                                            //
  //--------------------------------------------------------------------------------------------//
  
  /** \relates FieldContainer
  Outputs a formated stream with FieldContainer data. For debugging purposes.
  */
  template<class Scalar>
    std::ostream& operator << (std::ostream& os, const FieldContainer<Scalar>& container);
  
} // end namespace Intrepid

// include templated definitions
#include <Intrepid_FieldContainerDef.hpp>

#endif
