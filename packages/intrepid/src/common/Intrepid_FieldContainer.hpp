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
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
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
    
    /** \brief 1st dimension of the array */
    int dim0_;
    
    /** \brief 2nd dimension of the array */
    int dim1_;
    
    /** \brief 3rd dimension of the array */
    int dim2_;
    
    /** \brief 4th dimension of the array */
    int dim3_;
    
    /** \brief 5th dimension of the array */
    int dim4_;
    
  public:
    
    /** \brief Default destructor.
    */
    ~FieldContainer() {};
    
    
    /** \brief Default constructor.
    */
    FieldContainer() : dim0_(0), dim1_(0), dim2_(0), dim3_(0), dim4_(0) 
    {
      data_.resize(0);
      dimensions_.resize(0);
    } ;
    
    
    /** \brief Copy constructor.
      */
    FieldContainer(const FieldContainer& right);
 
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                              Constructors of FieldContainer class                          //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    /** \brief Creates a rank-1 FieldContainer with the specified dimension, initialized by 0.
      
      \param dim0    [in]      - dimension for the only index 
      */
    FieldContainer(const int dim0);
    
    
    /** \brief Creates a rank-2 FieldContainer with the specified dimensions, initialized by 0.
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      */
    FieldContainer(const int dim0,
                   const int dim1);

    
    /** \brief Creates a rank-3 FieldContainer with the specified dimensions, initialized by 0.
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
      */
    FieldContainer(const int dim0,
                   const int dim1,
                   const int dim2);

    
    /** \brief Creates a rank-4 FieldContainer with the specified dimensions, initialized by 0.
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
      \param dim3    [in]      - dimension for the 4th index 
      */
    FieldContainer(const int dim0,
                   const int dim1,
                   const int dim2,
                   const int dim3);
    
    
    /** \brief Creates a rank-5 FieldContainer with the specified dimensions, initialized by 0.
      
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
    
    
    /** \brief Creates a  FieldContainer of arbitrary rank,, initialized by 0, using dimensions 
      specified in an array. The size of the input array implicitely defines the rank of the 
      container and its capacity is defined by the specified dimensions. 
      
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
      multi-index, for rank-1 containers. 
      
      \param i0         [in]        - 1st index
      */
    int getEnumeration(const int i0) const;
    
    
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
    
    
    /** \brief  Returns the multi-index of a value, based on its enumeration, as a list, for
      rank-1 containers.
      
      \param i0         [out]        - 1st index
      \param valueEnum  [in]         - enumeration of the value (its order relative to the container)      
    */
    void getMultiIndex(int &     i0,
                       const int valueEnum) const;
    
    
    /** \brief  Returns the multi-index of a value, based on its enumeration, as a list, for
      rank-2 containers.
      
      \param i0         [out]        - 1st index
      \param i1         [out]        - 2nd index
      \param valueEnum  [in]         - enumeration of the value (its order relative to the container)      
    */
    void getMultiIndex(int &     i0,
                       int &     i1,
                       const int valueEnum) const;
    
    
    /** \brief  Returns the multi-index of a value, based on its enumeration, as a list, for
      rank-3 containers.
      
      \param i0         [out]        - 1st index
      \param i1         [out]        - 2nd index
      \param i2         [out]        - 3rd index
      \param valueEnum  [in]         - enumeration of the value (its order relative to the container)      
      */
    void getMultiIndex(int &     i0,
                       int &     i1,
                       int &     i2,
                       const int valueEnum) const;
    
    
    /** \brief  Returns the multi-index of a value, based on its enumeration, as a list, for
      rank-4 containers.
      
      \param i0         [out]        - 1st index
      \param i1         [out]        - 2nd index
      \param i2         [out]        - 3rd index
      \param i3         [out]        - 4th index
      \param valueEnum  [in]         - enumeration of the value (its order relative to the container)      
      */
    void getMultiIndex(int &     i0,
                       int &     i1,
                       int &     i2,
                       int &     i3,
                       const int valueEnum) const;
    
    
    /** \brief  Returns the multi-index of a value, based on its enumeration, as a list, for
      rank-5 containers.
      
      \param i0         [out]        - 1st index
      \param i1         [out]        - 2nd index
      \param i2         [out]        - 3rd index
      \param i3         [out]        - 4th index
      \param i4         [out]        - 5th index
      \param valueEnum  [in]         - enumeration of the value (its order relative to the container)      
      */
    void getMultiIndex(int &     i0,
                       int &     i1,
                       int &     i2,
                       int &     i3,
                       int &     i4,
                       const int valueEnum) const;
    
    
    /** \brief Returns the multi-index of a value, based on its enumeration, as an array, for
      containers of arbitrary rank.
      
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
    
    
    /** \brief Clears FieldContainer to trivial container (one with rank = 0 and size = 0)
    */
    void clear();
    
    
    /** \brief Resizes FieldContainer to a rank-1 container with the specified dimension, initialized by 0. 
      
      \param dim0    [in]      - dimension for the 1st index 
    */
    void resize(const int dim0);
    
    
    /** \brief Resizes FieldContainer to a rank-2 container with specified dimensions, initialized by 0. 

      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
    */
    void resize(const int dim0,
                const int dim1);
    
    
    /** \brief Resizes FieldContainer to a rank-3 container with specified dimensions, initialized by 0. 
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
    */
    void resize(const int dim0,
                const int dim1,
                const int dim2);

    
    /** \brief Resizes FieldContainer to a rank-4 container with specified dimensions, initialized by 0. 
      
      \param dim0    [in]      - dimension for the 1st index 
      \param dim1    [in]      - dimension for the 2nd index 
      \param dim2    [in]      - dimension for the 3rd index 
      \param dim3    [in]      - dimension for the 4th index 
    */
    void resize(const int dim0,
                const int dim1,
                const int dim2,
                const int dim3);
    
    
    /** \brief Resizes FieldContainer to a rank-5 container with specified dimensions, initialized by 0. 
      
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
    
    
    /** \brief Resizes FieldContainer to arbitrary rank container, initialized by 0, with dimensions 
      specified in the input array. The size of this array implicitely defined the rank of the FieldContainer.
      
      \param newDimensions[in]          - new upper values for index ranges
    */
    void resize(const Teuchos::Array<int>& newDimensions);
        
    
    /** \brief Resizes FieldContainer to have the same rank and dimensions as another FieldContainer,
      and initializes by 0.
      
      \param anotherContainer[in]          - a FieldContainer
      */
    void resize(const FieldContainer<Scalar>& anotherContainer);
    
    
    /** \brief Resizes FieldContainer to a container whose rank depends on the specified field and 
      operator types and the space dimension, initialized by 0. The admissible combinations of these 
      arguments, the rank of the resulitng container and its dimensions are summarized in the following table:
      \verbatim
      |--------------------|-------------------|-------------------|-------------------|
      |operator/field rank |       rank 0      | rank 1 2D/3D      | rank 2 2D/3D      |
      |--------------------|-------------------|-------------------|-------------------|
      |       VALUE        | (P,F)             | (P,F,D)           | (P,F,D,D)         |
      |--------------------|-------------------|-------------------|-------------------|
      |     GRAD, D1       | (P,F,D)           | (P,F,D,D)         | (P,F,D,D,D)       |
      |--------------------|-------------------|-------------------|-------------------|
      |        CURL        | (P,F,D) (undef3D) | (P,F)/(P,F,D)     | (P,F,D)/(P,F,D,D) |
      |--------------------|-------------------|-------------------|-------------------|
      |        DIV         | (P,F,D) (only 1D) | (P,F)             | (P,F,D)           |
      |--------------------|-------------------|-------------------|-------------------|
      |    D1,D2,..,D10    | (P,F,K)           | (P,F,D,K)         | (P,F,D,D,K)       |
      |--------------------|-------------------|-------------------|-------------------|
      
      |------|----------------------|---------------------------|
      |      |         Index        |         Dimension         |
      |------|----------------------|---------------------------|
      |   P  |         point        |  0 <= P < numPoints       |
      |   F  |         field        |  0 <= F < numFields       |
      |   D  |   field coordinate   |  0 <= D < spaceDim        |
      |   K  |   enumeration of Dk  |  0 <= K < DkCardinality   |
      |------|----------------------|---------------------------|
      \endverbatim
      \remarks 
      \li Enumeration of Dk (derivatives of total order k) follows the lexicographical order of 
      the partial derivatives; see getDkEnumeration() for details.
      
      
      \param numPoints       [in]        - number of evaluation points
      \param numFields       [in]        - number of fields that will be evaluated
      \param fieldType       [in]        - type of the field whose basis will be evaluated
      \param operatorType    [in]        - type of the operator that will be applied to the basis
      \param spaceDim        [in]        - dimension of the ambient space
    */
    void resize(const int       numPoints,
                const int       numFields,
                const EField    fieldType,
                const EOperator operatorType,
                const int       spaceDim);
    
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
    
    
    /**\brief Fills an existing FieldContainer with Scalars stored in a Teuchos::Array without changing 
      rank and dimensions of the container. Size of the input array must match the size of the container.
      
      \param dataArray[in]               - new values
    */
    void setValues(const Teuchos::Array<Scalar>& dataArray);
    
    
    /** \brief Fills an existing FieldContainer with Scalars referenced by <var>dataPtr</var> without
      changing rank and dimensions of the container. Number of data must match the size of the container.
      
      \param dataPtr  [in]               - new values
      \param numData  [in]               - number of values
    */
    void setValues(const Scalar* dataPtr, 
                   const int     numData); 
    
    
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

    template<class ArrayType>
    void contractScalar(ArrayType &                     outputValues,
                        const FieldContainer<Scalar> &  leftValues,
                        const ECompEngine               compEngine) const;
    
    
    template<class ArrayType>
    void contractVector(ArrayType &                     outputValues,
                        const FieldContainer<Scalar> &  leftValues,
                        const ECompEngine               compEngine) const;

    
    template<class ArrayType>
    void contractTensor(ArrayType &                     outputValues,
                        const FieldContainer<Scalar> &  leftValues,
                        const ECompEngine               compEngine) const;
    
    
    template<class ArrayType>
    void contractScalarData(ArrayType &        outputValues,
                            const ArrayType &  inputData,
                            const ECompEngine  compEngine) const;
 
    
    
    template<class ArrayType>
    void contractVectorData(ArrayType &        outputValues,
                            const ArrayType &  inputData,
                            const ECompEngine  compEngine) const;

    
    template<class ArrayType>
    void contractTensorData(ArrayType &        outputValues,
                            const ArrayType &  inputData,
                            const ECompEngine  compEngine) const;
    
    template<class ArrayType>
    void multiplyScalarData(const ArrayType &  inputData);
    
    
    template<class ArrayType>
    void multiplyVectorData(FieldContainer<Scalar> outputValues,
                            const ArrayType &  inputData);

    
    
    /** \brief Overloaded () operators for rank-1 containers. Data <strong>cannot</strong> be modified.
      
      \param i0         [in]        - 1st index
    */
    const Scalar& operator () (const int i0) const;
    
    /** \brief Overloaded () operators for rank-1 containers. Data <strong>can</strong> be modified.
      
      \param i0         [in]        - 1st index
    */
    Scalar&       operator () (const int i0);
    
    
    /** \brief Overloaded () operators for rank-2 containers. Data <strong>cannot</strong> be modified.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
    */
    const Scalar& operator () (const int i0,
                               const int i1) const;
    
    /** \brief Overloaded () operators for rank-2 containers. Data <strong>can</strong> be modified.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
    */
    Scalar&       operator () (const int i0,
                               const int i1);
    
    
    /** \brief Overloaded () operator for rank-3 containers. Data <strong>cannot</strong> be modified.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
    */
    const Scalar& operator () (const int i0,
                               const int i1,
                               const int i2) const;
    
    /** \brief Overloaded () operator for rank-3 containers. Data <strong>can</strong> be modified.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      */
    Scalar&       operator () (const int i0,
                               const int i1,
                               const int i2);
    
    
    /** \brief Overloaded () operator for rank-4 containers. Data <strong>cannot</strong> be modified.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      \param i3         [in]        - 4th index
    */
    const Scalar& operator () (const int i0,
                               const int i1,
                               const int i2,
                               const int i3) const;
    
    /** \brief Overloaded () operator for rank-4 containers. Data <strong>can</strong> be modified.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      \param i3         [in]        - 4th index
    */
    Scalar&       operator () (const int i0,
                               const int i1,
                               const int i2,
                               const int i3);
    
    
    /** \brief Overloaded () operator for rank-5 containers. Data <strong>cannot</strong> be modified.
      
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
    
    /** \brief Overloaded () operator for rank-5 containers. Data <strong>can</strong> be modified.
      
      \param i0         [in]        - 1st index
      \param i1         [in]        - 2nd index
      \param i2         [in]        - 3rd index
      \param i3         [in]        - 4th index
      \param i4         [in]        - 5th index
    */
    Scalar&       operator () (const int i0,
                               const int i1,
                               const int i2,
                               const int i3,
                               const int i4);
    
        
    /** \brief   Overloaded [] operator. Returns value based on its enumeration.
      Data <strong>cannot</strong> be modified.
    */
    const Scalar & operator [] (const int address) const;

    
    /** \brief   Overloaded [] operator. Returns value based on its enumeration.
      Data <strong>can</strong> be modified.
      */
    Scalar &       operator [] (const int address);
    
    
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
 
  
  
  
  
  /** \relates VarContainer
    Field rank is defined as the number of indices needed to specify the value. Intrepid provides
    three basic types of fields: scalar, vectors and tensors with ranks 0,1, and 2 respectively. 
    The scalar field types are FIELD_FORM_0 and FIELD _FORM_3. The vector field types are 
    FIELD_FORM_1, FIELD_FORM_2 and FIELD_VECTOR. The tensor field type is FIELD_TENSOR.
    
    \param fieldType         [in]     - field type whose rank we want to know
    */
  int getFieldRank(const EField fieldType); 
  
  
  
  /**\relates varContainer
    Returns rank of an operator. When an operator acts on a field with a certain rank, the result
    is a field with another rank. Operator rank is defined as the difference between the ranks of
    the output field and the input field:
    
    Rank(OPERATOR) = Rank(OPERATOR(FIELD)) - Rank(FIELD)
    
    This function is used by VarContainer class to determine how many indices are needed to store
    the values of OPERATOR(FIELD). By knowing the rank of the operator and the rank of the field we
    can figure out the rank of the result:
    
    Rank(OPERATOR(FIELD)) = Rank(OPERATOR) + Rank(FIELD)
    
    Operator ranks are defined in the following table: (X denotes undefined, below slash means 3D).
    By default, in 1D any operator other than VALUE has rank 1, i.e., GRAD, CURL and DIV default to
  d/dx and Dk are the higher-order derivatives d^k/dx^k. Only scalar functions are allowed in 1D.
    
    \verbatim
    |--------|------|--------------------------|----------|----------|----------|----------|
    | field  | rank |  field type ENUM         |   VALUE  | GRAD, Dk |   CURL   |    DIV   |
    |--------|------|--------------------------|----------|----------|----------|----------|
    | scalar |   0  |  FORM_0, FORM_3          |     0    |     1    | 3-dim/X  |     X    |
    | vector |   1  |  FORM_1, FORM_2, VECTOR  |     0    |     1    | dim - 3  |    -1    |
    | tensor |   2  |  TENSOR                  |     0    |     1    | dim - 3  |    -1    |
    |--------|------|--------------------------|----------|----------|----------|----------|
    |   1D   |   0  |  FORM_0, FORM_3 only     |     0    |     1    |     1    |     1    |
    |--------|------|--------------------------|----------|----------|----------|----------|
    \endverbatim
    
    \param operatorType     [in]    - the operator acting on the field
    \param fieldRank        [in]    - rank of the field being acted upon (use getFieldRank)
    \param spaceDim         [in]    - dimension of the space where all happens (1,2 or 3)
    */
  int getOperatorRank(const EOperator operatorType,
                      const int       fieldRank,
                      const int       spaceDim);
  
  
  
  /**\relates VarContainer 
    Returns order of an operator: ranges from 0 (for OPERATOR_VALUE) to 10 (OPERATOR_D10)
    
    \param operatorType       [in]    - type of the operator whose order we want to know
    */
  int getOperatorOrder(const EOperator operatorType);
  
  
  
  /** \relates VarContainer
    Returns enumeration of a partial derivative of order k based on multiplicities of the partial
    derivatives dx, dy and dz. Enumeration follows the lexicographical order of the partial 
    derivatives multiplicities. For example, the 10 derivatives of order 3 in 3D are enumerated as:
    
    D3 = {(3,0,0), (2,1,0), (2,0,1), (1,2,0), (1,1,1), (1,0,2), (0,3,0), (0,2,1), (0,1,2), (0,0,3)}
  
  Enumeration formula for lexicographically ordered partial derivatives of order k is given by
    
    <table>   
    <tr> <td>\f$i(xMult)               = 0\f$</td>              <td>in 1D (only 1 derivative)</td> </tr>
    <tr> <td>\f$i(xMult, yMult)        = yMult\f$</td>                              <td>in 2D</td> </tr>
    <tr> <td>\f$i(xMult, yMult, zMult) = zMult + \displaystyle\sum_{r = 0}^{k - xMult} r\f$</td> <td>in 3D</td> </tr>
    </table>
    
    where the order k of Dk is implicitly defined by xMult + yMult + zMult. Space dimension is
    implicitly defined by the default values of the multiplicities of y and z derivatives.
  
  \param xMult            [in]    - multiplicity of the partial derivative in x
    \param yMult            [in]    - multiplicity of the partial derivative in y (default = -1)
  \param zMult            [in]    - multiplicity of the partial derivative in z (default = -1)
  */
  int getDkEnumeration(const int xMult,
                       const int yMult = -1,
                       const int zMult = -1);
  
  
  
  /** \relates VarContainer
    Returns multiplicities of the partials in x, y, and z based on the enumeration of the partial
    derivative, its order and the space dimension. Inverse of the getDkEnumeration() method.
    
    \param partialMult      [out]    - array with partial derivative multiplicities
    \param derivativeEnum   [in]     - enumeration of the partial derivative
    \param derivativeOrder  [in]     - order of Dk
    \param spaceDim         [in]     - space dimension
    */
  void getDkMultiplicities(Teuchos::Array<int>& partialMult,
                           const int derivativeEnum,
                           const int derivativeOrder,
                           const int spaceDim);
  
  
  
  /** \relates VarContainer
    Returns cardinality of Dk, i.e., the number of all derivatives of order k. The set of all 
    partial derivatives of order k is isomorphic to the set of all multisets of cardinality k with 
    elements taken from the sets {x}, {x,y}, and {x,y,z} in 1D, 2D, and 3D respectively. 
    For example, the partial derivative
    
    \f$\displaystyle D\{1,2,1\}f = \frac{d^4 f}{dx dy^2 dz}\f$   maps to multiset
    \f$\{x, y, z\}\f$ with multiplicities \f$\{1,2,1\}\f$ 
    
    The number of all such multisets is given by the binomial coefficient
    
    \f$
    \begin{pmatrix}
  spaceDim + k - 1 \\
    spaceDim - 1
    \end{pmatrix}
  \f$
    
1D: cardinality = 1\n
2D: cardinality = k + 1\n
3D: cardinality = (k + 1)*(k + 2)/2
    
    \param derivativeOrder  [in]     - order of Dk
    \param spaceDim         [in]     - space dimension
    */
  int getDkCardinality(const int derivativeOrder,
                       const int spaceDim);
  
  //===========================================================================//
  //                                                                           //
  //              End Function declarations related to VarContainer            //
  //                                                                           //
  //===========================================================================//
  
} // end namespace Intrepid

// include templated definitions
#include <Intrepid_FieldContainerDef.hpp>

#endif
