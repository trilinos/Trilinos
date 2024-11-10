// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_FieldContainer.hpp
    \brief  Header file for utility class to provide multidimensional containers.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_FIELDCONTAINER_HPP
#define INTREPID_FIELDCONTAINER_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Shards_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_Assert.hpp"

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
  template<class Scalar, int ArrayTypeId=0>
  class FieldContainer {
  public:
    //! The template parameter of this class; the type of objects stored.
    typedef Scalar scalar_type;

  protected:
    
    /** \brief Array to store the multi-indexed quantity 
    */
    Teuchos::ArrayRCP<Scalar> data_;

    typedef typename Teuchos::ArrayRCP<Scalar>::iterator data_ptr_t;
    data_ptr_t data_ptr_;    

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
      data_ptr_ = Teuchos::NullIteratorTraits<data_ptr_t>::getNull();
      dimensions_.resize(0);
    } ;
    
    
    /** \brief Copy constructor.
    */
    FieldContainer(const FieldContainer& right);
 
    
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
    
    
    /** \brief Creates a FieldContainer of arbitrary rank, using dimensions specified in the
               <var><b>dimensions</b></var> array, and fills it by deep-copying data from a
               Teuchos::ArrayView array (which implicitly doubles as Teuchos::ArrayRCP or
               Teuchos::Array). If the input data array is a Teuchos::ArrayRCP, then '()' should
               be appended to it when calling this function. This forces direct conversion to a
               Teuchos::ArrayView, and prevents the call to the shallow-copy constructor that
               takes a Teuchos::ArrayRCP.
      
        \param dimensions[in]           - array with container dimensions
        \param data[in]                 - array with container values
    */
    FieldContainer(const Teuchos::Array<int>&        dimensions,
                   const Teuchos::ArrayView<Scalar>& data);


    /** \brief Creates a FieldContainer of arbitrary rank, using dimensions specified in the
               <var><b>dimensions</b></var> array, and wraps (shallow-copies) the data pointed
               to by the input Teuchos::ArrayRCP array. If a deep copy is desired instead,
               one can force the use of the constructor that takes Teuchos::ArrayView by
               appending () to the input Teuchos::ArrayRCP parameter. This forces direct
               conversion to a Teuchos::ArrayView.
      
        \param dimensions[in]           - array with container dimensions
        \param data[in]                 - array with container values
    */
    FieldContainer(const Teuchos::Array<int>&       dimensions,
                   const Teuchos::ArrayRCP<Scalar>& data);


    /** \brief Creates a FieldContainer of arbitrary rank, using dimensions specified in the
               <var><b>dimensions</b></var> array, and either wraps (shallow-copies) Scalar*
               <var><b>data</b></var>, or deep-copies it, based on the value of the parameter
               <var><b>deep_copy</b></var>. Memory management through FieldContainer, via
               its Teuchos::ArrayRCP data member, can be enabled.
      
        \param dimensions[in]           - array with container dimensions
        \param data[in]                 - array with container values
        \param deep_copy[in]            - if true, then deep-copy, otherwise shallow-copy; default: false
        \param owns_mem[in]             - if true, the field container will manage memory; default: false
    */
    FieldContainer(const Teuchos::Array<int>&    dimensions,
                   Scalar*                       data,
                   const bool                    deep_copy = false,
                   const bool                    owns_mem  = false);


    /** \brief Creates a FieldContainer either as a wrapper of the shards::Array<Scalar,shards::NaturalOrder>
               array <var><b>data</b></var>, or as its deep copy, based on the value of the parameter
               <var><b>deep_copy</b></var>. Memory management through FieldContainer, via
               its Teuchos::ArrayRCP data member, can be enabled.
      
        \param data[in]                 - array with container values
        \param deep_copy[in]            - if true, then deep-copy, otherwise shallow-copy; default: false
        \param owns_mem[in]             - if true, the field container will manage memory; default: false
    */
    FieldContainer(const shards::Array<Scalar,shards::NaturalOrder>&  data,
                   const bool                                         deep_copy = false,
                   const bool                                         owns_mem  = false);



    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                            Access methods of FieldContainer class                          //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    /** \brief Return rank of the FieldContainer = number of indices used to tag the multi-indexed value
      */
    int rank() const;
    
    
    /** \brief Returns size of the FieldContainer defined as the product of its dimensions.
      */
    int size() const;
    
    
    /** \brief Returns array with the dimensions of the container
      */
    template<class Vector>
    void dimensions(Vector& dimensions) const;
    
    
    /** \brief Returns the specified dimension
      
      \param whichDim     [in]      - order of the dimension we want to get
      */
    int dimension(const int whichDim) const;

    
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
    
    
    /** \brief Returns the multi-index of a value, based on its enumeration, as a vector, for
      containers of arbitrary rank. The template argument must support the following subset of
      std::vector interface:
      \li   Vector.size()
      \li   Vector.resize()
      \li   Vector[]
      
      \param multiIndex   [out]       - vector containg multi-index of the specified enumeration
      \param valueEnum    [in]        - enumeration of the value (its order relative to the container)
    */
    template<class Vector>
    void getMultiIndex(Vector&    multiIndex,
                       const int  valueEnum) const;
    
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
    void resize(const FieldContainer<Scalar, ArrayTypeId>& anotherContainer);
    
    
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
      \param spaceType       [in]        - type of the function space whose basis will be evaluated
      \param operatorType    [in]        - type of the operator that will be applied to the basis
      \param spaceDim        [in]        - dimension of the ambient space
    */
    void resize(const int             numPoints,
                const int             numFields,
                const EFunctionSpace  spaceType,
                const EOperator       operatorType,
                const int             spaceDim);
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                     Methods to read and write values to FieldContainer                     //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//

    
    /** \brief Initializes a field container by assigning <var>value</var> to all its elements.  
      */
    void initialize(const Scalar value = 0); 
        
    
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
    void setValues(const Teuchos::ArrayView<Scalar>& dataArray);
    
    
    /** \brief Fills an existing FieldContainer with Scalars referenced by <var>dataPtr</var> without
      changing rank and dimensions of the container. Number of data must match the size of the container.
      
      \param dataPtr  [in]               - new values
      \param numData  [in]               - number of values
    */
    void setValues(const Scalar* dataPtr, 
                   const int     numData); 
    
    
    /** \brief Exposes data of FieldContainer, data can be modified.
    */
    Teuchos::ArrayRCP<Scalar> getData() {
      return data_;
    }    


    /** \brief Exposes data of FieldContainer, data cannot be modified.
    */
    Teuchos::ArrayRCP<const Scalar> getData() const {
      return data_;
    }    


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
  template<class Scalar, int ArrayTypeId>
    std::ostream& operator << (std::ostream& os, const FieldContainer<Scalar, ArrayTypeId>& container);
 
  
} // end namespace Intrepid

/*
template<class FadType, class Scalar>
struct Return_Type< Intrepid::FieldContainer<FadType>, Scalar>{
      typedef FadType& return_type;
      typedef FadType  const_return_type;
};

template<class FadType, class Scalar>
struct Return_Type<const Intrepid::FieldContainer<FadType>, Scalar>{
      typedef FadType& return_type;
      typedef FadType  const_return_type;
};
*/

// include templated definitions
#include <Intrepid_FieldContainerDef.hpp>

#endif

/**
\page md_array_page               Multi-dimensional array (MD array) template arguments
 
 
 \section md_array_intro_sec       Introduction
 
 This page describes basic requirements for multi-dimensional array (MD array) template arguments in 
 Intrepid. MD array is a fundamental concept for managing a wide range of numerical data arising in
 PDE-based simulations. An MD array is a collection of members of a given scalar data type that are 
 identified through a multi-index. Therefore, allowing for an MD array template argument provides a
 convenient mechanism to share numerical data between Intrepid and user applications. 
 
 The scalar data type members of a multi-dimensional array are stored in a contiguous block of memory.  
 Contiguous memory storage is not necessary to the concept of an array; however, in order to achieve 
 the best performance from cache-based computer memory systems contiguity of member storage has proven 
 to be essential.
 
 
 \section md_array_def_sec        Definitions
 
 The following rules and definitions apply to all MD arrays used as template arguments in Intrepid.
 
 \li A scalar value that is uniquely identified by a set of interegs <var>{i0,i1,...,iN}</var>  is called 
 a multi-indexed value;
 \li The set of all multi-indexed values such that <var>0 <= ik < dim_k</var> is called MD array;
 \li The integer <var>ik</var> is the kth index of the MD array; the N-tuple <var>{i0,...,iN}</var>
 is the multi-index of the MD array;
 \li The integer <var>dim_k</var> is the kth dimension of the MD array; the N-tuple <var>{dim_0,...,dim_N}</var>
 is the multi-dimension of the MD array
 \li The integer <var>N+1</var> is the rank of the MD array;
 \li A map <var>{i0,...,iN} -> {0,1,2,...}</var> from the set of all multi-indices to the set of the
 natural numbers is called enumeration of the MD array;
 \li The numerical position of an element indexed by <var>{i0,...,iN}</var>, established by the
 enumeration is called ordinal number of that element or simply ordinal.
 
 Enumeration of all MD arrays passed as template arguments to Intrepid must follow the natural 
 lexicographical order: the leftmost index <var>i0</var> changes last and the rightmost index 
 <var>iN</var> changes first. In summary, an MD array pased to Intrepid should comply with the
 following rules:
 
 \li the indices are zero-based;
 \li dimensions are counted from 0; e.g., a rank-4 array has dimensions <var>{dim_0,dim_1,dim_2,dim_3}</var>
 \li the enumeration is induced by the natural lexicographical order;
 \li the MD array is not strongly type on its dimensions and rank.
 
 
 \section md_array_interface_sec  MD Array interface
 
 An MD array type passed as a template argument to Intrepid is expected to implement the following 
 minimal interface:
 
 \li int rank()                      - returns number of dimensions
 \li int dimension(dim_k)            - returns the kth dimension (dimensions are dim0, dim1, etc.)
 \li int size()                      - returns size, i.e., dim0*dim1*...*dim_k
 \li const Scalar& operator(i,j,...,k)  - const accessor using multi-index
 \li       Scalar& operator(i,j,...,k)  - non-const accessor using multi-index
 \li const Scalar& operator[i]          - const accessor using the ordinal of the array element
 \li       Scalar& operator[i]          - non-const accessor using the ordinal of the array element
 
 
 \section md_array_notation_sec   MD Array notation and examples
 
 In addition to the generic index and dimension notation <var>ik</var> and <var>dim_k</var> it is 
 convenient to introduce data-specific notation for indices and dimensions of MD arrays that recur
 in PDE-based simulation codes. 
 \verbatim
 |-------------------------------------------------------------------------------------------------|
 |   Index type              | Dimension |  Description                                            |
 |---------------------------|-----------|---------------------------------------------------------|
 |   point                   |     P     |  number of points stored in an MD array                 |
 |   vertex                  |     V     |  number of nodes stored in an MD aray                |
 |   field                   |     F     |  number of fields stored in an MD array                 |
 |   basis field             |     B     |  number of basis fields stored in an MD array           |
 |   cell                    |     C     |  number of cells stored in an MD array                  |
 |   field coordinate        |     D     |  space dimension                                        |
 |   derivative ordinal      |     K     |  cardinality of the set of kth derivatives              |
 |-------------------------------------------------------------------------------------------------|
 \endverbatim
 
 \remarks 
 \li The totality of all derivatives whose order equals k (OPERATOR_Dk in Intrepid) forms a multiset;
 see http://mathworld.wolfram.com/Multiset.html . In Intrepid this multiset is enumerated using the 
 lexicographical order of the partial derivatives; see getDkEnumeration() for details.
 
 
 This notation is used throughout Intrepid's documentation as a means to add further clarity to the 
 specifications of MD array passed to and returned from Intrepid methods. The following is a list of
 typical MD arrays that arise in PDE-based simulation codes:
 \verbatim
 |-------------------------------------------------------------------------------------------------|
 | Rank | Multi-dimension | Multi-index    | Description                                           |
 |-------------------------------------------------------------------------------------------------|
 |  1   | (P)             | (p)            | Scalar (rank 0) field evaluated at P points           |
 |  2   | (P,D)           | (p,d)          | Vector (rank 1) field evaluated at P points           |
 |  3   | (P,D,D)         | (p,d,d)        | Tensor (rank 2) field evaluated at P points           |           
 |-------------------------------------------------------------------------------------------------|
 |  2   | (P,F)           | (p,f)          | F scalar fields evaluated at P points                 |
 |  3   | (P,F,D)         | (p,f,d)        | F vector fields evaluated at P points                 |
 |  4   | (P,F,D,D)       | (p,f,d,d)      | F tensor fields evaluated at P points                 |
 |-------------------------------------------------------------------------------------------------|
 |  3   | (P,F,K)         | (p,f,k)        | kth deriv. of F scalar fields evaluated at P points   |       
 |  4   | (P,F,D,K)       | (p,f,d,k)      | kth deriv. of F vector fields evaluated at P points   |       
 |  5   | (P,F,D,D,K)     | (p,f,d,d,k)    | kth deriv. of F tensor fields evaluated at P points   |           
 |-------------------------------------------------------------------------------------------------|
 |  3   | (C,V,D)         | (c,v,d )       | Vertex coords. of C cells having V vertices each      |
 |  3   | (C,P,D)         | (c,p,d )       | Coords. of C*P points in C cells, P per cell          |
 |-------------------------------------------------------------------------------------------------| 
 \endverbatim
 
 
 \section md_array_intrepid_sec   MD Array implementaion in Intrepid
 
 The FieldContainer class provides an implementation of an MD array type that is used throughout Intrepid.
 A FieldContainer object is templated on a Scalar type. Its rank and dimensions are runtime parameters,
 i.e., a FieldContainer is not strongly typed on rank and dimension.
 
 */
 
 






#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

