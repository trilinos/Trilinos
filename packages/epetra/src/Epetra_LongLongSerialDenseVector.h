/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_LONGLONGSERIALDENSEVECTOR_H
#define EPETRA_LONGLONGSERIALDENSEVECTOR_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_LongLongSerialDenseMatrix.h"

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

//! Epetra_LongLongSerialDenseVector: A class for constructing and using dense vectors.

/*! The Epetra_LongLongSerialDenseVector class enables the construction and use of integer-valued,
    dense vectors.  It derives from the Epetra_LongLongSerialDenseMatrix class.

The Epetra_LongLongSerialDenseVector class is intended to provide convenient vector notation but derives all signficant
functionality from Epetra_LongLongSerialDenseMatrix.

<b>Constructing Epetra_LongLongSerialDenseVector Objects</b>

There are three Epetra_LongLongSerialDenseVector constructors.  The first constructs a zero-length object which should be made
to appropriate length using the Size() or Resize() functions and then filled with the [] or () operators.
The second constructs an object sized to the dimension specified, which should be filled with the [] or () operators.
The third is a constructor that accepts user
data as a 1D array, and the fourth is a copy constructor. The third constructor has
two data access modes (specified by the Epetra_DataAccess argument):
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the object.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and
only use the View mode in a secondary optimization phase.

<b>Extracting Data from Epetra_LongLongSerialDenseVector Objects</b>

Once a Epetra_LongLongSerialDenseVector is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.

*/


//=========================================================================
class EPETRA_LIB_DLL_EXPORT Epetra_LongLongSerialDenseVector : public Epetra_LongLongSerialDenseMatrix{

  public:

  //! Default constructor; defines a zero size object.
  /*!
    Epetra_LongLongSerialDenseVector objects defined by the default constructor should be sized with the
    Size() or Resize functions.
    Values should be defined by using the [] or () operators.
   */
  Epetra_LongLongSerialDenseVector();

  //! Sized constructor; defines a variable-sized object
  /*!
    \param In
           Length - Length of vector.

    Epetra_LongLongSerialDenseVector objects defined by the sized constructor are already sized to the
    dimension given as a parameter. All values are initialized to 0. Calling this constructor
    is equivalent to using the default constructor, and then calling the Size function on it.
    Values should be defined by using the [] or () operators.
   */
  Epetra_LongLongSerialDenseVector(int Length_in);

  //! Set object values from one-dimensional array.
  /*!
    \param In
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
           Values - Pointer to an array of integer numbers containing the values.
    \param In
           Length - Length of vector.

     See Detailed Description section for further discussion.
  */
  Epetra_LongLongSerialDenseVector(Epetra_DataAccess CV_in, long long* Values_in, int Length_in);

  //! Epetra_LongLongSerialDenseVector copy constructor.

  Epetra_LongLongSerialDenseVector(const Epetra_LongLongSerialDenseVector& Source);

  //! Set length of a Epetra_LongLongSerialDenseVector object; init values to zero.
  /*!
    \param In
           Length - Length of vector object.

     Allows user to define the dimension of a Epetra_LongLongSerialDenseVector. This function can
     be called at any point after construction.  Any values that were previously in this object are
     destroyed and the resized vector starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int Size(int Length_in) {return(Epetra_LongLongSerialDenseMatrix::Shape(Length_in, 1));};

  //! Resize a Epetra_LongLongSerialDenseVector object.
  /*!
    \param In
           Length - Length of vector object.

     Allows user to define the dimension of a Epetra_LongLongSerialDenseVector. This function can
     be called at any point after construction.  Any values that were previously in this object are
     copied into the new size.  If the new shape is smaller than the original, the first Length values
     are copied to the new vector.

    \return Integer error code, set to 0 if successful.
  */
  int Resize(int Length_in) {return(Epetra_LongLongSerialDenseMatrix::Reshape(Length_in, 1));};

  //! Epetra_LongLongSerialDenseVector destructor.
  virtual ~Epetra_LongLongSerialDenseVector ();

  //bring the base-class operator() into the current scope, in order to tell the
  //compiler that we intend to overload it, rather than hide it.
  using Epetra_LongLongSerialDenseMatrix::operator();

  //! Element access function.
  /*!
    Returns the specified element of the vector.
    \return Specified element in vector.

    \warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    long long& operator () (int Index);

  //! Element access function.
  /*!
    Returns the specified element of the vector.
    \return Specified element in vector.

    \warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    const long long& operator () (int Index) const;

  //! Element access function.
  /*!
    Returns the specified element of the vector.
    \return Specified element in vector.

    \warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    long long& operator [] (int Index);

  //! Element access function.
  /*!
    Returns the specified element of the vector.
    \return Specified element in vector.

    \warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    const long long& operator [] (int Index) const;

  //! Set vector values to random numbers.
  /*!
    LongLongSerialDenseVector uses the random number generator provided by Epetra_Util.
    The vector values will be set to random values on the interval (0, 2^31 - 1).

    \return Integer error code, set to 0 if successful.
  */
  int Random();

  //! Returns length of vector.
  int Length() const {return(M_);};

  //! Returns pointer to the values in vector.
  long long* Values() {return(A_);};

  //! Returns const pointer to the values in vector.
  const long long* Values() const {return(A_);};

  //! Returns the data access mode of the \e this vector.
  Epetra_DataAccess CV() const {return(CV_);};

  //! Copy from one vector to another.
  /*!
    The operator= allows one to copy the values from one existing LongLongSerialDenseVector to another.
    The left hand side vector will take on the data access mode of the right hand side vector.

    \return Values of the left hand side vector are modified by the values of the right hand side vector.
  */
    Epetra_LongLongSerialDenseVector& operator = (const Epetra_LongLongSerialDenseVector& Source);

    //! @name I/O methods
  //@{
  //! Print service methods; defines behavior of ostream << operator.
  virtual void Print(std::ostream& os) const;
  //@}

  //! @name Expert-only unsupported methods
  //@{

  //Bring the base-class MakeViewOf method into the current scope so that the
  //compiler knows we intend to overload it, rather than hide it.
  using Epetra_LongLongSerialDenseMatrix::MakeViewOf;

  //! Reset an existing LongLongSerialDenseVector to point to another Vector.
  /*! Allows an existing LongLongSerialDenseVector to become a View of another
    vector's data, regardless of the DataAccess mode of the Source vector.
    It is assumed that the Source vector is an independent vector, and
    no checking is done to verify this.

    This is used by Epetra_CrsGraph in the OptimizeStorage method. It is used
    so that an existing (Copy) vector can be converted to a View. This frees up
    memory that CrsGraph no longer needs.

    @param Source The LongLongSerialDenseVector this will become a view of.

    \return Integer error code, set to 0 if successful.

    \warning This method is extremely dangerous and should only be used by experts.
  */

  int MakeViewOf(const Epetra_LongLongSerialDenseVector& Source);
  //@}
};

// inlined definitions of op() and op[]
//=========================================================================
inline long long& Epetra_LongLongSerialDenseVector::operator() (int Index) {
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if(Index >= M_ || Index < 0)
    throw ReportError("Index = " + toString(Index) +
                      " Out of Range 0 - " + toString(M_-1),-1);
#endif
  return(A_[Index]);
}
//=========================================================================
inline const long long& Epetra_LongLongSerialDenseVector::operator() (int Index) const {
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if(Index >= M_ || Index < 0)
    throw ReportError("Index = " + toString(Index) +
                      " Out of Range 0 - " + toString(M_-1),-1);
#endif
   return(A_[Index]);
}
//=========================================================================
inline long long& Epetra_LongLongSerialDenseVector::operator [] (int Index) {
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if(Index >= M_ || Index < 0)
    throw ReportError("Index = " + toString(Index) +
                      " Out of Range 0 - " + toString(M_-1),-1);
#endif
   return(A_[Index]);
}
//=========================================================================
inline const long long& Epetra_LongLongSerialDenseVector::operator [] (int Index) const {
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if(Index >= M_ || Index < 0)
    throw ReportError("Index = " + toString(Index) +
                      " Out of Range 0 - " + toString(M_-1),-1);
#endif
   return(A_[Index]);
}
//=========================================================================

#endif // EPETRA_NO_64BIT_GLOBAL_INDICES

#endif /* EPETRA_LONGLONGSERIALDENSEVECTOR_H */
