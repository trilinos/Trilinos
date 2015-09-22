//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

//! Komplex_Vector: A class for constructing and using dense vectors on a parallel computer.

/*! The Komplex_Vector class enables the construction and use of real-valued, 
    double-precision dense vectors in a distributed memory environment.  The distribution of the dense
    vector is determined in part by a Epetra_Comm object and a Epetra_Map (or Epetra_LocalMap
    or Epetra_BlockMap).

    This class is derived from the Komplex_MultiVector class.  As such, it has full access
    to all of the functionality provided in the Komplex_MultiVector class.

<b> Distributed Global vs. Replicated Local</b>
<ul>
  <li> Distributed Global Vectors - In most instances, a multi-vector will be partitioned
       across multiple memory images associated with multiple processors.  In this case, there is 
       a unique copy of each element and elements are spread across all processors specified by 
       the Epetra_Comm communicator.
  <li> Replicated Local Vectors - Some algorithms use vectors that are too small to
       be distributed across all processors.  Replicated local vectors handle
       these types of situation.
</ul>

<b>Constructing Komplex_Vectors</b>

There are four Komplex_Vector constructors.  The first is a basic constructor that allocates
space and sets all values to zero, the second is a general constructor, the first is a copy
constructor, and the fourth works with user data.  This constructor has two data access modes:
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the vector.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and 
only use the View mode in a secondary optimization phase.

All Komplex_Vector constructors require a map argument that describes the layout of elements
on the parallel machine.  Specifically, 
\c map is a Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object describing the desired
memory layout for the vector.

There are four different Komplex_Vector constructors:
<ul>
  <li> Basic - All values are zero.
  <li> General - Create a new Komplex_Vector from two Epetra_Vectors.
  <li> Copy - Copy an existing Komplex_Vector.
  <li> Copy or make view of a vector from a Komplex_MultiVector object.
</ul>

<b>Vector and Utility Functions</b>

Once a Komplex_Vector is constructed, a variety of mathematical functions can be applied to
the vector.  Specifically:
<ul>
  <li> Vector Updates.
  <li> \e p Norms.
  </ul>

\warning A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is required for all 
  Komplex_Vector constructors.
*/

//=========================================================================
class Komplex_Vector : public Komplex_MultiVector {

  public:

  //@{ \name Constructors/destructors.
  //! Basic Komplex_Vector constuctor.
  /*! Creates a Komplex_Vector object and fills with zero values.  

    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In 
           zeroOut - If <tt>true</tt> then the allocated memory will be zeroed
                     out initially.  If <tt>false</tt> then this memory will not
                     be touched, which can be significantly faster.

	   \warning Note that, because Epetra_LocalMap
	   derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
	   for all three types of Epetra map classes.
    \return Pointer to a Komplex_Vector.
  */
  Komplex_Vector(const Epetra_BlockMap & Map, bool zeroOut = true);

  //! General Komplex_Vector constructor.
  /*!
	\param In
		 Map - A Epetra_LocalMap, Epetra_Map, or Epetra_BlockMap
	\param In
		 br - A Epetra_Vector containing the real parts of the complex vector
	\parma In
		 bi - A Epetra_Vector containing the imaginary parts of the complex vector
  */
  Komplex_Vector(const Epetra_BlockMap & Map, const Epetra_Vector & br,
		     const Epetra_Vector & bi);

  //! Komplex_Vector copy constructor.
  Komplex_Vector(const Komplex_Vector & Source);

  //! Set vector values from a vector in an existing Komplex_MultiVector.
  /*!
    \param In 
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
           Source - An existing fully constructed Komplex_MultiVector.
    \param In
           Index - Index of vector to access.  

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Komplex_Vector(Epetra_DataAccess CV, const Epetra_BlockMap & Map, 
		     const Komplex_MultiVector & Source, int Index);

  //! Komplex_Vector destructor.  
  virtual ~Komplex_Vector();
  //@}
  
  //@{ \name Post-construction modification routines.

  //! Replace values in a vector with a given indexed list of values, indices are in global index space.
  /*!
     Replace the Indices[i] entry in the \e this object with Values[i], for i=0; i<NumEntries.  The indices
     are in global index space.

    \param In
           NumEntries - Number of vector entries to modify.
    \param In
           Values - Values which will replace existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in global index space corresponding to Values.
    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int ReplaceGlobalValues(int NumEntries, double * Values, int * Indices);

  //! Replace values in a vector with a given indexed list of values, indices are in local index space.
  /*!
     Replace the Indices[i] entry in the \e this object with Values[i], for i=0; i<NumEntries.  The indices
     are in local index space.

    \param In
           NumEntries - Number of vector entries to modify.
    \param In
           Values - Values which will replace existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in local index space corresponding to Values.
    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int ReplaceMyValues(int NumEntries, double * Values, int * Indices);

  //! Sum values into a vector with a given indexed list of values, indices are in global index space.
  /*!
     Sum Values[i] into the Indices[i] entry in the \e this object, for i=0; i<NumEntries.  The indices
     are in global index space.

    \param In
           NumEntries - Number of vector entries to modify.
    \param In
           Values - Values which will be added to existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in global index space corresponding to Values.
    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int SumIntoGlobalValues(int NumEntries, double * Values, int * Indices);

  //! Sum values into a vector with a given indexed list of values, indices are in local index space.
  /*!
     Sum Values[i] into the Indices[i] entry in the \e this object, for i=0; i<NumEntries.  The indices
     are in local index space.

    \param In
           NumEntries - Number of vector entries to modify.
    \param In
           Values - Values which will be added to existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in local index space corresponding to Values.
    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int SumIntoMyValues(int NumEntries, double * Values, int * Indices);

  // Blockmap Versions

  //! Replace values in a vector with a given indexed list of values at the specified BlockOffset, indices are in global index space.
  /*!
     Replace the Indices[i] entry in the \e this object with Values[i], for i=0; i<NumEntries.  The indices
     are in global index space.  This method is intended for vector that are defined using block maps.  In this situation, 
     an index value is associated with one or more vector entries, depending on the element size of the given index.
     The BlockOffset argument indicates which vector entry to modify as an offset from the first vector entry associated with
     the given index.  The offset is used for each entry in the input list.

    \param In
           NumEntries - Number of vector entries to modify.
    \param In
           BlockOffset - Offset from the first vector entry associated with each of the given indices.
    \param In
           Values - Values which will replace existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in global index space corresponding to Values.
    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int ReplaceGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices);

  //! Replace values in a vector with a given indexed list of values at the specified BlockOffset, indices are in local index space.
  /*!
     Replace the (Indices[i], BlockOffset) entry in the \e this object with Values[i], for i=0; i<NumEntries.  The indices
     are in local index space.  This method is intended for vector that are defined using block maps.  In this situation, 
     an index value is associated with one or more vector entries, depending on the element size of the given index.
     The BlockOffset argument indicates which vector entry to modify as an offset from the first vector entry associated with
     the given index.  The offset is used for each entry in the input list.

    \param In
           NumEntries - Number of vector entries to modify.
    \param In
           BlockOffset - Offset from the first vector entry associated with each of the given indices.
    \param In
           Values - Values which will replace existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in local index space corresponding to Values.
    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int ReplaceMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices);

  //! Sum values into a vector with a given indexed list of values at the specified BlockOffset, indices are in global index space.
  /*!
     Sum Values[i] into the Indices[i] entry in the \e this object, for i=0; i<NumEntries.  The indices
     are in global index space.  This method is intended for vectors that are defined using block maps.  In this situation, 
     an index value is associated with one or more vector entries, depending on the element size of the given index.
     The BlockOffset argument indicates which vector entry to modify as an offset from the first vector entry associated with
     the given index.  The offset is used for each entry in the input list.

    \param In
           NumEntries - Number of vector entries to modify.
    \param In
           BlockOffset - Offset from the first vector entry associated with each of the given indices.
    \param In
           Values - Values which will replace existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in global index space corresponding to Values.
    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int SumIntoGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices);

  //! Sum values into a vector with a given indexed list of values at the specified BlockOffset, indices are in local index space.
  /*!
     Sum Values[i] into the Indices[i] entry in the \e this object, for i=0; i<NumEntries.  The indices
     are in local index space.  This method is intended for vector that are defined using block maps.  In this situation, 
     an index value is associated with one or more vector entries, depending on the element size of the given index.
     The BlockOffset argument indicates which vector entry to modify as an offset from the first vector entry associated with
     the given index.  The offset is used for each entry in the input list.

    \param In
           NumEntries - Number of vector entries to modify.
    \param In
           BlockOffset - Offset from the first vector entry associated with each of the given indices.
    \param In
           Values - Values which will be added to existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in local index space corresponding to Values.
    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int SumIntoMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices);
  //@}

  //@{ \name Mathematical methods.

    //! Scale the current values of the \e this vector, \e this = ScalarValue*\e this.
  /*!
    \param In
	     ScalarValue - Scale value.
    \param Out
	     \e This - Vector with scaled values.
    \return Integer error code, set to 0 if successful.
  */
  int Scale(double ScalarValue);

  //! Replace vector values with scaled values of A, \e this = ScalarA*A.
  /*!
    \param In
	     ScalarA - Scale value.
    \param In
	     A - Vector to copy.
    \param Out
	     \e This - Vector with values overwritten by scaled values of A.
    \return Integer error code, set to 0 if successful.
  */
  int Scale(double ScalarA, const Komplex_Vector & A);

  //! Compute 1-norm of the \e this vector.
  /*!
    \param Out
	     Result - Result contains 1-norm of the \e this vector.
    \return Integer error code, set to 0 if successful.
  */
  int Norm1(double & Result) const;

  //! Compute 2-norm of the \e this vector.
  /*!
    \param Out
	     Result - Result contains 2-norm of the \e this vector.
    \return Integer error code, set to 0 if successful.
  */
  int Norm2(double & Result) const;

  //! Compute Inf-norm of the \e this vector.
  /*!
    \param Out
	     Result - Result contains Inf-norm of the \e this vector.
    \return Integer error code, set to 0 if successful.
  */
  int NormInf(double & Result) const;
  //@}

  //@{ \name Overloaded operators

  //! = Operator.
  /*!
    \param In
  	     A - Komplex_Vector to copy.
    \return Komplex_Vector.
  */
  Komplex_Vector & operator = (const Komplex_Vector & Source);

  //! Element access function.
  /*!
    \return V[Index].
  */
    double & operator [] (int index);

  //! Element access function.
  /*!
    \return V[Index].
  */
  const double & operator [] (int index) const;
  //@}

  //@{ \name Attribute access functions

  //! Returns the length of the vector.
  int Length() const; 

  /*! Replace map, only if new map has same point-structure as current map.
      return 0 if map is replaced, -1 if not.
   */
  int ReplaceMap(const Epetra_BlockMap & map);
  //@}  

  //@{ \name I/O methods

  //! Print method
  void Print(ostream & os) const;
  //@}

  protected:

  private:
  int GlobalLength_;
  int MyLength_;
  int NumElements_;
  int IndexBase_;
 };