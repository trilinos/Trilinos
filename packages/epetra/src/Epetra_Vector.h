
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_VECTOR_H_
#define _EPETRA_VECTOR_H_

#include "Epetra_MultiVector.h"
class Epetra_Map;

//! Epetra_Vector: A class for constructing and using dense vectors on a parallel computer.

/*! The Epetra_Vector class enables the construction and use of real-valued, 
    double-precision dense vectors in a distributed memory environment.  The distribution of the dense
    vector is determined in part by a Epetra_Comm object and a Epetra_Map (or Epetra_LocalMap
    or Epetra_BlockMap).

    This class is derived from the Epetra_MultiVector class.  As such, it has full access
    to all of the functionality provided in the Epetra_MultiVector class.

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

<b>Constructing Epetra_Vectors</b>

There are four Epetra_Vector constructors.  The first is a basic constructor that allocates
space and sets all values to zero, the second is a 
copy constructor. The third and fourth constructors work with user data.  These constructors have
two data access modes:
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the vector.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and 
only use the View mode in a secondary optimization phase.

All Epetra_Vector constructors require a map argument that describes the layout of elements
on the parallel machine.  Specifically, 
\c map is a Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object describing the desired
memory layout for the vector.

There are four different Epetra_Vector constructors:
<ul>
  <li> Basic - All values are zero.
  <li> Copy - Copy an existing vector.
  <li> Copy from or make view of user double array.
  <li> Copy or make view of a vector from a Epetra_MultiVector object.
</ul>

<b>Extracting Data from Epetra_Vectors</b>

Once a Epetra_Vector is constructed, it is possible to extract a copy of the values or create
a view of them.

\warning ExtractView functions are \e extremely dangerous from a data hiding perspective.
For both ExtractView fuctions, there is a corresponding ExtractCopy function.  We
strongly encourage users to develop code using ExtractCopy functions first and 
only use the ExtractView functions in a secondary optimization phase.

There are two Extract functions:
<ul>
  <li> ExtractCopy - Copy values into a user-provided array.
  <li> ExtractView - Set user-provided array to point to Epetra_Vector data.
</ul>

<b>Vector and Utility Functions</b>

Once a Epetra_Vector is constructed, a variety of mathematical functions can be applied to
the vector.  Specifically:
<ul>
  <li> Dot Products.
  <li> Vector Updates.
  <li> \e p Norms.
  <li> Weighted Norms.
  <li> Minimum, Maximum and Average Values.
</ul>

The final useful function is Flops().  Each Epetra_Vector object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Epetra_Time class, one can get accurate parallel performance
numbers.

\warning A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is required for all 
  Epetra_Vector constructors.

*/

//=========================================================================
class Epetra_Vector : public Epetra_MultiVector {

  public:

  //@{ \name Constructors/destructors.
  //! Basic Epetra_Vector constuctor.
  /*! Creates a Epetra_Vector object and fills with zero values.  

    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

	   \warning Note that, because Epetra_LocalMap
	   derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
	   for all three types of Epetra map classes.

    \return Pointer to a Epetra_Vector.

  */
  Epetra_Vector(const Epetra_BlockMap& Map);

  //! Epetra_Vector copy constructor.
  
  Epetra_Vector(const Epetra_Vector& Source);
  
  //! Set vector values from user array.
  /*!
    \param In 
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
           V - Pointer to an array of double precision numbers..

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Epetra_Vector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *V);

  //! Set vector values from a vector in an existing Epetra_MultiVector.
  /*!
    \param In 
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
           Source - An existing fully constructed Epetra_MultiVector.
    \param In
           Index - Index of vector to access.  

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Epetra_Vector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int Index);

  //! Epetra_Vector destructor.  
    virtual ~Epetra_Vector ();
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
           Values - Values which will replace existing values in vector, of length NumEntries.
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
           Values - Values which will replace existing values in vector, of length NumEntries.
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
           Values - Values which will replace existing values in vector, of length NumEntries.
    \param In
           Indices - Indices in local index space corresponding to Values.

    \return Integer error code, set to 0 if successful, set to 1 if one or more indices are not associated with calling processor.
  */
  int SumIntoMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices);
  //@}

  //@{ \name Extraction methods


  //! Put vector values into user-provided array.
  /*!
    \param Out
           V - Pointer to memory space that will contain the vector values.  

    \return Integer error code, set to 0 if successful.
  */
  int ExtractCopy(double *V);
  
  //! Set user-provided address of V.
  /*!
    \param Out
           V - Address of a pointer to that will be set to point to the values of the vector.  

    \return Integer error code, set to 0 if successful.
  */
  int ExtractView(double **V);
  //@}

  //@{ \name Overloaded operators

  //! Element access function.
  /*!
    \return V[Index].
  */
    double& operator [] (int index) { return Values_[index]; }
  //! Element access function.
  /*!
    \return V[Index].
  */
    const double& operator [] (int index) const { return Values_[index]; }
    //@}
    
 private:

    int ChangeValues(int NumEntries, int BlockOffset, double * Values, int * Indices, bool IndicesGlobal, bool SumInto);

};

#endif /* _EPETRA_VECTOR_H_ */
