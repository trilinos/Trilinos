
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

#ifndef EPETRA_INTSERIALDENSEVECTOR_H
#define EPETRA_INTSERIALDENSEVECTOR_H

#include "Epetra_Object.h" 
#include "Epetra_IntSerialDenseMatrix.h"

//! Epetra_IntSerialDenseVector: A class for constructing and using dense vectors.

/*! The Epetra_IntSerialDenseVector class enables the construction and use of integer-valued, 
    dense vectors.  It derives from the Epetra_IntSerialDenseMatrix class.

The Epetra_IntSerialDenseVector class is intended to provide convenient vector notation but derives all signficant 
functionality from Epetra_IntSerialDenseMatrix.

<b>Constructing Epetra_IntSerialDenseVector Objects</b>

There are three Epetra_IntSerialDenseVector constructors.  The first constructs a zero-length object which should be made
to appropriate length using the Size() or Resize() functions and then filled with the [] or () operators. 
The second is a constructor that accepts user
data as a 1D array, the third is a copy constructor. The second constructor has
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

<b>Extracting Data from Epetra_IntSerialDenseVector Objects</b>

Once a Epetra_IntSerialDenseVector is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.

*/


//=========================================================================
class Epetra_IntSerialDenseVector : public Epetra_IntSerialDenseMatrix{

  public:
  
  //! Default constructor; defines a zero size object.
  /*!
    Epetra_IntSerialDenseVector objects defined by the default constructor should be sized with the 
    Size() or Resize functions.  
    Values should be defined by using the [] or () operators.
   */
  Epetra_IntSerialDenseVector(void);
  
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
  Epetra_IntSerialDenseVector(Epetra_DataAccess CV, int *Values, int Length);
  
  //! Epetra_IntSerialDenseVector copy constructor.
  
  Epetra_IntSerialDenseVector(const Epetra_IntSerialDenseVector& Source);
  
  //! Set length of a Epetra_IntSerialDenseVector object; init values to zero.
  /*!
    \param In 
           Length - Length of vector object.

	   Allows user to define the dimension of a Epetra_IntSerialDenseVector. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   destroyed and the resized vector starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int Size(int Length) {return(Epetra_IntSerialDenseMatrix::Shape(Length, 1));};
  
  //! Resize a Epetra_IntSerialDenseVector object.
  /*!
    \param In 
           Length - Length of vector object.

	   Allows user to define the dimension of a Epetra_IntSerialDenseVector. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   copied into the new size.  If the new shape is smaller than the original, the first Length values
	   are copied to the new vector.

    \return Integer error code, set to 0 if successful.
  */
  int Resize(int Length) {return(Epetra_IntSerialDenseMatrix::Reshape(Length, 1));};

  //! Epetra_IntSerialDenseVector destructor.  
  virtual ~Epetra_IntSerialDenseVector ();

  //! Element access function.
  /*!
    Returns the specified element of the vector.  Bounds checking is enforced.
    \return Specified element in vector.
  */
    int& operator () (int Index);

  //! Element access function.
  /*!
    Returns the specified element of the vector.  Bounds checking is enforced.
    \return Specified element in vector.
  */
    const int& operator () (int Index) const;

  //! Element access function.
  /*!
    Returns the specified element of the vector.
    \return Specified element in vector.

    \warning No bounds checking is done unless Epetra is compiled with EPETRA_ARRAY_BOUNDS_CHECK.
  */
    int& operator [] (int Index);

  //! Column access function.
  /*!
    Returns the specified element of the vector.
    \return Specified element in vector.

    \warning No bounds checking is done unless Epetra is compiled with EPETRA_ARRAY_BOUNDS_CHECK.
  */
    const int& operator [] (int Index) const;
    
  //! Returns length of vector.
  int Length()  const {return(M_);};

  //! Returns pointer to the values in vector.
  int * Values()  const {return(A_);};
};

#endif /* EPETRA_INTSERIALDENSEVECTOR_H */
