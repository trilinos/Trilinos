#ifndef _PETRA_RDP_DENSEVECTOR_H_
#define _PETRA_RDP_DENSEVECTOR_H_

//! Petra_RDP_DenseVector: A class for constructing and using dense vectors.

/*! The Petra_RDP_DenseVector class enables the construction and use of real-valued, 
    double-precision dense vectors.  It is built on the BLAS and LAPACK and derives from the Petra_RDP_DenseMatrix class.

The Petra_RDP_DenseVector class is intended to provide convenient vector notation but derives all signficant 
functionality from Petra_RDP_DenseMatrix.

<b>Constructing Petra_RDP_DenseVector Objects</b>

There are three Petra_RDP_DenseVector constructors.  The first constructs a zero-length object which should be made
to appropriate length using the Size() or Resize() functions and then filled with the [] or () operators. 
The second is a constructor that accepts user
data as a 1D array, the third is a copy constructor. The second constructor has
two data access modes (specified by the Petra_DataAccess argument):
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the object.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and 
only use the View mode in a secondary optimization phase.

<b>Extracting Data from Petra_RDP_DenseVector Objects</b>

Once a Petra_RDP_DenseVector is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.


The final useful function is Flops().  Each Petra_RDP_DenseVector object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Petra_Time class, one can get accurate parallel performance
numbers.


*/
#include "Petra_Petra.h" 

#include "Petra_RDP_DenseMatrix.h"


//=========================================================================
class Petra_RDP_DenseVector : public Petra_RDP_DenseMatrix{

  public:
  
  //! Default constructor; defines a zero size object.
  /*!
    Petra_RDP_DenseVector objects defined by the default constructor should be sized with the 
    Size() or Resize functions.  
    Values should be defined by using the [] or () operators.
   */
  Petra_RDP_DenseVector(void);
  
  //! Set object values from one-dimensional array.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In
           Values - Pointer to an array of double precision numbers containing the values.
    \param In 
           Length - Length of vector.

	   See Detailed Description section for further discussion.
  */
  Petra_RDP_DenseVector(Petra_DataAccess CV, double *Values, int Length);
  
  //! Petra_RDP_DenseVector copy constructor.
  
  Petra_RDP_DenseVector(const Petra_RDP_DenseVector& Source);
  
  //! Set length of a Petra_RDP_DenseVector object; init values to zero.
  /*!
    \param In 
           Length - Length of vector object.

	   Allows user to define the dimension of a Petra_RDP_DenseVector. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   destroyed and the resized vector starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int Size(int Length) {return(Petra_RDP_DenseMatrix::Shape(Length, 1));};
  
  //! Resize a Petra_RDP_DenseVector object.
  /*!
    \param In 
           Length - Length of vector object.

	   Allows user to define the dimension of a Petra_RDP_DenseVector. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   copied into the new size.  If the new shape is smaller than the original, the first Length values
	   are copied to the new vector.

    \return Integer error code, set to 0 if successful.
  */
  int Resize(int Length) {return(Petra_RDP_DenseMatrix::Reshape(Length, 1));};

  //! Petra_RDP_DenseVector destructor.  
  virtual ~Petra_RDP_DenseVector ();

  //! Element access function.
  /*!
    Returns the specified element of the vector.  Bounds checking is enforced.
    \return Specified element in vector.
  */
    double& operator () (int Index);

  //! Element access function.
  /*!
    Returns the specified element of the vector.  Bounds checking is enforced.
    \return Specified element in vector.
  */
    const double& operator () (int Index) const;

  //! Element access function.
  /*!
    Returns the specified element of the vector.
    \return Specified element in vector.

    \warning No bounds checking is done unless Petra is compiled with PETRA_ARRAY_BOUNDS_CHECK.
  */
    double& operator [] (int Index);

  //! Column access function.
  /*!
    Returns the specified element of the vector.
    \return Specified element in vector.

    \warning No bounds checking is done unless Petra is compiled with PETRA_ARRAY_BOUNDS_CHECK.
  */
    const double& operator [] (int Index) const;
    
  //! Returns length of vector.
  int Length()  const {return(M_);};

  //! Returns pointer to the values in vector.
  double * Values()  const {return(A_);};
};

#endif /* _PETRA_RDP_DENSEVECTOR_H_ */
