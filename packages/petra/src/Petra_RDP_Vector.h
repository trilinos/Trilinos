#ifndef _PETRA_RDP_VECTOR_H_
#define _PETRA_RDP_VECTOR_H_

//! Petra_RDP_Vector: A class for constructing and using dense vectors on a parallel computer.

/*! The Petra_RDP_Vector class enables the construction and use of real-valued, 
    double-precision dense vectors in a distributed memory environment.  The distribution of the dense
    vector is determined in part by a Petra_Comm object and a Petra_Map (or Petra_LocalMap
    or Petra_BlockMap).

    This class is derived from the Petra_RDP_MultiVector class.  As such, it has full access
    to all of the functionality provided in the Petra_RDP_MultiVector class.

<b> Distributed Global vs. Replicated Local</b>
<ul>
  <li> Distributed Global Vectors - In most instances, a multi-vector will be partitioned
       across multiple memory images associated with multiple processors.  In this case, there is 
       a unique copy of each element and elements are spread across all processors specified by 
       the Petra_Comm communicator.
  <li> Replicated Local Vectors - Some algorithms use vectors that are too small to
       be distributed across all processors.  Replicated local vectors handle
       these types of situation.
</ul>

<b>Constructing Petra_RDP_Vectors</b>

There are four Petra_RDP_Vector constructors.  The first is a basic constructor that allocates
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

All Petra_RDP_Vector constructors require a map argument that describes the layout of elements
on the parallel machine.  Specifically, 
\c map is a Petra_Map, Petra_LocalMap or Petra_BlockMap object describing the desired
memory layout for the vector.

There are four different Petra_RDP_Vector constructors:
<ul>
  <li> Basic - All values are zero.
  <li> Copy - Copy an existing vector.
  <li> Copy from or make view of user double array.
  <li> Copy or make view of a vector from a Petra_MultiVector object.
</ul>

<b>Extracting Data from Petra_RDP_Vectors</b>

Once a Petra_RDP_Vector is constructed, it is possible to extract a copy of the values or create
a view of them.

\warning ExtractView functions are \e extremely dangerous from a data hiding perspective.
For both ExtractView fuctions, there is a corresponding ExtractCopy function.  We
strongly encourage users to develop code using ExtractCopy functions first and 
only use the ExtractView functions in a secondary optimization phase.

There are two Extract functions:
<ul>
  <li> ExtractCopy - Copy values into a user-provided array.
  <li> ExtractView - Set user-provided array to point to Petra_RDP_Vector data.
</ul>

<b>Vector and Utility Functions</b>

Once a Petra_RDP_Vector is constructed, a variety of mathematical functions can be applied to
the vector.  Specifically:
<ul>
  <li> Dot Products.
  <li> Vector Updates.
  <li> \e p Norms.
  <li> Weighted Norms.
  <li> Minimum, Maximum and Average Values.
</ul>

The final useful function is Flops().  Each Petra_RDP_Vector object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Petra_Time class, one can get accurate parallel performance
numbers.

\warning A Petra_Map, Petra_LocalMap or Petra_BlockMap object is required for all 
  Petra_RDP_Vector constructors.

*/
#include "Petra_Petra.h" 
#include "Petra_Map.h"
#include "Petra_LocalMap.h"
#include "Petra_RDP_MultiVector.h"


//=========================================================================
class Petra_RDP_Vector : public Petra_RDP_MultiVector {

  // Give ostream << function some access to private and protected data/functions.

  friend ostream& operator << (ostream& os, const Petra_RDP_Vector& A);
  public:

  //! Basic Petra_RDP_Vector constuctor.
  /*! Creates a Petra_RDP_Vector object and fills with zero values.  

    \param In 
           Map - A Petra_LocalMap, Petra_Map or Petra_BlockMap.

	   \warning Note that, because Petra_LocalMap
	   derives from Petra_Map and Petra_Map derives from Petra_BlockMap, this constructor works
	   for all three types of Petra map classes.

    \return Pointer to a Petra_RDP_Vector.

  */
  Petra_RDP_Vector(const Petra_BlockMap& Map);

  //! Petra_RDP_Vector copy constructor.
  
  Petra_RDP_Vector(const Petra_RDP_Vector& Source);
  
  //! Set vector values from user array.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Petra_LocalMap, Petra_Map or Petra_BlockMap.
    \param In
           V - Pointer to an array of double precision numbers..

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Petra_RDP_Vector(Petra_DataAccess CV, const Petra_BlockMap& Map, double *V);

  //! Set vector values from a vector in an existing Petra_RDP_MultiVector.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Petra_LocalMap, Petra_Map or Petra_BlockMap.
    \param In
           Source - An existing fully constructed Petra_RDP_Vector.
    \param In
           Index - Index of vector to access.  

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Petra_RDP_Vector(Petra_DataAccess CV, const Petra_RDP_Vector& Source, int Index);

  //! Petra_RDP_Vector destructor.  
    virtual ~Petra_RDP_Vector ();

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

  //! Element access function.
  /*!
    \return V[Index].
  */
    double& operator [] (int index);
  //! Element access function.
  /*!
    \return V[Index].
  */
    const double& operator [] (int index) const;
    
};

//! << operator will work for Petra_RDP_Vectors.
ostream& operator << (ostream& os, const Petra_RDP_Vector& A);

#endif /* _PETRA_RDP_VECTOR_H_ */
