#ifndef _NLS_PETRAVECTOR_H
#define _NLS_PETRAVECTOR_H


#include "NLS_Vector.H"	// Base Class
#include "Petra_RDP_Vector.h"

//! Implementation of NLS_Vector pure virtual class for Petra_RDP_Vectors.

class NLS_PetraVector : public NLS_Vector {

 public:			

  //! Construct an NLS_PetraVector by copying map and/or elements of
  //! a Petra_RDP_Vector.  Creates a deep copy of the vector.
  /*! The input vector is used to specify the map.  The contents are
    only copied if doCopyEntries is true.  If false, the vector is 
    initialized to all zeros.  */
  NLS_PetraVector(const Petra_RDP_Vector& copyFrom, 
<<<<<<< NLS_PetraVector.H
		  bool doCopyEntries = false);
=======
		  bool doCopyEntries);

  //! Construct an NLS_PetraVector by pointing to an external Petra_RDP_Vector.
  //! Creates a shallow copy (view) of the original vector.
  /*! The vector may be modified, but it should not be deleted when
    this is deleted.  The difference between this constructor and the
    other constuctor is the CONSTness of the Petra_RDP_Vector.  */
  NLS_PetraVector(Petra_RDP_Vector& pointTo);
>>>>>>> 1.4

  //! Deconstruct NLS_PetraVector
  /*! The memory for the actual Petra_RDP_Vector should only be
      deleted if this object created that memory itself. */
  ~NLS_PetraVector();

  //! Deep copy operator. Copies both the vector and its map.
  NLS_PetraVector& operator=(const NLS_PetraVector& copyFrom);
<<<<<<< NLS_PetraVector.H
=======

>>>>>>> 1.3

<<<<<<< NLS_PetraVector.H
  // Derived
  NLS_PetraVector& scale(double alpha);

  // Derived
  double norm() const;
=======
  void random();
  
  void putScalar (double scalar);

  void update(double scalarA, const NLS_Vector& y, double scalar);
>>>>>>> 1.3

  //! See above for definition.
  void update(double scalarA, const NLS_PetraVector& y, double scalar);

  void scale(double alpha);

  double dot(const NLS_Vector& y) const;

  //! See above for definition.
  double dot(const NLS_PetraVector& y) const;

  double normInf() const;

  double norm1() const;
  
  double norm2() const;

  double norm2(const NLS_Vector& weights) const;

  //! Weighted 2-norm.  See above for definition.
  double norm2(const NLS_PetraVector& weights) const;

  void Abs(const NLS_Vector& y);

  //! See above for definition.
  void Abs(const NLS_PetraVector& y);

  //! Returns the minimum value of vector.
  double minValue() const;

  //! Returns the maximum value of vector.
  double maxValue() const;

  //! Returns the average value of vector.
  double meanValue() const;

 private:
  
  //! Pointer to petra vector owned by this object
  Petra_RDP_Vector* petraVec;		

  //! True if we should delete x when this object is deconstructed.
  bool doDeletePetraVec; 

};



#endif
