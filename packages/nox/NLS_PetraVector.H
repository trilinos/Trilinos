// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

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
		  bool doCopyEntries = false);

  //! Construct an NLS_PetraVector by pointing to an external Petra_RDP_Vector.
  //! Creates a shallow copy (view) of the original vector.
  /*! The vector may be modified, but it should not be deleted when
    this is deleted.  The difference between this constructor and the
    other constuctor is the CONSTness of the Petra_RDP_Vector.  */
  NLS_PetraVector(Petra_RDP_Vector& pointTo);

  //! Deconstruct NLS_PetraVector
  /*! The memory for the actual Petra_RDP_Vector should only be
      deleted if this object created that memory itself. */
  ~NLS_PetraVector();

  //! Deep copy operator. Copies both the vector and its map.
  NLS_PetraVector& operator=(const NLS_PetraVector& copyFrom);

 
  NLS_Vector& init(double value);

  virtual NLS_Vector& abs(const NLS_Vector& y) = 0;

  virtual NLS_Vector& abs(const NLS_PetraVector& y) = 0;

  virtual NLS_Vector& copy(const NLS_Vector& y, double scale = 1.0) = 0;

  virtual NLS_Vector& copy(const NLS_PetraVector& y, double scale = 1.0) = 0;

  NLS_Vector& update(double alpha, const NLS_Vector& y, double beta);

  //! See above for definition.
  NLS_Vector& update(double alpha, const NLS_PetraVector& y, double beta);

  NLS_Vector&  scale(double alpha);

  virtual NLS_Vector* newcopy() = 0;

  double infnorm() const;

  double onenorm() const;
  
  double norm() const;

  double norm(const NLS_Vector& weights) const;

  //! Weighted 2-norm.  See above for definition.
  double norm(const NLS_PetraVector& weights) const;

  double dot(const NLS_Vector& y) const;

  //! See above for definition.
  double dot(const NLS_PetraVector& y) const;

 private:
  
  //! Pointer to petra vector owned by this object
  Petra_RDP_Vector* petraVec;		

  //! True if we should delete x when this object is deconstructed.
  bool doDeletePetraVec; 

};



#endif
