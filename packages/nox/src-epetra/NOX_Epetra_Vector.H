// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_VECTOR_H
#define NOX_EPETRA_VECTOR_H

#include "NOX_Abstract_Vector.H" // base class

// Forward declaration
class Epetra_Vector;

namespace NOX {

//! %NOX %Epetra support.
namespace Epetra {

//! Implementation of NOX::Abstract::Vector for %Epetra vectors.
class Vector : public virtual NOX::Abstract::Vector {

 public:			

  //! Construct by copying map and/or elements of an Epetra_Vector.
  /*! Option for creating a "View": Construct a NOX::Epetra::Vector 
    using the SAME memory location used in the source Epetra_Vector.
    WARNING: A View can be dangerous since multiple objects can access
    the same memory locations.  
   */
  Vector(Epetra_Vector& source, NOX::CopyType type = NOX::DeepCopy,
	 bool createView = false);

  //! Construct by copying map and/or elements of an Epetra_Vector.
  Vector(const Epetra_Vector& source, NOX::CopyType type = NOX::DeepCopy);

  //! Copy constructor.
  Vector(const NOX::Epetra::Vector& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destruct Vector.
  ~Vector();

  //@{ \name Access to underlying Petra vector.

  //! Get reference to underlying Epetra vector.
  virtual Epetra_Vector& getEpetraVector();

  //! Get const reference to underlying Epetra vector.
  virtual const Epetra_Vector& getEpetraVector() const;
 
  //@}

  //@{ \name Initialization methods.

  // derived
  virtual NOX::Abstract::Vector& init(double value);

  // derived
  virtual NOX::Abstract::Vector& random(bool useSeed = false, int seed = 1);

  //! Copies source vector into "this".
  virtual NOX::Abstract::Vector& operator=(const Epetra_Vector& y);

  // derived
  virtual NOX::Abstract::Vector& operator=(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);
  
  // derived
  virtual NOX::Abstract::Vector& abs(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);

  // derived
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Epetra::Vector& y);
  virtual NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);

  //@}

  //@{ \name Update methods.

  // derived
  virtual NOX::Abstract::Vector& scale(double gamma);

  // derived
  virtual NOX::Abstract::Vector& scale(const NOX::Epetra::Vector& a);
  virtual NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);

  // derived
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Epetra::Vector& a, 
			     double gamma = 0.0);
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			     double gamma = 0.0);

  // derived
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Epetra::Vector& a, 
			     double beta, const NOX::Epetra::Vector& b,
			     double gamma = 0.0);
  virtual NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
			     double beta, const NOX::Abstract::Vector& b,
			     double gamma = 0.0);


  //@}

  //@{ \name Creating new Vectors. 

  // derived
  virtual NOX::Abstract::Vector* clone(CopyType type = DeepCopy) const;

#ifdef HAVE_NOX_MULTIVECS

  /*! 
   * \brief Create a MultiVector with \c numVecs+1 columns out of an array of 
   * Vectors.  The vector stored under \c this will be the first column with
   * the remaining \c numVecs columns given by \c vecs.
   *
   * The implementation here creates a NOX::Epetra::MultiVector with
   * either Shape or Deep copies of the supplied vectors.
   */
  virtual NOX::Abstract::MultiVector* 
  createMultiVector(const NOX::Abstract::Vector* const* vecs,
		    int numVecs, NOX::CopyType type = NOX::DeepCopy) const;

  /*! 
   * \brief Create a MultiVector with \c numVecs columns.  
   *
   * The implementation here creates a NOX::Epetra::MultiVector with
   * either Shape or Deep copies of the supplied vector.
   */
  virtual NOX::Abstract::MultiVector* 
  createMultiVector(int numVecs, NOX::CopyType type = NOX::DeepCopy) const;
#endif

  //@}

  //@{ \name Norms.

  // derived
  virtual double norm(NOX::Abstract::Vector::NormType type = TwoNorm) const;

  // derived
  virtual double norm(const NOX::Epetra::Vector& weights) const;
  virtual double norm(const NOX::Abstract::Vector& weights) const;

  //@}

  //@{ \name Dot products

  // derived
  virtual double dot(const NOX::Epetra::Vector& y) const;
  virtual double dot(const NOX::Abstract::Vector& y) const;

  //@}

  // derived
  virtual int length() const;

  // derived
  virtual void print() const;

 protected:
  
  //! Pointer to petra vector owned by this object
  Epetra_Vector* epetraVec;

  //! True if *this created (i.e. allocated the memory for) the underlying Epetra_Vector.
  bool ownsEpetraVector;

};
} // namespace Epetra
} // namespace NOX

#endif
