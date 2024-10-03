// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_VECTORSPACE_L2_H
#define NOX_EPETRA_VECTORSPACE_L2_H

#include "NOX_Epetra_VectorSpace.H" // base class

namespace NOX {

namespace Epetra {

  /*! \brief Concrete class for an L2 vector space.

    This class allows users to override the inner product and norm
    used by the NOX::Epetra::Vector class.  The most frequent use of
    this class is for introducing a weighted norm throughout NOX.

  */
  class VectorSpaceL2 : public NOX::Epetra::VectorSpace {

  public:

  //! Constructor
  VectorSpaceL2();

  //! Destructor
  virtual ~VectorSpaceL2();

  virtual double innerProduct(const Epetra_Vector& a,
                  const Epetra_Vector& b) const;

  virtual double norm(const Epetra_Vector& a,
              NOX::Abstract::Vector::NormType =
              NOX::Abstract::Vector::TwoNorm) const;

};
} // namespace Epetra
} // namespace NOX

#endif
