// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_VECTORSPACE_H
#define NOX_EPETRA_VECTORSPACE_H

#include "NOX_Abstract_Vector.H"  // For NormType

// Forward declaration
class Epetra_Vector;

namespace NOX {

namespace Epetra {

  /*! \brief Pure virtual base class for the vector space used by NOX::Epetra::Vectors.

    This class allows users to override the inner product and norm
    used by the NOX::Epetra::Vector class.  The most frequent use of
    this class is for introducing a weighted norm throughout NOX.

  */
  class VectorSpace {

  public:

    //! Constructor
    VectorSpace() {};

    //! Destructor
    virtual ~VectorSpace() {};

    //! Computes the inner product: <a,b>.
    virtual double innerProduct(const Epetra_Vector& a,
                const Epetra_Vector& b) const = 0;

    //! Computes the norm.
    /*! For an L2 norm, the computation is: sqrt( <a,a> ). */
    virtual double norm(const Epetra_Vector& a,
            NOX::Abstract::Vector::NormType =
            NOX::Abstract::Vector::TwoNorm) const = 0;

  };
} // namespace Epetra
} // namespace NOX

#endif
