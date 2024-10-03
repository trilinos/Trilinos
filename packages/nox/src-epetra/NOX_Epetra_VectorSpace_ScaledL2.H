// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_VECTORSPACE_SCALEDL2_H
#define NOX_EPETRA_VECTORSPACE_SCALEDL2_H

#include "NOX_Epetra_VectorSpace.H" // base class
#include "Teuchos_RCP.hpp"
#include "NOX_Epetra_Scaling.H"

// Forward declarations
class Epetra_Vector;

namespace NOX {

namespace Epetra {

  /*! \brief Concrete class for a weighted L2 vector space.

  This class allows users to override the inner product and norm
  used by the NOX::Epetra::Vector class.  The most frequent use of
  this class is for introducing a weighted norm throughout NOX.

  */
  class VectorSpaceScaledL2 : public NOX::Epetra::VectorSpace {

  public:

    //! Constructor
    VectorSpaceScaledL2(const Teuchos::RCP<NOX::Epetra::Scaling>& s,
            NOX::Epetra::Scaling::ScaleType st =
            NOX::Epetra::Scaling::Left);

    //! Destructor
    virtual ~VectorSpaceScaledL2();

    //! Computes a scaled inner product.
    /*! Computes a scaled inner product: \f$ <Da, Db> \f$ where
      \f$D\f$ is the set of scaling vectors associated with either
      left of right scaling.
    */
    virtual double innerProduct(const Epetra_Vector& a,
                const Epetra_Vector& b) const;

    //! Computes the scaled norm.
    /*! Computes the scaled norm using \f$ Da \f$ where \f$D\f$ is the
      set of scaling vectors associated with either left of right
      scaling.
    */
    virtual double norm(const Epetra_Vector& a,
            NOX::Abstract::Vector::NormType =
            NOX::Abstract::Vector::TwoNorm) const;

  protected:

    //! Scaling vector used in the inner product
    Teuchos::RCP<NOX::Epetra::Scaling> scalingPtr;

    //! Scaling type to apply to vector space
    NOX::Epetra::Scaling::ScaleType scaleType;

    //! Temporary vector used in scaling computations
    mutable Teuchos::RCP<Epetra_Vector> tmpVectorPtr;

  };
} // namespace Epetra
} // namespace NOX

#endif
