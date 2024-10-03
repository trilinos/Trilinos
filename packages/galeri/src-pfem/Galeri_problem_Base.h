// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_PROBLEM_BASE_H
#define GALERI_PROBLEM_BASE_H

class Epetra_RowMatrix;
class Epetra_FEVector;
class Epetra_MultiVector;
class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseVector;

namespace Galeri {
  namespace quadrature {
    class Element;
  }
  namespace grid {
    class Loadable;
  }
}

namespace Galeri {
namespace problem {
  
class Base : public Galeri::core::Object 
{
  public:
    virtual ~Base() {}

    // FIXME: Only MultiVector and not FEVector??
    virtual void integrate(Galeri::grid::Loadable& domain,
                         Epetra_RowMatrix& A,
                         Epetra_FEVector& RHS) = 0;

    virtual void computeNorms(Galeri::grid::Loadable& domain,
                              const Epetra_MultiVector& solution) = 0;

    virtual void computeNorms(Galeri::grid::Loadable& domain,
                              const Epetra_MultiVector& solution,
                              double& solNormL2, double& solSemiNormH1,
                              double& exaNormL2, double& exaSemiNormH1,
                              double& errNormL2, double& errSemiNormH1) = 0;

    virtual void integrateOverElement(Galeri::quadrature::Element& QE,
                                      Epetra_SerialDenseMatrix& ElementLHS, 
                                      Epetra_SerialDenseMatrix& ElementRHS) = 0;

    virtual void computeNormOverElement(Galeri::quadrature::Element& QE,
                                        Epetra_SerialDenseMatrix& elementSol,
                                        Epetra_SerialDenseMatrix& elementNorm) = 0;

    virtual void computeErrorOverElement(Galeri::quadrature::Element& QE,
                                         Epetra_SerialDenseMatrix& elementSol,
                                         Epetra_SerialDenseMatrix& elementNorm) = 0;

    virtual void computeNormOverElement(Galeri::quadrature::Element& QE,
                                        Epetra_SerialDenseMatrix& elementNorm) = 0;
};

} // namespace problem
} // namespace Galeri
#endif
