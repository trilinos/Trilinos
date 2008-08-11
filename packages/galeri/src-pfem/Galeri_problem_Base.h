// @HEADER
// ************************************************************************
//
//            Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
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
