// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
