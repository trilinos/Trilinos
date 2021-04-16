// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  constraint.hpp
    \brief Defines the SimOpt constraint for the 'poisson' example.
*/

#ifndef FEMDATA_HPP
#define FEMDATA_HPP

#include "../../../TOOLS/pde.hpp"
#include "../../../TOOLS/assembler.hpp"

// Do not instantiate the template in this translation unit.
extern template class Assembler<double>;

template<class Real>
class FEM_Data {
private:
  const ROL::Ptr<PDE<Real>> pde_;
  ROL::Ptr<Assembler<Real>> assembler_;

  ROL::Ptr<Tpetra::CrsMatrix<>>                                        matS0_;
  std::vector<std::vector<std::vector<ROL::Ptr<Tpetra::CrsMatrix<>>>>> matS_;

  ROL::Ptr<Tpetra::MultiVector<>> vecR_;

  int M_, N_, T_;

  void assemble(Teuchos::ParameterList &list) {
    ROL::Ptr<std::vector<Real>>     z = ROL::makePtr<std::vector<Real>>(M_*N_*T_,0);
    ROL::Ptr<Tpetra::MultiVector<>> u = assembler_->createStateVector();
    u->scale(static_cast<Real>(0));

    assembler_->assemblePDEResidual(vecR_,pde_,u,ROL::nullPtr,z);
    assembler_->assemblePDEJacobian1(matS0_,pde_,u,ROL::nullPtr,z);
    matS_.resize(M_);
    for (int i = 0; i < M_; ++i) {
      matS_[i].resize(N_);
      for (int j = 0; j < N_; ++j) {
        matS_[i][j].resize(T_);
        for (int k = 0; k < T_; ++k) {
          z->assign(M_*N_*T_,static_cast<Real>(0));
          (*z)[i + M_*(j + N_*k)] = static_cast<Real>(1);
          assembler_->assemblePDEJacobian1(matS_[i][j][k],pde_,u,ROL::nullPtr,z);
        }
      }
    }
  }

public:
  FEM_Data(const ROL::Ptr<PDE<Real>>                &pde,
           const ROL::Ptr<MeshManager<Real>>        &meshMgr,
           const ROL::Ptr<const Teuchos::Comm<int>> &comm,
           Teuchos::ParameterList                   &list,
           std::ostream                             &outStream = std::cout)
    : pde_(pde) {
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),meshMgr,comm,list,outStream);
    assembler_->setCellNodes(*pde_);

    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(list.sublist("Problem"), "Young's Modulus");
    M_ = list.sublist("Problem").get("Number of Horizontal Cells",10);
    N_ = list.sublist("Problem").get("Number of Vertical Cells",20);
    T_ = ym.size();

    assemble(list);
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  const ROL::Ptr<PDE<Real>> getPDE(void) const {
    return pde_;
  }

  const ROL::Ptr<const Tpetra::MultiVector<>> getForce(void) const {
    return vecR_;
  }

  const ROL::Ptr<const Tpetra::CrsMatrix<>> getS0(void) const {
    return matS0_;
  }

  const ROL::Ptr<const Tpetra::CrsMatrix<>> getS(const unsigned i, const unsigned j, const unsigned k) const {
    return matS_[i][j][k];
  }

  void applyS0(Tpetra::MultiVector<> &Su, const Tpetra::MultiVector<> &u, bool trans = false) const {
    if (trans) {
     matS0_->apply(u,Su,Teuchos::TRANS);
    }
    else {
      matS0_->apply(u,Su);
    }
  }

  void applyS(const unsigned i, const unsigned j, const unsigned k,
                 Tpetra::MultiVector<> &Su, const Tpetra::MultiVector<> &u, bool trans = false) const {
    if (trans) {
     matS_[i][j][k]->apply(u,Su,Teuchos::TRANS);
    }
    else {
      matS_[i][j][k]->apply(u,Su);
    }
  }

  /***************************************************************************/
  /* Output routines.                                                        */
  /***************************************************************************/
  void printMeshData(std::ostream &outStream) const {
    assembler_->printMeshData(outStream);
  }

  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<>> &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<>> vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }
  /***************************************************************************/
  /* End of output routines.                                                 */
  /***************************************************************************/
};

#endif
