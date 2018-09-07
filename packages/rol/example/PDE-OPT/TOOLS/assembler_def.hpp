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

/*! \file  assembler_def.hpp
    \brief Finite element assembly class.
*/

#ifndef ROL_PDEOPT_ASSEMBLER_DEF_H
#define ROL_PDEOPT_ASSEMBLER_DEF_H

#include "assembler.hpp"

/*****************************************************************************/
/******************* PUBLIC MEMBER FUNCTIONS *********************************/
/*****************************************************************************/
// Constuctor: Discretization set from ParameterList
template<class Real>
Assembler<Real>::Assembler(
        const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs,
        const ROL::Ptr<const Teuchos::Comm<int>> &comm,
        Teuchos::ParameterList &parlist,
        std::ostream &outStream)
  : isJ1Transposed_(false), isJ2Transposed_(false),
    isJuoTransposed_(false), isJunTransposed_(false), isJzfTransposed_(false) {
  setCommunicator(comm,parlist,outStream);
  setBasis(basisPtrs,parlist,outStream);
  setDiscretization(parlist,ROL::nullPtr,outStream);
  setParallelStructure(parlist,outStream);
  setCellNodes(outStream);
}

// Constructor: Discretization set from MeshManager input
template<class Real>
Assembler<Real>::Assembler(
        const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs,
        const ROL::Ptr<MeshManager<Real>> &meshMgr,
        const ROL::Ptr<const Teuchos::Comm<int>> &comm,
        Teuchos::ParameterList &parlist,
        std::ostream &outStream)
  : isJ1Transposed_(false), isJ2Transposed_(false),
    isJuoTransposed_(false), isJunTransposed_(false), isJzfTransposed_(false) {
  setCommunicator(comm,parlist,outStream);
  setBasis(basisPtrs,parlist,outStream);
  setDiscretization(parlist,meshMgr,outStream);
  setParallelStructure(parlist,outStream);
  setCellNodes(outStream);
}

template<class Real>
void Assembler<Real>::setCellNodes(PDE<Real> &pde) const {
  // Set PDE cell nodes
  pde.setFieldPattern(dofMgr_->getFieldPattern());
  pde.setCellNodes(volCellNodes_, bdryCellNodes_, bdryCellLocIds_);
}

template<class Real>
void Assembler<Real>::setCellNodes(DynamicPDE<Real> &pde) const {
  // Set PDE cell nodes
  pde.setFieldPattern(dofMgr_->getFieldPattern());
  pde.setCellNodes(volCellNodes_, bdryCellNodes_, bdryCellLocIds_);
}

/***************************************************************************/
/* PDE assembly routines                                                   */
/***************************************************************************/
template<class Real>
void Assembler<Real>::assemblePDEResidual(ROL::Ptr<Tpetra::MultiVector<>> &r,
                                    const ROL::Ptr<PDE<Real>> &pde,
                                    const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                    const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                    const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEResidual);
  #endif
  // Initialize residual vectors if not done so
  if ( r == ROL::nullPtr ) { // Unique components of residual vector
    r = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueResidualMap_, 1, true);
  }
  if ( pde_vecR_overlap_ == ROL::nullPtr ) { // Overlapping components of residual vector
    pde_vecR_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapResidualMap_, 1, true);
  }
  // Get u_coeff from u and z_coeff from z
  ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
  getCoeffFromStateVector(u_coeff,u);
  getCoeffFromControlVector(z_coeff,z);
  // Assemble
  ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
  pde->residual(localVal,u_coeff,z_coeff,z_param);
  assembleFieldVector( r, localVal, pde_vecR_overlap_ );
}

template<class Real>
void Assembler<Real>::assemblePDEJacobian1(ROL::Ptr<Tpetra::CrsMatrix<>> &J1,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian1);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize Jacobian matrices
    if ( J1 == ROL::nullPtr ) {
      J1 = ROL::makePtr<Tpetra::CrsMatrix<>>(matJ1Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Jacobian_1(localVal,u_coeff,z_coeff,z_param);
    assembleFieldMatrix( J1, localVal );
    isJ1Transposed_ = false;
  }
  catch ( Exception::Zero & ez ) {
    throw Exception::Zero(">>> (Assembler::assemblePDEJacobian1): Jacobian is zero.");
  }
  catch ( Exception::NotImplemented & eni ) {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian1): Jacobian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEJacobian2(ROL::Ptr<Tpetra::CrsMatrix<>> &J2,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian2);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize Jacobian matrices
    if ( J2 == ROL::nullPtr ) {
      J2 = ROL::makePtr<Tpetra::CrsMatrix<>>(matJ2Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Jacobian_2(localVal,u_coeff,z_coeff,z_param);
    assembleFieldMatrix( J2, localVal );
    isJ2Transposed_ = false;
  }
  catch ( Exception::Zero & ez ) {
    throw Exception::Zero(">>> (Assembler::assemblePDEJacobian2): Jacobian is zero.");
  }
  catch ( Exception::NotImplemented & eni ) {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian2): Jacobian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEJacobian3(ROL::Ptr<Tpetra::MultiVector<>> &J3,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian3);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Initialize Jacobian storage if not done so already
      if (J3 == ROL::nullPtr) {
        J3 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueResidualMap_, size, true);
      }
      if ( pde_vecJ3_overlap_ == ROL::nullPtr) {
        pde_vecJ3_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapResidualMap_, size, true);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      pde->Jacobian_3(localVal,u_coeff,z_coeff,z_param);
      assembleParamFieldMatrix( J3, localVal, pde_vecJ3_overlap_ );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assemblePDEJacobian3): Jacobian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian3): Jacobian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian3): Jacobian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEHessian11(ROL::Ptr<Tpetra::CrsMatrix<>> &H11,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian11);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H11 == ROL::nullPtr ) {
      H11 = ROL::makePtr<Tpetra::CrsMatrix<>>(matH11Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_11(localVal,l_coeff,u_coeff,z_coeff,z_param);
    assembleFieldMatrix( H11, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assemblePDEHessian11): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian11): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEHessian12(ROL::Ptr<Tpetra::CrsMatrix<>> &H12,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian12);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H12 == ROL::nullPtr ) {
      H12 = ROL::makePtr<Tpetra::CrsMatrix<>>(matH12Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_12(localVal,l_coeff,u_coeff,z_coeff,z_param);
    assembleFieldMatrix( H12, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assemblePDEHessian12): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian12): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEHessian13(ROL::Ptr<Tpetra::MultiVector<>> &H13,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian13);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u, z_coeff from z and l_coeff from l
      ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff, l_coeff;
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Initialize Jacobian storage if not done so already
      if (H13 == ROL::nullPtr) {
        H13 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, size, true);
      }
      if ( pde_vecH13_overlap_ == ROL::nullPtr) {
        pde_vecH13_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, size, true);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      pde->Hessian_13(localVal,l_coeff,u_coeff,z_coeff,z_param);
      assembleParamFieldMatrix( H13, localVal, pde_vecH13_overlap_ );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian13): Hessian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian13): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian13): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEHessian21(ROL::Ptr<Tpetra::CrsMatrix<>> &H21,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian21);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H21 == ROL::nullPtr ) {
      H21 = ROL::makePtr<Tpetra::CrsMatrix<>>(matH21Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_21(localVal,l_coeff,u_coeff,z_coeff,z_param);
    assembleFieldMatrix( H21, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assemblePDEHessian21): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian21): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEHessian22(ROL::Ptr<Tpetra::CrsMatrix<>> &H22,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian22);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H22 == ROL::nullPtr ) {
      H22 = ROL::makePtr<Tpetra::CrsMatrix<>>(matH22Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_22(localVal,l_coeff,u_coeff,z_coeff,z_param);
    assembleFieldMatrix( H22, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assemblePDEHessian22): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian22): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEHessian23(ROL::Ptr<Tpetra::MultiVector<>> &H23,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian23);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u, z_coeff from z and l_coeff from l
      ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff, l_coeff;
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Initialize Jacobian storage if not done so already
      if (H23 == ROL::nullPtr) {
        H23 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, size, true);
      }
      if ( pde_vecH23_overlap_ == ROL::nullPtr) {
        pde_vecH23_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, size, true);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      pde->Hessian_23(localVal,l_coeff,u_coeff,z_coeff,z_param);
      assembleParamFieldMatrix( H23, localVal, pde_vecH23_overlap_ );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian23): Hessian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian23): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian23): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assemblePDEHessian31(ROL::Ptr<Tpetra::MultiVector<>> &H31,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian31);
  #endif
  assemblePDEHessian13(H31,pde,l,u,z,z_param);
}

template<class Real>
void Assembler<Real>::assemblePDEHessian32(ROL::Ptr<Tpetra::MultiVector<>> &H32,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian32);
  #endif
  assemblePDEHessian23(H32,pde,l,u,z,z_param);
}

template<class Real>
void Assembler<Real>::assemblePDEHessian33(ROL::Ptr<std::vector<std::vector<Real>>> &H33,
                                     const ROL::Ptr<PDE<Real>> &pde,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian33);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u, z_coeff from z and l_coeff from l
      ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff, l_coeff;
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Initialize Jacobian storage if not done so already
      if (H33 == ROL::nullPtr) {
        std::vector<Real> col(size,static_cast<Real>(0));
        H33 = ROL::makePtr<std::vector<std::vector<Real>>>(size,col);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> tmp(size,ROL::nullPtr);
      std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> localVal(size,tmp);
      pde->Hessian_33(localVal,l_coeff,u_coeff,z_coeff,z_param);
      assembleParamMatrix( H33, localVal );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian33): Hessian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian33): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian33): Hessian not implemented.");
  }
}
/***************************************************************************/
/* End of PDE assembly routines.                                           */
/***************************************************************************/

/***************************************************************************/
/* DynamicPDE assembly routines                                            */
/***************************************************************************/
template<class Real>
void Assembler<Real>::assembleDynPDEResidual(ROL::Ptr<Tpetra::MultiVector<>> &r,
                                       const ROL::Ptr<DynamicPDE<Real>> &pde,
                                       const ROL::TimeStamp<Real> &ts,
                                       const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                       const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                       const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                       const ROL::Ptr<const std::vector<Real>> &z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEResidual);
  #endif
  // Initialize residual vectors if not done so
  if ( r == ROL::nullPtr ) { // Unique components of residual vector
    r = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueResidualMap_, 1, true);
  }
  if ( pde_vecR_overlap_ == ROL::nullPtr ) { // Overlapping components of residual vector
    pde_vecR_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapResidualMap_, 1, true);
  }
  // Get u_coeff from u and z_coeff from z
  ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff;
  getCoeffFromStateVector(uo_coeff,uo);
  getCoeffFromStateVector(un_coeff,un);
  getCoeffFromControlVector(z_coeff,z);
  // Assemble
  ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
  pde->residual(localVal,ts,uo_coeff,un_coeff,z_coeff,z_param);
  assembleFieldVector( r, localVal, pde_vecR_overlap_ );
}

template<class Real>
void Assembler<Real>::assembleDynPDEJacobian_uo(ROL::Ptr<Tpetra::CrsMatrix<>> &J,
                                          const ROL::Ptr<DynamicPDE<Real>> &pde,
                                          const ROL::TimeStamp<Real> &ts,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                          const ROL::Ptr<const std::vector<Real>> &z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEJacobian_uo);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize Jacobian matrices
    if ( J == ROL::nullPtr ) {
      J = ROL::makePtr<Tpetra::CrsMatrix<>>(matJ1Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Jacobian_uo(localVal,ts,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( J, localVal );
    isJuoTransposed_ = false;
  }
  catch ( Exception::Zero & ez ) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEJacobian_uo): Jacobian is zero.");
  }
  catch ( Exception::NotImplemented & eni ) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEJacobian_uo): Jacobian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEJacobian_un(ROL::Ptr<Tpetra::CrsMatrix<>> &J,
                                          const ROL::Ptr<DynamicPDE<Real>> &pde,
                                          const ROL::TimeStamp<Real> &ts,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                          const ROL::Ptr<const std::vector<Real>> &z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEJacobian_un);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize Jacobian matrices
    if ( J == ROL::nullPtr ) {
      J = ROL::makePtr<Tpetra::CrsMatrix<>>(matJ1Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Jacobian_un(localVal,ts,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( J, localVal );
    isJunTransposed_ = false;
  }
  catch ( Exception::Zero & ez ) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEJacobian_un): Jacobian is zero.");
  }
  catch ( Exception::NotImplemented & eni ) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEJacobian_un): Jacobian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEJacobian_zf(ROL::Ptr<Tpetra::CrsMatrix<>> &J,
                                          const ROL::Ptr<DynamicPDE<Real>> &pde,
                                          const ROL::TimeStamp<Real> &ts,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                          const ROL::Ptr<const std::vector<Real>> &z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEJacobian_zf);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize Jacobian matrices
    if ( J == ROL::nullPtr ) {
      J = ROL::makePtr<Tpetra::CrsMatrix<>>(matJ2Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Jacobian_zf(localVal,ts,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( J, localVal );
    isJzfTransposed_ = false;
  }
  catch ( Exception::Zero & ez ) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEJacobian_zf): Jacobian is zero.");
  }
  catch ( Exception::NotImplemented & eni ) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEJacobian_zf): Jacobian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEJacobian_zp(ROL::Ptr<Tpetra::MultiVector<>> &J,
                                          const ROL::Ptr<DynamicPDE<Real>> &pde,
                                          const ROL::TimeStamp<Real> &ts,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                          const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEJacobian_zp);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff;
      getCoeffFromStateVector(uo_coeff,uo);
      getCoeffFromStateVector(un_coeff,un);
      getCoeffFromControlVector(z_coeff,z);
      // Initialize Jacobian storage if not done so already
      if (J == ROL::nullPtr) {
        J = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueResidualMap_, size, true);
      }
      if ( pde_vecJ3_overlap_ == ROL::nullPtr) {
        pde_vecJ3_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapResidualMap_, size, true);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      pde->Jacobian_zp(localVal,ts,uo_coeff,un_coeff,z_coeff,z_param);
      assembleParamFieldMatrix( J, localVal, pde_vecJ3_overlap_ );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleDynPDEJacobian_zp): Jacobian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEJacobian_zp): Jacobian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEJacobian_zp): Jacobian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_uo_uo(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_uo_uo);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH11Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_uo_uo(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_uo_uo): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_uo_uo): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_uo_un(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_uo_un);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH11Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_uo_un(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_uo_un): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_uo_un): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_uo_zf(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_uo_zf);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH12Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_uo_zf(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_uo_zf): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_uo_zf): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_uo_zp(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_uo_zp);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u, z_coeff from z and l_coeff from l
      ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
      getCoeffFromStateVector(uo_coeff,uo);
      getCoeffFromStateVector(un_coeff,un);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Initialize Jacobian storage if not done so already
      if (H == ROL::nullPtr) {
        H = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, size, true);
      }
      if ( pde_vecH13_overlap_ == ROL::nullPtr) {
        pde_vecH13_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, size, true);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      pde->Hessian_uo_zp(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
      assembleParamFieldMatrix( H, localVal, pde_vecH13_overlap_ );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_uo_zp): Hessian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_uo_zp): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_uo_zp): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_un_uo(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_un_uo);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH11Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_un_uo(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_un_uo): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_un_uo): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_un_un(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_un_un);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH11Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_un_un(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_un_un): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_un_un): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_un_zf(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_un_zf);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH12Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_un_zf(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_un_zf): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_un_zf): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_un_zp(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_un_zp);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u, z_coeff from z and l_coeff from l
      ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
      getCoeffFromStateVector(uo_coeff,uo);
      getCoeffFromStateVector(un_coeff,un);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Initialize Jacobian storage if not done so already
      if (H == ROL::nullPtr) {
        H = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, size, true);
      }
      if ( pde_vecH13_overlap_ == ROL::nullPtr) {
        pde_vecH13_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, size, true);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      pde->Hessian_un_zp(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
      assembleParamFieldMatrix( H, localVal, pde_vecH13_overlap_ );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_un_zp): Hessian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_un_zp): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_un_zp): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_zf_uo(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_zf_uo);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH21Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_zf_uo(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_zf_uo): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_zf_uo): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_zf_un(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_zf_un);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH21Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_zf_un(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_zf_un): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_zf_un): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_zf_zf(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_zf_zf);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
    getCoeffFromStateVector(uo_coeff,uo);
    getCoeffFromStateVector(un_coeff,un);
    getCoeffFromControlVector(z_coeff,z);
    getCoeffFromStateVector(l_coeff,l);
    // Initialize Hessian storage if not done so already
    if ( H == ROL::nullPtr ) {
      H = ROL::makePtr<Tpetra::CrsMatrix<>>(matH22Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->Hessian_zf_zf(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
    assembleFieldMatrix( H, localVal );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_zf_zf): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_zf_zf): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_zf_zp(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_zf_zp);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u, z_coeff from z and l_coeff from l
      ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
      getCoeffFromStateVector(uo_coeff,uo);
      getCoeffFromStateVector(un_coeff,un);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Initialize Jacobian storage if not done so already
      if (H == ROL::nullPtr) {
        H = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, size, true);
      }
      if ( pde_vecH23_overlap_ == ROL::nullPtr) {
        pde_vecH23_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, size, true);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      pde->Hessian_zf_zp(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
      assembleParamFieldMatrix( H, localVal, pde_vecH23_overlap_ );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_zf_zp): Hessian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_zf_zp): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_zf_zp): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_zp_uo(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_zp_uo);
  #endif
  assembleDynPDEHessian_uo_zp(H,pde,ts,l,uo,un,z,z_param);
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_zp_un(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_zp_un);
  #endif
  assembleDynPDEHessian_un_zp(H,pde,ts,l,uo,un,z,z_param);
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_zp_zf(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_zp_zf);
  #endif
  assembleDynPDEHessian_zf_zp(H,pde,ts,l,uo,un,z,z_param);
}

template<class Real>
void Assembler<Real>::assembleDynPDEHessian_zp_zp(ROL::Ptr<std::vector<std::vector<Real>>> &H,
                                            const ROL::Ptr<DynamicPDE<Real>> &pde,
                                            const ROL::TimeStamp<Real> &ts,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                            const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                            const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleDynPDEHessian_zp_zp);
  #endif
  if ( z_param != ROL::nullPtr ) {
    try {
      int size = static_cast<int>(z_param->size());
      // Get u_coeff from u, z_coeff from z and l_coeff from l
      ROL::Ptr<Intrepid::FieldContainer<Real>> uo_coeff, un_coeff, z_coeff, l_coeff;
      getCoeffFromStateVector(uo_coeff,uo);
      getCoeffFromStateVector(un_coeff,un);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Initialize Jacobian storage if not done so already
      if (H == ROL::nullPtr) {
        std::vector<Real> col(size,static_cast<Real>(0));
        H = ROL::makePtr<std::vector<std::vector<Real>>>(size,col);
      }
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> tmp(size,ROL::nullPtr);
      std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> localVal(size,tmp);
      pde->Hessian_zp_zp(localVal,ts,l_coeff,uo_coeff,un_coeff,z_coeff,z_param);
      assembleParamMatrix( H, localVal );
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleDynPDEHessian_zp_zp): Hessian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_zp_zp): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleDynPDEHessian_zp_zp): Hessian not implemented.");
  }
}
/***************************************************************************/
/* End of PDE assembly routines.                                           */
/***************************************************************************/

/***************************************************************************/
/* QoI assembly routines                                                   */
/***************************************************************************/
template<class Real>
Real Assembler<Real>::assembleQoIValue(const ROL::Ptr<QoI<Real>> &qoi,
                                       const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                       const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                       const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIValue);
  #endif
  Real val(0);
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    val  = qoi->value(localVal,u_coeff,z_coeff,z_param);
    val += assembleScalar( localVal );
  }
  catch ( Exception::Zero & ez ) {
    val = static_cast<Real>(0);
  }
  catch ( Exception::NotImplemented & eni ) {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIValue): Value not implemented.");
  }
  return val;
}

template<class Real>
void Assembler<Real>::assembleQoIGradient1(ROL::Ptr<Tpetra::MultiVector<>> &g1,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient1);
  #endif
  try {
    // Initialize state QoI gradient vectors
    if ( g1 == ROL::nullPtr ) {
      g1 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
    }
    if ( qoi_vecG1_overlap_ == ROL::nullPtr ) {
      qoi_vecG1_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
    }
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    qoi->gradient_1(localVal,u_coeff,z_coeff,z_param);
    assembleFieldVector( g1, localVal, qoi_vecG1_overlap_ );
  }
  catch ( Exception::Zero & ez ) {
    throw Exception::Zero(">>> (Assembler::assembleQoIGradient1): Gradient is zero.");
  }
  catch ( Exception::NotImplemented & eni ) {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient1): Gradient not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIGradient2(ROL::Ptr<Tpetra::MultiVector<>> &g2,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient2);
  #endif
  try {
    // Initialize control gradient vectors
    if ( g2 == ROL::nullPtr ) {
      g2 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
    }
    if ( qoi_vecG2_overlap_ == ROL::nullPtr ) {
      qoi_vecG2_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
    }
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    qoi->gradient_2(localVal,u_coeff,z_coeff,z_param);
    assembleFieldVector( g2, localVal, qoi_vecG2_overlap_ );
  }
  catch ( Exception::Zero & ez ) {
    throw Exception::Zero(">>> (Assembler::assembleQoIGradient2): Gradient is zero.");
  }
  catch ( Exception::NotImplemented & eni ) {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient2): Gradient not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIGradient3(ROL::Ptr<std::vector<Real>> &g3,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient3);
  #endif
  if ( z_param != ROL::nullPtr ) {
    const int size = z_param->size();
    if ( g3 == ROL::nullPtr ) {
      g3 = ROL::makePtr<std::vector<Real>>(size,0);
    }
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      std::vector<Real> g = qoi->gradient_3(localVal,u_coeff,z_coeff,z_param);
      assembleParamVector( g3, localVal );
      for (int i = 0; i < size; ++i) {
        (*g3)[i] += g[i];
      }
    }
    catch ( Exception::Zero & ez ) {
      g3->assign(size,0);
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient3): Gradient not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient3): Gradient not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIHessVec11(ROL::Ptr<Tpetra::MultiVector<>> &H11,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec11);
  #endif
	  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff, u_coeff, z_coeff;
    getCoeffFromStateVector(v_coeff,v);
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize state-state HessVec vectors
    if ( H11 == ROL::nullPtr ) {
      H11 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
    }
    if ( qoi_vecH11_overlap_ == ROL::nullPtr ) {
      qoi_vecH11_overlap_  = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    qoi->HessVec_11(localVal,v_coeff,u_coeff,z_coeff,z_param);
    assembleFieldVector( H11, localVal, qoi_vecH11_overlap_ );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleQoIHessVec11): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec11): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIHessVec12(ROL::Ptr<Tpetra::MultiVector<>> &H12,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec12);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff, u_coeff, z_coeff;
    getCoeffFromControlVector(v_coeff,v);
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize state-control HessVec vectors
    if ( H12 == ROL::nullPtr ) {
      H12 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
    }
    if ( qoi_vecH12_overlap_ == ROL::nullPtr ) {
      qoi_vecH12_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    qoi->HessVec_12(localVal,v_coeff,u_coeff,z_coeff,z_param);
    assembleFieldVector( H12, localVal, qoi_vecH12_overlap_ );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleQoIHessVec12): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec12): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIHessVec13(ROL::Ptr<Tpetra::MultiVector<>> &H13,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const std::vector<Real>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> &z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec13);
  #endif
  if (z_param != ROL::nullPtr) {
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Initialize state-control HessVec vectors
      if ( H13 == ROL::nullPtr ) {
        H13 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
      }
      if ( qoi_vecH13_overlap_ == ROL::nullPtr ) {
        qoi_vecH13_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
      }
      // Assemble
      ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
      qoi->HessVec_13(localVal,v,u_coeff,z_coeff,z_param);
      assembleFieldVector( H13, localVal, qoi_vecH13_overlap_ );
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec13): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec13): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec13): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIHessVec21(ROL::Ptr<Tpetra::MultiVector<>> &H21,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec21);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff, u_coeff, z_coeff;
    getCoeffFromStateVector(v_coeff,v);
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize control-state HessVec vectors
    if ( H21 == ROL::nullPtr ) {
      H21 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
    }
    if ( qoi_vecH21_overlap_ == ROL::nullPtr ) {
      qoi_vecH21_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    qoi->HessVec_21(localVal,v_coeff,u_coeff,z_coeff,z_param);
    assembleFieldVector( H21, localVal, qoi_vecH21_overlap_ );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleQoIHessVec21): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec21): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIHessVec22(ROL::Ptr<Tpetra::MultiVector<>> &H22,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec22);
  #endif
  try {
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff, u_coeff, z_coeff;
    getCoeffFromControlVector(v_coeff,v);
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Initialize control-control HessVec vectors
    if ( H22 == ROL::nullPtr ) {
      H22 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
    }
    if ( qoi_vecH22_overlap_ == ROL::nullPtr ) {
      qoi_vecH22_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    qoi->HessVec_22(localVal,v_coeff,u_coeff,z_coeff,z_param);
    assembleFieldVector( H22, localVal, qoi_vecH22_overlap_ );
  }
  catch (Exception::Zero &ez) {
    throw Exception::Zero(">>> (Assembler::assembleQoIHessVec22): Hessian is zero.");
  }
  catch (Exception::NotImplemented &eni) {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec22): Hessian not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIHessVec23(ROL::Ptr<Tpetra::MultiVector<>> &H23,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const std::vector<Real>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec23);
  #endif
  if (z_param != ROL::nullPtr) {
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Initialize control-control HessVec vectors
      if ( H23 == ROL::nullPtr ) {
        H23 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
      }
      if ( qoi_vecH23_overlap_ == ROL::nullPtr ) {
        qoi_vecH23_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
      }
      // Assemble
      ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
      qoi->HessVec_23(localVal,v,u_coeff,z_coeff,z_param);
      assembleFieldVector( H23, localVal, qoi_vecH23_overlap_ );
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec23): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec23): Hessian not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec23): Hessian not implemented.");
  }

}

template<class Real>
void Assembler<Real>::assembleQoIHessVec31(ROL::Ptr<std::vector<Real>> &H31,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec31);
  #endif
  if ( z_param != ROL::nullPtr ) {
    const int size = z_param->size();
    if ( H31 == ROL::nullPtr ) {
      H31 = ROL::makePtr<std::vector<Real>>(size,0);
    }
    try {
      H31->assign(size,0);
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff, u_coeff, z_coeff;
      getCoeffFromStateVector(v_coeff,v);
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      std::vector<Real> H = qoi->HessVec_31(localVal,v_coeff,u_coeff,z_coeff,z_param);
      assembleParamVector( H31, localVal );
      for (int i = 0; i < size; ++i) {
        (*H31)[i] += H[i];
      }
    }
    catch ( Exception::Zero & ez ) {
      H31->assign(size,0);
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec31): HessVec not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec31): HessVec not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIHessVec32(ROL::Ptr<std::vector<Real>> &H32,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> & z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec32);
  #endif
  if ( z_param != ROL::nullPtr ) {
    const int size = z_param->size();
    if ( H32 == ROL::nullPtr ) {
      H32 = ROL::makePtr<std::vector<Real>>(size,0);
    }
    try {
      H32->assign(size,0);
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real>> v_coeff, u_coeff, z_coeff;
      getCoeffFromControlVector(v_coeff,v);
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      std::vector<Real> H = qoi->HessVec_32(localVal,v_coeff,u_coeff,z_coeff,z_param);
      assembleParamVector( H32, localVal );
      for (int i = 0; i < size; ++i) {
        (*H32)[i] += H[i];
      }
    }
    catch ( Exception::Zero & ez ) {
      H32->assign(size,0);
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec32): HessVec not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec32): HessVec not implemented.");
  }
}

template<class Real>
void Assembler<Real>::assembleQoIHessVec33(ROL::Ptr<std::vector<Real>> &H33,
                                     const ROL::Ptr<QoI<Real>> &qoi,
                                     const ROL::Ptr<const std::vector<Real>> &v,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                     const ROL::Ptr<const Tpetra::MultiVector<>> &z,
                                     const ROL::Ptr<const std::vector<Real>> &z_param) {
  #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec33);
  #endif
  if ( z_param != ROL::nullPtr ) {
    const int size = z_param->size();
    if ( H33 == ROL::nullPtr ) {
      H33 = ROL::makePtr<std::vector<Real>>(size,0);
    }
    try {
      H33->assign(size,0);
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff, z_coeff;
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Assemble
      std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> localVal(size,ROL::nullPtr);
      std::vector<Real> H = qoi->HessVec_33(localVal,v,u_coeff,z_coeff,z_param);
      assembleParamVector( H33, localVal );
      for (int i = 0; i < size; ++i) {
        (*H33)[i] += H[i];
      }
    }
    catch ( Exception::Zero & ez ) {
      H33->assign(size,0);
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec33): HessVec not implemented.");
    }
  }
  else {
    throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec33): HessVec not implemented.");
  }
}
/***************************************************************************/
/* End QoI assembly routines                                               */
/***************************************************************************/

/***************************************************************************/
/* Assemble and apply Riesz operator corresponding to simulation variables */
/***************************************************************************/
template<class Real>
void Assembler<Real>::assemblePDERieszMap1(ROL::Ptr<Tpetra::CrsMatrix<>> &R1,
                                     const ROL::Ptr<PDE<Real>> &pde) {
  try {
    // Initialize Riesz matrix if not done so already
    if ( R1 == ROL::nullPtr ) {
    R1 = ROL::makePtr<Tpetra::CrsMatrix<>>(matR1Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->RieszMap_1(localVal);
    assembleFieldMatrix( R1, localVal );
  }
  catch ( Exception::NotImplemented & eni ) {
    //throw Exception::NotImplemented(">>> (Assembler::assemblePDERieszMap1): Riesz map not implemented!");
  }
  catch ( Exception::Zero & ez ) {
    //throw Exception::Zero(">>> (Assembler::assemblePDERieszMap1): Riesz map is zero!");
  }
}
template<class Real>
void Assembler<Real>::assembleDynPDERieszMap1(ROL::Ptr<Tpetra::CrsMatrix<>> &R1,
                                        const ROL::Ptr<DynamicPDE<Real>> &pde) {
  try {
    // Initialize Riesz matrix if not done so already
    if ( R1 == ROL::nullPtr ) {
    R1 = ROL::makePtr<Tpetra::CrsMatrix<>>(matR1Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->RieszMap_1(localVal);
    assembleFieldMatrix( R1, localVal );
  }
  catch ( Exception::NotImplemented & eni ) {
    //throw Exception::NotImplemented(">>> (Assembler::assemblePDERieszMap1): Riesz map not implemented!");
  }
  catch ( Exception::Zero & ez ) {
    //throw Exception::Zero(">>> (Assembler::assemblePDERieszMap1): Riesz map is zero!");
  }
}
/***************************************************************************/
/* End of functions for Riesz operator of simulation variables.            */
/***************************************************************************/

/***************************************************************************/
/* Assemble and apply Riesz operator corresponding to optimization         */
/* variables                                                               */
/***************************************************************************/
template<class Real>
void Assembler<Real>::assemblePDERieszMap2(ROL::Ptr<Tpetra::CrsMatrix<>> &R2,
                                     const ROL::Ptr<PDE<Real>> &pde) {
  try {
    // Initialize Riesz matrix if not done so already
    if ( R2 == ROL::nullPtr ) {
      R2 = ROL::makePtr<Tpetra::CrsMatrix<>>(matR2Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->RieszMap_2(localVal);
    assembleFieldMatrix( R2, localVal );
  }
  catch ( Exception::NotImplemented & eni ) {
    //throw Exception::NotImplemented(">>> (Assembler::assemblePDERieszMap2): Riesz map not implemented!");
  }
  catch ( Exception::Zero & ez ) {
    //throw Exception::Zero(">>> (Assembler::assemblePDERieszMap2): Riesz map is zero!");
  }
}
template<class Real>
void Assembler<Real>::assembleDynPDERieszMap2(ROL::Ptr<Tpetra::CrsMatrix<>> &R2,
                                        const ROL::Ptr<DynamicPDE<Real>> &pde) {
  try {
    // Initialize Riesz matrix if not done so already
    if ( R2 == ROL::nullPtr ) {
      R2 = ROL::makePtr<Tpetra::CrsMatrix<>>(matR2Graph_);
    }
    // Assemble
    ROL::Ptr<Intrepid::FieldContainer<Real>> localVal;
    pde->RieszMap_2(localVal);
    assembleFieldMatrix( R2, localVal );
  }
  catch ( Exception::NotImplemented & eni ) {
    //throw Exception::NotImplemented(">>> (Assembler::assemblePDERieszMap2): Riesz map not implemented!");
  }
  catch ( Exception::Zero & ez ) {
    //throw Exception::Zero(">>> (Assembler::assemblePDERieszMap2): Riesz map is zero!");
  }
}
/***************************************************************************/
/* End of functions for Riesz operator of optimization variables.          */
/***************************************************************************/

/***************************************************************************/
/* Compute error routines.                                                 */
/***************************************************************************/
template<class Real>
Real Assembler<Real>::computeStateError(const ROL::Ptr<const Tpetra::MultiVector<>> &soln,
                                       const ROL::Ptr<Solution<Real>> &trueSoln,
                                       const int cubDeg,
                                       const ROL::Ptr<FieldHelper<Real>> &fieldHelper) const {
  Real totalError(0);
  // populate inCoeffs
  ROL::Ptr<Intrepid::FieldContainer<Real>> inCoeffs0;
  getCoeffFromStateVector(inCoeffs0, soln);
  // split fields
  int numFields = 1;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> inCoeffs;
  if (fieldHelper != ROL::nullPtr) {
    numFields = fieldHelper->numFields();
    fieldHelper->splitFieldCoeff(inCoeffs,inCoeffs0);
  }
  else {
    inCoeffs.push_back(inCoeffs0);
  }
  // compute error
  for (int fn = 0; fn < numFields; ++fn) {
    // create fe object for error computation
    Intrepid::DefaultCubatureFactory<Real> cubFactory;
    shards::CellTopology cellType = basisPtrs_[fn]->getBaseCellTopology();
    ROL::Ptr<Intrepid::Cubature<Real>> cellCub = cubFactory.create(cellType, cubDeg);
    ROL::Ptr<FE<Real>> fe
      = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrs_[fn],cellCub);

    // get dimensions
    int c = fe->gradN()->dimension(0);
    int p = fe->gradN()->dimension(2);
    int d = fe->gradN()->dimension(3);

    // evaluate input coefficients on fe basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> funcVals
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe->evaluateValue(funcVals, inCoeffs[fn]);

    // subtract off true solution
    std::vector<Real> x(d);
    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        for (int k=0; k<d; ++k) {
          x[k] = (*fe->cubPts())(i, j, k);
        }
        (*funcVals)(i, j) -= trueSoln->evaluate(x,fn);
      }
    }

    // compute norm squared of local error
    ROL::Ptr<Intrepid::FieldContainer<Real>> normSquaredError
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    fe->computeIntegral(normSquaredError, funcVals, funcVals, false);

    Real localErrorSum(0);
    Real globalErrorSum(0);
    for (int i=0; i<numCells_; ++i) {
      localErrorSum += (*normSquaredError)(i);
    }
    Teuchos::reduceAll<int, Real>(*comm_, Teuchos::REDUCE_SUM, 1, &localErrorSum, &globalErrorSum);
    totalError += globalErrorSum;
  }

  return std::sqrt(totalError);
}

template<class Real>
Real Assembler<Real>::computeControlError(const ROL::Ptr<const Tpetra::MultiVector<>> &soln,
                                          const ROL::Ptr<Solution<Real>> &trueSoln,
                                          const int cubDeg,
                                          const ROL::Ptr<FieldHelper<Real>> &fieldHelper) const {
  Real totalError(0);
  // populate inCoeffs
  ROL::Ptr<Intrepid::FieldContainer<Real>> inCoeffs0;
  getCoeffFromControlVector(inCoeffs0, soln);
  // split fields
  int numFields = 1;
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> inCoeffs;
  if (fieldHelper != ROL::nullPtr) {
    numFields = fieldHelper->numFields();
    fieldHelper->splitFieldCoeff(inCoeffs,inCoeffs0);
  }
  else {
    inCoeffs.push_back(inCoeffs0);
  }
  // compute error
  for (int fn = 0; fn < numFields; ++fn) {
    // create fe object for error computation
    Intrepid::DefaultCubatureFactory<Real> cubFactory;
    shards::CellTopology cellType = basisPtrs_[fn]->getBaseCellTopology();
    ROL::Ptr<Intrepid::Cubature<Real>> cellCub = cubFactory.create(cellType, cubDeg);
    ROL::Ptr<FE<Real>> fe
      = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrs_[fn],cellCub);

    // get dimensions
    int c = fe->gradN()->dimension(0);
    int p = fe->gradN()->dimension(2);
    int d = fe->gradN()->dimension(3);

    // evaluate input coefficients on fe basis
    ROL::Ptr<Intrepid::FieldContainer<Real>> funcVals
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
    fe->evaluateValue(funcVals, inCoeffs[fn]);

    // subtract off true solution
    std::vector<Real> x(d);
    for (int i=0; i<c; ++i) {
      for (int j=0; j<p; ++j) {
        for (int k=0; k<d; ++k) {
          x[k] = (*fe->cubPts())(i, j, k);
        }
        (*funcVals)(i, j) -= trueSoln->evaluate(x,fn);
      }
    }

    // compute norm squared of local error
    ROL::Ptr<Intrepid::FieldContainer<Real>> normSquaredError
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
    fe->computeIntegral(normSquaredError, funcVals, funcVals, false);

    Real localErrorSum(0);
    Real globalErrorSum(0);
    for (int i=0; i<numCells_; ++i) {
      localErrorSum += (*normSquaredError)(i);
    }
    Teuchos::reduceAll<int, Real>(*comm_, Teuchos::REDUCE_SUM, 1, &localErrorSum, &globalErrorSum);
    totalError += globalErrorSum;
  }

  return std::sqrt(totalError);
}
/***************************************************************************/
/* End of compute solution routines.                                       */
/***************************************************************************/

/***************************************************************************/
/* Output routines.                                                        */
/***************************************************************************/
template<class Real>
void Assembler<Real>::printMeshData(std::ostream &outStream) const {
  ROL::Ptr<Intrepid::FieldContainer<Real>> nodesPtr = meshMgr_->getNodes();
  ROL::Ptr<Intrepid::FieldContainer<int>>  cellToNodeMapPtr = meshMgr_->getCellToNodeMap();
  Intrepid::FieldContainer<Real>  &nodes = *nodesPtr;
  Intrepid::FieldContainer<int>   &cellToNodeMap = *cellToNodeMapPtr;
  if ( verbose_ && myRank_ == 0) {
    outStream << std::endl;
    outStream << "Number of nodes = " << meshMgr_->getNumNodes() << std::endl;
    outStream << "Number of cells = " << meshMgr_->getNumCells() << std::endl;
    outStream << "Number of edges = " << meshMgr_->getNumEdges() << std::endl;
  }
  // Print mesh to file.
  if ((myRank_ == 0)) {
    std::ofstream meshfile;
    meshfile.open("cell_to_node_quad.txt");
    for (int i=0; i<cellToNodeMap.dimension(0); ++i) {
      for (int j=0; j<cellToNodeMap.dimension(1); ++j) {
        meshfile << cellToNodeMap(i,j) << "  ";
      }
      meshfile << std::endl;
    }
    meshfile.close();
    
    meshfile.open("cell_to_node_tri.txt");
    for (int i=0; i<cellToNodeMap.dimension(0); ++i) {
      for (int j=0; j<3; ++j) {
        meshfile << cellToNodeMap(i,j) << "  ";
      }
      meshfile << std::endl;
      for (int j=2; j<5; ++j) {
        meshfile << cellToNodeMap(i,j%4) << "  ";
      }
      meshfile << std::endl;
    }
    meshfile.close();
   
    meshfile.open("nodes.txt");
    meshfile.precision(16);
    for (int i=0; i<nodes.dimension(0); ++i) {
      for (int j=0; j<nodes.dimension(1); ++j) {
        meshfile << std::scientific << nodes(i,j) << "  ";
      }
      meshfile << std::endl;
    }
    meshfile.close();
  }
} // prinf function end

template<class Real>
void Assembler<Real>::outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<>> &vec,
                                         const std::string &filename) const {
  Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<>> vecWriter;
  vecWriter.writeDenseFile(filename, vec);
  std::string mapfile = "map_" + filename;
  vecWriter.writeMapFile(mapfile, *(vec->getMap()));
}

template<class Real>
void Assembler<Real>::inputTpetraVector(ROL::Ptr<Tpetra::MultiVector<>> &vec,
                                        const std::string &filename) const {
  Tpetra::MatrixMarket::Reader< Tpetra::CrsMatrix<>> vecReader;
  auto map = vec->getMap();
  auto comm = map->getComm();
  auto vec_rcp = vecReader.readDenseFile(filename, comm, map);
  vec->scale(1.0, *vec_rcp);
}

template<class Real>
void Assembler<Real>::serialPrintStateEdgeField(const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                                const ROL::Ptr<FieldHelper<Real>> &fieldHelper,
                                                const std::string &filename,
                                                const ROL::Ptr<FE_CURL<Real>> &fe) const {
  const int c = fe->curlN()->dimension(0);
  const int f = fe->curlN()->dimension(1);
  const int p = 1, d = 3;

  ROL::Ptr<Intrepid::FieldContainer<Real>> u_coeff;
  getCoeffFromStateVector(u_coeff,u);

  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> U;
  fieldHelper->splitFieldCoeff(U, u_coeff);
  int numFields = U.size();

  // Transform cell center to physical
  ROL::Ptr<Intrepid::FieldContainer<Real>> rx
    = ROL::makePtr<Intrepid::FieldContainer<Real>>(p,d);
  ROL::Ptr<Intrepid::FieldContainer<Real>> px
    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
  fe->mapRefPointsToPhysical(px,rx);
  // Transform reference values into physical space.
  ROL::Ptr<Intrepid::FieldContainer<Real>> cellJac
    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d,d);
  ROL::Ptr<Intrepid::FieldContainer<Real>> cellJacInv
    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d,d);
  ROL::Ptr<shards::CellTopology> cellTopo
    = ROL::makePtr<shards::CellTopology>(basisPtrs_[0]->getBaseCellTopology());
  ROL::Ptr<Intrepid::FieldContainer<Real>> valReference
    = ROL::makePtr<Intrepid::FieldContainer<Real>>(f,p,d);
  basisPtrs_[0]->getValues(*valReference,*rx,Intrepid::OPERATOR_VALUE);
  ROL::Ptr<Intrepid::FieldContainer<Real>> valPhysical
    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p,d);
  Intrepid::CellTools<Real>::setJacobian(*cellJac,*rx,*volCellNodes_,*cellTopo);
  Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv, *cellJac);
  Intrepid::FunctionSpaceTools::HCURLtransformVALUE<Real>(*valPhysical,
                                                          *cellJacInv,
                                                          *valReference);

  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> uval(numFields);
  for (int k = 0; k < numFields; ++k) {
    uval[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    Intrepid::FunctionSpaceTools::evaluate<Real>(*uval[k], *U[k], *valPhysical);
  }
  // Print
  std::ofstream file;
  file.open(filename);
  for (int i = 0; i < c; ++i) {
    file << std::scientific << std::setprecision(15);
    for (int j = 0; j < d; ++j) {
      file << std::setw(25) << (*px)(i,0,j);
    }
    for (int k = 0; k < numFields; ++k) {
      for (int j = 0; j < d; ++j) {
        file << std::setw(25) << (*uval[k])(i,0,j);
      }
    }
    file << std::endl;
  }
  file.close();
}
/***************************************************************************/
/* End of output routines.                                                 */
/***************************************************************************/

/***************************************************************************/
/* Vector generation routines.                                             */
/***************************************************************************/
template<class Real>
const ROL::Ptr<const Tpetra::Map<>> Assembler<Real>::getStateMap(void) const {
  return myUniqueStateMap_;
}

template<class Real>
const ROL::Ptr<const Tpetra::Map<>> Assembler<Real>::getControlMap(void) const {
  return myUniqueControlMap_;
}

template<class Real>
const ROL::Ptr<const Tpetra::Map<>> Assembler<Real>::getResidualMap(void) const {
  return myUniqueResidualMap_;
}

template<class Real>
ROL::Ptr<Tpetra::MultiVector<>> Assembler<Real>::createStateVector(void) const {
  return ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
}

template<class Real>
ROL::Ptr<Tpetra::MultiVector<>> Assembler<Real>::createControlVector(void) const {
  return ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
}

template<class Real>
ROL::Ptr<Tpetra::MultiVector<>> Assembler<Real>::createResidualVector(void) const {
  return ROL::makePtr<Tpetra::MultiVector<>>(myUniqueResidualMap_, 1, true);
}
/***************************************************************************/
/* End of vector generation routines.                                      */
/***************************************************************************/

/***************************************************************************/
/* Accessor routines.                                                      */
/***************************************************************************/
template<class Real>
const ROL::Ptr<DofManager<Real>> Assembler<Real>::getDofManager(void) const {
  return dofMgr_;
}

template<class Real>
Teuchos::Array<int> Assembler<Real>::getCellIds(void) const {
  return myCellIds_;
}
/***************************************************************************/
/* End of accessor routines.                                               */
/***************************************************************************/

/*****************************************************************************/
/******************* PRIVATE MEMBER FUNCTIONS ********************************/
/*****************************************************************************/
template<class Real>
void Assembler<Real>::setCommunicator(const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                                      Teuchos::ParameterList &parlist,
                                      std::ostream &outStream) {
  comm_ = comm;
  // Get number of processors and my rank
  myRank_    = comm->getRank();
  numProcs_  = comm->getSize();
  // Parse parameter list
  verbose_ = parlist.sublist("PDE FEM").get("Verbose Output",false);
  if (verbose_ && myRank_==0 ) {
    outStream << "Initialized communicator. " << std::endl;
  }
  if (verbose_ && myRank_==0 ) {
    outStream << "Total number of processors: " << numProcs_ << std::endl;
  }
}

template<class Real>
void Assembler<Real>::setBasis(
       const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs,
       Teuchos::ParameterList &parlist,
       std::ostream &outStream) {
  basisPtrs_ = basisPtrs;
  if (verbose_ && myRank_==0) {
    outStream << "Initialized PDE." << std::endl;
  }
}

template<class Real>
void Assembler<Real>::setDiscretization(Teuchos::ParameterList &parlist,
                                  const ROL::Ptr<MeshManager<Real>> &meshMgr,
                                        std::ostream &outStream) {
  if ( meshMgr != ROL::nullPtr ) {
    // Use MeshManager object if supplied
    meshMgr_ = meshMgr;
  }
  else {
    // Otherwise construct MeshManager objective from parameter list
  }
  dofMgr_ = ROL::makePtr<DofManager<Real>>(meshMgr_,basisPtrs_);
  if (verbose_ && myRank_==0) {
    outStream << "Initialized discretization (MeshManager and DofManager)." << std::endl;
  }
}

template<class Real>
void Assembler<Real>::setParallelStructure(Teuchos::ParameterList &parlist,
                                           std::ostream &outStream) {
  int cellSplit = parlist.sublist("Geometry").get<int>("Partition type");
  /****************************************************/
  /*** Build parallel communication infrastructure. ***/
  /****************************************************/
  // Partition the cells in the mesh.  We use a basic quasi-equinumerous partitioning,
  // where the remainder, if any, is assigned to the last processor.
  Teuchos::Array<int> myGlobalIds;
  cellOffsets_.assign(numProcs_, 0);
  int totalNumCells = meshMgr_->getNumCells();
  int cellsPerProc  = totalNumCells / numProcs_;
  numCells_         = cellsPerProc;
  switch(cellSplit) {
    case 0:
      if (myRank_ == 0) {
        // remainder in the first
        numCells_ += totalNumCells % numProcs_;
      }
      for (int i=1; i<numProcs_; ++i) {
        cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc
                          + (static_cast<int>(i==1))*(totalNumCells % numProcs_);
      }
      break;
    case 1:
      if (myRank_ == numProcs_-1) {
        // remainder in the last
        numCells_ += totalNumCells % numProcs_;
      }
      for (int i=1; i<numProcs_; ++i) {
        cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc;
      }
      break;
    case 2:
      if (myRank_ < (totalNumCells%numProcs_)) {
        // spread remainder, starting from the first
        numCells_++;
      }
      for (int i=1; i<numProcs_; ++i) {
        cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc
                          + (static_cast<int>(i-1<(totalNumCells%numProcs_)));
      }
      break;
  }

  Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
  int numLocalDofs = cellDofs.dimension(1);
  if (verbose_) {
    outStream << "Cell offsets across processors: " << cellOffsets_
              << std::endl;
  }
  for (int i=0; i<numCells_; ++i) {
    myCellIds_.push_back(cellOffsets_[myRank_]+i);
    for (int j=0; j<numLocalDofs; ++j) {
      myGlobalIds.push_back( cellDofs(cellOffsets_[myRank_]+i,j) );
    }
  }
  std::sort(myGlobalIds.begin(), myGlobalIds.end());
  myGlobalIds.erase( std::unique(myGlobalIds.begin(),myGlobalIds.end()),myGlobalIds.end() );

  // Build maps.
  myOverlapStateMap_ = ROL::makePtr<Tpetra::Map<>>(
                       Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                       myGlobalIds, 0, comm_);
  //std::cout << std::endl << myOverlapMap_->getNodeElementList()<<std::endl;
  /** One can also use the non-member function:
      myOverlapMap_ = Tpetra::createNonContigMap<int,int>(myGlobalIds_, comm_);
      to build the overlap map.
  **/
  myUniqueStateMap_ = Tpetra::createOneToOne<int,int>(myOverlapStateMap_);
  //std::cout << std::endl << myUniqueMap_->getNodeElementList() << std::endl;
  myOverlapControlMap_  = myOverlapStateMap_;
  myUniqueControlMap_   = myUniqueStateMap_;
  myOverlapResidualMap_ = myOverlapStateMap_;
  myUniqueResidualMap_  = myUniqueStateMap_;
  //  myCellMap_ = ROL::makePtr<Tpetra::Map<>>(
  //               Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
  //               myCellIds_, 0, comm_);

  /****************************************/
  /*** Assemble global graph structure. ***/
  /****************************************/
  matJ1Graph_ = ROL::makePtr<Tpetra::CrsGraph<>>(myUniqueStateMap_, 0);
  Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
  for (int i=0; i<numCells_; ++i) {
    for (int j=0; j<numLocalDofs; ++j) {
      matJ1Graph_->insertGlobalIndices(cellDofs(myCellIds_[i],j),
        cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs));
    }
  }
  matJ1Graph_->fillComplete();
  matJ2Graph_  = matJ1Graph_;
  matR1Graph_  = matJ1Graph_;
  matR2Graph_  = matJ2Graph_;
  matH11Graph_ = matJ1Graph_;
  matH12Graph_ = matJ1Graph_;
  matH21Graph_ = matJ2Graph_;
  matH22Graph_ = matJ2Graph_;

  if (verbose_ && myRank_==0) {
    outStream << "Initialized parallel structures." << std::endl;
  }
}

template<class Real>
void Assembler<Real>::setCellNodes(std::ostream &outStream) {
  // Build volume cell nodes
  shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology();
  int spaceDim = cellType.getDimension();
  int numNodesPerCell = cellType.getNodeCount();
  volCellNodes_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numNodesPerCell, spaceDim);
  Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
  Intrepid::FieldContainer<int>  &ctn   = *meshMgr_->getCellToNodeMap();
  for (int i=0; i<numCells_; ++i) {
    for (int j=0; j<numNodesPerCell; ++j) {
      for (int k=0; k<spaceDim; ++k) {
        (*volCellNodes_)(i, j, k) = nodes(ctn(myCellIds_[i],j), k);
      }
    }
  }
  // Build boundary cell nodes
  bdryCellIds_    = meshMgr_->getSideSets((verbose_ && myRank_==0), outStream);
  int numSideSets = bdryCellIds_->size();
  if (numSideSets > 0) {
    bdryCellNodes_.resize(numSideSets);
    bdryCellLocIds_.resize(numSideSets);
    for (int i=0; i<numSideSets; ++i) {
      int numLocSides = (*bdryCellIds_)[i].size();
      bdryCellNodes_[i].resize(numLocSides);
      bdryCellLocIds_[i].resize(numLocSides);
      for (int j=0; j<numLocSides; ++j) {
        if ((*bdryCellIds_)[i][j].size() > 0) {
          int numCellsSide = (*bdryCellIds_)[i][j].size();
          for (int k=0; k<numCellsSide; ++k) {
            int idx = mapGlobalToLocalCellId((*bdryCellIds_)[i][j][k]);
            if (idx > -1) {
              bdryCellLocIds_[i][j].push_back(idx);
              //if (myRank_==1) {std::cout << "\nrank " << myRank_ << "   bcid " << i << "  " << j << "  " << k << "  " << myCellIds_[idx] << "  " << idx;}
            }
          }
          int myNumCellsSide = bdryCellLocIds_[i][j].size();
          if (myNumCellsSide > 0) {
            bdryCellNodes_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(myNumCellsSide, numNodesPerCell, spaceDim);
          }
          for (int k=0; k<myNumCellsSide; ++k) {
            for (int l=0; l<numNodesPerCell; ++l) {
              for (int m=0; m<spaceDim; ++m) {
                (*bdryCellNodes_[i][j])(k, l, m) = nodes(ctn(myCellIds_[bdryCellLocIds_[i][j][k]],l), m);
              }
            }
          }
        }
        //if ((myRank_==1) && (myNumCellsSide > 0)) {std::cout << (*bdryCellNodes_[i][j]);}
      }
    }
  }
  bdryCellNodes_.resize(numSideSets);
}

template<class Real>
int Assembler<Real>::mapGlobalToLocalCellId(const int & gid) {
  int minId = cellOffsets_[myRank_];
  int maxId = cellOffsets_[myRank_]+numCells_-1;
  if ((gid >= minId) && (gid <= maxId)) {
    return (gid - cellOffsets_[myRank_]);
  }
  else {
    return -1;
  }
}

template<class Real>
void Assembler<Real>::getCoeffFromStateVector(ROL::Ptr<Intrepid::FieldContainer<Real>> &xcoeff,
                                        const ROL::Ptr<const Tpetra::MultiVector<>> &x) const {
  if ( x != ROL::nullPtr ) {
    // Perform import onto myOverlapMap
    ROL::Ptr<Tpetra::MultiVector<>> xshared =
      ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
    Tpetra::Import<> importer(myUniqueStateMap_, myOverlapStateMap_);
    xshared->doImport(*x,importer,Tpetra::REPLACE);
    // Populate xcoeff
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int lfs = dofMgr_->getLocalFieldSize();
    xcoeff = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);
    Teuchos::ArrayRCP<const Real> xdata = xshared->get1dView();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<lfs; ++j) {
        (*xcoeff)(i,j) = xdata[xshared->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
      }
    }
    dofMgr_->transformToIntrepidPattern(xcoeff);
  }
  else {
    xcoeff = ROL::nullPtr;
  }
}

template<class Real>
void Assembler<Real>::getCoeffFromControlVector(ROL::Ptr<Intrepid::FieldContainer<Real>> &xcoeff,
                                          const ROL::Ptr<const Tpetra::MultiVector<>> &x) const {
  if ( x != ROL::nullPtr ) {
    // Perform import onto myOverlapMap
    ROL::Ptr<Tpetra::MultiVector<>> xshared =
      ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
    Tpetra::Import<> importer(myUniqueControlMap_, myOverlapControlMap_);
    xshared->doImport(*x,importer,Tpetra::REPLACE);
    // Populate xcoeff
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int lfs = dofMgr_->getLocalFieldSize();
    xcoeff = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);
    Teuchos::ArrayRCP<const Real> xdata = xshared->get1dView();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<lfs; ++j) {
        (*xcoeff)(i,j) = xdata[xshared->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
      }
    }
    dofMgr_->transformToIntrepidPattern(xcoeff);
  }
  else {
    xcoeff = ROL::nullPtr;
  }
}

/***************************************************************************/
/****** GENERIC ASSEMBLY ROUTINES ******************************************/
/***************************************************************************/
template<class Real>
Real Assembler<Real>::assembleScalar(ROL::Ptr<Intrepid::FieldContainer<Real>> &val) {
  // Assembly
  if ( val != ROL::nullPtr ) {
    Real myval(0), gval(0);
    for (int i=0; i<numCells_; ++i) {
      myval += (*val)(i);
    }
    Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,1,&myval,&gval);
    return gval;
  }
  return static_cast<Real>(0);
}

template<class Real>
void Assembler<Real>::assembleFieldVector(ROL::Ptr<Tpetra::MultiVector<>> &v,
                                          ROL::Ptr<Intrepid::FieldContainer<Real>> &val,
                                          ROL::Ptr<Tpetra::MultiVector<>> &vecOverlap) {
  // Set residual vectors to zero
  v->scale(static_cast<Real>(0));
  vecOverlap->scale(static_cast<Real>(0));
  // Get degrees of freedom
  Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
  int numLocalDofs = cellDofs.dimension(1);
  // Transform values
  dofMgr_->transformToFieldPattern(val);
  // assembly on the overlap map
  for (int i=0; i<numCells_; ++i) {
    for (int j=0; j<numLocalDofs; ++j) {
      vecOverlap->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),0,(*val)(i,j));
    }
  }
  // change map
  Tpetra::Export<> exporter(vecOverlap->getMap(), v->getMap()); // redistribution
  v->doExport(*vecOverlap, exporter, Tpetra::ADD);              // from the overlap map to the unique map
}

template<class Real>
void Assembler<Real>::assembleParamVector(ROL::Ptr<std::vector<Real>> &v,
                                          std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &val) {
  int size = v->size();
  v->assign(size,0);
  for (int i = 0; i < size; ++i) {
    dofMgr_->transformToFieldPattern(val[i]);
  }
  // Assembly
  std::vector<Real> myVal(size,0);
  for (int j = 0; j < size; ++j) {
    if ( val[j] != ROL::nullPtr ) {
      for (int i=0; i<numCells_; ++i) {
        myVal[j] += (*val[j])(i);
      }
    }
  }
  Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,size,&myVal[0],&(*v)[0]);
}

template<class Real>
void Assembler<Real>::assembleFieldMatrix(ROL::Ptr<Tpetra::CrsMatrix<>> &M,
                                          ROL::Ptr<Intrepid::FieldContainer<Real>> &val) {
  // Transform data
  dofMgr_->transformToFieldPattern(val);
  // Zero PDE Jacobian
  M->resumeFill(); M->setAllToScalar(static_cast<Real>(0));
  // Assemble PDE Jacobian
  Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
  int numLocalDofs = cellDofs.dimension(1);
  int numLocalMatEntries = numLocalDofs * numLocalDofs;
  Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
  Teuchos::ArrayRCP<const Real> valArrayRCP = val->getData();
  for (int i=0; i<numCells_; ++i) {
    for (int j=0; j<numLocalDofs; ++j) {
      M->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                             cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                             valArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
    }
  }
  M->fillComplete();
}

template<class Real>
void Assembler<Real>::assembleParamFieldMatrix(ROL::Ptr<Tpetra::MultiVector<>> &M,
                                               std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &val,
                                               ROL::Ptr<Tpetra::MultiVector<>> &matOverlap) {
  // Initialize res
  int size = M->getNumVectors();
  // Compute PDE local Jacobian wrt parametric controls
  for (int i = 0; i < size; ++i) {
    dofMgr_->transformToFieldPattern(val[i]);
  }
  // Assemble PDE Jacobian wrt parametric controls
  M->scale(static_cast<Real>(0));
  matOverlap->scale(static_cast<Real>(0));
  for (int k = 0; k < size; ++k) {
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    // assembly on the overlap map
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matOverlap->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                        k,(*val[k])[i*numLocalDofs+j]);
      }
    }
    // change map
    Tpetra::Export<> exporter(matOverlap->getMap(), M->getMap()); // redistribution
    M->doExport(*matOverlap, exporter, Tpetra::ADD);              // from the overlap map to the unique map
  }
}

template<class Real>
void Assembler<Real>::assembleParamMatrix(ROL::Ptr<std::vector<std::vector<Real>>> &M,
                                          std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &val) {
  // Initialize local matrix
  int size = M->size();
  std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> tmp(size,ROL::nullPtr);
  // Compute local matrix
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      dofMgr_->transformToFieldPattern(val[i][j]);
    }
  }
  // Assemble PDE Jacobian wrt parametric controls
  int cnt = 0, matSize = static_cast<int>(0.5*static_cast<Real>((size+1)*size));
  std::vector<Real> myMat(matSize,static_cast<Real>(0));
  std::vector<Real> globMat(matSize,static_cast<Real>(0));
  for (int k = 0; k < size; ++k) {
    for (int j = k; j < size; ++j) {
      for (int i=0; i<numCells_; ++i) {
        myMat[cnt] += (*(val[k][j]))(i);
      }
      cnt++;
    }
  }
  Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,matSize,&myMat[0],&globMat[0]);
  cnt = 0;
  for (int k = 0; k < size; ++k) {
    for (int j = k; j < size; ++j) {
      (*M)[k][j] += globMat[cnt];
      if ( j != k ) { 
        (*M)[j][k] += globMat[cnt];
      }
      cnt++;
    }
  }
}

#endif
