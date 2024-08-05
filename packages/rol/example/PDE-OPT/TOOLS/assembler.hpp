// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  assembler.hpp
    \brief Finite element assembly class.
*/

#ifndef ROL_PDEOPT_ASSEMBLER_H
#define ROL_PDEOPT_ASSEMBLER_H

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "Intrepid_DefaultCubatureFactory.hpp"

#include "fe.hpp"
#include "fe_curl.hpp"
#include "pde.hpp"
#include "dynpde.hpp"
#include "qoi.hpp"
#include "dofmanager.hpp"
#include "meshmanager.hpp"
#include "fieldhelper.hpp"

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    ROL::Ptr<Teuchos::Time> AssemblePDEResidual         = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Residual");
    ROL::Ptr<Teuchos::Time> AssemblePDEJacobian1        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian1");
    ROL::Ptr<Teuchos::Time> AssemblePDEJacobian2        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian2");
    ROL::Ptr<Teuchos::Time> AssemblePDEJacobian3        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian3");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian11        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian11");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian12        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian12");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian13        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian13");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian21        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian21");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian22        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian22");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian23        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian23");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian31        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian31");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian32        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian32");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian33        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian33");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEResidual      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Residual");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEJacobian_uo   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Jacobian_uo");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEJacobian_un   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Jacobian_un");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEJacobian_zf   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Jacobian_zf");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEJacobian_zp   = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Jacobian_zp");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_uo_uo = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_uo_uo");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_uo_un = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_uo_un");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_uo_zf = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_uo_zf");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_uo_zp = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_uo_zp");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_un_uo = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_un_uo");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_un_un = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_un_un");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_un_zf = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_un_zf");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_un_zp = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_un_zp");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_zf_uo = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_zf_uo");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_zf_un = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_zf_un");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_zf_zf = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_zf_zf");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_zf_zp = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_zf_zp");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_zp_uo = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_zp_uo");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_zp_un = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_zp_un");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_zp_zf = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_zp_zf");
    ROL::Ptr<Teuchos::Time> AssembleDynPDEHessian_zp_zp = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble DynamicPDE Hessian_zp_zp");
    ROL::Ptr<Teuchos::Time> AssembleQOIValue            = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Value");
    ROL::Ptr<Teuchos::Time> AssembleQOIGradient1        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient1");
    ROL::Ptr<Teuchos::Time> AssembleQOIGradient2        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient2");
    ROL::Ptr<Teuchos::Time> AssembleQOIGradient3        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient3");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec11        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec11");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec12        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec12");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec13        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec13");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec21        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec21");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec22        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec22");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec23        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec23");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec31        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec31");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec32        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec32");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec33        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec33");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessian11        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Hessian11");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessian12        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Hessian12");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessian21        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Hessian21");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessian22        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Hessian22");
  }
}
#endif

template<class Real>
class Solution {
public:
  virtual ~Solution() {}
  Solution(void) {}
  virtual Real evaluate(const std::vector<Real> &x, const int fieldNumber) const = 0;
};

template<class Real>
class Assembler {

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NO;
  typedef Tpetra::MultiVector<Real,LO,GO,NO> MV;
  typedef Tpetra::Operator<Real,LO,GO,NO> OP;

private:
  // Timers
//  ROL::Ptr<Teuchos::Time::Time> timerSolverFactorization_;
//  ROL::Ptr<Teuchos::Time::Time> timerSolverSubstitution_;
//  ROL::Ptr<Teuchos::Time::Time> timerAssemblyNonlinear_;
//  ROL::Ptr<Teuchos::Time::Time> timerSolverUpdate_;

  // Set in Constructor.
  bool verbose_;
  bool isJ1Transposed_, isJ2Transposed_, isJuoTransposed_, isJunTransposed_, isJzfTransposed_;

  // Set in SetCommunicator.
  ROL::Ptr<const Teuchos::Comm<int>> comm_;
  int myRank_, numProcs_;

  // Set in SetBasis.
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> basisPtrs1_, basisPtrs2_;

  // Set in SetDiscretization.
  ROL::Ptr<MeshManager<Real>> meshMgr_;
  ROL::Ptr<DofManager<Real>>  dofMgr1_, dofMgr2_;

  // Set in SetParallelStructure.
  int numCells_;
  Teuchos::Array<GO> myCellIds_;
  Teuchos::Array<GO> cellOffsets_;
  ROL::Ptr<const Tpetra::Map<>> myOverlapStateMap_;
  ROL::Ptr<const Tpetra::Map<>> myUniqueStateMap_;
  ROL::Ptr<const Tpetra::Map<>> myOverlapControlMap_;
  ROL::Ptr<const Tpetra::Map<>> myUniqueControlMap_;
  ROL::Ptr<const Tpetra::Map<>> myOverlapResidualMap_;
  ROL::Ptr<const Tpetra::Map<>> myUniqueResidualMap_;
  ROL::Ptr<Tpetra::CrsGraph<>>  matJ1Graph_;
  ROL::Ptr<Tpetra::CrsGraph<>>  matJ2Graph_;
  ROL::Ptr<Tpetra::CrsGraph<>>  matR1Graph_;
  ROL::Ptr<Tpetra::CrsGraph<>>  matR2Graph_;
  ROL::Ptr<Tpetra::CrsGraph<>>  matH11Graph_;
  ROL::Ptr<Tpetra::CrsGraph<>>  matH12Graph_;
  ROL::Ptr<Tpetra::CrsGraph<>>  matH21Graph_;
  ROL::Ptr<Tpetra::CrsGraph<>>  matH22Graph_;

  // Set in SetCellNodes.
  ROL::Ptr<Intrepid::FieldContainer<Real>> volCellNodes_;
  ROL::Ptr<std::vector<std::vector<std::vector<int>>>>  bdryCellIds_;
  //std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<int>>>>  bdryCellLocIds_;
  std::vector<std::vector<std::vector<int>>>  bdryCellLocIds_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> bdryCellNodes_;

  // Finite element vectors and matrices for PDE.
  ROL::Ptr<Tpetra::MultiVector<>> pde_vecR_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> pde_vecJ2_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> pde_vecJ3_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> pde_vecH13_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> pde_vecH23_overlap_;

  // Finite element vectors and matrices for QoI.
  ROL::Ptr<Tpetra::MultiVector<>> qoi_vecG1_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> qoi_vecG2_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> qoi_vecH11_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> qoi_vecH12_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> qoi_vecH13_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> qoi_vecH21_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> qoi_vecH22_overlap_;
  ROL::Ptr<Tpetra::MultiVector<>> qoi_vecH23_overlap_;

  bool store_;

private:

  void setCommunicator(const ROL::Ptr<const Teuchos::Comm<int>> &comm,
                       Teuchos::ParameterList &parlist,
                       std::ostream &outStream = std::cout);
  void setBasis(
         const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs1,
         const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs2,
         Teuchos::ParameterList &parlist,
         std::ostream &outStream = std::cout);
  void setDiscretization(Teuchos::ParameterList &parlist,
                         const ROL::Ptr<MeshManager<Real>> &meshMgr = ROL::nullPtr,
                         std::ostream &outStream = std::cout);
  void setParallelStructure(Teuchos::ParameterList &parlist,
                            std::ostream &outStream = std::cout);
  void setCellNodes(std::ostream &outStream = std::cout);
  int mapGlobalToLocalCellId(const int & gid);
  void getCoeffFromStateVector(ROL::Ptr<Intrepid::FieldContainer<Real>> &xcoeff,
                               const ROL::Ptr<const Tpetra::MultiVector<>> &x) const;
  void getCoeffFromControlVector(ROL::Ptr<Intrepid::FieldContainer<Real>> &xcoeff,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &x) const;

  /***************************************************************************/
  /****** GENERIC ASSEMBLY ROUTINES ******************************************/
  /***************************************************************************/
  Real assembleScalar(ROL::Ptr<Intrepid::FieldContainer<Real>> &val);
  void assembleFieldVector(ROL::Ptr<Tpetra::MultiVector<>> &v,
                           ROL::Ptr<Intrepid::FieldContainer<Real>> &val,
                           ROL::Ptr<Tpetra::MultiVector<>> &vecOverlap,
                           const ROL::Ptr<DofManager<Real>> &dofMgr);
  void assembleParamVector(ROL::Ptr<std::vector<Real>> &v,
                           std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &val);
  void assembleFieldMatrix(ROL::Ptr<Tpetra::CrsMatrix<>> &M,
                           ROL::Ptr<Intrepid::FieldContainer<Real>> &val,
                           const ROL::Ptr<DofManager<Real>> &dofMgr1,
                           const ROL::Ptr<DofManager<Real>> &dofMgr2);
  void assembleParamFieldMatrix(ROL::Ptr<Tpetra::MultiVector<>> &M,
                                std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>> &val,
                                ROL::Ptr<Tpetra::MultiVector<>> &matOverlap,
                                const ROL::Ptr<DofManager<Real>> &dofMgr);
  void assembleParamMatrix(ROL::Ptr<std::vector<std::vector<Real>>> &M,
                           std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real>>>> &val,
                           const ROL::Ptr<DofManager<Real>> &dofMgr);
  void transformToFieldPattern(const ROL::Ptr<Intrepid::FieldContainer<Real>> &array,
                               const ROL::Ptr<DofManager<Real>> &dofMgr1,
                               const ROL::Ptr<DofManager<Real>> &dofMgr2 = ROL::nullPtr) const;

public:
  // destructor
  virtual ~Assembler() {}
  // Constuctor: Discretization set from ParameterList
  Assembler(const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs,
          const ROL::Ptr<const Teuchos::Comm<int>> &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout);
  // Constructor: Discretization set from MeshManager input
  Assembler(const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs,
          const ROL::Ptr<MeshManager<Real>> &meshMgr,
          const ROL::Ptr<const Teuchos::Comm<int>> &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout);
  // Constuctor: Discretization set from ParameterList
  Assembler(const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs1,
          const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs2,
          const ROL::Ptr<const Teuchos::Comm<int>> &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout);
  // Constructor: Discretization set from MeshManager input
  Assembler(const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs1,
          const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real>>>> &basisPtrs2,
          const ROL::Ptr<MeshManager<Real>> &meshMgr,
          const ROL::Ptr<const Teuchos::Comm<int>> &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout);
  void setCellNodes(PDE<Real> &pde) const;
  void setCellNodes(DynamicPDE<Real> &pde) const;

  /***************************************************************************/
  /* PDE assembly routines                                                   */
  /***************************************************************************/
  void assemblePDEResidual(ROL::Ptr<Tpetra::MultiVector<>> &r,
                           const ROL::Ptr<PDE<Real>> &pde,
                           const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                           const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                           const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEJacobian1(ROL::Ptr<Tpetra::CrsMatrix<>> &J1,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEJacobian2(ROL::Ptr<Tpetra::CrsMatrix<>> &J2,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEJacobian3(ROL::Ptr<Tpetra::MultiVector<>> &J3,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEapplyJacobian1(ROL::Ptr<Tpetra::MultiVector<>> &Jv1,
                                 const ROL::Ptr<PDE<Real>> &pde,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEapplyAdjointJacobian1(ROL::Ptr<Tpetra::MultiVector<>> &Jv1,
                                        const ROL::Ptr<PDE<Real>> &pde,
                                        const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                        const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                        const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                        const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEapplyJacobian2(ROL::Ptr<Tpetra::MultiVector<>> &Jv2,
                                 const ROL::Ptr<PDE<Real>> &pde,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEapplyAdjointJacobian2(ROL::Ptr<Tpetra::MultiVector<>> &Jv2,
                                        const ROL::Ptr<PDE<Real>> &pde,
                                        const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                        const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                        const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                        const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian11(ROL::Ptr<Tpetra::CrsMatrix<>> &H11,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian12(ROL::Ptr<Tpetra::CrsMatrix<>> &H12,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian13(ROL::Ptr<Tpetra::MultiVector<>> &H13,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian21(ROL::Ptr<Tpetra::CrsMatrix<>> &H21,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian22(ROL::Ptr<Tpetra::CrsMatrix<>> &H22,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian23(ROL::Ptr<Tpetra::MultiVector<>> &H23,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian31(ROL::Ptr<Tpetra::MultiVector<>> &H31,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian32(ROL::Ptr<Tpetra::MultiVector<>> &H32,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEHessian33(ROL::Ptr<std::vector<std::vector<Real>>> &H33,
                            const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEapplyHessian11(ROL::Ptr<Tpetra::MultiVector<>> &Hv,
                                 const ROL::Ptr<PDE<Real>> &pde,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEapplyHessian12(ROL::Ptr<Tpetra::MultiVector<>> &Hv,
                                 const ROL::Ptr<PDE<Real>> &pde,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEapplyHessian21(ROL::Ptr<Tpetra::MultiVector<>> &Hv,
                                 const ROL::Ptr<PDE<Real>> &pde,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assemblePDEapplyHessian22(ROL::Ptr<Tpetra::MultiVector<>> &Hv,
                                 const ROL::Ptr<PDE<Real>> &pde,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  /***************************************************************************/
  /* End of PDE assembly routines.                                           */
  /***************************************************************************/

  /***************************************************************************/
  /* Dynamic PDE assembly routines                                           */
  /***************************************************************************/
  void assembleDynPDEResidual(ROL::Ptr<Tpetra::MultiVector<>> &r,
                              const ROL::Ptr<DynamicPDE<Real>> &pde,
                              const ROL::TimeStamp<Real> &ts,
                              const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                              const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                              const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                              const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEJacobian_uo(ROL::Ptr<Tpetra::CrsMatrix<>> &J,
                                 const ROL::Ptr<DynamicPDE<Real>> &pde,
                                 const ROL::TimeStamp<Real> &ts,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEJacobian_un(ROL::Ptr<Tpetra::CrsMatrix<>> &J,
                                 const ROL::Ptr<DynamicPDE<Real>> &pde,
                                 const ROL::TimeStamp<Real> &ts,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEJacobian_zf(ROL::Ptr<Tpetra::CrsMatrix<>> &J,
                                 const ROL::Ptr<DynamicPDE<Real>> &pde,
                                 const ROL::TimeStamp<Real> &ts,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEJacobian_zp(ROL::Ptr<Tpetra::MultiVector<>> &J,
                                 const ROL::Ptr<DynamicPDE<Real>> &pde,
                                 const ROL::TimeStamp<Real> &ts,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                 const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                 const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_uo_uo(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_uo_un(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_uo_zf(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_uo_zp(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_un_uo(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_un_un(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_un_zf(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_un_zp(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_zf_uo(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_zf_un(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_zf_zf(ROL::Ptr<Tpetra::CrsMatrix<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_zf_zp(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_zp_uo(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_zp_un(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_zp_zf(ROL::Ptr<Tpetra::MultiVector<>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleDynPDEHessian_zp_zp(ROL::Ptr<std::vector<std::vector<Real>>> &H,
                                   const ROL::Ptr<DynamicPDE<Real>> &pde,
                                   const ROL::TimeStamp<Real> &ts,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &l,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &uo,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &un,
                                   const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                                   const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  /***************************************************************************/
  /* End of DynamicPDE assembly routines.                                    */
  /***************************************************************************/

  /***************************************************************************/
  /* QoI assembly routines                                                   */
  /***************************************************************************/
  Real assembleQoIValue(const ROL::Ptr<QoI<Real>> &qoi,
                        const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                        const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                        const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIGradient1(ROL::Ptr<Tpetra::MultiVector<>> &g1,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIGradient2(ROL::Ptr<Tpetra::MultiVector<>> &g2,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIGradient3(ROL::Ptr<std::vector<Real>> &g3,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessVec11(ROL::Ptr<Tpetra::MultiVector<>> &H11,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessVec12(ROL::Ptr<Tpetra::MultiVector<>> &H12,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessVec13(ROL::Ptr<Tpetra::MultiVector<>> &H13,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const std::vector<Real>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleQoIHessVec21(ROL::Ptr<Tpetra::MultiVector<>> &H21,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessVec22(ROL::Ptr<Tpetra::MultiVector<>> &H22,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessVec23(ROL::Ptr<Tpetra::MultiVector<>> &H23,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const std::vector<Real>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessVec31(ROL::Ptr<std::vector<Real>> &H31,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessVec32(ROL::Ptr<std::vector<Real>> &H32,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessVec33(ROL::Ptr<std::vector<Real>> &H33,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const std::vector<Real>> &v,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr);
  void assembleQoIHessian11(ROL::Ptr<Tpetra::CrsMatrix<>> &H11,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessian12(ROL::Ptr<Tpetra::CrsMatrix<>> &H12,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessian21(ROL::Ptr<Tpetra::CrsMatrix<>> &H21,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  void assembleQoIHessian22(ROL::Ptr<Tpetra::CrsMatrix<>> &H22,
                            const ROL::Ptr<QoI<Real>> &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> & z_param = ROL::nullPtr);
  /***************************************************************************/
  /* End QoI assembly routines                                               */
  /***************************************************************************/

  /***************************************************************************/
  /* Assemble and apply Riesz operator corresponding to simulation variables */
  /***************************************************************************/
  void assemblePDERieszMap1(ROL::Ptr<Tpetra::CrsMatrix<>> &R1,
                            const ROL::Ptr<PDE<Real>> &pde);
  void assembleDynPDERieszMap1(ROL::Ptr<Tpetra::CrsMatrix<>> &R1,
                               const ROL::Ptr<DynamicPDE<Real>> &pde);
  /***************************************************************************/
  /* End of functions for Riesz operator of simulation variables.            */
  /***************************************************************************/

  /***************************************************************************/
  /* Assemble and apply Riesz operator corresponding to optimization         */
  /* variables                                                               */
  /***************************************************************************/
  void assemblePDERieszMap2(ROL::Ptr<Tpetra::CrsMatrix<>> &R2,
                            const ROL::Ptr<PDE<Real>> &pde);
  void assembleDynPDERieszMap2(ROL::Ptr<Tpetra::CrsMatrix<>> &R2,
                               const ROL::Ptr<DynamicPDE<Real>> &pde);
  /***************************************************************************/
  /* End of functions for Riesz operator of optimization variables.          */
  /***************************************************************************/

  /***************************************************************************/
  /* Compute error routines.                                                 */
  /***************************************************************************/
  Real computeStateError(const ROL::Ptr<const Tpetra::MultiVector<>> &soln,
                         const ROL::Ptr<Solution<Real>> &trueSoln,
                         const int cubDeg = 6,
                         const ROL::Ptr<FieldHelper<Real>> &fieldHelper = ROL::nullPtr) const;
  Real computeControlError(const ROL::Ptr<const Tpetra::MultiVector<>> &soln,
                           const ROL::Ptr<Solution<Real>> &trueSoln,
                           const int cubDeg = 6,
                           const ROL::Ptr<FieldHelper<Real>> &fieldHelper = ROL::nullPtr) const;
  /***************************************************************************/
  /* End of compute solution routines.                                       */
  /***************************************************************************/

  /***************************************************************************/
  /* Output routines.                                                        */
  /***************************************************************************/
  void printMeshData(std::ostream &outStream) const;
  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<>> &vec,
                          const std::string &filename) const;
  void inputTpetraVector(ROL::Ptr<Tpetra::MultiVector<>> &vec,
                         const std::string &filename) const;
  void printDataPDE(const ROL::Ptr<PDE<Real>> &pde,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                    const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                    const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) const;
  void printCellAveragesPDE(const ROL::Ptr<PDE<Real>> &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                            const ROL::Ptr<const Tpetra::MultiVector<>> &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real>> &z_param = ROL::nullPtr) const;
  void serialPrintStateEdgeField(const ROL::Ptr<const Tpetra::MultiVector<>> &u,
                                 const ROL::Ptr<FieldHelper<Real>> &fieldHelper,
                                 const std::string &filename,
                                 const ROL::Ptr<FE_CURL<Real>> &fe) const;
  /***************************************************************************/
  /* End of output routines.                                                 */
  /***************************************************************************/

  /***************************************************************************/
  /* Vector generation routines.                                             */
  /***************************************************************************/
  const ROL::Ptr<const Tpetra::Map<>> getStateMap(void) const;
  const ROL::Ptr<const Tpetra::Map<>> getControlMap(void) const;
  const ROL::Ptr<const Tpetra::Map<>> getResidualMap(void) const;
  ROL::Ptr<Tpetra::MultiVector<>> createStateVector(void) const;
  ROL::Ptr<Tpetra::MultiVector<>> createControlVector(void) const;
  ROL::Ptr<Tpetra::MultiVector<>> createResidualVector(void) const;
  /***************************************************************************/
  /* End of vector generation routines.                                      */
  /***************************************************************************/

  /***************************************************************************/
  /* Accessor routines.                                                      */
  /***************************************************************************/
  const ROL::Ptr<DofManager<Real>> getDofManager(void) const;
  const ROL::Ptr<DofManager<Real>> getDofManager2(void) const;
  Teuchos::Array<GO> getCellIds(void) const;
  /***************************************************************************/
  /* End of accessor routines.                                               */
  /***************************************************************************/
}; // class Assembler

#include "assembler_def.hpp"

#endif
