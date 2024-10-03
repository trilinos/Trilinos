// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MiniEM_FullDarcyPreconditionerFactory.hpp"

#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

#include "Teko_SolveInverseFactory.hpp"

#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"

#include "Teuchos_as.hpp"
#include "Teuchos_Time.hpp"

#include "Teko_TpetraHelpers.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Panzer_NodeType.hpp"
#include "PanzerDiscFE_config.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraThyraWrappers.hpp"
#include "ml_epetra_utils.h"
#include "EpetraExt_MatrixMatrix.h"
#endif
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"

#include "Panzer_String_Utilities.hpp"

#include "MiniEM_Utils.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {


///////////////////////////////////////
// FullDarcyPreconditionerFactory  //
///////////////////////////////////////

Teko::LinearOp FullDarcyPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & /* state */) const
{
   typedef double Scalar;
   typedef int LocalOrdinal;
   typedef panzer::GlobalOrdinal GlobalOrdinal;
   typedef panzer::TpetraNodeType Node;

   Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("DarcyPreconditioner::build")));

   // Output stream for debug information
   RCP<Teuchos::FancyOStream> debug = Teuchos::null;
   if (doDebug)
     debug = Teko::getOutputStream();

   // Check that system is right size
   int rows = Teko::blockRowCount(blo);
   int cols = Teko::blockColCount(blo);
   TEUCHOS_ASSERT(rows==cols);
   TEUCHOS_ASSERT(rows==2);

   // Notation:
   // 0 - Hgrad
   // 1 - Hcurl
   // 2 - HDiv
   // 3 - HVol

   // M_k(a) - mass matrix on space k with weight a
   // D_k - derivative from space k to k+1

   // The block matrix is
   //
   // | Q_u  K   |
   // | Kt   Q_sigma |
   //
   // where
   // Q_u = 1/dt * M_3(1)
   // K   = -M_3(1) * D_2
   // Kt  = D_2^T * M_3(1)
   // Q_sigma = M_2(1/kappa)

   // S_sigma = Q_sigma - Kt * Q_u^-1 * K
   //     = M_2(1/kappa) + dt * D_2^T * M_3(1) * D_2
   //
   // -> grad-div term scales like dt
   //
   // for refmaxwell: Q_rho = M_1(1/dt) so that the addon is:
   // M_2(1) * D_1 * M_1(1/dt)^-1 * D_1^T * M_2(1)

   // Modify the system
   if (simplifyFaraday_) {
     RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();
     *out << std::endl;
     *out << "*** WARNING ***" << std::endl;
     *out << "We are modifying the linear system. That's not a friendly thing to do." << std::endl;
     *out << std::endl;

     Teko::LinearOp Q_u  = Teko::getBlock(0,0,blo);
     Teko::LinearOp id_u = getIdentityMatrix(Q_u, 1/dt);
     Teko::LinearOp hoDiv  = Teko::scale(-1.0, getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Div")));
     Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);
     Teko::LinearOp Q_sigma   = Teko::getBlock(1,1,blo);
     blo->beginBlockFill(2,2);
     Teko::setBlock(0,0,blo,id_u);
     Teko::setBlock(0,1,blo,hoDiv);
     Teko::setBlock(1,0,blo,Kt);
     Teko::setBlock(1,1,blo,Q_sigma);
     blo->endBlockFill();
   }

   // Extract the blocks
   Teko::LinearOp Q_u   = Teko::getBlock(0,0,blo);
   Teko::LinearOp K     = Teko::getBlock(0,1,blo);
   Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);
   Teko::LinearOp Q_sigma   = Teko::getBlock(1,1,blo);

   // discrete curl and its transpose
   Teko::LinearOp Div, DivT;
   if (use_discrete_div_) {
     Div = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Div"));
     DivT = Teko::explicitTranspose(Div);
   }

   // Set up the Schur complement
   Teko::LinearOp S_sigma;
   {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Schur complement"));
     S_sigma = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("DarcySchurComplement AUXILIARY_FACE"));
   }

   // Check whether we are using Tpetra or Epetra
   RCP<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > checkTpetra = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(Q_sigma);
   bool useTpetra = nonnull(checkTpetra);

   /////////////////////////////////////////////////
   // Debug and matrix dumps                      //
   /////////////////////////////////////////////////

   if (dump) {
     if (Div != Teuchos::null) {
       Teko::LinearOp K2 = Teko::explicitMultiply(Q_u, Div);
       Teko::LinearOp diffK;

       if (useTpetra) {
         RCP<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(K2,true);
         RCP<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp2 = Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(tOp);
         RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tOp2->getTpetraOperator(),true);
         crsOp->scale(dt);
         diffK = Teko::explicitAdd(K, Teko::scale(1.0,K2));

       } else {
         diffK = Teko::explicitAdd(K, Teko::scale(dt,K2));
       }

       TEUCHOS_ASSERT(Teko::infNorm(diffK) < 1.0e-8 * Teko::infNorm(K));
     }
   }
   describeAndWriteMatrix("Q_u",*Q_u,debug,dump);
   describeAndWriteMatrix("K",*K,debug,dump);
   describeAndWriteMatrix("Kt",*Kt,debug,dump);
   describeAndWriteMatrix("Q_sigma",*Q_sigma,debug,dump);
   if (Div != Teuchos::null)
     describeAndWriteMatrix("Div",*Div,debug,dump);
   describeAndWriteMatrix("S_sigma",*S_sigma,debug,dump);


   /////////////////////////////////////////////////
   // Set up inverses for sub-blocks              //
   /////////////////////////////////////////////////

   // Inverse of B mass matrix
   Teko::LinearOp invQ_u;
   if (!simplifyFaraday_) {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Inverse Q_u"));
     // Are we building a solver or a preconditioner?
     if (useAsPreconditioner) {
       invQ_u = Teko::buildInverse(*invLib.getInverseFactory("Q_u Preconditioner"),Q_u);
     } else {
       Teko::LinearOp invDiagQ_u = Teko::buildInverse(*invLib.getInverseFactory("Q_u Preconditioner"),Q_u);
       describeMatrix("invDiagQ_u",*invDiagQ_u,debug);
       invQ_u = Teko::buildInverse(*invLib.getInverseFactory("Q_u Solve"),Q_u, invDiagQ_u);
     }
   }

   // Solver for S_sigma
   Teko::LinearOp invS_sigma;
   {
     Teuchos::TimeMonitor tm1(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Solver S_sigma"));

     if ((S_sigma_prec_type_ == "MueLuRefDarcy") || (S_sigma_prec_type_ == "MueLuRefMaxwell") || (S_sigma_prec_type_ == "MueLu")) {// refDarcy

       RCP<Teko::InverseFactory> S_sigma_prec_factory;
       Teuchos::ParameterList S_sigma_prec_pl;
       S_sigma_prec_factory = invLib.getInverseFactory("S_sigma Preconditioner");
       S_sigma_prec_pl = *S_sigma_prec_factory->getParameterList();
       Teuchos::ParameterList& muelulist = S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_);

       // interpolations between HCurl spaces of different orders
       std::vector<Teko::LinearOp> interpolationsHCurl;
       // interpolations between HDiv spaces of different orders
       std::vector<Teko::LinearOp> interpolationsHDiv;
       // discrete curls HCurl -> HDiv
       std::vector<Teko::LinearOp> discreteCurls;
       // discrete curls HCurl -> HDiv
       std::vector<Teko::LinearOp> discreteGradients;
       // Schur complements
       std::vector<Teko::LinearOp> schurComplements;
       // projection of Schur complements into HCurl using discrete curls
       std::vector<Teko::LinearOp> projectedSchurComplements;
       if (pCoarsenSchedule_.size() > 0) {
         {
           using assembleType = std::pair<std::string, std::vector<Teko::LinearOp>& >;

           std::vector<assembleType> assemble = {{"Discrete Curl",                                discreteCurls},
                                                 {"DarcySchurComplement AUXILIARY_FACE",          schurComplements},
                                                 {"ProjectedDarcySchurComplement AUXILIARY_EDGE", projectedSchurComplements}};

           for (auto ait = assemble.begin(); ait != assemble.end(); ++ait) {
             std::string operatorName = ait->first;
             std::vector<Teko::LinearOp>& operatorVector = ait->second;
             operatorVector.push_back(getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg(operatorName)));
             describeAndWriteMatrix(operatorName, *operatorVector.back(), debug, dump);
           }

           auto it = pCoarsenSchedule_.begin();
           int p = *it;
           ++it;
           while (it != pCoarsenSchedule_.end()) {
             int q = *it;

             assemble = {{"Discrete Curl "+std::to_string(q),                                 discreteCurls},
                         {"DarcySchurComplement AUXILIARY_FACE_"+std::to_string(q),           schurComplements},
                         {"ProjectedDarcySchurComplement AUXILIARY_EDGE_"+std::to_string(q),  projectedSchurComplements},
                         {"Interpolation Hdiv "+std::to_string(q) + "->" + std::to_string(p), interpolationsHDiv}};

             for (auto ait = assemble.begin(); ait != assemble.end(); ++ait) {
               std::string operatorName = ait->first;
               std::vector<Teko::LinearOp>& operatorVector = ait->second;
               operatorVector.push_back(getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg(operatorName)));
               describeAndWriteMatrix(operatorName, *operatorVector.back(), debug, dump);
             }

             p = q;
             ++it;
           }
         }

       }

       RCP<Teuchos::ParameterList> refDarcyList;
       if (pCoarsenSchedule_.size() > 0) {
         int maxLevels = interpolationsHDiv.size()+1;

         // Make sure MueLu only creates levels for the operators that we pass in
         muelulist.set("max levels", maxLevels);

         bool implicitTranspose = muelulist.get("transpose: use implicit", false);

         // Cannot explicitly transpose matrix-free operators
         if (!implicitTranspose && isMatrixFreeOperator(interpolationsHDiv[0])) {
           implicitTranspose = true;
           muelulist.set("transpose: use implicit", true);
         }
         std::string smootherType = muelulist.get("smoother: type", "HIPTMAIR");
         std::transform(smootherType.begin(), smootherType.end(), smootherType.begin(), ::toupper);
         if (smootherType == "HIPTMAIR") {
           if (isMatrixFreeOperator(interpolationsHDiv[0]))
             muelulist.sublist("smoother: params").set("hiptmair: implicit transpose", true);
           TEUCHOS_ASSERT_EQUALITY(muelulist.get<bool>("transpose: use implicit"),
                                   muelulist.sublist("smoother: params").get<bool>("hiptmair: implicit transpose"));
         }

         std::vector<Teko::LinearOp> interpolationsHDivT;
         if (!implicitTranspose)
           // Get restriction operators
           for (int lvl = 1; lvl < maxLevels; ++lvl) {
             Teko::LinearOp interpT = Teko::explicitTranspose(interpolationsHDiv[lvl-1]);
             interpolationsHDivT.push_back(interpT);
           }

         for (int lvl = 0; lvl < maxLevels; ++lvl) {
           Teuchos::ParameterList& lvlList = muelulist.sublist("level " + std::to_string(lvl) + " user data");
           if (lvl > 0) {
             lvlList.set("A",schurComplements[lvl]);
             lvlList.set("P",interpolationsHDiv[lvl-1]);
             if (!implicitTranspose)
               lvlList.set("R",interpolationsHDivT[lvl-1]);
           }
           // Operators for Hiptmair smoothing
           lvlList.set("NodeMatrix",projectedSchurComplements[lvl]);
           lvlList.set("D0",discreteCurls[lvl]);
         }

         refDarcyList = Teuchos::rcpFromRef(muelulist.sublist("level " + std::to_string(maxLevels-1) + " user data"));
       } else {
         refDarcyList = Teuchos::rcpFromRef(muelulist);
       }

       std::string postFix = "";
       std::string postFixStrong = "";
       if (pCoarsenSchedule_.size() > 0) {
         postFix = "_1";
         postFixStrong = " 1";
       }

       // get auxiliary matrices for RefDarcy
       Teko::LinearOp M1_beta = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix dt weighted AUXILIARY_EDGE"+postFix));
       Teko::LinearOp M1_alpha = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix 1/kappa weighted AUXILIARY_EDGE"+postFix));
       Teko::LinearOp Mk_one = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_FACE"+postFix));
       Teko::LinearOp Mk_1_one = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_EDGE"+postFix));
       Teko::LinearOp Curl = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"+postFixStrong));
       Teko::LinearOp Grad = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"+postFixStrong));
       Teko::LinearOp Mk_1_invBeta = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix 1/dt weighted AUXILIARY_EDGE"+postFix));
       Teko::LinearOp Mk_2_invAlpha = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix kappa weighted AUXILIARY_NODE"+postFix));

       describeAndWriteMatrix("DiscreteGradient",*Grad,debug,dump);
       describeAndWriteMatrix("DiscreteCurl",*Curl,debug,dump);
       describeAndWriteMatrix("Mk_1_invBeta",*Mk_1_invBeta,debug,dump);
       describeAndWriteMatrix("Mk_2_invAlpha",*Mk_2_invAlpha,debug,dump);

       RCP<const Thyra::DiagonalLinearOpBase<Scalar> > invMk_1_invBeta;
       {
         Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Lumped diagonal Mk_1_invBeta"));

         // Get inverse of lumped Mk_1_invBeta
         RCP<Thyra::VectorBase<Scalar> > ones = Thyra::createMember(Mk_1_invBeta->domain());
         RCP<Thyra::VectorBase<Scalar> > diagonal = Thyra::createMember(Mk_1_invBeta->range());
         Thyra::assign(ones.ptr(),1.0);
         // compute lumped diagonal
         Thyra::apply(*Mk_1_invBeta,Thyra::NOTRANS,*ones,diagonal.ptr());
         Thyra::reciprocal(*diagonal,diagonal.ptr());
         invMk_1_invBeta = rcp(new Thyra::DefaultDiagonalLinearOp<Scalar>(diagonal));
       }

       RCP<const Thyra::DiagonalLinearOpBase<Scalar> > invMk_2_invAlpha;
       {
         Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Lumped diagonal Mk_2_invAlpha"));

         // Get inverse of lumped Mk_2_invAlpha
         RCP<Thyra::VectorBase<Scalar> > ones = Thyra::createMember(Mk_2_invAlpha->domain());
         RCP<Thyra::VectorBase<Scalar> > diagonal = Thyra::createMember(Mk_2_invAlpha->range());
         Thyra::assign(ones.ptr(),1.0);
         // compute lumped diagonal
         Thyra::apply(*Mk_2_invAlpha,Thyra::NOTRANS,*ones,diagonal.ptr());
         Thyra::reciprocal(*diagonal,diagonal.ptr());
         invMk_2_invAlpha = rcp(new Thyra::DefaultDiagonalLinearOp<Scalar>(diagonal));
       }

       refDarcyList->set("Dk_1",Curl);
       refDarcyList->set("Dk_2",Grad);
       refDarcyList->set("D0",Grad);

       refDarcyList->set("M1_beta",M1_beta);
       refDarcyList->set("M1_alpha",M1_alpha);

       refDarcyList->set("Mk_one",Mk_one);
       refDarcyList->set("Mk_1_one",Mk_1_one);

       refDarcyList->set("invMk_1_invBeta",invMk_1_invBeta);
       refDarcyList->set("invMk_2_invAlpha",invMk_2_invAlpha);

       // Get coordinates
       {
         if (useTpetra) {
           RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Coordinates = S_sigma_prec_pl.get<RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >("Coordinates");
           refDarcyList->set("Coordinates",Coordinates);
         }
#ifdef PANZER_HAVE_EPETRA_STACK
         else {
           RCP<Epetra_MultiVector> Coordinates = S_sigma_prec_pl.get<RCP<Epetra_MultiVector> >("Coordinates");
           refDarcyList->set("Coordinates",Coordinates);
         }
#endif // PANZER_HAVE_EPETRA_STACK
       }

       {
         Teko::InverseLibrary myInvLib = invLib;

         muelulist.set("Type",S_sigma_prec_type_);
         myInvLib.addInverse("S_sigma Preconditioner",muelulist);
         S_sigma_prec_factory = myInvLib.getInverseFactory("S_sigma Preconditioner");
       }

       // Are we building a solver or a preconditioner?
       {
         Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Build S_sigma preconditioner"));

         if (useAsPreconditioner)
           invS_sigma = Teko::buildInverse(*S_sigma_prec_factory,S_sigma);
         else {
           if (S_sigma_prec_.is_null())
             S_sigma_prec_ = Teko::buildInverse(*S_sigma_prec_factory,S_sigma);
           else
             Teko::rebuildInverse(*S_sigma_prec_factory,S_sigma, S_sigma_prec_);
           invS_sigma = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Solve"),S_sigma,S_sigma_prec_);
         }
       }
     }
#ifdef PANZER_HAVE_EPETRA_STACK
     else if (S_sigma_prec_type_ == "ML") {
       RCP<Teko::InverseFactory> S_sigma_prec_factory;
       Teuchos::ParameterList S_sigma_prec_pl;
       S_sigma_prec_factory = invLib.getInverseFactory("S_sigma Preconditioner");
       S_sigma_prec_pl = *S_sigma_prec_factory->getParameterList();

       double* x_coordinates = S_sigma_prec_pl.sublist("ML Settings").get<double*>("x-coordinates");
       S_sigma_prec_pl.sublist("ML Settings").sublist("graddiv: 11list").set("x-coordinates",x_coordinates);
       S_sigma_prec_pl.sublist("ML Settings").sublist("graddiv: 22list").set("x-coordinates",x_coordinates);
       double* y_coordinates = S_sigma_prec_pl.sublist("ML Settings").get<double*>("y-coordinates");
       S_sigma_prec_pl.sublist("ML Settings").sublist("graddiv: 11list").set("y-coordinates",y_coordinates);
       S_sigma_prec_pl.sublist("ML Settings").sublist("graddiv: 22list").set("y-coordinates",y_coordinates);
       double* z_coordinates = S_sigma_prec_pl.sublist("ML Settings").get<double*>("z-coordinates");
       S_sigma_prec_pl.sublist("ML Settings").sublist("graddiv: 11list").set("z-coordinates",z_coordinates);
       S_sigma_prec_pl.sublist("ML Settings").sublist("graddiv: 22list").set("z-coordinates",z_coordinates);

       // add discrete curl and face mass matrix
       Teko::LinearOp Curl = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"));
       Teko::LinearOp Grad = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"));

       RCP<const Epetra_CrsMatrix> D1 = get_Epetra_CrsMatrix(*Curl);
       RCP<const Epetra_CrsMatrix> D0 = get_Epetra_CrsMatrix(*Grad);

       S_sigma_prec_pl.sublist("ML Settings").set("D1",D1);
       S_sigma_prec_pl.sublist("ML Settings").set("D0",D0);

       describeAndWriteMatrix("DiscreteGradient",*Grad,debug,dump);
       describeAndWriteMatrix("DiscreteCurl",*Curl,debug,dump);


       // We might have zero entries in the matrices, so instead of
       // setting all entries to one, we only modify the nonzero ones.
       Epetra_CrsMatrix D0_one(*D0);
       {
         int *rowptr;
         int *indices;
         double *values;
         D0_one.ExtractCrsDataPointers(rowptr, indices, values);
         for (int jj = 0; jj<D0_one.NumMyNonzeros(); jj++ )
           if (std::abs(values[jj])>1e-10)
             values[jj] = 1.0;
           else
             values[jj] = 0.0;
       }

       Epetra_CrsMatrix D1_one(*D1);
       {
         int *rowptr;
         int *indices;
         double *values;
         D1_one.ExtractCrsDataPointers(rowptr, indices, values);
         for (int jj = 0; jj<D1_one.NumMyNonzeros(); jj++ )
           if (std::abs(values[jj])>1e-10)
             values[jj] = 1.0;
           else
             values[jj] = 0.0;
       }

       RCP<Epetra_CrsMatrix> FaceNode;
       {
        RCP<Epetra_CrsMatrix> FaceNodeTemp = Teuchos::rcp(new Epetra_CrsMatrix(Copy,D1->RowMap(),0));
        EpetraExt::MatrixMatrix::Multiply(D1_one,false,D0_one,false, *FaceNodeTemp);
        FaceNode = Teuchos::rcp(new Epetra_CrsMatrix(Copy,D1->RowMap(),FaceNodeTemp->ColMap(),4));
        int nnzRow;
        int *indices;
        double *values;
        int *newIndices;
        double *newValues;
        newIndices = new int[4];
        newValues = new double[4];
        for (int jj = 0; jj < 4; jj++)
          newValues[jj] = 1.0;
        for (int row = 0; row<FaceNodeTemp->NumMyRows(); row++) {
          FaceNodeTemp->ExtractMyRowView(row, nnzRow, values, indices);
          for (int jj = 0; jj < nnzRow; jj++) {
            int nnzRowNew = 0;
            if (std::abs(values[jj])>1e-10) {
              newIndices[nnzRowNew] = indices[jj];
              nnzRowNew += 1;
            }
            if (nnzRowNew > 0) {
              int ret = FaceNode->InsertMyValues(row, nnzRowNew, newValues, newIndices);
              TEUCHOS_ASSERT(ret == 0);
            }
          }
        }
        FaceNode->FillComplete(FaceNodeTemp->DomainMap(), FaceNodeTemp->RangeMap());
       }

       RCP<const Epetra_CrsMatrix> FaceNodeConst = FaceNode;
       S_sigma_prec_pl.sublist("ML Settings").set("FaceNode",FaceNodeConst);

       if (dump) {
         EpetraExt::RowMatrixToMatrixMarketFile("EdgeNode.dat", D0_one);
         EpetraExt::BlockMapToMatrixMarketFile("rowmap_EdgeNode.dat", D0_one.RowMap());
         EpetraExt::BlockMapToMatrixMarketFile("colmap_EdgeNode.dat", D0_one.ColMap());
         EpetraExt::BlockMapToMatrixMarketFile("domainmap_EdgeNode.dat", D0_one.DomainMap());
         EpetraExt::BlockMapToMatrixMarketFile("rangemap_EdgeNode.dat", D0_one.RangeMap());
         EpetraExt::RowMatrixToMatrixMarketFile("FaceEdge.dat", D1_one);
         EpetraExt::BlockMapToMatrixMarketFile("rowmap_FaceEdge.dat", D1_one.RowMap());
         EpetraExt::BlockMapToMatrixMarketFile("colmap_FaceEdge.dat", D1_one.ColMap());
         EpetraExt::BlockMapToMatrixMarketFile("domainmap_FaceEdge.dat", D1_one.DomainMap());
         EpetraExt::BlockMapToMatrixMarketFile("rangemap_FaceEdge.dat", D1_one.RangeMap());
         EpetraExt::RowMatrixToMatrixMarketFile("FaceNode.dat", *FaceNodeConst);
         EpetraExt::BlockMapToMatrixMarketFile("rowmap_FaceNode.dat", FaceNodeConst->RowMap());
         EpetraExt::BlockMapToMatrixMarketFile("colmap_FaceNode.dat", FaceNodeConst->ColMap());
         EpetraExt::BlockMapToMatrixMarketFile("domainmap_FaceNode.dat", FaceNodeConst->DomainMap());
         EpetraExt::BlockMapToMatrixMarketFile("rangemap_FaceNode.dat", FaceNodeConst->RangeMap());
       }

       Teko::LinearOp M1_alpha = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix 1/kappa weighted AUXILIARY_EDGE"));
       RCP<const Epetra_CrsMatrix> M1 = get_Epetra_CrsMatrix(*M1_alpha);
       Epetra_CrsMatrix * TMT_Agg_Matrix;
       ML_Epetra::ML_Epetra_PtAP(*M1, *D0, TMT_Agg_Matrix,false);
       RCP<const Epetra_CrsMatrix> TMT_Agg_MatrixConst = Teuchos::rcp(TMT_Agg_Matrix);
       S_sigma_prec_pl.sublist("ML Settings").set("K0",TMT_Agg_MatrixConst);

       if (dump) {
         EpetraExt::RowMatrixToMatrixMarketFile("TMT.dat", *TMT_Agg_MatrixConst);
         EpetraExt::BlockMapToMatrixMarketFile("rowmap_TMT.dat", TMT_Agg_MatrixConst->RowMap());
         EpetraExt::BlockMapToMatrixMarketFile("colmap_TMT.dat", TMT_Agg_MatrixConst->ColMap());
         EpetraExt::BlockMapToMatrixMarketFile("domainmap_TMT.dat", TMT_Agg_MatrixConst->DomainMap());
         EpetraExt::BlockMapToMatrixMarketFile("rangemap_TMT.dat", TMT_Agg_MatrixConst->RangeMap());
       }

       {
         Teko::InverseLibrary myInvLib = invLib;
         S_sigma_prec_pl.set("Type",S_sigma_prec_type_);
         myInvLib.addInverse("S_sigma Preconditioner",S_sigma_prec_pl);
         S_sigma_prec_factory = myInvLib.getInverseFactory("S_sigma Preconditioner");
       }

       // Are we building a solver or a preconditioner?
       {
         Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Build S_sigma preconditioner"));

         if (useAsPreconditioner)
           invS_sigma = Teko::buildInverse(*S_sigma_prec_factory,S_sigma);
         else {
           if (S_sigma_prec_.is_null())
             S_sigma_prec_ = Teko::buildInverse(*S_sigma_prec_factory,S_sigma);
           else
             Teko::rebuildInverse(*S_sigma_prec_factory,S_sigma, S_sigma_prec_);
           invS_sigma = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Solve"),S_sigma,S_sigma_prec_);
         }
       }
     }
#endif
     else {
       if (useAsPreconditioner)
         invS_sigma = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Preconditioner"),S_sigma);
       else {
         if (S_sigma_prec_.is_null())
           S_sigma_prec_ = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Preconditioner"),S_sigma);
         else
           Teko::rebuildInverse(*invLib.getInverseFactory("S_sigma Preconditioner"),S_sigma, S_sigma_prec_);
         invS_sigma = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Solve"),S_sigma,S_sigma_prec_);
       }
     }
   }


   /////////////////////////////////////////////////
   // Build block  inverse matrices               //
   /////////////////////////////////////////////////

   if (!use_discrete_div_) {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Block preconditioner"));

     // Inverse blocks
     std::vector<Teko::LinearOp> diag(2);
     diag[0] = invQ_u;
     diag[1] = invS_sigma;

     // Upper tri blocks
     Teko::BlockedLinearOp U = Teko::createBlockedOp();
     Teko::beginBlockFill(U,rows,rows);
     Teko::setBlock(0,0,U,Q_u);
     Teko::setBlock(1,1,U,S_sigma);
     Teko::setBlock(0,1,U,K);
     Teko::endBlockFill(U);

     Teko::LinearOp invU = Teko::createBlockUpperTriInverseOp(U,diag);

     if (!useAsPreconditioner) {
       Teko::BlockedLinearOp invL = Teko::createBlockedOp();
       Teko::LinearOp id_u = Teko::identity(Teko::rangeSpace(Q_u));
       Teko::LinearOp id_sigma = Teko::identity(Teko::rangeSpace(Q_sigma));
       Teko::beginBlockFill(invL,rows,rows);
       Teko::setBlock(0,0,invL,id_u);
       Teko::setBlock(1,0,invL,Teko::multiply(Thyra::scale(-1.0, Kt), invQ_u));
       Teko::setBlock(1,1,invL,id_sigma);
       Teko::endBlockFill(invL);

       return Teko::multiply(invU, Teko::toLinearOp(invL));
     } else
       return invU;
   } else {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Block preconditioner"));

     Teko::LinearOp id_u = Teko::identity(Teko::rangeSpace(Q_u));

     // Inverse blocks
     std::vector<Teko::LinearOp> diag(2);
     diag[0] = id_u;
     diag[1] = invS_sigma;

     // Upper tri blocks
     Teko::BlockedLinearOp U = Teko::createBlockedOp();
     Teko::beginBlockFill(U,rows,rows);
     Teko::setBlock(0,0,U,id_u);
     Teko::setBlock(1,1,U,S_sigma);
     Teko::setBlock(0,1,U,Thyra::scale(-dt, Div));
     Teko::endBlockFill(U);

     Teko::LinearOp invU = Teko::createBlockUpperTriInverseOp(U,diag);

     if (!useAsPreconditioner) {
       Teko::BlockedLinearOp invL = Teko::createBlockedOp();
       Teko::LinearOp id_sigma = Teko::identity(Teko::rangeSpace(Q_sigma));
       Teko::beginBlockFill(invL,rows,rows);
       Teko::setBlock(0,0,invL,id_u);
       Teko::setBlock(1,0,invL,Thyra::scale(-1.0, Kt));
       Teko::setBlock(1,1,invL,id_sigma);
       Teko::endBlockFill(invL);

       if (!simplifyFaraday_) {
         Teko::BlockedLinearOp invDiag = Teko::createBlockedOp();
         Teko::beginBlockFill(invDiag,rows,rows);
         Teko::setBlock(0,0,invDiag,invQ_u);
         Teko::setBlock(1,1,invDiag,id_sigma);
         Teko::endBlockFill(invDiag);

         return Teko::multiply(invU, Teko::multiply(Teko::toLinearOp(invL), Teko::toLinearOp(invDiag)));
       } else
         return Teko::multiply(invU, Teko::toLinearOp(invL));
     } else
       return invU;
   }
}

//! Initialize from a parameter list
void FullDarcyPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   /////////////////////
   // Solver options  //
   // //////////////////

   params = pl;

   // The discrete divergence is wonky.
   // use_discrete_div_     = params.get("Use discrete div",false);
   use_discrete_div_ = false;
   dump                   = params.get("Dump",false);
   doDebug                = params.get("Debug",false);
   useAsPreconditioner    = params.get("Use as preconditioner",false);
   simplifyFaraday_       = params.get("Simplify Faraday",false) && use_discrete_div_;

   if(pl.isSublist("S_sigma Preconditioner") && pl.sublist("S_sigma Preconditioner").isParameter("Type"))
     S_sigma_prec_type_ = pl.sublist("S_sigma Preconditioner").get<std::string>("Type");
   else
     S_sigma_prec_type_ = "";

   // Output stream for debug information
   RCP<Teuchos::FancyOStream> debug = Teuchos::null;
   if (doDebug)
     debug = Teko::getOutputStream();

   //////////////////////////////////
   // Set up sub-solve factories   //
   //////////////////////////////////

   // New inverse lib to add inverse factories to
   invLib = *getInverseLibrary();

   // Q_u solve
   if (pl.isParameter("Q_u Solve")) {
     Teuchos::ParameterList cg_pl = pl.sublist("Q_u Solve");
     invLib.addInverse("Q_u Solve",cg_pl);
   }

   // Q_u preconditioner
   Teuchos::ParameterList Q_u_prec_pl = pl.sublist("Q_u Preconditioner");
   invLib.addStratPrecond("Q_u Preconditioner",
                          Q_u_prec_pl.get<std::string>("Prec Type"),
                          Q_u_prec_pl.sublist("Prec Types"));

   dt = params.get<double>("dt");

   if ((S_sigma_prec_type_ == "MueLuRefDarcy") || (S_sigma_prec_type_ == "MueLuRefMaxwell") || (S_sigma_prec_type_ == "MueLu")) { // RefDarcy based solve

     if (params.isType<std::string>("p coarsen schedule")) {
       std::string pCoarsenScheduleStr = params.get<std::string>("p coarsen schedule");
       std::vector<std::string> pCoarsenScheduleVecStr;
       panzer::StringTokenizer(pCoarsenScheduleVecStr, pCoarsenScheduleStr, ",");
       panzer::TokensToInts(pCoarsenSchedule_, pCoarsenScheduleVecStr);
       TEUCHOS_ASSERT(pCoarsenSchedule_.size() > 0);
       if ((pCoarsenSchedule_.size() == 1) and (pCoarsenSchedule_[0] == 1))
         pCoarsenSchedule_.clear();
     }

     // S_sigma solve
     Teuchos::ParameterList ml_pl = pl.sublist("S_sigma Solve");
     invLib.addInverse("S_sigma Solve",ml_pl);

     // S_sigma preconditioner
     Teuchos::ParameterList S_sigma_prec_pl = pl.sublist("S_sigma Preconditioner");



     invLib.addInverse("S_sigma Preconditioner",S_sigma_prec_pl);
   }
#ifdef PANZER_HAVE_EPETRA_STACK
   else if (S_sigma_prec_type_ == "ML") {
     // S_sigma solve
     Teuchos::ParameterList ml_pl = pl.sublist("S_sigma Solve");
     invLib.addInverse("S_sigma Solve",ml_pl);

     // S_sigma preconditioner
     Teuchos::ParameterList S_sigma_prec_pl = pl.sublist("S_sigma Preconditioner");

     invLib.addInverse("S_sigma Preconditioner",S_sigma_prec_pl);
   }
#endif
   else {
     // S_sigma solve
     if (pl.isParameter("S_sigma Solve")) {
       Teuchos::ParameterList cg_pl = pl.sublist("S_sigma Solve");
       invLib.addInverse("S_sigma Solve",cg_pl);
     }

     // S_sigma preconditioner
     Teuchos::ParameterList S_sigma_prec_pl = pl.sublist("S_sigma Preconditioner");
     invLib.addStratPrecond("S_sigma Preconditioner",
                            S_sigma_prec_pl.get<std::string>("Prec Type"),
                            S_sigma_prec_pl.sublist("Prec Types"));
   }

}

} // namespace mini_em
