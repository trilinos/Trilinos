// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2013) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Glen Hansen (gahanse@sandia.gov), Irina Kalashnikova (ikalash@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_NullSpaceUtils.hpp"

#include "Piro_StratimikosUtils.hpp"

namespace Piro {

  MLRigidBodyModes::MLRigidBodyModes(int numPDEs_)
    : numPDEs(numPDEs_),
    numElasticityDim(0),
    nullSpaceDim(0),
    numSpaceDim(0),
    numScalar(0),
    mlUsed(false),
    mueLuUsed(false)
  {
  }

  void
    MLRigidBodyModes::setPiroPL(const Teuchos::RCP<Teuchos::ParameterList>& piroParams){

      const Teuchos::RCP<Teuchos::ParameterList> stratList = extractStratimikosParams(piroParams);

      if (Teuchos::nonnull(stratList) && stratList->isParameter("Preconditioner Type")) {
        if ("ML" == stratList->get<std::string>("Preconditioner Type")) {
          // ML preconditioner is used, get nodal coordinates from application
          mlList =
            sublist(sublist(sublist(stratList, "Preconditioner Types"), "ML"), "ML Settings");
          mlUsed = true;
        }
        else if ("MueLu" == stratList->get<std::string>("Preconditioner Type")) {
          // MueLu preconditioner is used, get nodal coordinates from application
          mueLuList =
            sublist(sublist(stratList, "Preconditioner Types"), "MueLu");
          mueLuUsed = true;
        }
        else if ("MueLu-Tpetra" == stratList->get<std::string>("Preconditioner Type")) {
          // MueLu preconditioner is used, get nodal coordinates from application
          mueLuList =
            sublist(sublist(stratList, "Preconditioner Types"), "MueLu-Tpetra");
          mueLuUsed = true;
        }
      }

    }

  void
    MLRigidBodyModes::updateMLPL(const Teuchos::RCP<Teuchos::ParameterList>& mlParams){

      mlList = mlParams;

      informML();

    }

  void
    MLRigidBodyModes::resize(const int numSpaceDim_, const int numNodes){

      numSpaceDim = numSpaceDim_;

      if ( (numNodes == 0) && (numSpaceDim_ > 0) ) x.resize(1);
      else x.resize(numNodes);
      if ( (numNodes == 0) && (numSpaceDim_ > 1) ) y.resize(1);
      else y.resize(numNodes);
      if ( (numNodes == 0) && (numSpaceDim_ > 2) ) z.resize(1);
      else z.resize(numNodes);
      
      //For MueLu: allocate vector holding the x, y and z coordinates concatenated
      if ( (numNodes == 0)) xyz.resize(1);
      else xyz.resize((numSpaceDim+1)*numNodes);

      //std<vector> rr is used for both MueLu and ML 
      if(nullSpaceDim > 0) rr.resize((nullSpaceDim + numScalar) * numPDEs * numNodes + 1);

    }

  void
    MLRigidBodyModes::getCoordArrays(double **xx, double **yy, double **zz){

      *xx = &x[0];
      *yy = &y[0];
      *zz = &z[0];

    }
  
  void
    MLRigidBodyModes::getCoordArraysMueLu(double **xxyyzz){

      *xxyyzz = &xyz[0];

    }

  void
    MLRigidBodyModes::setParameters(const int numPDEs_, const int numElasticityDim_,
        const int numScalar_, const int nullSpaceDim_){

      numPDEs = numPDEs_;
      numElasticityDim = numElasticityDim_;
      numScalar = numScalar_;
      nullSpaceDim = nullSpaceDim_;

    }

  void
    MLRigidBodyModes::informML(){

      //numPDEs = # PDEs
      //numElasticityDim = # elasticity dofs
      //nullSpaceDim = dimension of elasticity nullspace
      //numScalar = # scalar dofs coupled to elasticity

      mlList->set("x-coordinates", &x[0]);
      mlList->set("y-coordinates", &y[0]);
      mlList->set("z-coordinates", &z[0]);

      mlList->set("PDE equations", numPDEs);

      if (numElasticityDim > 0 ) {

        //    std::cout << "\nEEEEE setting ML Null Space for Elasticity-type problem of Dimension: "
        //          << numElasticityDim <<  " nodes  " << x.size() << " nullspace  " << nullSpaceDim << std::endl;
        //    std::cout << "\nIKIKIK number scalar dofs: " <<numScalar <<  ", number PDEs  " << numPDEs << std::endl;

        (void) Piro_ML_Coord2RBM(x.size(), &x[0], &y[0], &z[0], &rr[0], numPDEs, numScalar, nullSpaceDim);

        //const Epetra_Comm &comm = app->getMap()->Comm();
        //Epetra_Map map(nNodes*numPDEs, 0, comm);
        //Epetra_MultiVector rbm_mv(Copy, map, rbm, nNodes*numPDEs, nullSpaceDim + numScalar);
        //std::cout << "rbm: " << rbm_mv << std::endl;
        //for (int i = 0; i<nNodes*numPDEs*(nullSpaceDim+numScalar); i++)
        //   std::cout << rbm[i] << std::endl;
        //EpetraExt::MultiVectorToMatrixMarketFile("rbm.mm", rbm_mv);

        mlList->set("null space: type", "pre-computed");
        mlList->set("null space: dimension", nullSpaceDim);
        mlList->set("null space: vectors", &rr[0]);
        mlList->set("null space: add default vectors", false);

      }
    }

  //The following function returns the rigid body modes for elasticity problems.
  //It is a modification of the ML function ml_rbm.c, extended to the case that
  //NscalarDof scalar PDEs are coupled to an elasticity problem
  //Extended by IK, Feb. 2012

  void
    MLRigidBodyModes::Piro_ML_Coord2RBM(int Nnodes,
        double x[], double y[], double z[], double rbm[], int Ndof, int NscalarDof, int NSdim)
    {


      int vec_leng, ii, jj, offset, node, dof;

      vec_leng = Nnodes*Ndof;
      for (int i = 0; i < Nnodes*Ndof*(NSdim + NscalarDof); i++)
        rbm[i] = 0.0;


      for( node = 0 ; node < Nnodes; node++ )
      {
        dof = node*Ndof;
        switch( Ndof - NscalarDof )
        {
          case 6:
            for(ii=3;ii<6+NscalarDof;ii++){ /* lower half = [ 0 I ] */
              for(jj=0;jj<6+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;
              }
            }
            /* There is no break here and that is on purpose */
          case 3:
            for(ii=0;ii<3+NscalarDof;ii++){ /* upper left = [ I ] */
              for(jj=0;jj<3+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;
              }
            }
            for(ii=0;ii<3;ii++){ /* upper right = [ Q ] */
              for(jj=3+NscalarDof;jj<6+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                // std::cout <<"jj " << jj << " " << ii + jj << std::endl;
                if(ii == jj-3-NscalarDof) rbm[offset] = 0.0;
                else {
                  if (ii+jj == 4+NscalarDof) rbm[offset] = z[node];
                  else if ( ii+jj == 5+NscalarDof ) rbm[offset] = y[node];
                  else if ( ii+jj == 6+NscalarDof ) rbm[offset] = x[node];
                  else rbm[offset] = 0.0;
                }
              }
            }
            ii = 0; jj = 5+NscalarDof; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            ii = 1; jj = 3+NscalarDof; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            ii = 2; jj = 4+NscalarDof; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            break;

          case 2:
            for(ii=0;ii<2+NscalarDof;ii++){ /* upper left = [ I ] */
              for(jj=0;jj<2+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;
              }
            }
            for(ii=0;ii<2+NscalarDof;ii++){ /* upper right = [ Q ] */
              for(jj=2+NscalarDof;jj<3+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                if (ii == 0) rbm[offset] = -y[node];
                else {
                  if (ii == 1){  rbm[offset] =  x[node];}
                  else rbm[offset] = 0.0;
                }
              }
            }
            break;

          case 1:
            for (ii = 0; ii<1+NscalarDof; ii++) {
              for (jj=0; jj<1+NscalarDof; jj++) {
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii == jj) ? 1.0 : 0.0;
              }
            }
            break;

          default:
            TEUCHOS_TEST_FOR_EXCEPTION(
                true,
                std::logic_error,
                "Piro_ML_Coord2RBM: Ndof = " << Ndof << " not implemented\n"
                );
        } /*switch*/

      } /*for( node = 0 ; node < Nnodes; node++ )*/

      return;


    } /*Piro_ML_Coord2RBM*/




} // namespace Piro
