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
// Questions? Contact Glen Hansen (gahanse@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PIRO_NULLSPACEUTILS_HPP
#define PIRO_NULLSPACEUTILS_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#ifdef PIRO_HAS_TPETRA
#include "Tpetra_MultiVector.hpp"
#include "Kokkos_SerialNode.hpp"
#else
#include "Epetra_MultiVector.h"
#endif

namespace Piro {

class MLRigidBodyModes {

public:

   //! Construct RBM object
   MLRigidBodyModes(int numPDEs);

   //! Update the number of PDEs present
   void setNumPDEs(int numPDEs_){ numPDEs = numPDEs_; }

   //! Resize object as mesh changes
   void resize(const int numSpaceDim, const int numnodes);

   //! Set sizes of nullspace etc
   void setParameters(const int numPDEs, const int numElasticityDim, 
          const int numScalar, const int nullSpaceDim);

   //! Set Piro solver parameter list
   void setPiroPL(const Teuchos::RCP<Teuchos::ParameterList>& piroParams);

   //! Update the ML parameter list addresses and pass coordinare array addresses
   void updateMLPL(const Teuchos::RCP<Teuchos::ParameterList>& mlParams);

   //! Access the arrays to store the coordinates
   void getCoordArrays(double **x, double **y, double **z);
   
   //! Access the arrays to store the coordinates -- same as x, y and z but concatenated
   void getCoordArraysMueLu(double **xxyyzz);

   //! Is ML used on this problem?
   bool isMLUsed(){ return mlUsed; }
   
   //! Is MueLu used on this problem?
   bool isMueLuUsed(){ return mueLuUsed; }

   //! Pass coordinate arrays to ML
   void informML();
   
   //! Pass coordinate arrays to MueLu
#ifdef PIRO_HAS_TPETRA
   template<class ST, class LO, class GO, class Node>
   void informMueLu(Teuchos::RCP<Tpetra::MultiVector<ST,LO,GO,Node> > Coordinates); 
#else
   void informMueLu(Teuchos::RCP<Epetra_MultiVector> Coordinates); 
#endif

private:

    void Piro_ML_Coord2RBM(int Nnodes, double x[], double y[], double z[], double rbm[], int Ndof, int NscalarDof, int NSdim);

    int numPDEs;
    int numElasticityDim;
    int numScalar;
    int nullSpaceDim;
    int numSpaceDim;
    bool mlUsed;
    bool mueLuUsed; 
    Teuchos::RCP<Teuchos::ParameterList> mlList;
    Teuchos::RCP<Teuchos::ParameterList> mueLuList;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> xyz;
    std::vector<double> rr;

};

#ifdef PIRO_HAS_TPETRA
  template<class ST, class LO, class GO, class Node>
  void
    MLRigidBodyModes::informMueLu(Teuchos::RCP<Tpetra::MultiVector<ST,LO,GO,Node> > Coordinates) {

      std::cout << "in informMueLu!" << std::endl;
      //numPDEs = # PDEs
      mueLuList->set("Coordinates", Coordinates);
      mueLuList->set("number of equations", numPDEs);
      //TO DO: give MueLu rigid body modes!
  }
#endif

} // namespace Piro

#endif /* PIRO_NULLSPACEUTILS_HPP */
