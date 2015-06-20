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
#include "Teuchos_VerboseObject.hpp"

#ifdef PIRO_HAS_TPETRA
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Teuchos_VerboseObject.hpp"
#else
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
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
   void informMueLu(Teuchos::RCP<Tpetra::MultiVector<ST,LO,GO,Node> > Coordinates,
                    Teuchos::RCP<const Tpetra::Map<LO,GO,Node> > soln_map);
#else
   void informMueLu(Teuchos::RCP<Epetra_MultiVector> Coordinates,
                    Teuchos::RCP<const Epetra_Map> soln_map);
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
    MLRigidBodyModes::informMueLu(Teuchos::RCP<Tpetra::MultiVector<ST,LO,GO,Node> > Coordinates,
                                  Teuchos::RCP<const Tpetra::Map<LO,GO,Node> > soln_map) {

      //numPDEs = # PDEs
      mueLuList->set("Coordinates", Coordinates);
      mueLuList->set("number of equations", numPDEs);


      if (numElasticityDim > 0 ) {
        Teuchos::RCP<Teuchos::FancyOStream> out(Teuchos::VerboseObjectBase::getDefaultOStream());

        *out << "Setting rbms in informMueLu..." << std::endl;

       Teuchos::ArrayRCP<ST> xArray, yArray, zArray; //ArrayRCPs to x, y, and z coordinates
       xArray = Coordinates->getDataNonConst(0); //x
       if (Coordinates->getNumVectors() > 1)
         yArray = Coordinates->getDataNonConst(1); //y
       if (Coordinates->getNumVectors() > 2)
         zArray = Coordinates->getDataNonConst(2); //z

        (void) Piro_ML_Coord2RBM(x.size(), xArray.getRawPtr(), yArray.getRawPtr(),
                                 zArray.getRawPtr(), &rr[0], numPDEs, numScalar, nullSpaceDim);

        //get arrayview of rr
        Teuchos::ArrayView<ST> rrAV = Teuchos::arrayView(&rr[0], rr.size());

        //Create Tpetra_MultiVector of the rbms, to pass to MueLu.
        Teuchos::RCP<Tpetra::MultiVector<ST,LO,GO,Node> >Rbm =
          Teuchos::rcp(new Tpetra::MultiVector<ST,LO,GO,Node>(soln_map, rrAV, soln_map->getNodeNumElements(),
                           nullSpaceDim + numScalar));

        //IK, 8/18/14: it looks like  MueLu doesn't know the following 3 things
        //mueLuList->set("null space: type", "pre-computed");
        //mueLuList->set("null space: dimension", nullSpaceDim);
        //mueLuList->set("null space: add default vectors", false);
        mueLuList->set("Nullspace", Rbm);
        *out << "...done setting rbms!" << std::endl;
        //*out << "Rbm: ";
        //Rbm->describe(*out, Teuchos::VERB_EXTREME);


      }
  }

#endif

} // namespace Piro

#endif /* PIRO_NULLSPACEUTILS_HPP */
