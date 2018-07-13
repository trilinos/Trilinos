// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2013 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_COUPLEDRBMFACTORY_DEF_HPP
#define MUELU_COUPLEDRBMFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_CoupledRBMFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CoupledRBMFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~CoupledRBMFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoupledRBMFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    if (currentLevel.IsAvailable(nspName_, NoFactory::get()) == false && currentLevel.GetLevelID() == 0) {
      Input(currentLevel, "A");
      //Input(currentLevel,"Coordinates");
    }
    if (currentLevel.GetLevelID() !=0) {
      currentLevel.DeclareInput("Nullspace", GetFactory(nspName_).get(), this); /* ! "Nullspace" and nspName_ mismatch possible here */
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoupledRBMFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Structural acoustics nullspace factory", currentLevel);
    RCP<MultiVector> nullspace;
    if (currentLevel.GetLevelID() == 0) {
      if (currentLevel.IsAvailable(nspName_, NoFactory::get())) {
        nullspace = currentLevel.Get< RCP<MultiVector> >(nspName_, NoFactory::get());
        GetOStream(Runtime1) << "Use user-given rigid body modes " << nspName_ << ": nullspace dimension=" << nullspace->getNumVectors() << " nullspace length=" << nullspace->getGlobalLength() << std::endl;
      }
      else {
        RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
        RCP<MultiVector> Coords = Get< RCP<MultiVector> >(currentLevel,"Coordinates");
        GetOStream(Runtime1) << "Generating nullspace for structural acoustics: dimension = " << numPDEs_ << std::endl;
        RCP<const Map> xmap=A->getDomainMap();
        nullspace = MultiVectorFactory::Build(xmap, 6);
        Scalar zero = (Scalar) 0.0;
        nullspace -> putScalar(zero);
        ArrayRCP<Scalar> xnodes, ynodes, znodes;
        Scalar cx, cy, cz;
        ArrayRCP<Scalar> nsValues0, nsValues1, nsValues2, nsValues3, nsValues4, nsValues5;
        int nDOFs=xmap->getNodeNumElements();
        xnodes = Coords->getDataNonConst(0);
        ynodes = Coords->getDataNonConst(1);
        znodes = Coords->getDataNonConst(2);
        cx = Coords->getVector(0)->meanValue();
        cy = Coords->getVector(1)->meanValue();
        cz = Coords->getVector(2)->meanValue();
        nsValues0 = nullspace->getDataNonConst(0);
        nsValues1 = nullspace->getDataNonConst(1);
        nsValues2 = nullspace->getDataNonConst(2);
        nsValues3 = nullspace->getDataNonConst(3);
        nsValues4 = nullspace->getDataNonConst(4);
        nsValues5 = nullspace->getDataNonConst(5);
        for (int j=0; j<nDOFs; j+=numPDEs_) {
          Scalar one = (Scalar) 1.0;
          if( xmap->getGlobalElement(j) >= lastAcousticDOF_ ) {
            Scalar xdiff = xnodes[j]-cx;
            Scalar ydiff = ynodes[j]-cy;
            Scalar zdiff = znodes[j]-cz;
            // translation
            nsValues0[j+0] = one;
            nsValues1[j+1] = one;
            nsValues2[j+2] = one;
            // rotate around z-axis (x-y plane)
            nsValues3[j+0] = -ydiff;
            nsValues3[j+1] =  xdiff;
            // rotate around x-axis (y-z plane)
            nsValues4[j+1] = -zdiff;
            nsValues4[j+2] =  ydiff;
            // rotate around y-axis (x-z plane)
            nsValues5[j+0] =  zdiff;
            nsValues5[j+2] = -xdiff;
          }
          else {
            // translation
            nsValues0[j+0] = one;
            // insert random values and keep the top row for this node empty
            nsValues1[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues1[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues2[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues2[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues3[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues3[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues4[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues4[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues5[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues5[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
          }
        }
      } // end if "Nullspace" not available
    }
    else {
      nullspace = currentLevel.Get< RCP<MultiVector> >("Nullspace", GetFactory(nspName_).get());
    }
    Set(currentLevel, "Nullspace", nullspace);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
      void CoupledRBMFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildRBM(RCP<Matrix>& A, RCP<MultiVector>& Coords, RCP<MultiVector>& nullspace) const {
        GetOStream(Runtime1) << "Generating nullspace for structural acoustics: dimension = " << numPDEs_ << std::endl;
        RCP<const Map> xmap=A->getDomainMap();
        nullspace = MultiVectorFactory::Build(xmap, 6);
        Scalar zero = (Scalar) 0.0;
        nullspace -> putScalar(zero);
        ArrayRCP<Scalar> xnodes, ynodes, znodes;
        Scalar cx, cy, cz;
        ArrayRCP<Scalar> nsValues0, nsValues1, nsValues2, nsValues3, nsValues4, nsValues5;
        int nDOFs=xmap->getNodeNumElements();
        xnodes = Coords->getDataNonConst(0);
        ynodes = Coords->getDataNonConst(1);
        znodes = Coords->getDataNonConst(2);
        cx = Coords->getVector(0)->meanValue();
        cy = Coords->getVector(1)->meanValue();
        cz = Coords->getVector(2)->meanValue();
        nsValues0 = nullspace->getDataNonConst(0);
        nsValues1 = nullspace->getDataNonConst(1);
        nsValues2 = nullspace->getDataNonConst(2);
        nsValues3 = nullspace->getDataNonConst(3);
        nsValues4 = nullspace->getDataNonConst(4);
        nsValues5 = nullspace->getDataNonConst(5);
        for (int j=0; j<nDOFs; j+=numPDEs_) {
          Scalar one = (Scalar) 1.0;
          if( xmap->getGlobalElement(j) >= lastAcousticDOF_ ) {
            Scalar xdiff = xnodes[j]-cx;
            Scalar ydiff = ynodes[j]-cy;
            Scalar zdiff = znodes[j]-cz;
            // translation
            nsValues0[j+0] = one;
            nsValues1[j+1] = one;
            nsValues2[j+2] = one;
            // rotate around z-axis (x-y plane)
            nsValues3[j+0] = -ydiff;
            nsValues3[j+1] =  xdiff;
            // rotate around x-axis (y-z plane)
            nsValues4[j+1] = -zdiff;
            nsValues4[j+2] =  ydiff;
            // rotate around y-axis (x-z plane)
            nsValues5[j+0] =  zdiff;
            nsValues5[j+2] = -xdiff;
          }
          else {
            // translation
            nsValues0[j+0] = one;
            // insert random values and keep the top row for this node empty
            nsValues1[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues1[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues2[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues2[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues3[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues3[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues4[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues4[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues5[j+1] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
            nsValues5[j+2] = (Scalar) (((double) rand()) / ((double) RAND_MAX));
      }
    }
  }

} // namespace MueLu

#define MUELU_COUPLEDRBMFACTORY_SHORT
#endif // MUELU_COUPLEDRBMFACTORY_DEF_HPP
