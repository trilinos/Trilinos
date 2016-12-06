// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
#ifndef MUELU_GEOMETRICCOORDINATESTRANSFER_FACTORY_DEF_HPP
#define MUELU_GEOMETRICCOORDINATESTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_IO.hpp"

#include "MueLu_GeometricCoordinatesTransferFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> GeometricCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("coarseCoordinates",    Teuchos::null, "Factory for coarse coordinates generation");
    validParamList->set< int >                   ("write start",    -1, "first level at which coordinates should be written to file");
    validParamList->set< int >                   ("write end",       -1, "last level at which coordinates should be written to file");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeometricCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(coarseLevel, "coarseCoordinates");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void GeometricCoordinatesTransferFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

    using xdMV = Xpetra::MultiVector<double,LO,GO,NO>;

    GetOStream(Runtime0) << "Transferring coordinates" << std::endl;

    RCP<xdMV> coarseCoords=Get< RCP<xdMV> >(coarseLevel, "coarseCoordinates");
    Set(coarseLevel, "Coordinates", coarseCoords);

    // get user-provided printing parameters
    const ParameterList& pL = GetParameterList();
    int writeStart = pL.get<int>("write start"), writeEnd = pL.get<int>("write end");
    if (writeStart == 0 && fineLevel.GetLevelID() == 0 && writeStart <= writeEnd) {
      RCP<xdMV>    fineCoords = Get< RCP<xdMV> >(fineLevel, "Coordinates");
      std::ostringstream buf;
      buf << fineLevel.GetLevelID();
      std::string fileName = "coordinates_before_rebalance_level_" + buf.str() + ".m";
      Xpetra::IO<double,LO,GO,NO>::Write(fileName,*fineCoords);
    }
    if (writeStart <= coarseLevel.GetLevelID() && coarseLevel.GetLevelID() <= writeEnd) {
      std::ostringstream buf;
      buf << coarseLevel.GetLevelID();
      std::string fileName = "coordinates_before_rebalance_level_" + buf.str() + ".m";
      Xpetra::IO<double,LO,GO,NO>::Write(fileName,*coarseCoords);
    }
  }

} // namespace MueLu

#endif // MUELU_GEOMETRICCOORDINATESTRANSFER_FACTORY_DEF_HPP
