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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_PARAMETERLISTINTERPRETER_DEF_HPP
#define MUELU_PARAMETERLISTINTERPRETER_DEF_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_ParameterListInterpreter_decl.hpp"
#include "MueLu_FactoryFactory.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ParameterListInterpreter(Teuchos::ParameterList & paramList) {
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ParameterListInterpreter(const std::string & xmlFileName) {
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);
    SetParameterList(*paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameterList(const Teuchos::ParameterList & paramList) {

    /*
    std::cout << "Parameter List:" << std::endl
              << paramList
              << std::endl;
    */

    // Parameter List Parsing:
    // ---------
    //   <ParameterList name="MueLu">
    //    <ParameterList name="Matrix">
    //   </ParameterList>
    if (paramList.isSublist("Matrix")) {
      operatorList_ = paramList.sublist("Matrix");
      //TODO: should be validate here.
    }

    // Parameter List Parsing:
    // ---------
    //   <ParameterList name="MueLu">
    //    <ParameterList name="Factories"> <== call BuildFactoryMap() on this parameter list
    //    ...
    //    </ParameterList>
    //   </ParameterList>
    FactoryMap factoryMap;
    if (paramList.isSublist("Factories")) {
      this->BuildFactoryMap(paramList.sublist("Factories"), factoryMap, factoryMap);
    }

    // Parameter List Parsing:
    // ---------
    //   <ParameterList name="MueLu">
    //    <ParameterList name="Hierarchy">
    //     <Parameter name="verbose"  type="string" value="Warnings"/> <== get
    //     <Parameter name="numDesiredLevel" type="int" value="10"/>   <== get
    //
    //     <ParameterList name="firstLevel">                           <== parse first args and call BuildFactoryMap() on the rest of this parameter list
    //      ...
    //     </ParameterList>
    //    </ParameterList>
    //   </ParameterList>
    if (paramList.isSublist("Hierarchy")) {

      Teuchos::ParameterList hieraList = paramList.sublist("Hierarchy"); // copy because list temporally modified (remove 'id')

      // Get hierarchy options
      this->numDesiredLevel_ = 10; /* default should be provided by the Hierarchy class */;
      if(hieraList.isParameter("numDesiredLevel")) { this->numDesiredLevel_ = hieraList.get<int>("numDesiredLevel"); hieraList.remove("numDesiredLevel"); }

      this->maxCoarseSize_ = 50; /* default should be provided by the Hierarchy class */;
      if(hieraList.isParameter("maxCoarseSize")) { this->maxCoarseSize_ = hieraList.get<int>("maxCoarseSize"); hieraList.remove("maxCoarseSize"); }

      //TODO Move this its own class or MueLu::Utils?
      std::map<std::string,MsgType> verbMap;
      //for developers
      verbMap["Errors"] = Errors;
      verbMap["Warnings0"] = Warnings0;
      verbMap["Warnings00"] = Warnings00;
      verbMap["Warnings1"] = Warnings1;
      verbMap["PerfWarnings"] = PerfWarnings;
      verbMap["Runtime0"] = Runtime0;
      verbMap["Runtime1"] = Runtime1;
      verbMap["RuntimeTimings"] = RuntimeTimings;
      verbMap["NoTimeReport"] = NoTimeReport;
      verbMap["Parameters0"] = Parameters0;
      verbMap["Parameters1"] = Parameters1;
      verbMap["Statistics0"] = Statistics0;
      verbMap["Statistics1"] = Statistics1;
      verbMap["Timings0"] = Timings0;
      verbMap["Timings1"] = Timings1;
      verbMap["TimingsByLevel"] = TimingsByLevel;
      verbMap["External"] = External;
      verbMap["Debug"] = Debug;
      //for users and developers
      verbMap["None"] = None;
      verbMap["Low"] = Low;
      verbMap["Medium"] = Medium;
      verbMap["High"] = High;
      verbMap["Extreme"] = Extreme;
      if(hieraList.isParameter("verbosity")) {
        std::string vl = hieraList.get<std::string>("verbosity");
        hieraList.remove("verbosity");
        //TODO Move this to its own class or MueLu::Utils?
        if (verbMap.find(vl) != verbMap.end())
          this->verbosity_ = verbMap[vl];
        else
          TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::ParameterListInterpreter():: invalid verbosity level");
      }

      if (hieraList.isParameter("dependencyOutputLevel"))
        this->graphOutputLevel_ = hieraList.get<int>("dependencyOutputLevel");


      // Get level configuration
      for (Teuchos::ParameterList::ConstIterator param = hieraList.begin(); param != hieraList.end(); ++param) {
        const std::string & paramName  = hieraList.name(param);

        if (hieraList.isSublist(paramName)) {
          Teuchos::ParameterList levelList = hieraList.sublist(paramName); // copy because list temporally modified (remove 'id')

          int startLevel = 0;       if(levelList.isParameter("startLevel"))      { startLevel      = levelList.get<int>("startLevel");      levelList.remove("startLevel"); }
          int numDesiredLevel = 1;  if(levelList.isParameter("numDesiredLevel")) { numDesiredLevel = levelList.get<int>("numDesiredLevel"); levelList.remove("numDesiredLevel"); }

          // Parameter List Parsing:
          // ---------
          //   <ParameterList name="firstLevel">
          //      <Parameter name="startLevel"       type="int" value="0"/>
          //      <Parameter name="numDesiredLevel"  type="int" value="1"/>
          //      <Parameter name="verbose"          type="string" value="Warnings"/>
          //
          //      [] <== call BuildFactoryMap() on the rest of the parameter list
          //
          //  </ParameterList>
          FactoryMap levelFactoryMap;
          BuildFactoryMap(levelList, factoryMap, levelFactoryMap);
          RCP<FactoryManagerBase> m = rcp(new FactoryManager(levelFactoryMap));

          if (startLevel >= 0) {
            this->AddFactoryManager(startLevel, numDesiredLevel, m);
          } else if (startLevel == -1) { // -1 == coarsest level
            this->SetFactoryManagerCoarsestLevel(m);
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::ParameterListInterpreter():: invalid level id");
          }
        } /* TODO: else { } */

      }
    }
  }

  // Parameter List Parsing:
  // Create an entry in factoryMap for each parameter of the list paramList
  // ---------
  //   <ParameterList name="...">
  //     <Parameter name="smootherFact0" type="string" value="TrilinosSmoother"/>
  //
  //     <ParameterList name="smootherFact1">
  //       <Parameter name="type" type="string" value="TrilinosSmoother"/>
  //       ...
  //     </ParameterList>
  //    </ParameterList>
  //
  //TODO: static?
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildFactoryMap(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn, FactoryMap & factoryMapOut) const {
    for (Teuchos::ParameterList::ConstIterator param = paramList.begin(); param != paramList.end(); ++param) {
      const std::string             & paramName  = paramList.name(param);
      const Teuchos::ParameterEntry & paramValue = paramList.entry(param);

      //TODO: do not allow name of existing MueLu classes (can be tested using FactoryFactory)

      factoryMapOut[paramName] = FactoryFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>().BuildFactory(paramValue, factoryMapIn);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupMatrix(Matrix & Op) const {
    if(operatorList_.isParameter("PDE equations")) {
      int nPDE = operatorList_.get<int>("PDE equations");
      if (Op.GetFixedBlockSize() != nPDE)
        this->GetOStream(Warnings0,  0) << "Warning: setting matrix block size to " << nPDE << " (value of \"PDE equations\" parameter in the list) "
            << "instead of " << Op.GetFixedBlockSize() << " (provided matrix)." << std::endl;
      Op.SetFixedBlockSize(nPDE);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupExtra(Hierarchy & H) const {

    // Transfert data from the parameter list to the fine level

#ifdef MUELU_OLD_
    Array< ArrayView<const double> > arrayOfPtrs;

    if(operatorList_.isParameter("xcoords")) {
      const Array<double> & coords = operatorList_.get<Array<double> >("xcoords");
      arrayOfPtrs.push_back(coords());
      // std::cout << coords << std::endl;
    }
    if(operatorList_.isParameter("ycoords")) {
      TEUCHOS_TEST_FOR_EXCEPTION(!operatorList_.isParameter("xcoords"), Exceptions::RuntimeError, "ycoords specified but no xcoords");
      const Array<double> & coords = operatorList_.get<Array<double> >("ycoords");
      arrayOfPtrs.push_back(coords());
    }
    if(operatorList_.isParameter("zcoords")) {
      TEUCHOS_TEST_FOR_EXCEPTION(!operatorList_.isParameter("ycoords"), Exceptions::RuntimeError, "zcoords specified but no ycoords");
      const Array<double> & coords = operatorList_.get<Array<double> >("zcoords");
      arrayOfPtrs.push_back(coords());
    }

    if (arrayOfPtrs.size() != 0) {
      std::cout << "ParameterListInterpreter:: Coordinates found! (dim=" << arrayOfPtrs.size() << ")" << std::endl;
      RCP<Level>     lvl0 = H.GetLevel(0);
      RCP<const Map> map  = lvl0->Get<RCP<Matrix> >("A")->getRowMap();
      RCP<MultiVector> coordinates = MultiVectorFactory::Build(map, arrayOfPtrs, arrayOfPtrs.size());

      lvl0->Set("Coordinates", coordinates);
    }
#endif

    // TODO: should I remove the coordinates from the ParameterList to release the RCP?
    // TODO: const coordinates...

    RCP<Level>     lvl0 = H.GetLevel(0);
    //int dim=-1; // DEBUG ONLY

    if(operatorList_.isParameter("XCoordinates")) {
      ArrayRCP<SC> coords = operatorList_.get<ArrayRCP<SC> >("XCoordinates");
      lvl0->Set("XCoordinates", coords);
      //dim=1;
    }
    if(operatorList_.isParameter("YCoordinates")) {
      TEUCHOS_TEST_FOR_EXCEPTION(!operatorList_.isParameter("XCoordinates"), Exceptions::RuntimeError, "YCoordinates specified but no XCoordinates");
      ArrayRCP<SC> coords = operatorList_.get<ArrayRCP<SC> >("YCoordinates");
      lvl0->Set("YCoordinates", coords);
      //dim=2;
    }
    if(operatorList_.isParameter("ZCoordinates")) {
      TEUCHOS_TEST_FOR_EXCEPTION(!operatorList_.isParameter("YCoordinates"), Exceptions::RuntimeError, "ZCoordinates specified but no YCoordinates");
      ArrayRCP<SC> coords = operatorList_.get<ArrayRCP<SC> >("ZCoordinates");
      lvl0->Set("ZCoordinates", coords);
      //dim=3;
    }

    //    GetOStream(Runtime1, 0) << "MueLu::ParameterListInterpreter: Coordinates found! (dim=" << dim << ")" << std::endl;
  }

} // namespace MueLu

#endif // MUELU_PARAMETERLISTINTERPRETER_DEF_HPP
