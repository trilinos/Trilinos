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
    //    <ParameterList name="Operator">
    //   </ParameterList>
    if (paramList.isSublist("Operator")) {
      operatorList_ = paramList.sublist("Operator");
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
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupOperator(Operator & Op) const {
    if(operatorList_.isParameter("PDE equations")) {
      int nPDE = operatorList_.get<int>("PDE equations");
      Op.SetFixedBlockSize(nPDE);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupExtra(Hierarchy & H) const {
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
      RCP<const Map> map  = lvl0->Get<RCP<Operator> >("A")->getRowMap();
      RCP<MultiVector> coordinates = MultiVectorFactory::Build(map, arrayOfPtrs, arrayOfPtrs.size());
      
      lvl0->Set("Coordinates", coordinates);
    }
  }

} // namespace MueLu

#endif // MUELU_PARAMETERLISTINTERPRETER_DEF_HPP
