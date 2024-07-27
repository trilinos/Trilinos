// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Export.hpp>

//
//  This tests vector export to a target vector with a map that has
//  no elements on some processors.
//

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  using Teuchos::rcp;

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Details::DefaultTypes::node_type NO;
  typedef Tpetra::Map<LO, GO, NO> mapType;
  typedef Tpetra::Vector<GO, LO, GO, NO> vectorType;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(-1);

  Array<GO> eltList;
  eltList.push_back(0);
  size_t num_proc = comm->getSize();
  RCP<const mapType> sourceMap = rcp(new mapType(num_proc,eltList(),0,comm));
  eltList.clear();
  if (comm->getRank()==0)
    eltList.push_back(0);
  RCP<const mapType> targetMap = rcp(new mapType(1,eltList(),0,comm));

  // uncomment the following to get wrong stats for the target map
  //size_t nge=1;
  //size_t nle=0;
  //if (comm->getRank()==0) nle=1;
  //RCP<const mapType> targetMap = rcp(new mapType(nge,nle,0,comm));

  fos->setOutputToRootOnly(0);
  *fos << "=======\nsourceMap\n=======" << std::endl;
  fos->setOutputToRootOnly(-1);
  sourceMap->describe(*fos,Teuchos::VERB_EXTREME);
  fos->setOutputToRootOnly(0);
  *fos << "=======\ntargetMap\n=======" << std::endl;
  fos->setOutputToRootOnly(-1);
  targetMap->describe(*fos,Teuchos::VERB_EXTREME);

  comm->barrier();

  RCP<vectorType> sourceVec = rcp(new vectorType(sourceMap) );
  ArrayRCP<GO> data = sourceVec->getDataNonConst();
  data[0] = 102;
  data = Teuchos::null;
  fos->setOutputToRootOnly(0);
  *fos << "=======\nsourceVec\n=======" << std::endl;
  fos->setOutputToRootOnly(-1);
  sourceVec->describe(*fos,Teuchos::VERB_EXTREME);

  comm->barrier();

  RCP<vectorType> targetVec = rcp(new vectorType(targetMap) );
  RCP<Tpetra::Export<LO,GO,NO> > exporter = rcp(new Tpetra::Export<LO,GO,NO>(sourceMap, targetMap));
  targetVec->doExport(*sourceVec,*exporter,Tpetra::ADD);
  fos->setOutputToRootOnly(0);
  *fos << "=======\ntargetVec\n=======" << std::endl;
  fos->setOutputToRootOnly(-1);
  targetVec->describe(*fos,Teuchos::VERB_EXTREME);

  comm->barrier();

  // check entry 0 for targetVec on proc 0 == 102*num_proc's
  ArrayRCP<const GO> target_data = targetVec->getData();
  if (comm->getRank() == 0) {
    GO val = target_data[0];
    GO val_expected = 102 * num_proc;
    *fos << "\nchecking " << val << " == " << val_expected << " : ";
    if (val == val_expected)
      *fos << "passed!" << std::endl;
    else
      *fos << "failed!" << std::endl;
  }

  return EXIT_SUCCESS;
}
