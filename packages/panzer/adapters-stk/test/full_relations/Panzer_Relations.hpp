#ifndef __Panzer_Relations_hpp__
#define __Panzer_Relations_hpp__

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>

// #include "Kokkos_DynRankView.hpp"
#include "Intrepid2_FieldContainer.hpp"

#include "PanzerCore_config.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

#include <vector>

namespace panzer {
class Relations{
public:
  Relations(Teuchos::RCP<panzer::ConnManager<int,int> > conn);


protected:
  Teuchos::RCP<panzer::ConnManager<int,int> > conn_;

  int dimension_;
  int num_blocks_;

  // THis is blocks, element, dimension, gids.
  std::vector<std::vector<std::vector<std::vector<unsigned> > > > base_element_mapping_;

  std::vector<shards::CellTopology> element_block_topologies_;
};

}
#endif
