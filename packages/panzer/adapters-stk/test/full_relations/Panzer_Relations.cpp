#include "Panzer_Relations.hpp"
#include "Panzer_NodalFieldPattern.hpp"
namespace panzer {

Relations::Relations(Teuchos::RCP<panzer::ConnManager<int,int> > conn) :
  conn_(conn){

  std::vector<std::string> block_ids;
  conn_->getElementBlockIds(block_ids);
  num_blocks_ = block_ids.size();
  TEUCHOS_ASSERT(num_blocks_ > 0);

  conn_->getElementBlockTopologies(element_block_topologies_);
  dimension_ = element_block_topologies_[0].getDimension();

  for (int i=0;i<num_blocks_; ++i) {
    TEUCHOS_ASSERT(element_block_topologies_[i].getDimension() == dimension_);
    std::cout << block_ids[i]<< " had block type "<<element_block_topologies_[i].getName()<<std::endl;
  }

  base_element_mapping_.resize(num_blocks_);
  for (int i=0; i< num_blocks_; ++i) {
    const std::vector<int> &elements = conn_->getElementBlock(block_ids[i]);
    base_element_mapping_[i].resize(dimension_);
    for (int idim=dimension_-1; idim >= 0; ++idim) {

    }
  }
}

}
