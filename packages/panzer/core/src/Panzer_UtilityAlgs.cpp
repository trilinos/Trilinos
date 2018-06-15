#include "Panzer_UtilityAlgs.hpp"

namespace panzer{

void reorder(std::vector<int> & order,std::function<void(int,int)> swapper) 
{
  // each entry has to be sorted
  for(std::size_t i=0;i<order.size();i++) { 

    // a, b, c
    // 2, 0, 1

    // here we are following a linked list until
    // the entry is correct
    while(order[i]!=i) {
      int nearIndex = order[i];         // 2 : 1
      int farIndex  = order[nearIndex]; // 1 : 0

      // handle the user defined swap of indices
      swapper(nearIndex,farIndex);       // a, c, b : c, a, b

      order[order[i]] = nearIndex;      // 
      order[i]        = farIndex;       // 1, 0, 2 : 0, 1, 2
    }
  }

  // at the end of this, order vector will be sorted
}

} // end namespace panzer
