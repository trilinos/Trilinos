#ifndef STK_MIDDLE_MESH_BOUNDING_BOX_SEARCH_OPTS_H
#define STK_MIDDLE_MESH_BOUNDING_BOX_SEARCH_OPTS_H

namespace stk {
namespace middle_mesh {
namespace search {

struct BoundingBoxSearchOpts
{
  double normalDirectionFactor  = 1.0;  // when forming bounding boxes using the normal vectors, multiplier for how
                                        // far in the normal direction to go
  double initialExpansionFactor = 0.0;  // multiplier for how much to expand bounding boxes before doing the coarse
                                        // search the first time
  double initialExpansionSum    = 0.0;  // additive factor for how much to expand the bounding boxes before doing the
                                        // coarse search the first time
  double repeatExpansionFactor = 1.5;   // multiplier for how much to expand the bounding boxes if the first coarse search
                                        // fails to find a match for all entities
  double repeatExpansionSum    = 0.0;   // additive factor for how much to expand the bounding boxes if the first coarse search
                                        // fails to find a match for all entities
};

}
}
}

#endif