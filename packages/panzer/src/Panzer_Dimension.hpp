
#ifndef PANZER_DIMENSION_HPP
#define PANZER_DIMENSION_HPP

#include "Shards_Array.hpp"

namespace panzer {

  //! Spatial Dimension Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Dim)
  //! Integration Point Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(IP)
  //! Basis Point Tag (generalization of the NODE)
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(BASIS)
  //! Node Point Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(NODE)
  //! Generic Point Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Point)
  //! Cell Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
  //! Dummy Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Dummy)

}

#endif
