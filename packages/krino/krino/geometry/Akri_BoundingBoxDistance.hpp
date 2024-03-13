#ifndef KRINO_KRINO_GEOMETRY_AKRI_BOUNDINGBOXDISTANCE_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_BOUNDINGBOXDISTANCE_HPP_
#include <stk_math/StkVector.hpp>
#include <Akri_BoundingBox.hpp>

namespace krino {

template<class REAL, unsigned DIM>
REAL min_possible_closest_squared_distance(const BoundingBox_T<REAL,DIM> & bbox, const stk::math::Vec<REAL,3> & queryPt);

template<class REAL, unsigned DIM>
REAL max_possible_closest_squared_distance(const BoundingBox_T<REAL,DIM> & bbox, const stk::math::Vec<REAL,3> & queryPt);

template<class REAL, unsigned DIM>
REAL max_possible_closest_squared_distance_between_contained_points(const BoundingBox_T<REAL,DIM> & box1, const BoundingBox_T<REAL,DIM> & box2);

}



#endif /* KRINO_KRINO_GEOMETRY_AKRI_BOUNDINGBOXDISTANCE_HPP_ */
