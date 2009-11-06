#ifndef stk_util_search_PrintBoundingBox_hpp
#define stk_util_search_PrintBoundingBox_hpp

#include <vector>
#include <ostream>

#include <stk_util/diag/Writer.hpp>
#include <stk_search/BoundingBox.hpp>

namespace stk {
namespace search {
namespace box {

template <class K, class T, int DIMENSION>
stk::diag::Writer &
operator<<(
  stk::diag::Writer &                   dout,
  const PointBoundingBox<K, T, DIMENSION> &  point_)
{
  if (dout.shouldPrint()) {
    dout << "(" << point_.key << ": ";
    for (int i=0; i<DIMENSION; ++i) {
      dout << point_.center[i];
      if (i != DIMENSION-1)
        dout << ",";
      else
        dout << ")";
    }
  }
  return dout;
}

template <class K, class T, int DIMENSION>
std::ostream &
operator<<(
    std::ostream &                   dout,
    const PointBoundingBox<K, T, DIMENSION> &  point_)
{
  dout << "(" << point_.key << ": ";
  for (int i=0; i<DIMENSION; ++i) {
    dout << point_.center[i];
    if (i != DIMENSION-1)
      dout << ",";
    else
      dout << ")";
  }
  return dout;
}

template <class K, class T, int DIMENSION>
stk::diag::Writer &
operator<<(
  stk::diag::Writer &                   dout,
  const AxisAlignedBoundingBox<K,T, DIMENSION> &  box_)
{
  if (dout.shouldPrint()) {
    dout << "(" << box_.key << ": ";
    for (int i=0; i<2*DIMENSION; ++i) {
      dout << box_.box[i];
      if (i != DIMENSION-1)
        dout << ",";
      else
        dout << ")";
    }
  }
  return dout;
}

template <class K, class T, int DIMENSION>
std::ostream &
operator<<(
    std::ostream &                   dout,
    const AxisAlignedBoundingBox<K, T, DIMENSION> &  box_)
{
  dout << "(" << box_.key << ": ";
  for (int i=0; i<2*DIMENSION; ++i) {
    dout << box_.box[i];
    if (i != DIMENSION-1)
      dout << ",";
    else
      dout << ")";
  }
  return dout;
}


template <class K, class T, int DIMENSION>
stk::diag::Writer &
operator<<(
  stk::diag::Writer &                   dout,
  const SphereBoundingBox<K, T, DIMENSION> &  sphere_)
{
  if (dout.shouldPrint()) {
    dout << "(" << sphere_.key << ", center: ";
    for (int i=0; i<DIMENSION; ++i) {
      dout << sphere_.center[i] << "," ;
    }
    dout << " radius: " << sphere_.radius << ")";
  }
  return dout;
}

template <class K, class T, int DIMENSION>
std::ostream &
operator<<(
    std::ostream &                   dout,
    const SphereBoundingBox<K, T, DIMENSION> &  sphere_)
{
  dout << "(" << sphere_.key << ", center: ";
  for (int i=0; i<DIMENSION; ++i) {
    dout << sphere_.center[i] << "," ;
  }
  dout << " radius: " << sphere_.radius << ")";
  return dout;
}


template <class K, class T, int DIMENSION>
std::ostream &
operator<<(
    std::ostream &                   dout,
    const std::vector<PointBoundingBox<K, T, DIMENSION> > &  points)
{
  dout << "[Size: " << points.size() << "]\n";
  for (unsigned i = 0; i < points.size(); ++i)
    dout << points[i] <<std::endl;
  return dout;
}

template <class K, class T, int DIMENSION>
std::ostream &
operator<<(
    std::ostream &                   dout,
    const std::vector<AxisAlignedBoundingBox<K, T, DIMENSION> > &  boxes)
{
  dout << "[Size: " << boxes.size() << "]\n";
  for (unsigned i = 0; i < boxes.size(); ++i)
    dout << boxes[i] <<std::endl;
  return dout;
}

template <class K, class T, int DIMENSION>
std::ostream &
operator<<(
    std::ostream &                   dout,
    const std::vector<SphereBoundingBox<K, T, DIMENSION> > &  spheres)
{
  dout << "[Size: " << spheres.size() << "]\n";
  for (unsigned i = 0; i < spheres.size(); ++i)
    dout << spheres[i] <<std::endl;
  return dout;
}

} // namespace box
} // namespace search
} // namespace stk

#endif //  stk_util_search_PrintBoundingBox_hpp
