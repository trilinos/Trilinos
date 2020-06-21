#ifndef SIDEGEOMETRY_HPP
#define SIDEGEOMETRY_HPP

#include "stk_math/StkVector.hpp"

namespace stk { namespace balance { namespace internal {

class SideGeometry
{
public:
  explicit SideGeometry(size_t numNodes);
  virtual ~SideGeometry() = default;

  virtual const stk::math::Vector3d & node(int index) const = 0;
  virtual stk::math::Vector3d centroid() const = 0;
  virtual stk::math::Vector3d closest_proj_on_face(const stk::math::Vector3d & point) const = 0;

  double min_distance_to_point(const stk::math::Vector3d & point) const;
  bool are_nodes_close_to_side(const SideGeometry & otherSide, double tolerance);

protected:
  size_t m_numNodes;
};

class PointGeometry : public SideGeometry
{
public:
  PointGeometry(const stk::math::Vector3d & n);
  ~PointGeometry() override = default;

  const stk::math::Vector3d & node(int index) const override;
  stk::math::Vector3d centroid() const override;
  stk::math::Vector3d closest_proj_on_face(const stk::math::Vector3d & point) const override;

private:
  stk::math::Vector3d m_nodeData;
};

class LineGeometry : public SideGeometry
{
public:
  LineGeometry(const stk::math::Vector3d & n0,
               const stk::math::Vector3d & n1);
  ~LineGeometry() override = default;

  const stk::math::Vector3d & node(int index) const override;
  stk::math::Vector3d centroid() const override;
  stk::math::Vector3d closest_proj_on_face(const stk::math::Vector3d & point) const override;

private:
  stk::math::Vector3d m_nodeData[2];
};

class TriGeometry : public SideGeometry
{
public:
  TriGeometry(const stk::math::Vector3d & n0,
              const stk::math::Vector3d & n1,
              const stk::math::Vector3d & n2);
  ~TriGeometry() override = default;

  const stk::math::Vector3d & node(int index) const override;
  stk::math::Vector3d centroid() const override;
  stk::math::Vector3d closest_proj_on_face(const stk::math::Vector3d & point) const override;

private:
  stk::math::Vector3d m_nodeData[3];
};

class QuadGeometry : public SideGeometry
{
public:
  QuadGeometry(const stk::math::Vector3d & n0,
               const stk::math::Vector3d & n1,
               const stk::math::Vector3d & n2,
               const stk::math::Vector3d & n3);
  ~QuadGeometry() override = default;

  const stk::math::Vector3d & node(int index) const override;
  stk::math::Vector3d centroid() const override;
  stk::math::Vector3d closest_proj_on_face(const stk::math::Vector3d & point) const override;

private:
  stk::math::Vector3d normal_dir() const;

  stk::math::Vector3d m_nodeData[4];
};


}}}

#endif
