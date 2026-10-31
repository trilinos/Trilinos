// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <PySurface.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>
#include <stdexcept>

#include <Akri_BoundingBox.hpp>
#include <Akri_Surface.hpp>
#include <Akri_MeshSurface.hpp>
#include <Akri_AnalyticSurf.hpp>

namespace pykrino {

namespace py = pybind11;

static krino::BoundingBox
compute_point_bbox(std::vector<stk::math::Vector3d> points)
{
  krino::BoundingBox bbox;

  for ( auto && x : points )
    bbox.accommodate( x );

  return bbox;
}

static std::vector<double>
compute_surface_signed_distance_values(krino::Surface & surface, const std::vector<stk::math::Vector3d> & queryPoints)
{
  const krino::BoundingBox pointBbox = compute_point_bbox(queryPoints);
  surface.prepare_to_compute(0., pointBbox, 0.);

  std::vector<double> signedDistances;
  for (auto & x : queryPoints)
  {
    const double signedDist = surface.possibly_truncated_point_signed_distance(x, 0., 0.);
    signedDistances.push_back(signedDist);
  }
  return signedDistances;
}

PySurface::PySurface(std::shared_ptr<krino::Surface> surface)
  : mySurface(std::move(surface))
{
    if (!mySurface) {
        throw std::runtime_error("Cannot construct PySurface from null krino::Surface");
    }
}

py::array_t<double> PySurface::bounding_box(const double scale) const
{
  krino::BoundingBox box;
  mySurface->insert_into(box);

  if (scale != 1.0) {
      box.scale(scale);
  }

  const auto lower = box.get_min();
  const auto upper = box.get_max();

    py::array_t<double> result({2, 3});
    auto r = result.mutable_unchecked<2>();

    r(0, 0) = lower[0];
    r(0, 1) = lower[1];
    r(0, 2) = lower[2];

    r(1, 0) = upper[0];
    r(1, 1) = upper[1];
    r(1, 2) = upper[2];

    return result;
}

py::array_t<double> PySurface::signed_distance(py::array_t<double, py::array::c_style | py::array::forcecast> points) const
{
    py::buffer_info info = points.request();

    if (info.ndim != 2 || info.shape[1] != 3)
    {
        throw std::runtime_error("points must be a NumPy array with shape (N, 3)");
    }

    const ssize_t num_points = info.shape[0];
    const double * data = static_cast<const double *>(info.ptr);

    std::vector<stk::math::Vector3d> query_points;
    query_points.reserve(static_cast<std::size_t>(num_points));

    for (ssize_t i = 0; i < num_points; ++i)
    {
        query_points.emplace_back(
            data[3 * i + 0],
            data[3 * i + 1],
            data[3 * i + 2]
        );
    }

    std::vector<double> signed_distances;

    {
        py::gil_scoped_release release;
        signed_distances =
            compute_surface_signed_distance_values(*mySurface, query_points);
    }

    if (signed_distances.size() != static_cast<std::size_t>(num_points))
    {
        throw std::runtime_error("signed distance result has unexpected length");
    }

    py::array_t<double> result(num_points);
    auto r = result.mutable_unchecked<1>();

    for (ssize_t i = 0; i < num_points; ++i)
    {
        r(i) = signed_distances[static_cast<std::size_t>(i)];
    }

    return result;
}

PySphere::PySphere(
        const std::array<double, 3> & center,
        const double radius,
        const int sign)
: PySurface(std::make_shared<krino::Sphere>(
    stk::math::Vector3d(center[0], center[1], center[2]),
    radius,
    sign)) {}

PyCylinder::PyCylinder(
    const std::array<double, 3> & endpoint0,
    const std::array<double, 3> & endpoint1,
    const double radius,
    const int sign)
: PySurface(std::make_shared<krino::Cylinder>(
    stk::math::Vector3d(endpoint0[0], endpoint0[1], endpoint0[2]),
    stk::math::Vector3d(endpoint1[0], endpoint1[1], endpoint1[2]),
    radius,
    sign)) {}

PySTLSurface::PySTLSurface(const std::string & filename)
: PySurface(std::make_shared<krino::STLSurface>(filename)) {}

void bind_surfaces(py::module_ &m)
{
    m.doc() = "Python bindings for generic Krino surface signed-distance queries";

    py::class_<PySurface, std::shared_ptr<PySurface>>(m, "Surface")
        .def(
            "bounding_box",
            &PySurface::bounding_box,
            py::arg("scale") = 1.0,
            "Return bounding box as a NumPy array with shape (2, 3)"
        )
        .def(
            "signed_distance",
            &PySurface::signed_distance,
            py::arg("points"),
            "Compute signed distances for points with shape (N, 3)"
        );

    py::class_<PySTLSurface, PySurface, std::shared_ptr<PySTLSurface>>(m, "STLSurface")
        .def(
            py::init<const std::string &>(),
            py::arg("filename")
        );

    py::class_<PySphere, PySurface, std::shared_ptr<PySphere>>(m, "Sphere")
        .def(
            py::init<const std::array<double, 3> &, double>(),
            py::arg("center"),
            py::arg("radius")
        );

    py::class_<PyCylinder, PySurface, std::shared_ptr<PyCylinder>>(m, "Cylinder")
        .def(
            py::init<
                const std::array<double, 3> &,
                const std::array<double, 3> &,
                double
            >(),
            py::arg("point"),
            py::arg("axis"),
            py::arg("radius")
        );
}

}

