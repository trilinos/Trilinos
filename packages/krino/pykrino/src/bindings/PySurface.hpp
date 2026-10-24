// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_PYKRINO_SRC_PYSURFACE_HPP_
#define KRINO_PYKRINO_SRC_PYSURFACE_HPP_

#include <pybind11/numpy.h>
#include <Akri_Surface.hpp>


namespace pykrino {

class PySurface
{
public:
    explicit PySurface(std::shared_ptr<krino::Surface> surface);
    pybind11::array_t<double> bounding_box(const double scale = 1.0) const;
    pybind11::array_t<double> signed_distance(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> points) const;
protected:
    std::shared_ptr<krino::Surface> mySurface;
};

class PySphere : public PySurface
{
public:
  PySphere(
        const std::array<double, 3> & center,
        const double radius,
        const int sign=1);
};

class PyCylinder : public PySurface
{
public:
  PyCylinder(
      const std::array<double, 3> & endpoint0,
      const std::array<double, 3> & endpoint1,
      const double radius,
      const int sign=1);
};

class PySTLSurface : public PySurface
{
public:
  PySTLSurface(const std::string & filename);
};

void bind_surfaces(pybind11::module_ &m);

}

#endif /* KRINO_PYKRINO_SRC_PYSURFACE_HPP_ */
