// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_PYKRINO_SRC_PYLOGREDIRECTER_HPP_
#define KRINO_PYKRINO_SRC_PYLOGREDIRECTER_HPP_

#include <pybind11/pybind11.h>
#include <Akri_LogRedirecter.hpp>

namespace pykrino {

class PyLogRedirecter {
public:
  PyLogRedirecter();
  ~PyLogRedirecter();

  void clear();
  std::string get_log() const;

private:
  krino::LogRedirecter myRedirecter;
};

void bind_log_redirecter(pybind11::module_ &m);

} // namespace pykrino

#endif /* KRINO_PYKRINO_SRC_PYLOGREDIRECTER_HPP_ */
