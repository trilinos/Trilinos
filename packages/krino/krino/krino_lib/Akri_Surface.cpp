// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Surface.hpp>
#include <Akri_DiagWriter.hpp>

namespace krino{

void
Surface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{
  ThrowErrorMsgIf(NULL != my_transformation,
      "This surface with type (" << type()
      << ") has motion specified, but the prepare_to_compute() method has not been implemented yet to support motion.");
}

BoundingBox
Surface::get_bounding_box()
{
  ThrowRuntimeError("This surface with type (" << type() << ") has not implemented get_bounding_box().");
}

} // namespace krino
