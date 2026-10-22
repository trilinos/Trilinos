// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_LEVELSET_AKRI_BOUNDINGSURFACE_HPP_
#define KRINO_KRINO_LEVELSET_AKRI_BOUNDINGSURFACE_HPP_
#include <memory>
#include <string>

#include <stk_util/util/ReportHandler.hpp>
#include <Akri_Surface.hpp>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class BoundingSurface {
public:
  static BoundingSurface & build( stk::mesh::MetaData & meta, const std::string & name, std::shared_ptr<Surface> surf);
  const std::string & name() const { return myName; }
  const Surface & surface() const { STK_ThrowAssert(mySurface); return *mySurface; }
  Surface & surface() { STK_ThrowAssert(mySurface); return *mySurface; }
private:
  BoundingSurface(const std::string & surfName, std::shared_ptr<Surface> surf);
  std::string myName;
  std::shared_ptr<Surface> mySurface;
};

}


#endif /* KRINO_KRINO_LEVELSET_AKRI_BOUNDINGSURFACE_HPP_ */
