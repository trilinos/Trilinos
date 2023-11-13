/*
 * Akri_BoundingSurface.hpp
 *
 *  Created on: Jan 24, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_BOUNDINGSURFACE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_BOUNDINGSURFACE_HPP_
#include <memory>
#include <string>
#include <vector>

#include <stk_util/util/ReportHandler.hpp>
#include <Akri_Surface.hpp>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class BoundingSurface {
public:
  static BoundingSurface & build( stk::mesh::MetaData & meta, const std::string & name, Surface * surf );
  const std::string & name() const { return myName; }
  const Surface & surface() const { STK_ThrowAssert(mySurface); return *mySurface; }
  Surface & surface() { STK_ThrowAssert(mySurface); return *mySurface; }
private:
  BoundingSurface(const std::string & surfName, Surface * surf);
  std::string myName;
  std::unique_ptr<Surface> mySurface;
};

}


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_BOUNDINGSURFACE_HPP_ */
