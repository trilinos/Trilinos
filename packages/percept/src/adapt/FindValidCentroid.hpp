// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef FindValidCentroid_hpp
#define FindValidCentroid_hpp


#include <percept/PerceptMesh.hpp>
#include <percept/mesh/geometry/volume/VolumeUtil.hpp>

namespace percept {

  class FindValidCentroid {
  public:
    PerceptMesh& m_eMesh;
    const int ndiv;
    const bool m_debug;
    bool m_use_finite_volume;
    FindValidCentroid(PerceptMesh& eM, int nd=5, bool deb=false, bool use_finite_volume=false) :
        m_eMesh(eM), ndiv(nd), m_debug(deb), m_use_finite_volume(use_finite_volume)
    {
      if (eM.getProperty("FindValidCentroid_use_finite_volume") == "true")
        m_use_finite_volume = true;
    }

    double getVolumes(std::vector<double>& volumes, stk::mesh::Entity element);
    double metric(std::vector<double>& volumes, bool& foundBad);
    // return if changed
    bool findCentroid(stk::mesh::Entity element, double *c_p, std::vector<stk::mesh::Entity>& nodes, stk::mesh::Entity c_node);
  };

}
#endif
