// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GeometryVerifier_hpp
#define GeometryVerifier_hpp

#include <iostream>


#include "Shards_CellTopology.hpp"
//#include "Teuchos_GlobalMPISession.hpp"

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

using namespace shards;


  namespace percept
  {

    template<typename T> class Histogram;

    class GeometryVerifier
    {
      //typedef std::set<std::pair<char > invalid_edge_set_type;
      int m_dump;
      double m_badJacobian;  // default and settable value for bad jacobian
      bool m_checkLocalJacobians;
      double getEquiVol(CellTopology& cell_topo);
      bool m_use_finite_volume;
    public:
      GeometryVerifier(int dump=0, double badJac=0.0, bool checkLocalJacobians = true, bool use_finite_volume= false);
      bool isGeometryBad(stk::mesh::BulkData& bulk, bool printTable=false, std::vector<double> *volume_histogram=0);
    };

  }//namespace percept

#endif
