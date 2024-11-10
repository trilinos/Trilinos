// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/FindValidCentroid.hpp>

#if defined(STK_BUILT_FOR_SIERRA) && !STK_PERCEPT_LITE
#include <percept/mesh/geometry/volume/sierra_only/FiniteVolumeMesh.hpp>
#endif

namespace percept {

  double FindValidCentroid::getVolumes(std::vector<double>& volumes, stk::mesh::Entity element)
  {
    volumes.resize(0);
    const CellTopologyData *elem_cell_topo_data = m_eMesh.get_cell_topology(element);
    double element_volume = m_eMesh.volume(element, m_eMesh.get_coordinates_field(), elem_cell_topo_data);
    if (m_eMesh.numChildren(element) == 0)
      return element_volume;
    std::vector<stk::mesh::Entity> children;
    m_eMesh.getChildren(element, children, false, false);
    volumes.resize(children.size(), 0.0);
    for (unsigned ii=0; ii < children.size(); ++ii)
      {
        const CellTopologyData *cell_topo_data = m_eMesh.get_cell_topology(children[ii]);
        volumes[ii] = m_eMesh.volume(children[ii], m_eMesh.get_coordinates_field(), cell_topo_data);
#if defined(STK_BUILT_FOR_SIERRA) && !STK_PERCEPT_LITE
        if (m_use_finite_volume)
          {
            volumes[ii] = std::numeric_limits<double>::max();
            double sc_volume[8];
            FiniteVolumeMesh3D fvm(*m_eMesh.get_bulk_data());

            shards::CellTopology cell_topo(cell_topo_data);
            unsigned numNodes = cell_topo.getNodeCount();
            fvm.elementVolume(children[ii], sc_volume);
            for (unsigned in=0; in < numNodes; ++in)
              {
                volumes[ii] = std::min(volumes[ii], sc_volume[in]);
              }
          }
#endif
      }
    return element_volume;
  }

  double FindValidCentroid::metric(std::vector<double>& volumes, bool& foundBad)
  {
    double met = 0.0;
    foundBad = false;
    for (unsigned ii=0; ii < volumes.size(); ++ii)
      {
        if (volumes[ii] <= 0.0)
          {
            met += -volumes[ii];
            foundBad = true;
            break;
          }
      }
    if (foundBad)
      return met;

    met = 0.0;
    for (unsigned ii=0; ii < volumes.size(); ++ii)
      {
        met += 1.0/volumes[ii];
      }
    return met;
  }

  // return if changed
  bool FindValidCentroid::findCentroid(stk::mesh::Entity element, double *c_p, std::vector<stk::mesh::Entity>& nodes, stk::mesh::Entity c_node)
  {
    // only for coordinate field
    const int fieldDim = 3;
    stk::mesh::FieldBase *field = m_eMesh.get_coordinates_field();
    unsigned nsz = nodes.size();
    double *c_node_p = (double *)stk::mesh::field_data(*field, c_node);

    // check volumes
    stk::topology topo = m_eMesh.topology(element);
    if (topo == stk::topology::PYRAMID_5 || topo == stk::topology::WEDGE_6 || topo == stk::topology::HEX_8)
      {
        std::vector<double> volumes;
        bool foundBad = false;
        getVolumes(volumes, element);
        double met = metric(volumes, foundBad);
        if (!foundBad)
            return false;

        double metMin = std::numeric_limits<double>::max();
        double c_p_min[3] = {0,0,0};
        bool foundGood = false;
        for (int ix = 1; ix < ndiv; ++ix)
          {
            double xi = double(ix)/double(ndiv);
            for (int iy = 1; iy < ndiv; ++iy)
              {
                double eta = double(iy)/double(ndiv);
                for (int iz = 1; iz < ndiv; ++iz)
                  {
                    double zeta = double(iz)/double(ndiv);

                    double hexBases[] = {
                      (1-xi)*(1-eta)*(1-zeta),
                      (  xi)*(1-eta)*(1-zeta),
                      (  xi)*(  eta)*(1-zeta),
                      (1-xi)*(  eta)*(1-zeta),
                      (1-xi)*(1-eta)*(  zeta),
                      (  xi)*(1-eta)*(  zeta),
                      (  xi)*(  eta)*(  zeta),
                      (1-xi)*(  eta)*(  zeta)
                    };
                    double pyrBases[] = {
                      (1-xi)*(1-eta)*(1-zeta),
                      (  xi)*(1-eta)*(1-zeta),
                      (  xi)*(  eta)*(1-zeta),
                      (1-xi)*(  eta)*(1-zeta),
                      zeta
                    };
                    double wedgeBases[] = {
                      (1-xi-eta)*(1-zeta),
                      (  xi    )*(1-zeta),
                      (     eta)*(1-zeta),
                      (1-xi-eta)*(  zeta),
                      (  xi    )*(  zeta),
                      (     eta)*(  zeta)
                    };

                    double *bases = 0;
                    switch(topo.value()) {
                    case stk::topology::PYRAMID_5:
                      bases = &pyrBases[0];
                      break;
                    case stk::topology::WEDGE_6:
                      bases = &wedgeBases[0];
                      break;
                    case stk::topology::HEX_8:
                      bases = &hexBases[0];
                      break;
                    default:
                      VERIFY_MSG("bad topo");
                    }

                    for (int isp = 0; isp < fieldDim; isp++)
                      {
                        c_node_p[isp] = 0.0;
                      }

                    for (unsigned ipts=0; ipts < nsz; ipts++)
                      {
                        stk::mesh::Entity node = nodes[ipts];
                        double *  field_data = m_eMesh.field_data_inlined(field, node);
                        if (field_data)
                          {
                            for (int isp = 0; isp < fieldDim; isp++)
                              {
                                c_node_p[isp] += field_data[isp]*bases[ipts];
                              }
                          }
                      }

                    foundBad = false;
                    getVolumes(volumes, element);
                    met = metric(volumes, foundBad);
                    if (!foundBad)
                      {
                        foundGood = true;
                        if (met < metMin)
                          {
                            metMin = met;
                            for (int isp = 0; isp < fieldDim; isp++)
                              {
                                c_p_min[isp] = c_node_p[isp];
                              }
                          }
                      }

                  }
              }
          }

        if (foundGood)
          {
            for (int isp = 0; isp < fieldDim; isp++)
              {
                c_node_p[isp] = c_p_min[isp];
                c_p[isp] = c_p_min[isp];
              }
          }
        else // output_debug_data function here
          {
            std::ostringstream str;
            getVolumes(volumes, element);
            str << "negative volumes in refined mesh - couldn't find centroid, m_use_finite_volume= " << m_use_finite_volume << " parent topo= " << topo << " element= " << m_eMesh.print_entity_compact(element);

            if (1)
              {
                std::set<stk::mesh::Entity> ll;
                ll.insert(element);
                std::vector<stk::mesh::Entity> children;
                m_eMesh.getChildren(element, children, false, false);
                for (unsigned ii=0; ii < children.size(); ++ii)
                  {
                    ll.insert(children[ii]);
                  }
                std::string ff = "bad-fvc-"+toString(m_eMesh.get_rank())+".vtk";
                m_eMesh.dump_vtk(ff, false, &ll);
              }
            VERIFY_MSG(str.str());
          }

      }
    return true;
  }

}
