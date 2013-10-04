#include <stk_adapt/NodeRegistry.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_percept/mesh/mod/smoother/SpacingFieldUtil.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <set>
#include <typeinfo>

namespace stk {
  namespace adapt {

    bool s_compare_using_entity_impl = false;
    static bool s_allow_empty_sub_dims = true; // for uniform refine, this should be false for testing

    static int s_nsz_parent = 0;
    static bool s_element_is_ghost = false;
    static stk::mesh::Selector *s_oldPartSelector = 0;

    // FIXME
    //static double m_min_spacing_factor = 0.5;  // reproduce (mostly) old behavior (centroids will be slightly diff)
    static double m_min_spacing_factor = 0.05;

    /// fill
    ///    @param subDimEntity with the stk::mesh::EntityId's of
    ///    the ordinal @param iSubDimOrd sub-dimensional entity of
    ///    @param element of rank
    ///    @param needed_entity_rank
    ///
    void NodeRegistry::
    noInline_getSubDimEntity(SubDimCell_SDCEntityType& subDimEntity, const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
    {
      subDimEntity.clear();
      // in the case of elements, we don't share any nodes so we just make a map of element id to node
      if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
        {
          subDimEntity.insert( element );
          //!!subDimEntity.insert(m_eMesh.identifier(element));
          return;
        }

      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      //CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::mesh::MetaData::NODE_RANK);

      const unsigned *  inodes = 0;
      unsigned nSubDimNodes = 0;
      static const unsigned edge_nodes_2[2] = {0,1};
      static const unsigned face_nodes_3[3] = {0,1,2};
      static const unsigned face_nodes_4[4] = {0,1,2,3};

      // special case for faces in 3D
      if (needed_entity_rank == m_eMesh.face_rank() && needed_entity_rank == m_eMesh.entity_rank(element))
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          if (nSubDimNodes ==3 )
            inodes = face_nodes_3;
          else
            inodes = face_nodes_4;

        }
      // special case for edges in 2D
      else if (needed_entity_rank == m_eMesh.edge_rank() && needed_entity_rank == m_eMesh.entity_rank(element))
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          if (nSubDimNodes == 2 )
            {
              inodes = edge_nodes_2;
            }
          else
            {
              throw std::runtime_error("NodeRegistry bad for edges");
            }
        }
      else if (needed_entity_rank == m_eMesh.edge_rank())
        {
          inodes = cell_topo_data->edge[iSubDimOrd].node;
          nSubDimNodes = 2;
        }
      else if (needed_entity_rank == m_eMesh.face_rank())
        {
          nSubDimNodes = cell_topo_data->side[iSubDimOrd].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = cell_topo_data->side[iSubDimOrd].node;
        }

      //subDimEntity.reserve(nSubDimNodes);
      for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
        {
          subDimEntity.insert( elem_nodes[inodes[jnode]].entity() );
          //!!subDimEntity.insert(elem_nodes[inodes[jnode]].entity()->identifier());
        }

    }

    double NodeRegistry::spacing_edge(std::vector<stk::mesh::Entity>& nodes,
                                      unsigned iv0, unsigned iv1, unsigned nsz, unsigned nsp,  double lspc[8][3], double den_xyz[3], double *coord[8])
    {
      VERIFY_OP_ON(nsz, ==, 2, "hmmm");
      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      subDimEntity.clear();
      subDimEntity.insert(nodes[iv0]);
      subDimEntity.insert(nodes[iv1]);
      bool swapped=false;
      if (nodes[iv0] != subDimEntity[0])
        swapped = true;

      static SubDimCellData new_SubDimCellData;
      static SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;
      //bool is_not_empty_but_data_cleared = (!is_empty && nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().size() == 0);
      if (is_empty) {
        if (0)
          {
            std::cout << "tmp srk P[" << m_eMesh.get_rank() << "] s_nsz_parent = " << s_nsz_parent << " s_element_is_ghost = " << s_element_is_ghost
                      << " iv0= " << iv0 << " iv1= " << iv1
                      << " n0 =  " << m_eMesh.identifier(nodes[iv0]) << " n1= " << m_eMesh.identifier(nodes[iv1])
                      << std::endl;
          }
        return 0.5;
      }
      if (0 && is_empty)
        {
          m_eMesh.dump_vtk(nodes[iv0], "node-iv0.vtk", s_oldPartSelector);
          m_eMesh.dump_vtk(nodes[iv1], "node-iv1.vtk", s_oldPartSelector);
        }
      VERIFY_OP_ON(is_empty, ==, false, "hmmm");

      //unsigned iv[2]={iv0,iv1};
      double alp[2]={0.0,0.0};
      double alps[2]={0.0,0.0};
      double alp1[2]={0.0,0.0};

      double alpsum=0.0;
      for (unsigned ipts=0; ipts < nsz; ipts++)
        {
          alp[ipts] = nodeId_elementOwnderId.get<SDC_DATA_SPACING>()[ipts];
          alps[ipts] = alp[ipts];
          VERIFY_OP_ON(alp[ipts], >=, 0.0, "hmmm33");
          alpsum += alp[ipts];
        }
      for (unsigned ipts=0; ipts < nsz; ipts++)
        {
          alp[ipts] /= alpsum;
        }
      const double fac=3.0, facden = 4.0;
      for (unsigned ipts=0; ipts < nsz; ipts++)
        {
          double lsum=0.0;
          double lsum1=0.0;
          for (unsigned jpts=0; jpts < nsz; jpts++)
            {
              lsum1 += (jpts==ipts?0:alp[jpts]);
              lsum += alp[jpts];
            }
          alp1[ipts] = (alp[ipts] + lsum1*fac)/(facden*lsum);

        }
      double sum=0.0;
      for (unsigned ipts=0; ipts < nsz; ipts++)
        {
          sum += alp1[ipts];
        }
      for (unsigned ipts=0; ipts < nsz; ipts++)
        {
          alp1[ipts] /= sum;
        }
      if (!(alp1[0] <= 1.0)) {
        std::cout << "alp1[0] = " << alp1[0] << " sum= " << sum << " alpsum= " << alpsum
                  << " alp= " << alp[0] << " " << alp[1]
                  << " alps= " << alps[0] << " " << alps[1] << std::endl;
      }
      VERIFY_OP_ON(alp1[0], <=, 1.0, "hmmm35");
      if (swapped) alp1[0] = 1.0-alp1[0];
      double candidate_alpha = 1.0-alp1[0];
      if (candidate_alpha < m_min_spacing_factor) candidate_alpha=m_min_spacing_factor;
      return candidate_alpha;
    }

    static void normalize_spacing_0(unsigned nsz, unsigned nsp, double spc[8][3], double den_xyz[3])
    {
      for (unsigned isp = 0; isp < nsp; isp++)
        {
          double den = 0.0;
          unsigned ipts=0;
          for (ipts=0; ipts < nsz; ipts++)
            {
              spc[ipts][isp] = 1.0/spc[ipts][isp];
              den += spc[ipts][isp];
            }
          for (ipts=0; ipts < nsz; ipts++)
            {
              spc[ipts][isp] /= den;
            }
          // now it's a fraction [0,1], check if it's too big
          for (ipts=0; ipts < nsz; ipts++)
            {
              if ( spc[ipts][isp] > 1.0 - m_min_spacing_factor)
                {
                  spc[ipts][isp] = 1.0 - m_min_spacing_factor;
                  for (unsigned jpts=0; jpts < nsz; jpts++)
                    {
                      if (ipts != jpts)
                        spc[ipts][isp] = m_min_spacing_factor/((double)(nsz-1));
                    }

                  break;
                }
            }
          // now renormalize it
          den = 0.0;
          for (ipts=0; ipts < nsz; ipts++)
            {
              den += spc[ipts][isp];
            }
          for (ipts=0; ipts < nsz; ipts++)
            {
              spc[ipts][isp] /= den;
            }

        }
    }


    static void check_for_min_spacing(unsigned nsz, unsigned nsp, double weights[8][3])
    {
      VERIFY_OP_ON(nsz, ==, 2, "bad nsz");
      for (unsigned isp = 0; isp < nsp; isp++)
        {
          int ifnd = -1;
          for (unsigned ipts=0; ipts < nsz; ipts++)
            {
              if (weights[ipts][isp] > 1.0 - m_min_spacing_factor)
                {
                  //weights[ipts][isp] = 1.0 - m_min_spacing_factor;
                  ifnd = ipts;
                }
            }
          if (ifnd >= 0)
            {
              weights[ifnd][isp] = 1.0 - m_min_spacing_factor;
              weights[(ifnd+1)%2][isp] = m_min_spacing_factor;
            }
        }
    }


    void NodeRegistry::normalize_spacing(stk::mesh::Entity element, std::vector<stk::mesh::Entity> &nodes,
                                         unsigned nsz, unsigned nsp, double spc[8][3], double den_xyz[3], double *coord[8])
    {
      s_nsz_parent = nsz;

      double fac = 0.0, facden=0.0;
      switch(nsz) {
      case 2:
        fac = 3.0; facden=4.0;
        break;
      case 4:
        fac = 41./3./7.; facden=16./7.; // heuristic
        break;
      case 8:
        fac = 411./259.; facden=64./37.; // heuristic
        break;
      default:
        normalize_spacing_0(nsz,nsp,spc,den_xyz);
        return;
      }

      //facden = 1+(double(nsz)-1)*fac;
      double lspc[8][3];
      for (unsigned isp = 0; isp < nsp; isp++)
        {
          for (unsigned ipts=0; ipts < nsz; ipts++)
            {
              lspc[ipts][isp] = spc[ipts][isp];
            }
        }

      for (unsigned isp = 0; isp < nsp; isp++)
        {
          if (nsz == 4)
            {
              double alp01 = spacing_edge(nodes, 0, 1, 2, nsp,  lspc, den_xyz, coord);
              double alp32 = spacing_edge(nodes, 3, 2, 2, nsp,  lspc, den_xyz, coord);
              double x = 0.5*(alp01+alp32);
              double alp12 = spacing_edge(nodes, 1, 2, 2, nsp,  lspc, den_xyz, coord);
              double alp03 = spacing_edge(nodes, 0, 3, 2, nsp,  lspc, den_xyz, coord);
              double y = 0.5*(alp12+alp03);
              if (isp == 0)
                {
                  if (0 && m_eMesh.identifier(element) == 6659)
                    {
                      PerceptMesh::get_static_instance()->print(element, false);
                      std::cout
                        << " alp01= " << alp01
                        << " alp32= " << alp32
                        << " alp12= " << alp12
                        << " alp03= " << alp03
                        << std::endl;
                    }

                  spc[0][0] = (1-x)*(1-y);
                  spc[1][0] = x*(1-y);
                  spc[2][0] = x*y;
                  spc[3][0] = (1-x)*y;
                }
              else
                {
                  for (unsigned ipts=0; ipts < nsz; ipts++)
                    {
                      spc[ipts][isp] = spc[ipts][0];
                    }
                }
            }
          else if (nsz == 8)
            {
              double alp01 = spacing_edge(nodes, 0, 1, 2, nsp,  lspc, den_xyz, coord);
              double alp32 = spacing_edge(nodes, 3, 2, 2, nsp,  lspc, den_xyz, coord);
              double alp45 = spacing_edge(nodes, 4, 5, 2, nsp,  lspc, den_xyz, coord);
              double alp76 = spacing_edge(nodes, 7, 6, 2, nsp,  lspc, den_xyz, coord);
              double x = 0.25*(alp01+alp32+alp45+alp76);
              double alp12 = spacing_edge(nodes, 1, 2, 2, nsp,  lspc, den_xyz, coord);
              double alp03 = spacing_edge(nodes, 0, 3, 2, nsp,  lspc, den_xyz, coord);
              double alp56 = spacing_edge(nodes, 5, 6, 2, nsp,  lspc, den_xyz, coord);
              double alp47 = spacing_edge(nodes, 4, 7, 2, nsp,  lspc, den_xyz, coord);
              double y = 0.25*(alp12+alp03+alp56+alp47);
              double alp04 = spacing_edge(nodes, 0, 4, 2, nsp,  lspc, den_xyz, coord);
              double alp15 = spacing_edge(nodes, 1, 5, 2, nsp,  lspc, den_xyz, coord);
              double alp26 = spacing_edge(nodes, 2, 6, 2, nsp,  lspc, den_xyz, coord);
              double alp37 = spacing_edge(nodes, 3, 7, 2, nsp,  lspc, den_xyz, coord);
              double z = 0.25*(alp04+alp15+alp26+alp37);
              if (isp == 0)
                {
                  spc[0][0] = (1-x)*(1-y)*(1-z);
                  spc[1][0] = x*(1-y)*(1-z);
                  spc[2][0] = x*y*(1-z);
                  spc[3][0] = (1-x)*y*(1-z);
                  spc[4][0] = (1-x)*(1-y)*z;
                  spc[5][0] = x*(1-y)*z;
                  spc[6][0] = x*y*z;
                  spc[7][0] = (1-x)*y*z;
                }
              else
                {
                  for (unsigned ipts=0; ipts < nsz; ipts++)
                    {
                      spc[ipts][isp] = spc[ipts][0];
                    }
                }
            }
          else if (nsz == 2)
            {
              for (unsigned ipts=0; ipts < nsz; ipts++)
                {
                  double lsum=0.0;
                  double lsum1=0.0;
                  for (unsigned jpts=0; jpts < nsz; jpts++)
                    {
                      lsum1 += (jpts==ipts?0:lspc[jpts][isp]);
                      lsum += lspc[jpts][isp];
                    }
                  spc[ipts][isp] = (lspc[ipts][isp] + lsum1*fac)/(facden*lsum);
                }
              check_for_min_spacing(nsz, nsp, spc);
            }
          else
            {
              throw std::logic_error("nsz wrong - logic error");
            }
          // these values end up as weights on the coordinates
          double sum=0.0;
          for (unsigned ipts=0; ipts < nsz; ipts++)
            {
              sum += spc[ipts][isp];
            }
          for (unsigned ipts=0; ipts < nsz; ipts++)
            {
              spc[ipts][isp] /= sum;
              if (0 && nsz == 8 && isp == 1)
                std::cout << "tmp srk spc[" << ipts << "]= " << spc[ipts][isp] << std::endl;
            }
        }
    }

    /// makes coordinates of this new node be the centroid of its sub entity - this version does it for all new nodes
    void NodeRegistry::prolongate(stk::mesh::FieldBase *field, unsigned *subDimSize_in)
    {
      EXCEPTWATCH;
      bool do_respect_spacing = m_eMesh.get_respect_spacing();
      // called from main code
      if (do_respect_spacing && !subDimSize_in)
        {
          // recurse to specialize to compute edges first
          unsigned subDimSize_2 = 2;
          prolongate(field, &subDimSize_2);
          // then signal compute all non-edges
          subDimSize_2 = 0;
          prolongate(field, &subDimSize_2);
          return;
        }

      //unsigned *null_u = 0;
      stk::mesh::FieldBase *spacing_field    = m_eMesh.get_field("ref_spacing_field");
      const mesh::Part *oldPart = m_eMesh.getPart(UniformRefinerPatternBase::getOldElementsPartName()+toString(m_eMesh.element_rank()));
      VERIFY_OP_ON(oldPart, !=, 0, "hmmm");
      stk::mesh::Selector oldPartSelector (*oldPart);
      s_oldPartSelector = &oldPartSelector;

      int spatialDim = m_eMesh.get_spatial_dim();
      int fieldDim = spatialDim;
      stk::mesh::EntityRank field_rank = stk::mesh::MetaData::NODE_RANK;
      {
        EXCEPTWATCH;
        unsigned nfr = field->restrictions().size();
        //if (print_info) std::cout << "P[" << p_rank << "] info>    number of field restrictions= " << nfr << std::endl;
        for (unsigned ifr = 0; ifr < nfr; ifr++)
          {
            const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
            //mesh::Part& frpart = metaData.get_part(fr.ordinal());
            field_rank = fr . entity_rank();
            fieldDim = fr.dimension() ;
          }
      }
      // FIXME for interpolation of element fields
      if (field_rank != stk::mesh::MetaData::NODE_RANK)
        {
          if (field_rank == stk::mesh::MetaData::ELEMENT_RANK)
            {

            }
          return;
        }

      SubDimCellToDataMap::iterator iter;

      for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
          if (do_respect_spacing && *subDimSize_in == 2 && subDimEntity.size() != 2)
            continue;

          SubDimCellData& nodeId_elementOwnderId = (*iter).second;

          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          unsigned owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          unsigned owning_elementRank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

          stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);
          VERIFY_OP_ON(owning_element, !=, stk::mesh::Entity(), "bad entity");

          unsigned char owning_elementSubDimOrd = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_ORDINAL>();
          VERIFY_OP_ON(owning_elementSubDimOrd, >, 0, "hmm 2");
          --owning_elementSubDimOrd ;

          Double2& nodeId_spacing = nodeId_elementOwnderId.get<SDC_DATA_SPACING>();

          static const SubDimCellData empty_SubDimCellData;

          bool is_empty = (nodeId_elementOwnderId == empty_SubDimCellData);

          if (is_empty)
            {
              throw std::runtime_error("prolongate(field) empty cell found");
            }

          if (nodeIds_onSE.size() != 1)
            {
              continue;
            }

          stk::mesh::EntityRank needed_entity_rank = stk::mesh::MetaData::NODE_RANK;
          // SPECIAL CASE
          // SPECIAL CASE
          // SPECIAL CASE
          if (subDimEntity.size() == 1)
            {
              needed_entity_rank = stk::mesh::MetaData::ELEMENT_RANK;
            }

          if (!m_eMesh.is_valid(nodeIds_onSE[0]))
            {
              continue;
            }

          stk::mesh::Entity c_node = nodeIds_onSE[0];

          if (!m_eMesh.is_valid(c_node))
            {
              throw std::runtime_error("prolongate(field): bad node found 0.0");
            }

          std::vector<double> c_p(fieldDim,0);
          bool doPrint = false;
          std::vector<stk::mesh::Entity> nodes(8, stk::mesh::Entity());
          unsigned nsz = 0;
          stk::mesh::Entity element_p = stk::mesh::Entity();

          if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
            {
              EXCEPTWATCH;
              {
                SDCEntityType elementId = *subDimEntity.begin();
                //!!element_p = get_entity_element(*m_eMesh.get_bulk_data(), stk::mesh::MetaData::ELEMENT_RANK, elementId);
                element_p = elementId;
                if (!m_eMesh.is_valid(element_p))
                  {
                    throw std::runtime_error("prolongate(field): bad elem found 2");
                  }
              }

              // replaced with a check for parent/child 04/09/13 srk
              //if (!oldPartSelector(element_p))
              //  continue;
              //                 if (m_eMesh.isLeafElement(element_p))
              //                   {
              //                     continue;
              //                   }
              if (!m_eMesh.isParentElement(element_p, false))
                {
                  continue;
                }

              stk::mesh::Entity element = element_p;
              bool element_is_ghost = m_eMesh.isGhostElement(element);
              s_element_is_ghost = element_is_ghost;
              if (element_is_ghost)
                {
                  //std::cout << "tmp found ghost" << std::endl;
                }
              else
                {
                  const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::mesh::MetaData::NODE_RANK);
                  unsigned npts = elem_nodes.size();
                  nsz = npts;
                  nodes.resize(nsz, stk::mesh::Entity());
                  c_p.resize(fieldDim,0);
                  //if (npts == 2) doPrint=true;
                  //double dnpts = elem_nodes.size();
                  for (unsigned ipts = 0; ipts < npts; ipts++)
                    {
                      stk::mesh::Entity node = elem_nodes[ipts].entity();
                      if (!m_eMesh.is_valid(node))
                        {
                          throw std::runtime_error("prolongate(field): bad node found 1.0");
                        }
                      nodes[ipts] = node;
                    }
                }
            }
          else
            {
              nsz = subDimEntity.size();
              nodes.resize(nsz, stk::mesh::Entity());
              c_p.resize(fieldDim,0);

              if (do_respect_spacing && nsz == 4 &&  spatialDim == 3 && owning_elementRank == m_eMesh.element_rank())
                {
                  std::vector<stk::mesh::Entity> side_node_entities;
                  m_eMesh.element_side_nodes(owning_element, owning_elementSubDimOrd, m_eMesh.face_rank(), side_node_entities);
                  VERIFY_OP_ON(side_node_entities.size(), ==, 4, "hmmm 3");
                  for (unsigned ipts=0; ipts < side_node_entities.size(); ipts++)
                    {
                      nodes[ipts] = side_node_entities[ipts];
                    }
                }
              else
                {
                  unsigned ipts=0;
                  for (SubDimCell_SDCEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids, ++ipts)
                    {
                      SDCEntityType nodeId = *ids;
                      stk::mesh::Entity node = nodeId;
                      nodes[ipts]=node;
                    }
                }
            }

          {
            //if ( (spacing_field && (spacing_field != field) && subDimEntity.size() == 2))
            // FIXME for quadratic elements
            if (do_respect_spacing && (nsz <= 8 && spacing_field && (spacing_field != field) ) )
              {
#if defined(__IBMCPP__)
                throw std::runtime_error("\nERROR: respect spacing and smoothing is not supported on IBM CPP platforms.");
#else
                EXCEPTWATCH;
                unsigned ipts=0;
                SpacingFieldUtil sfu(m_eMesh);

                double * coord[8] = {0,0,0,0,0,0,0,0};
                double * field_data[8] = {0,0,0,0,0,0,0,0};
                double * spacing[8] = {0,0,0,0,0,0,0,0};

                for (ipts=0; ipts < nsz; ipts++)
                  {
                    coord[ipts] = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), nodes[ipts]);
                    field_data[ipts] = m_eMesh.field_data_inlined(field, nodes[ipts]);
                    spacing[ipts] = m_eMesh.field_data_inlined(spacing_field, nodes[ipts]);
                  }

                double spc[8][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
                double den = 0.0;
                // cast to suppress "set but not used warnings" from gcc 4.6.3
                (void)den;
                double den_xyz[3] = {0,0,0};
                double unit_edge_vec[3] = {0,0,0};
                if (nsz == 2)
                  {
                    den = 0.0;
                    for (ipts=0; ipts < nsz; ipts++)
                      {
                        double len = 0.0;
                        for (int isp = 0; isp < spatialDim; isp++)
                          {
                            len += (coord[1][isp] - coord[0][isp])*(coord[1][isp] - coord[0][isp]);
                          }
                        VERIFY_OP_ON(len, >, 1.e-20, "bad len in prolongate");
                        for (int isp = 0; isp < spatialDim; isp++)
                          {
                            unit_edge_vec[isp] = (coord[1][isp] - coord[0][isp]) / std::sqrt(len);
                          }
                        spc[ipts][0] = sfu.spacing_at_node_in_direction(unit_edge_vec, nodes[ipts], &oldPartSelector);

                        spc[ipts][0] = std::fabs(spc[ipts][0])/std::sqrt(len);
                      }
                    for (ipts=0; ipts < nsz; ipts++)
                      {
                        nodeId_spacing[ipts] = spc[ipts][0];
                        for (int isp = 1; isp < spatialDim; isp++)
                          {
                            spc[ipts][isp] = spc[ipts][0];
                          }
                      }

                    if (0)
                      for (ipts=0; ipts < nsz; ipts++)
                        for (int isp = 0; isp < spatialDim; isp++)
                          {
                            std::cout << "y = " << coord[ipts][1] << " new spc[" << ipts << "]= " << spc[ipts][0] << std::endl;
                          }

                  }
                else
                  {
                    for (ipts=0; ipts < nsz; ipts++)
                      {
                        for (int isp = 0; isp < spatialDim; isp++)
                          {
                            spc[ipts][isp] = spacing[ipts][isp];
                          }
                      }

                  }
                normalize_spacing(element_p, nodes, nsz, spatialDim, spc, den_xyz, coord);
                if (0 && nsz==2 && (coord[1][0] < 1.e-3 && coord[0][0] < 1.e-3))
                  for (ipts=0; ipts < nsz; ipts++)
                    for (int isp = 0; isp < spatialDim; isp++)
                      {
                        std::cout << "y = " << coord[ipts][1] << " new spc[" << ipts << "]= " << spc[ipts][1] << std::endl;
                      }


                for (ipts=0; ipts < nsz; ipts++)
                  {
                    if (field_data[ipts])
                      {
                        for (int isp = 0; isp < fieldDim; isp++)
                          {
                            c_p[isp] += field_data[ipts][isp]*spc[ipts][isp];
                          }
                      }
                  }
#endif // !defined(__IBMCPP__)
              }
            else
              {
                EXCEPTWATCH;
                double dnpts = nsz;
                unsigned ipts=0;
                for (ipts=0; ipts < nsz; ipts++)
                  {
                    stk::mesh::Entity node = nodes[ipts];
                    if (!m_eMesh.is_valid(node))
                      {
                        throw std::runtime_error("prolongate(field): bad node found 2.0");
                      }
                    //double *  coord = m_eMesh.field_data(field, *node, null_u);
                    double *  field_data = m_eMesh.field_data_inlined(field, node);

                    if (doPrint && field_data)
                      {
                        //const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
                        //CellTopology cell_topo(cell_topo_data);

                        std::cout << "tmp NodeRegistry::prolongate(field) npts= " << subDimEntity.size() << " ipts= " << ipts
                                  << " field_data= " << field_data[0] << " " << field_data[1] << " " << field_data[2] << std::endl;
                      }

                    if (field_data)
                      {
                        for (int isp = 0; isp < fieldDim; isp++)
                          {
                            c_p[isp] += field_data[isp]/dnpts;
                          }
                      }
                  }
              }
          }

          // set coords
          {
            EXCEPTWATCH;

            //double *  c_coord = m_eMesh.field_data(field, *c_node, null_u);
            double *  c_coord = m_eMesh.field_data_inlined(field, c_node);

            if (c_coord)
              {
                for (int isp = 0; isp < fieldDim; isp++)
                  {
                    c_coord[isp] = c_p[isp];
                  }

                if (doPrint)
                  std::cout << "tmp NodeRegistry::prolongate(field) c_coord= " << c_coord[0] << " " << c_coord[1] << " " << c_coord[2] << std::endl;

              }
          }
        }
    } // prolongate(stk::mesh::FieldBase *)


    static void get_rbar_parts(stk::percept::PerceptMesh& eMesh, std::vector<std::string>& block_names_include, std::set<stk::mesh::Part *>& rbar_parts)
    {
      rbar_parts.clear();

      stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();
      bool found_include_only_block = false;
      for (unsigned ib = 0; ib < block_names_include.size(); ib++)
        {
          bool foundPart = false;
          for (mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
            {
              stk::mesh::Part * part = *i_part ;

              std::string bname = block_names_include[ib];
              if ('+' == bname[0])
                found_include_only_block = true;
              bname = bname.substr(1, bname.length()-1);
              if (part->name() == bname)
                {
                  foundPart = true;
                  break;
                }
            }
          if (!foundPart)
            {
              std::string msg = "UniformRefinerPattern::setNeededParts unknown block name: " + block_names_include[ib];
              throw std::runtime_error(msg.c_str());
            }
        }

      VERIFY_OP_ON(found_include_only_block, ==, true, "have to specify with +");
      for (mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
        {
          stk::mesh::Part *  part = *i_part ;
          if ( stk::mesh::is_auto_declared_part(*part) )
            continue;

          //bool doThisPart = (block_names_ranks[stk::mesh::MetaData::ELEMENT_RANK].size() == 0);
          bool doThisPart = false;

          if (!doThisPart)
            {
              // we found one block with a "+", so this means include only the actual specified list of blocks, except for those excluded with "-"
              if (found_include_only_block)
                {
                  doThisPart = false;
                  for (unsigned ib = 0; ib < block_names_include.size(); ib++)
                    {
                      std::string bname = block_names_include[ib];
                      if ('+' == bname[0])
                        {
                          bname = bname.substr(1, bname.length()-1);
                          if (part->name() == bname)
                            {
                              doThisPart = true;
                              break;
                            }
                        }
                    }
                }
              else
                // do them all, except for excludes
                {
                  doThisPart = true;
                }

              // check for excludes
              if (doThisPart)
                {
                  for (unsigned ib = 0; ib < block_names_include.size(); ib++)
                    {
                      std::string bname = block_names_include[ib];
                      if ('-' == bname[0])
                        {
                          bname = bname.substr(1, bname.length()-1);
                          if (part->name() == bname)
                            {
                              doThisPart = false;
                              break;
                            }
                        }
                    }
                }
            }
          if (doThisPart)
            {
              //std::cout << "tmp srk rbar part = " << part->name() << std::endl;
              rbar_parts.insert(part);
            }
        }
    }

    /// Check for adding new rbars - these are used for joint modeling in Salinas
    /// This version does it in bulk and thus avoids repeats on shared sub-dim entities.

#define DEBUG_ADD_RBARS 0
    void NodeRegistry::add_rbars(std::vector<std::vector<std::string> >& rbar_types )
    {
      static std::vector<stk::mesh::Part*> add_parts(1, static_cast<stk::mesh::Part*>(0));
      static std::vector<stk::mesh::Part*> remove_parts;

      std::vector<std::string>& vstr = rbar_types[m_eMesh.element_rank()];
      std::set< stk::mesh::Part * > parts_set;
      get_rbar_parts(m_eMesh, vstr, parts_set);
      stk::mesh::PartVector parts(parts_set.begin(), parts_set.end());

      unsigned nparts = parts.size();
      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];
          typedef std::pair<stk::mesh::Entity, stk::mesh::Entity> NewBarType;
          typedef std::vector<stk::mesh::Entity> EntityVector;
          std::vector<NewBarType> new_elems;
          std::vector<EntityVector> new_elems_attached_rbars;

          if (DEBUG_ADD_RBARS && !m_eMesh.get_rank())
            std::cout << "P[" << m_eMesh.get_rank() << "] NodeRegistry::add_rbars Part[" << ipart << "]= " << part.name() << std::endl;
          if ( !m_eMesh.get_rank())
            std::cout << "P[" << m_eMesh.get_rank() << "] Info: Adding rbar elements as requested by user for block[" << ipart << "]= " << part.name()
                      << "\n  NOTE:  This block is automatically ignored during refinement."
                      << std::endl;
          //std::string part_name = part.name();

          if (stk::mesh::is_auto_declared_part(part))
            continue;

          const CellTopologyData *const topology = m_eMesh.get_cell_topology(part);
          (void)topology;
          stk::mesh::Selector selector(part);

          add_parts[0] = &part;

          SubDimCellToDataMap::iterator iter;

          // Step 1. loop over all pseudo-face/edges in the NodeRegistry

          for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
            {
              const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
              SubDimCellData& nodeId_elementOwnderId = (*iter).second;

              NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
              //std::cout << "P[" << p_rank << "] info>     Part[" << ipart << "]= " << part.name()
              //              << " topology = " << (topology?shards::CellTopology(topology).getName():"null")
              //              << std::endl;

              bool found = false;
              //stk::mesh::EntityRank needed_entity_rank = stk::mesh::MetaData::NODE_RANK;
              //
              // SPECIAL CASE
              // SPECIAL CASE
              // SPECIAL CASE
              //
              if( subDimEntity.size() == 1)
                {
                  //needed_entity_rank = stk::mesh::MetaData::ELEMENT_RANK;
                  continue;
                }

              //if (subDimEntity.size() == 2) // skip beams  FIXME - need topology
              //  continue;

              // Step 2. for all nodes of this pseudo face/edge...
              stk::mesh::Entity common_node;
              bool common_node_is_valid = false;

              std::vector<stk::mesh::Entity> attached_rbars;

              for (SubDimCell_SDCEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                {
                  SDCEntityType nodeId = *ids;
                  stk::mesh::Entity node = nodeId;
                  found = false;
                  percept::MyPairIterRelation beams (m_eMesh, node, m_eMesh.element_rank());

                  if (DEBUG_ADD_RBARS > 1 && !m_eMesh.get_rank())
                    {
                      for (unsigned ii=0; ii < beams.size(); ii++)
                        {
                          if (selector(m_eMesh.bucket(beams[ii].entity())))
                            {
                              percept::MyPairIterRelation beam_nodes (m_eMesh, beams[ii].entity(), m_eMesh.node_rank());
                              VERIFY_OP_ON(beam_nodes.size(), ==, 2, "rbar issue");
                              std::cout << "node= " << m_eMesh.identifier(node) << " beam_nodes[" << m_eMesh.identifier(beams[ii].entity()) << "]= { "
                                        << std::setw(20) << m_eMesh.identifier(beam_nodes[0].entity())
                                        << std::setw(20) << m_eMesh.identifier(beam_nodes[1].entity())
                                        << std::endl;
                            }
                        }
                    }

                  // Step 3. for all beams attached to this node that belong to the current rbar block...
                  for (unsigned ii=0; ii < beams.size(); ii++)
                    {
                      if (selector(m_eMesh.bucket(beams[ii].entity())))
                        {
                          // Step 3a.  we found a beam in the current rbar block attached to current node
                          found = true;

                          attached_rbars.push_back(beams[ii].entity());

                          percept::MyPairIterRelation beam_nodes ( m_eMesh, beams[ii].entity(), m_eMesh.node_rank());
                          VERIFY_OP_ON(beam_nodes.size(), ==, 2, "rbar issue");

                          // Step 4. for beam nodes, find the node that is common to all rbars, the other node is current one from Step 2
                          for (unsigned jj=0; jj < beam_nodes.size(); jj++)
                            {
                              if (beam_nodes[jj].entity() != node)
                                {
                                  if (common_node_is_valid)
                                    {
                                      VERIFY_OP_ON(m_eMesh.identifier(common_node), ==, m_eMesh.identifier(beam_nodes[jj].entity()), "rbar issue2: please rerun and exclude rbar blocks from refinement using --block_name option ");
                                    }
                                  else
                                    {
                                      common_node = beam_nodes[jj].entity();
                                      common_node_is_valid = true;
                                    }
                                }
                            }
                          break;
                        }
                    }
                  if (!found)
                    break;
                } // loop over nodes of pseudo face/edge

              if (found)
                {
                  // create new beam element, add to part
                  unsigned nidsz = nodeIds_onSE.size();
                  for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
                    {
                      stk::mesh::Entity c_node = nodeIds_onSE[i_nid];

                      if (!m_eMesh.is_valid(c_node))
                        {
                          continue;
                        }

                      // only try to add element if I am the owner
                      if (m_eMesh.owner_rank(c_node) == m_eMesh.get_parallel_rank())
                        {
                          NewBarType new_elem(c_node, common_node);
                          new_elems.push_back(new_elem);
                          new_elems_attached_rbars.push_back(attached_rbars);
                        }
                    }
                }
            } // m_cell_2_data_map iter

          // create new beam elements, add to part
          if (new_elems.size())
            {
              if (DEBUG_ADD_RBARS && !m_eMesh.get_rank())
                std::cout << "for Part[" << ipart << "] = " << part.name() << " creating " << new_elems.size() << " new rbars" << std::endl;
              if (1 && !m_eMesh.get_rank())
                std::cout << "P[0] Info: ... for block[" << ipart << "] = " << part.name() << " creating " << new_elems.size() << " new rbars" << std::endl;
              vector<stk::mesh::Entity> new_elements;
              m_eMesh.createEntities( m_eMesh.element_rank(), new_elems.size(), new_elements);
              for (unsigned i=0; i < new_elems.size(); i++)
                {
                  std::pair<stk::mesh::Entity, stk::mesh::Entity>& new_elem = new_elems[i];
                  stk::mesh::Entity newElement = new_elements[i];
                  m_eMesh.get_bulk_data()->declare_relation(newElement, new_elem.first, 0);
                  m_eMesh.get_bulk_data()->declare_relation(newElement, new_elem.second, 1);

                  m_eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );

                  UniformRefinerPatternBase::prolongateElementFields(m_eMesh, new_elems_attached_rbars[i], newElement);

                  if (DEBUG_ADD_RBARS > 1 && !m_eMesh.get_rank())
                    {
                      std::cout << "P[" << m_eMesh.get_rank() << "] adding rbar<" << new_elem.first << ", " << new_elem.second << "> to   Part[" << ipart << "]= " << part.name()
                        //<< " topology = " << (topology ? shards::CellTopology(topology).getName() : "null")
                                << std::endl;
                    }
                }
            }
        }
      if (DEBUG_ADD_RBARS) std::cout << "tmp add_rbars " << std::endl;

    }

    void NodeRegistry::init_comm_all()
    {
      if (m_comm_all)
        delete m_comm_all;
      m_comm_all = new stk::CommAll(m_eMesh.parallel());
    }

    void NodeRegistry::init_entity_repo()
    {
      for (unsigned i = 0; i < stk::percept::EntityRankEnd; i++) m_entity_repo[i].clear();
    }

    void NodeRegistry::clear_dangling_nodes(SetOfEntities* nodes_to_be_deleted)
    {
      bool debug = false;
      if (debug) std::cout <<  "tmp srk NodeRegistry::clear_dangling_nodes start" << std::endl;
      double cpu_0 = stk::cpu_time();

      stk::mesh::BulkData &bulk_data = *m_eMesh.get_bulk_data();

      SubDimCellToDataMap::iterator iter;
      SubDimCellToDataMap& map = getMap();

      std::vector<const SubDimCell_SDCEntityType *> to_erase;
      int num_delete=0;

      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;
          stk::mesh::EntityId owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          stk::mesh::EntityRank      owning_element_rank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          VERIFY_OP_ON(nodeIds_onSE.size(), ==, nodeIds_onSE.m_entity_id_vector.size(), "NodeRegistry::clear_dangling_nodes id vector/size mismatch");
          unsigned nnodes = nodeIds_onSE.size();
          NodeIdsOnSubDimEntityType node_to_keep(0);
          //std::vector<stk::mesh::Entity> node_to_keep;
          std::vector<stk::mesh::EntityId> node_id_to_keep(0);

          if (debug)
            {
              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_element_rank, owning_elementId);
              if (!m_eMesh.is_valid(owning_element))
                {
                  std::cout << "owning_element is invalid in clear_dangling_nodes  owning_elementId= " << owning_elementId << " owning_element_rank= " << owning_element_rank
                            << " subDimCellData.size=" << subDimEntity.size()
                            << std::endl;
                  throw std::logic_error("logic: clear_dangling_nodes hmmm #111.0");
                }

            }

          for (unsigned inode=0; inode < nnodes; inode++)
            {
              stk::mesh::EntityId id = nodeIds_onSE.m_entity_id_vector[inode];
              if (!m_eMesh.is_valid(nodeIds_onSE[inode])) {
                if (debug) std::cout << " invalid deleting " << id << " sd=" << m_eMesh.identifier(subDimEntity[0]) << " " << (nnodes>=2?(int)m_eMesh.identifier(subDimEntity[1]):-1)
                                     << " owning_elementId= " << owning_elementId
                                     <<  std::endl;
                continue;
              }
              stk::mesh::EntityId id_check = m_eMesh.identifier(nodeIds_onSE[inode]);
              VERIFY_OP_ON(id_check, ==, id, "NodeRegistry::clear_dangling_nodes id");

              if (nodes_to_be_deleted && nodes_to_be_deleted->find(nodeIds_onSE[inode]) != nodes_to_be_deleted->end())
                {
                  if (debug) std::cout << " deleting " << id << " sd=" << m_eMesh.identifier(subDimEntity[0]) << " " << (nnodes>=2?(int)m_eMesh.identifier(subDimEntity[1]):-1)
                                       << " owning_elementId= " << owning_elementId
                                       <<  std::endl;
                  ++num_delete;
                }
              else if (!nodes_to_be_deleted && stk::mesh::Deleted == bulk_data.state(nodeIds_onSE[inode]) )
                {
                  if (debug) std::cout << " 2 deleting " << id << " sd=" << m_eMesh.identifier(subDimEntity[0]) << " " <<  (nnodes>=2?(int)m_eMesh.identifier(subDimEntity[1]):-1)
                                       << " owning_elementId= " << owning_elementId
                                       <<  std::endl;
                  ++num_delete;
                }
              else
                {
                  node_to_keep.push_back(nodeIds_onSE[inode]);
                  node_id_to_keep.push_back(id);
                }
            }
          nodeIds_onSE = node_to_keep;
          nodeIds_onSE.m_entity_id_vector = node_id_to_keep;
          if (nodeIds_onSE.size() != nodeIds_onSE.m_entity_id_vector.size())
            {
              std::cout << "NodeRegistry::clear_dangling_nodes id vector/size mismatch 1 size= " << nodeIds_onSE.size() << " id.size= " << nodeIds_onSE.m_entity_id_vector.size() << std::endl;
            }
          VERIFY_OP_ON(nodeIds_onSE.size(), ==, nodeIds_onSE.m_entity_id_vector.size(), "NodeRegistry::clear_dangling_nodes id vector/size mismatch 1");

          if (nodeIds_onSE.size() == 0)
            {
              stk::mesh::EntityId id = nodeIds_onSE.m_entity_id_vector[0];
              if (debug) std::cout << " to_erase deleting " << id << " sd=" << m_eMesh.identifier(subDimEntity[0]) << " " << (nnodes>=2?(int)m_eMesh.identifier(subDimEntity[1]):-1)
                                   << " owning_elementId= " << owning_elementId
                                   <<  std::endl;
              to_erase.push_back(&(iter->first));
            }

        }
      if (debug) std::cout << "tmp srk NodeRegistry::clear_dangling_nodes num_delete= " << num_delete <<  std::endl;
      if (to_erase.size())
        {
          s_compare_using_entity_impl = true;
          if (debug) std::cout << "tmp srk NodeRegistry::clear_dangling_nodes nodeIds_onSE.size() != node_to_keep.size()), to_erase= " << to_erase.size() <<  std::endl;
          for (unsigned i=0; i < to_erase.size(); i++)
            {
              if ( to_erase[i]->size() &&  map.find(*to_erase[i]) != map.end())
                {
                  map.erase(*to_erase[i]);
                }
            }
          s_compare_using_entity_impl = false;
        }

      // check
      if (1)
        {
          for (iter = map.begin(); iter != map.end(); ++iter)
            {
              SubDimCellData& nodeId_elementOwnderId = (*iter).second;
              NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
              VERIFY_OP_ON(nodeIds_onSE.size(), ==, nodeIds_onSE.m_entity_id_vector.size(), "NodeRegistry::clear_dangling_nodes id vector/size mismatch after erase");
            }
        }
      double cpu_1 = stk::cpu_time();
      if (debug) std::cout <<  "tmp srk NodeRegistry::clear_dangling_nodes end, time= " << (cpu_1-cpu_0) << std::endl;
    }

    void NodeRegistry::clear_elements_to_be_deleted(SetOfEntities* elements_to_be_deleted)
    {
      bool debug = false;
      if (debug) std::cout <<  "tmp srk NodeRegistry::clear_elements_to_be_deleted start" << std::endl;

      SubDimCellToDataMap::iterator iter;
      SubDimCellToDataMap& map = getMap();

      std::vector<const SubDimCell_SDCEntityType *> to_erase;
      int num_delete=0;

      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;
          stk::mesh::EntityId owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          stk::mesh::EntityRank      owning_element_rank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          VERIFY_OP_ON(nodeIds_onSE.size(), ==, nodeIds_onSE.m_entity_id_vector.size(), "NodeRegistry::clear_elements_to_be_deleted id vector/size mismatch");

          stk::mesh::Entity owning_element = stk::mesh::Entity();
          if (owning_elementId)
            {
              owning_element = m_eMesh.get_bulk_data()->get_entity(owning_element_rank, owning_elementId);
              if (!m_eMesh.is_valid(owning_element))
                {
                  std::cout << "owning_element is invalid in clear_elements_to_be_deleted  owning_elementId= " << owning_elementId << " owning_element_rank= " << owning_element_rank
                            << " subDimCellData.size=" << subDimEntity.size()
                            << std::endl;
                  throw std::logic_error("logic: clear_elements_to_be_deleted hmmm #211.0");
                }
            }

          if (!owning_elementId || (elements_to_be_deleted && elements_to_be_deleted->find(owning_element) != elements_to_be_deleted->end()))
            {
              if (debug) std::cout << " deleting "
                                   << " owning_elementId= " << owning_elementId
                                   << " subDimEntity size= " << subDimEntity.size()
                                   <<  std::endl;
              ++num_delete;
              to_erase.push_back(&(iter->first));
            }
        }
      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;
          stk::mesh::EntityId owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          if (!owning_elementId) continue;
          stk::mesh::EntityRank      owning_element_rank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
          stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_element_rank, owning_elementId);
          if (!m_eMesh.is_valid(owning_element))
            {
              std::cout << "owning_element is invalid in clear_elements_to_be_deleted  owning_elementId= " << owning_elementId << " owning_element_rank= " << owning_element_rank
                        << " subDimCellData.size=" << subDimEntity.size()
                        << std::endl;
              throw std::logic_error("logic: clear_elements_to_be_deleted hmmm #211.1");
            }
        }
      if (debug) std::cout << "tmp srk NodeRegistry::clear_elements_to_be_deleted num_delete= " << num_delete <<  std::endl;
      if (to_erase.size())
        {
          s_compare_using_entity_impl = true;
          for (unsigned i=0; i < to_erase.size(); i++)
            {
              if ( to_erase[i]->size() &&  map.find(*to_erase[i]) != map.end())
                {
                  map.erase(*to_erase[i]);
                }
            }
          s_compare_using_entity_impl = false;
        }

    }

    void NodeRegistry::initialize()
    {
      m_cell_2_data_map.clear();
      init_entity_repo();
    }

    void NodeRegistry::
    beginRegistration()
    {
      if (CHECK_DB) checkDB("beginRegistration");

      m_nodes_to_ghost.resize(0);
      m_pseudo_entities.clear();
      m_state = NRS_START_REGISTER_NODE;
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::beginRegistration" << std::endl;
    }

    void NodeRegistry::
    endRegistration()
    {
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endRegistration start" << std::endl;

      if (CHECK_DB) checkDB("endRegistration");

      removeUnmarkedSubDimEntities();

      if (CHECK_DB) checkDB("endRegistration - after removeUnmarkedSubDimEntities");

      m_eMesh.get_bulk_data()->modification_begin();
      this->createNewNodesInParallel();
      m_nodes_to_ghost.resize(0);

      if (CHECK_DB) checkDB("endRegistration - after createNewNodesInParallel");

#if STK_ADAPT_NODEREGISTRY_DO_REHASH
      m_cell_2_data_map.rehash(m_cell_2_data_map.size());
#endif
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endRegistration end" << std::endl;

      if (CHECK_DB) checkDB("endRegistration - after rehash");

      m_state = NRS_END_REGISTER_NODE;

    }

    void NodeRegistry::
    beginLocalMeshMods()
    {
    }

    void NodeRegistry::
    endLocalMeshMods()
    {
    }

    void NodeRegistry::
    beginCheckForRemote()
    {
      m_state = NRS_START_CHECK_FOR_REMOTE;
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::beginCheckForRemote " << std::endl;
    }

    void NodeRegistry::
    endCheckForRemote()
    {
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endCheckForRemote start " << std::endl;
      stk::ParallelMachine pm = m_eMesh.parallel();
      int failed = 0;
      stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

      this->allocateBuffers();

#if STK_ADAPT_NODEREGISTRY_DO_REHASH
      m_cell_2_data_map.rehash(m_cell_2_data_map.size());
#endif

      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endCheckForRemote end " << std::endl;

      m_state = NRS_END_CHECK_FOR_REMOTE;

    }

    void NodeRegistry::
    beginGetFromRemote()
    {
      m_state = NRS_START_GET_FROM_REMOTE;
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::beginGetFromRemote  " << std::endl;

    }
    void NodeRegistry::
    endGetFromRemote()
    {
      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endGetFromRemote start " << std::endl;
      stk::ParallelMachine pm = m_eMesh.parallel();
      int failed = 0;
      stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

      if (CHECK_DB) checkDB("endGetFromRemote before comm");

      this->communicate();

      failed = 0;
      stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

      if (m_useCustomGhosting)
        {
          stk::mesh::Ghosting & ghosting = m_eMesh.get_bulk_data()->create_ghosting( std::string("new_nodes") );
          vector<stk::mesh::EntityKey> receive;
          ghosting.receive_list( receive );
          m_eMesh.get_bulk_data()->change_ghosting( ghosting, m_nodes_to_ghost, receive);
        }

      failed = 0;
      stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

      m_eMesh.get_bulk_data()->modification_end();

      if (CHECK_DB) checkDB("endGetFromRemote after mod end");

      if (m_useCustomGhosting) setAllReceivedNodeData();

      if (CHECK_DB) checkDB("endGetFromRemote after setAllReceivedNodeData");

      if (m_debug)
        std::cout << "P[" << m_eMesh.get_rank() << "] tmp NodeRegistry::endGetFromRemote end " << std::endl;

      m_state = NRS_END_GET_FROM_REMOTE;
    }

    void NodeRegistry::setAllReceivedNodeData()
    {
      EXCEPTWATCH;
      SubDimCellToDataMap::iterator iter;
      stk::mesh::PartVector empty_parts;

      for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
        {
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;

          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

          if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
            {
              std::cout << " nodeIds_onSE.m_entity_id_vector.size() = " << nodeIds_onSE.m_entity_id_vector.size()
                        << " nodeIds_onSE.size() = " <<  nodeIds_onSE.size() << std::endl;

              throw std::logic_error("NodeRegistry:: setAllReceivedNodeData logic err #0");
            }

          for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
            {
              if (!nodeIds_onSE.m_entity_id_vector[ii])
                {
                  nodeIds_onSE[ii] = stk::mesh::Entity();
                }
              else
                {
                  stk::mesh::Entity node = get_entity_node_I(*m_eMesh.get_bulk_data(), stk::mesh::MetaData::NODE_RANK, nodeIds_onSE.m_entity_id_vector[ii]);  // FIXME
                  if (0) std::cout << "tmp P[" << m_eMesh.get_rank() << "] NodeRegistry::setAllReceivedNodeData id= " << nodeIds_onSE.m_entity_id_vector[ii] << std::endl;
                  if (!m_eMesh.is_valid(node))
                    {
                      if (m_useCustomGhosting)
                        {
                          throw std::logic_error("NodeRegistry:: setAllReceivedNodeData logic err #3");
                        }
                      else
                        {
                          //std::cout << "tmp P[" << m_eMesh.get_rank() << "] NodeRegistry::setAllReceivedNodeData id= " << nodeIds_onSE.m_entity_id_vector[ii] << std::endl;
                          node = m_eMesh.get_bulk_data()->declare_entity(m_eMesh.node_rank(), nodeIds_onSE.m_entity_id_vector[ii], empty_parts);
                          if (!m_eMesh.is_valid(node)) throw std::logic_error("NodeRegistry:: setAllReceivedNodeData logic err #3.1");
                        }
                    }
                  if (!m_eMesh.is_valid(node)) throw std::logic_error("NodeRegistry:: setAllReceivedNodeData logic err #3.2");
                  nodeIds_onSE[ii] = node;
                }
            }
        }
    }

    /// when a sub-dim entity is visited during node registration but is flagged as not being marked, and thus not requiring
    ///   any new nodes, we flag it with NR_MARK_NONE, then remove it here
    void NodeRegistry::removeUnmarkedSubDimEntities()
    {
      EXCEPTWATCH;
      SubDimCellToDataMap::iterator iter;
      bool debug = false;

      for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
        {
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;

          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          if (nodeIds_onSE.size())
            {
              unsigned mark = nodeIds_onSE.m_mark;
              unsigned is_marked = mark & NR_MARK;
              unsigned is_not_marked = mark & NR_MARK_NONE;

              if (!is_marked && is_not_marked)
                {
                  if (debug) std::cout << "tmp SRK FOUND mark = " << mark << " NR_MARK =" << NR_MARK << " NR_MARK_NONE= " << NR_MARK_NONE << std::endl;
                  // check if the node is a "hanging node" in which case it has relations
                  bool found = false;
                  for (unsigned in=0; in < nodeIds_onSE.size(); ++in)
                    {
                      if(!m_eMesh.is_valid(nodeIds_onSE[in])) continue;
                      if (debug) std::cout << "is valid = " << m_eMesh.is_valid(nodeIds_onSE[in]) << std::endl;
                      size_t num_rels = m_eMesh.get_bulk_data()->count_relations(nodeIds_onSE[in]);
                      if (num_rels)
                        {
                          if (debug) std::cout << "tmp SRK num_rels is non-zero in removeUnmarkedSubDimEntities, id= " << m_eMesh.identifier(nodeIds_onSE[in]) <<  std::endl;
                          found = true;
                          break;
                        }
                    }

                  if (!found)
                    {
                      nodeIds_onSE.resize(0);
                    }
                }
            }
        }
    }

    //bool is_empty( const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd);

    /// Register the need for a new node on the sub-dimensional entity @param subDimEntity on element @param element.
    /// If the element is a ghost element, the entity is still registered: the locality/ownership of the new entity
    /// can be determined by the locality of the element (ghost or not).
    bool NodeRegistry::registerNeedNewNode(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes)
    {
      bool ret_val = false;
      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);

      static SubDimCellData new_SubDimCellData;
      static SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;
      bool is_not_empty_but_data_cleared = (!is_empty && nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().size() == 0);

      // if empty or if my id is the smallest, make this element the owner
      stk::mesh::EntityId db_id = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
      stk::mesh::EntityId db_rank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

      stk::mesh::EntityId element_id = m_eMesh.identifier(element);
      stk::mesh::EntityId element_rank = m_eMesh.entity_rank(element);
      bool should_put_in_id = (element_id < db_id);
      bool should_put_in_rank_gt = (element_rank > db_rank);
      bool should_put_in_rank_gte = (element_rank >= db_rank);
      bool should_put_in = should_put_in_rank_gt || (should_put_in_id && should_put_in_rank_gte);

      unsigned smark=0;
      if (!is_empty)
        {
          unsigned& mark = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().m_mark;
          if (needNodes)
            mark |= NR_MARK;
          else
            mark |= NR_MARK_NONE;
          smark = mark;
        }

      /// once it's in, the assertion should be:
      ///   owning_elementId < non_owning_elementId && owning_elementRank >= non_owning_elementRank
      ///
      if (is_empty || is_not_empty_but_data_cleared || should_put_in)
        {
          // new SubDimCellData SDC_DATA_OWNING_ELEMENT_KEY
          // CHECK

          unsigned numNewNodes = needed_entity_rank.second;

          NodeIdsOnSubDimEntityType nid_new(numNewNodes, stk::mesh::Entity(), smark);
          //NodeIdsOnSubDimEntityType nid_new(numNewNodes, stk::mesh::Entity(), 0u);
          if (needNodes)
            nid_new.m_mark |= NR_MARK;
          else
            nid_new.m_mark |= NR_MARK_NONE;

          smark = nid_new.m_mark;

          // create SubDimCellData for the map rhs
          // add one to iSubDimOrd for error checks later
          VERIFY_OP_ON(m_eMesh.identifier(element), !=, 0, "hmmm registerNeedNewNode #1");
          VERIFY_OP_ON(m_eMesh.is_valid(element), ==, true, "hmmm registerNeedNewNode #2");
          if (is_empty || is_not_empty_but_data_cleared)
            {
              SubDimCellData data(nid_new, stk::mesh::EntityKey(m_eMesh.entity_rank(element), m_eMesh.identifier(element)), iSubDimOrd+1 );
              putInMap(subDimEntity,  data);
            }
          else
            {
              NodeIdsOnSubDimEntityType& nid = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
              SubDimCellData data(nid, stk::mesh::EntityKey(m_eMesh.entity_rank(element), m_eMesh.identifier(element)), iSubDimOrd+1 );
              putInMap(subDimEntity,  data);
            }
          ret_val = true;
        }

      if (DEBUG_NR_UNREF)
        {
          SubDimCellData* a_nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
          SubDimCellData& a_nodeId_elementOwnderId = (a_nodeId_elementOwnderId_ptr ? *a_nodeId_elementOwnderId_ptr : empty_SubDimCellData);
          bool a_is_empty = a_nodeId_elementOwnderId_ptr == 0;
          unsigned& gotMark = a_nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().m_mark;
          unsigned nidSize = a_nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().size();

          std::ostringstream sout;
          sout << "P[" << m_eMesh.get_rank() << "] registerNeedNewNode:: element= " << m_eMesh.identifier(element) << " nidSize= " << nidSize
               << " nid= " << (nidSize ? (int)m_eMesh.identifier(nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>()[0]) : -1);

          //m_eMesh.print(sout, element,false,true);
          sout << " smark= " << smark << " gotMark= " << gotMark << " needNodes= " << needNodes << " isG= " << m_eMesh.isGhostElement(element)
               << " is_empty= " << is_empty
               << " a_is_empty= " << a_is_empty
               << " should_put_in= " << should_put_in
               << " should_put_in_id= " << should_put_in_id
               << " should_put_in_rank_gt= " << should_put_in_rank_gt
               << " should_put_in_rank_gte= " << should_put_in_rank_gte
               << " needed_entity_rank= "
               << needed_entity_rank.first << " subDimEntity= ";

          for (unsigned k=0; k < subDimEntity.size(); k++)
            {
              sout << " " << m_eMesh.identifier(subDimEntity[k]) << " ";
            }
          std::cout << sout.str() << std::endl;
        }

      return ret_val;
    }


    void NodeRegistry::query(stk::mesh::EntityId elementId, unsigned rank, unsigned iSubDimOrd, std::string msg)
    {
      stk::mesh::Entity element = m_eMesh.get_entity(m_eMesh.element_rank(), elementId);
      if (!m_eMesh.is_valid(element))
        {
          std::cout  << "P[" << m_eMesh.get_rank() << "] query:: " << msg << " element is not valid= " << elementId << std::endl;
          return;
        }
      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      getSubDimEntity(subDimEntity, element, rank, iSubDimOrd);

      static SubDimCellData new_SubDimCellData;
      static SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;

      unsigned smark=0;
      if (!is_empty)
        {
          unsigned& mark = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().m_mark;
          smark = mark;
        }

      if (DEBUG_NR_DEEP) // && (m_eMesh.identifier(subDimEntity[0]) == 19 && m_eMesh.identifier(subDimEntity[1]) == 20))
        {
          unsigned& gotMark = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().m_mark;
          unsigned nidSize = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().size();

          std::ostringstream sout;
          sout << "P[" << m_eMesh.get_rank() << "] query:: " << msg << " element= " << m_eMesh.identifier(element) << " nidSize= " << nidSize
               << " nid= " << (nidSize ? (int)m_eMesh.identifier(nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>()[0]) : -1);
          sout << " smark= " << smark << " gotMark= " << gotMark << " isG= " << m_eMesh.isGhostElement(element)
               << " is_empty= " << is_empty << " subDimEntity= ";

          for (unsigned k=0; k < subDimEntity.size(); k++)
            {
              sout << " " << m_eMesh.identifier(subDimEntity[k]) << " ";
            }
#if 0
          for (unsigned k=0; k < subDimEntity.size(); k++)
            {
              std::cout << " " << (stk::mesh::EntityId) subDimEntity[k] << " ";
            }
#endif
          std::cout << sout.str() << std::endl;
        }
    }

    /// Replace element ownership - FIXME for parallel
    /// When remeshing during unrefinement, replace ownership of sub-dim entities by non-deleted elements
    bool NodeRegistry::replaceElementOwnership(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes)
    {
      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      noInline_getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);

      static SubDimCellData new_SubDimCellData;
      static SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;
      bool is_not_empty_but_data_cleared = (!is_empty && nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>().size() == 0);
      (void)is_not_empty_but_data_cleared;

      // if empty or if my id is the smallest, make this element the owner
      stk::mesh::EntityId db_id = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
      stk::mesh::EntityId db_rank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

      bool should_put_in_id = (m_eMesh.identifier(element)  < db_id);
      bool should_put_in_rank_gt = (m_eMesh.entity_rank(element) > db_rank);
      bool should_put_in_rank_gte = (m_eMesh.entity_rank(element) >= db_rank);
      bool should_put_in = should_put_in_rank_gt || (should_put_in_id && should_put_in_rank_gte);
      if (db_id == 0u)
        should_put_in = true;

      if (!is_empty && should_put_in)
        {
          nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = stk::mesh::EntityKey(m_eMesh.entity_rank(element), m_eMesh.identifier(element));
          return true;
        }

      return false;
    }

    /// check the newly registered node from the registry, which does one of three things, depending on what mode we are in:
    ///   1. counts buffer in prep for sending (just does a pack)
    ///   2. packs the buffer (after buffers are alloc'd)
    ///   3. returns the new node after all communications are done
    bool NodeRegistry::checkForRemote(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes_notUsed)
    {
      EXCEPTWATCH;
      static SubDimCellData empty_SubDimCellData;
      static CommDataType buffer_entry;

      bool isGhost = m_eMesh.isGhostElement(element);
      if (!isGhost) return true;

      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank.first, iSubDimOrd);

      stk::CommAll& comm_all = *m_comm_all;
      unsigned proc_rank = comm_all.parallel_rank();
      unsigned owner_proc_rank = m_eMesh.owner_rank(element);

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;

      if (is_empty)
        {
          std::cout << "element= " << element
                    << " needed_entity_rank= " << needed_entity_rank.first<< " " << needed_entity_rank.second << std::endl;
          std::cout << "subDimEntity= " << subDimEntity << std::endl;
          std::cout << "nodeId_elementOwnderId= " << nodeId_elementOwnderId << std::endl;
          std::cout << "empty_SubDimCellData= " << empty_SubDimCellData << std::endl;
          throw std::logic_error("NodeRegistry::checkForRemote no data (is_empty=true) - logic error.");
          return false;
        }
      else
        {
          stk::mesh::EntityId owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          stk::mesh::EntityRank owning_elementRank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          unsigned nidsz = nodeIds_onSE.size();

          // error check
          bool isNotOK = (m_eMesh.identifier(element) < owning_elementId) && (m_eMesh.entity_rank(element) > owning_elementRank) ;

          if ( isNotOK )
            {
              std::cout << "P[" << proc_rank << "] elem id = " << m_eMesh.identifier(element)
                        << " nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = "
                        << owning_elementId
                        << std::endl;
              throw std::logic_error("NodeRegistry::checkForRemote logic: owning element info is wrong");
            }

          stk::mesh::EntityRank erank = owning_elementRank;
          stk::mesh::Entity owning_element = get_entity_element(*m_eMesh.get_bulk_data(), erank, owning_elementId);

          if (!m_eMesh.is_valid(owning_element))
            {
              std::cout << "owning_elementId = " << owning_elementId << " erank = " << erank << std::endl;
              throw std::logic_error("NodeRegistry::checkForRemote logic: owning_element is null");
            }

          bool owning_element_is_ghost = m_eMesh.isGhostElement(owning_element);

          // if this element is a ghost, and the owning element of the node is not a ghost, send info
          //   to ghost element's owner proc
          if (!owning_element_is_ghost && isGhost)
            {
              buffer_entry = CommDataType(
                                          needed_entity_rank.first,
                                          iSubDimOrd,
                                          m_eMesh.key(element)
                                          );

              if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
                {
                  throw std::logic_error("NodeRegistry::checkForRemote logic err #0.1");
                }

              for (unsigned iid = 0; iid < nidsz; iid++)
                {
                  if (!m_eMesh.is_valid(nodeIds_onSE[iid]))
                    {
                      throw std::logic_error("logic: hmmm #5.0");
                    }

                  if (!nodeIds_onSE.m_entity_id_vector[iid])
                    {
                      throw std::logic_error("NodeRegistry::checkForRemote logic err #0.2");
                    }

                  //stk::mesh::Entity new_node = get_entity_node_Ia(*m_eMesh.get_bulk_data(), Node, nodeIds_onSE, iid);
                  stk::mesh::Entity new_node = nodeIds_onSE[iid];
                  if (!m_eMesh.is_valid(new_node))
                    {
                      throw std::logic_error("NodeRegistry::checkForRemote logic: new_node is null");
                    }

                  if (m_eMesh.owner_rank(new_node) == m_eMesh.get_rank())
                    m_nodes_to_ghost.push_back( stk::mesh::EntityProc(new_node, owner_proc_rank) );
                }

              {
                m_comm_all->send_buffer( owner_proc_rank ).pack< CommDataType > (buffer_entry);
                NodeIdsOnSubDimEntityType& nids = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                nids.pack(m_eMesh, m_comm_all->send_buffer( owner_proc_rank ));
              }
            }
        }
      return true;
    }

    bool NodeRegistry::getFromRemote(const stk::mesh::Entity element, NeededEntityType& needed_entity_rank, unsigned iSubDimOrd, bool needNodes_notUsed)
    {
      return checkForRemote(element, needed_entity_rank, iSubDimOrd, needNodes_notUsed);
    }


    /// makes coordinates of this new node be the centroid of its sub entity
    void NodeRegistry::prolongateCoords(const stk::mesh::Entity element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
    {
      prolongateField(element, needed_entity_rank, iSubDimOrd, m_eMesh.get_coordinates_field());
    }

    void NodeRegistry::prolongateField(const stk::mesh::Entity element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd, stk::mesh::FieldBase *field)
    {
      //EXCEPTWATCH;

      int spatialDim = m_eMesh.get_spatial_dim();
      stk::mesh::EntityRank field_rank = stk::mesh::MetaData::NODE_RANK;
      {
        unsigned nfr = field->restrictions().size();
        for (unsigned ifr = 0; ifr < nfr; ifr++)
          {
            const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
            field_rank = fr . entity_rank();
            spatialDim = fr.dimension() ;
          }
      }

      if (field_rank != stk::mesh::MetaData::NODE_RANK)
        {
          return;
        }

      unsigned *null_u = 0;
      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
      static SubDimCellData empty_SubDimCellData;

      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;

      if (s_allow_empty_sub_dims && is_empty)
        {
          return;
        }
      if (is_empty)
        {
          const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(*m_eMesh.get_bulk_data(), element);
          shards::CellTopology cell_topo(cell_topo_data);

          std::cout << "NodeRegistry::prolongateField: no node found, cell_topo = " << cell_topo.getName()
                    << "\n subDimEntity= " << subDimEntity
                    << "\n element= " << element
                    << "\n m_eMesh.entity_rank(element) = " << m_eMesh.entity_rank(element)
                    << "\n needed_entity_rank= " << needed_entity_rank
                    << "\n iSubDimOrd= " << iSubDimOrd << std::endl;
          throw std::runtime_error("prolongateField: no node found");
        }
      NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
      if (nodeIds_onSE.size() != 1)
        {
          throw std::logic_error("logic error: prolongateField not ready for multiple nodes, or there are 0 nodes on the marked quantity");
        }
      stk::mesh::Entity c_node = nodeIds_onSE[0];

      if (!m_eMesh.is_valid(c_node))
        {
          throw std::runtime_error("prolongateField: bad node found 0");
        }

      double c_p[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      bool doPrint = false;

      if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
        {
          const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::mesh::MetaData::NODE_RANK);
          unsigned npts = elem_nodes.size();
          double dnpts = elem_nodes.size();
          for (unsigned ipts = 0; ipts < npts; ipts++)
            {
              stk::mesh::Entity node = elem_nodes[ipts].entity();
              if (!m_eMesh.is_valid(node))
                {
                  throw std::runtime_error("prolongateField: bad node found 1");
                }
              double *  coord = m_eMesh.field_data(field, node, null_u);

              if (doPrint && coord)
                {
                  const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
                  shards::CellTopology cell_topo(cell_topo_data);

                  std::cout << "tmp NodeRegistry::prolongateField cell_topo = " << cell_topo.getName() << " ipts= " << ipts
                            << " coord= " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
                }

              if (coord)
                {
                  for (int isp = 0; isp < spatialDim; isp++)
                    {
                      c_p[isp] += coord[isp]/dnpts;
                    }
                }
            }

        }
      else
        {
          double dnpts = subDimEntity.size();
          for (SubDimCell_SDCEntityType::iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ids++)
            {
              SDCEntityType nodeId = *ids;

              stk::mesh::Entity node = nodeId;
              if (!m_eMesh.is_valid(node))
                {
                  throw std::runtime_error("prolongateField: bad node found 2");
                }
              double *  coord = m_eMesh.field_data(field, node, null_u);
              if (coord)
                {
                  for (int isp = 0; isp < spatialDim; isp++)
                    {
                      c_p[isp] += coord[isp]/dnpts;
                    }
                }
            }
        }
      double *  c_coord = m_eMesh.field_data(field, c_node, null_u);
      if (c_coord)
        {
          for (int isp = 0; isp < spatialDim; isp++)
            {
              c_coord[isp] = c_p[isp];
            }
        }
    }

    /// do interpolation for all fields
    void NodeRegistry::prolongateFields(const stk::mesh::Entity element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
    {
      const stk::mesh::FieldVector & fields = m_eMesh.get_fem_meta_data()->get_fields();
      unsigned nfields = fields.size();
      for (unsigned ifld = 0; ifld < nfields; ifld++)
        {
          stk::mesh::FieldBase *field = fields[ifld];
          const stk::mesh::DataTraits & data_traits = field->data_traits();
          if (data_traits.is_floating_point)
            prolongateField(element, needed_entity_rank, iSubDimOrd, field);
          else
            {
              std::cout << "Not prolonging possible integer Field: " << field->name() << " field= " << typeid(*field).name()
                        <<"\n is_void = " << data_traits.is_void
                        <<"\n is_integral = " << data_traits.is_integral
                        <<"\n is_floating_point = " << data_traits.is_floating_point
                        <<"\n is_array = " << data_traits.is_array
                        <<"\n is_pointer = " << data_traits.is_pointer
                        <<"\n is_enum = " << data_traits.is_enum
                        <<"\n is_class = " << data_traits.is_class
                        <<"\n is_pod = " << data_traits.is_pod
                        <<"\n is_signed  = " << data_traits.is_signed
                        <<"\n is_unsigned = " << data_traits.is_unsigned
                        <<"\n alignment_of = " << data_traits.alignment_of
                        <<"\n stride_of = " << data_traits.stride_of
                        <<"\n name = " << data_traits.name
                        << std::endl;
            }

        }
    }

    /// do interpolation for all fields
    void NodeRegistry::prolongateFields()
    {
      const stk::mesh::FieldVector & fields = m_eMesh.get_fem_meta_data()->get_fields();
      unsigned nfields = fields.size();
      for (unsigned ifld = 0; ifld < nfields; ifld++)
        {
          stk::mesh::FieldBase *field = fields[ifld];
          // only do it for reals - integers are not really defined as being able to be prolonged
          const stk::mesh::DataTraits & data_traits = field->data_traits();
          if (data_traits.is_floating_point)
            prolongate(field);
          else
            {
              std::cout << "Not prolonging possible integer Field: " << field->name() << " field= " << typeid(*field).name()
                        <<"\n is_void = " << data_traits.is_void
                        <<"\n is_integral = " << data_traits.is_integral
                        <<"\n is_floating_point = " << data_traits.is_floating_point
                        <<"\n is_array = " << data_traits.is_array
                        <<"\n is_pointer = " << data_traits.is_pointer
                        <<"\n is_enum = " << data_traits.is_enum
                        <<"\n is_class = " << data_traits.is_class
                        <<"\n is_pod = " << data_traits.is_pod
                        <<"\n is_signed  = " << data_traits.is_signed
                        <<"\n is_unsigned = " << data_traits.is_unsigned
                        <<"\n alignment_of = " << data_traits.alignment_of
                        <<"\n stride_of = " << data_traits.stride_of
                        <<"\n name = " << data_traits.name
                        << std::endl;
            }
        }
    }

    /// check for adding new nodes to existing parts based on sub-entity part ownership
    void NodeRegistry::addToExistingParts(const stk::mesh::Entity element,  stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
    {
      const std::vector< stk::mesh::Part * > & parts = m_eMesh.get_fem_meta_data()->get_parts();

      unsigned nparts = parts.size();

      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
      static  SubDimCellData empty_SubDimCellData;
      SubDimCellData* nodeId_elementOwnderId_ptr = getFromMapPtr(subDimEntity);
      SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
      bool is_empty = nodeId_elementOwnderId_ptr == 0;

      if (s_allow_empty_sub_dims && is_empty)
        {
          return;
        }

      if (is_empty)
        {
          throw std::runtime_error("addToExistingParts: no node found");
        }
      NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
      unsigned nidsz = nodeIds_onSE.size();

      for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
        {
          stk::mesh::Entity c_node = nodeIds_onSE[i_nid];

          if (!m_eMesh.is_valid(c_node))
            {
              std::cout << "addToExistingParts: " <<  nodeIds_onSE[i_nid] << " i_nid= " << i_nid << " nidsz= " << nidsz
                        << " needed_entity_rank= " << needed_entity_rank << " iSubDimOrd= " << iSubDimOrd << std::endl;
              throw std::runtime_error("addToExistingParts: bad node found 0.1");
            }

          for (unsigned ipart=0; ipart < nparts; ipart++)
            {
              stk::mesh::Part& part = *parts[ipart];
              stk::mesh::Selector selector(part);

              if (stk::mesh::is_auto_declared_part(part))
                continue;

              bool found = true;
              if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
                {
                  const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::mesh::MetaData::NODE_RANK);
                  unsigned npts = elem_nodes.size();
                  for (unsigned ipts = 0; ipts < npts; ipts++)
                    {
                      stk::mesh::Entity node = elem_nodes[ipts].entity();
                      if (!m_eMesh.is_valid(node))
                        {
                          throw std::runtime_error("addToExistingParts: bad node found 1.1");
                        }
                      if (!selector(m_eMesh.bucket(node)))
                        {
                          found = false;
                          break;
                        }
                    }
                }
              else
                {
                  for (SubDimCell_SDCEntityType::iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                    {
                      SDCEntityType nodeId = *ids;
                      stk::mesh::Entity node = nodeId;
                      if (!m_eMesh.is_valid(node))
                        {
                          throw std::runtime_error("addToExistingParts: bad node found 2.1");
                        }
                      if (!selector(m_eMesh.bucket(node)))
                        {
                          found = false;
                          break;
                        }
                    }
                }
              if (found)
                {
                  // add to part
                  std::vector<stk::mesh::Part*> add_parts(1, &part);
                  std::vector<stk::mesh::Part*> remove_parts;
                  const unsigned part_rank = part.primary_entity_rank();

                  if (part_rank == stk::mesh::MetaData::NODE_RANK)
                    {
                      if (m_eMesh.bucket(c_node).owned())
                        m_eMesh.get_bulk_data()->change_entity_parts( c_node, add_parts, remove_parts );
                    }
                }
            }
        }
    }

    /// Check for adding new nodes to existing parts based on sub-entity part ownership.
    /// This version does it in bulk and thus avoids repeats on shared sub-dim entities.
    void NodeRegistry::addToExistingPartsNew()
    {
      static std::vector<stk::mesh::Part*> add_parts(1, static_cast<stk::mesh::Part*>(0));
      static std::vector<stk::mesh::Part*> remove_parts;

      const std::vector< stk::mesh::Part * > & parts = m_eMesh.get_fem_meta_data()->get_parts();

      unsigned nparts = parts.size();
      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];

          if (stk::mesh::is_auto_declared_part(part))
            continue;

          const unsigned part_rank = part.primary_entity_rank();

          if (part_rank == stk::mesh::MetaData::NODE_RANK)
            {
              stk::mesh::Selector selector(part);

              add_parts[0] = &part;

              SubDimCellToDataMap::iterator iter;

              for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
                {
                  const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
                  SubDimCellData& nodeId_elementOwnderId = (*iter).second;

                  NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

                  bool found = true;
                  stk::mesh::EntityRank needed_entity_rank = stk::mesh::MetaData::NODE_RANK;

                  // SPECIAL CASE
                  if( subDimEntity.size() == 1)
                    {
                      needed_entity_rank = stk::mesh::MetaData::ELEMENT_RANK;
                    }

                  if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
                    {
                      stk::mesh::Entity element_p = stk::mesh::Entity();
                      {
                        SDCEntityType elementId = subDimEntity[0];
                        element_p = elementId;
                        if (!m_eMesh.is_valid(element_p))
                          {
                            stk::mesh::EntityId        owning_elementId    = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
                            stk::mesh::EntityRank      owning_element_rank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
                            stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_element_rank, owning_elementId);
                            if (!m_eMesh.isGhostElement(owning_element))
                              {
                                std::cout << "subDimEntity size= " << subDimEntity.size() << std::endl;
                                VERIFY_OP_ON(m_eMesh.is_valid(element_p), ==, true, "Bad elem found 2 in addToExistingPartsNew 2");
                              }
                            continue;
                          }
                      }

                      stk::mesh::Entity element = element_p;

                      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::mesh::MetaData::NODE_RANK);
                      unsigned npts = elem_nodes.size();
                      for (unsigned ipts = 0; ipts < npts; ipts++)
                        {
                          stk::mesh::Entity node = elem_nodes[ipts].entity();
                          if (!m_eMesh.is_valid(node))
                            {
                              std::cout << "subDimEntity size= " << subDimEntity.size() << std::endl;
                            }
                          VERIFY_OP_ON(m_eMesh.is_valid(node), ==, true, "Bad node in addToExistingPartsNew 1");
                          if (!selector(m_eMesh.bucket(node)))
                            {
                              found = false;
                              break;
                            }
                        }
                    }
                  else
                    {
                      for (SubDimCell_SDCEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                        {
                          SDCEntityType nodeId = *ids;
                          stk::mesh::Entity node = nodeId;
                          if (!m_eMesh.is_valid(node))
                            {
                              stk::mesh::EntityId        owning_elementId    = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
                              stk::mesh::EntityRank      owning_element_rank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
                              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_element_rank, owning_elementId);
                              if (!m_eMesh.isGhostElement(owning_element))
                                {
                                  std::cout << "subDimEntity size= " << subDimEntity.size() << std::endl;
                                  VERIFY_OP_ON(m_eMesh.is_valid(node), ==, true, "Bad node in addToExistingPartsNew 2");
                                }
                              continue;
                            }
                          if (!selector(m_eMesh.bucket(node)))
                            {
                              found = false;
                              break;
                            }
                        }
                    }
                  if (found)
                    {
                      // add to part
                      unsigned nidsz = nodeIds_onSE.size();

                      for (unsigned i_nid = 0; i_nid < nidsz; i_nid++)
                        {
                          stk::mesh::Entity c_node = nodeIds_onSE[i_nid];

                          if (!m_eMesh.is_valid(c_node))
                            {
                              // note, this is ok - a null node here can come from a ghost element
                              continue;
                            }

                          // only try to add to part if I am the owner
                          if (m_eMesh.owner_rank(c_node) == m_eMesh.get_parallel_rank())
                            m_eMesh.get_bulk_data()->change_entity_parts( c_node, add_parts, remove_parts );
                        }
                    }
                }
            }
        }
      //std::cout << "tmp addToExistingPartsNew...done " << std::endl;
    }

    void NodeRegistry::
    doForAllSubEntities(ElementFunctionPrototype function, const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      shards::CellTopology cell_topo(cell_topo_data);

      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;
          stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_rank == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
            }
          else if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              (this ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true);

            } // iSubDimOrd
        } // ineed_ent
    }

    void NodeRegistry::
    getSubDimEntity(SubDimCell_SDCEntityType& subDimEntity, const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
    {
      subDimEntity.clear();
      // in the case of elements, we don't share any nodes so we just make a map of element id to node
      if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
        {
          subDimEntity.insert( element );
          return;
        }

      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      stk::mesh::Entity const * const elem_nodes = m_eMesh.get_bulk_data()->begin_nodes(element);

      const unsigned *  inodes = 0;
      unsigned nSubDimNodes = 0;
      static const unsigned edge_nodes_2[2] = {0,1};
      static const unsigned face_nodes_3[3] = {0,1,2};
      static const unsigned face_nodes_4[4] = {0,1,2,3};

      // special case for faces in 3D
      if (needed_entity_rank == m_eMesh.face_rank() && needed_entity_rank == m_eMesh.entity_rank(element))
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          if (nSubDimNodes ==3 )
            inodes = face_nodes_3;
          else
            inodes = face_nodes_4;

        }
      // special case for edges in 2D
      else if (needed_entity_rank == m_eMesh.edge_rank() && needed_entity_rank == m_eMesh.entity_rank(element))
        {
          nSubDimNodes = cell_topo_data->vertex_count;

          if (nSubDimNodes == 2 )
            {
              inodes = edge_nodes_2;
            }
          else
            {
              throw std::runtime_error("NodeRegistry bad for edges");
            }
        }
      else if (needed_entity_rank == m_eMesh.edge_rank())
        {
          inodes = cell_topo_data->edge[iSubDimOrd].node;
          nSubDimNodes = 2;
        }
      else if (needed_entity_rank == m_eMesh.face_rank())
        {
          nSubDimNodes = cell_topo_data->side[iSubDimOrd].topology->vertex_count;
          // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
          inodes = cell_topo_data->side[iSubDimOrd].node;
        }

      //subDimEntity.reserve(nSubDimNodes);
      for (unsigned jnode = 0; jnode < nSubDimNodes; jnode++)
        {
          subDimEntity.insert( elem_nodes[inodes[jnode]] );
        }
    }

    unsigned NodeRegistry::total_size()
    {
      unsigned sz=0;

      for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
        {
          SubDimCellData& data = (*cell_iter).second;
          NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();

          sz += nodeIds_onSE.size();
        }
      return sz;
    }

    unsigned NodeRegistry::local_size()
    {
      unsigned sz=0;
      for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
        {
          SubDimCellData& data = (*cell_iter).second;

          stk::mesh::EntityId owning_elementId = data.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();

          NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();
          if (nodeIds_onSE.size())
            {
              stk::mesh::EntityRank erank = data.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
              stk::mesh::Entity owning_element = get_entity_element(*m_eMesh.get_bulk_data(), erank, owning_elementId);

              if (!m_eMesh.is_valid(owning_element))
                {
                  std::cout << "tmp owning_element = null, owning_elementId= " << owning_elementId
                            << std::endl;
                  throw std::logic_error("logic: hmmm #5.2");
                }
              if (!m_eMesh.isGhostElement(owning_element))
                {
                  sz += nodeIds_onSE.size();
                }
            }
        }
      return sz;
    }

    //========================================================================================================================
    // low-level interface
    // FIXME

    void NodeRegistry::checkDB(std::string msg)
    {
      if (DEBUG_NR_DEEP)
        {
          //query(265, 1, 2, msg);
          //query(253, 1, 0, msg);
        }

      if (0)
        {
          unsigned sz=0;
          for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
            {
              SubDimCellData& data = (*cell_iter).second;
              stk::mesh::EntityId owning_elementId = data.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
              NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();

              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(stk::mesh::MetaData::ELEMENT_RANK, owning_elementId);
              if (!m_eMesh.is_valid(owning_element))
                throw std::logic_error("logic: hmmm #5.3");
              bool isGhost = m_eMesh.isGhostElement(owning_element);
              if (!m_eMesh.isGhostElement(owning_element))
                {
                  ++sz;
                }
              std::ostringstream str;
              const SubDimCell_SDCEntityType& subDimEntity = (*cell_iter).first;
              str << "subDimEntity.size() = " << subDimEntity.size();
              for (unsigned i=0; i < subDimEntity.size(); i++)
                {
                  str << " s[" << i << "] = " << m_eMesh.identifier(subDimEntity[i]);
                }

              std::cout << "P[" << m_eMesh.get_rank() << "] owning_elementId = "  << owning_elementId << " isGhostElement = " << isGhost
                        << " nodeId = " << nodeIds_onSE << " subDimEntity= " << str.str() << std::endl;
            }
        }
      if (0)
        {
          std::cout << "P[" << m_eMesh.get_rank() << "] NodeRegistry::checkDB start msg= " << msg << " m_cell_2_data_map.size= " << m_cell_2_data_map.size() << std::endl;
          for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
            {
              const SubDimCell_SDCEntityType& subDimEntity = (*cell_iter).first;
              SubDimCellData&            subDimCellData      = (*cell_iter).second;
              stk::mesh::EntityId        owning_elementId    = subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
              stk::mesh::EntityRank      owning_element_rank = subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_element_rank, owning_elementId);

              NodeIdsOnSubDimEntityType& nodeIds_onSE        = subDimCellData.get<SDC_DATA_GLOBAL_NODE_IDS>();

              int count=0;
              for (SubDimCell_SDCEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids, ++count)
                {
                  SDCEntityType nodeId = *ids;
                  stk::mesh::Entity node = nodeId;
                  if (!m_eMesh.is_valid(node))
                    {
                      std::cout << "subDimEntity size= " << subDimEntity.size() << " count= " << count << std::endl;
                    }
                  VERIFY_OP_ON(m_eMesh.is_valid(node), ==, true, "Bad node in checkDB 1");
                }

              for (unsigned i=0; i < nodeIds_onSE.size(); i++)
                {
                  stk::mesh::Entity node = nodeIds_onSE[i];
                  stk::mesh::EntityId nodeId = nodeIds_onSE.m_entity_id_vector[i];
                  if (!m_eMesh.is_valid(node)) {
                    std::cout << "i= " << i << " size= " << nodeIds_onSE.size()
                              << " node= " << node << " nodeId= " << nodeId
                              << " subDimEntity.size= " << subDimEntity.size()
                              << " owning_elementId = " << owning_elementId
                              << " owning_element_rank = " << owning_element_rank
                              << " owning_element = " << owning_element
                              << " owning_element.isGhostElement= " << m_eMesh.isGhostElement(owning_element)
                              << std::endl;
                  }
                  VERIFY_OP_ON(node, !=, stk::mesh::Entity(), "checkDB #11.1");
                  VERIFY_OP_ON(nodeId, !=, 0, "checkDB #11.1.1");
                  VERIFY_OP_ON(m_eMesh.identifier(node), ==, nodeId, "checkDB #11.2");
                  stk::mesh::Entity node_0 = m_eMesh.get_bulk_data()->get_entity(0, nodeId);

                  VERIFY_OP_ON(node, ==, node_0, "checkDB #11.3");
                  VERIFY_OP_ON(m_eMesh.identifier(node_0), ==, nodeId, "checkDB #11.4");
                }
            }
          std::cout << "P[" << m_eMesh.get_rank() << "] NodeRegistry::checkDB end msg= " << msg << std::endl;
        }

      if (1)
        {
          std::cout << "P[" << m_eMesh.get_rank() << "] NodeRegistry::checkDB start 1 msg= " << msg << std::endl;
          for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
            {
              const SubDimCell_SDCEntityType& subDimEntity = (*cell_iter).first;
              SubDimCellData&            subDimCellData      = (*cell_iter).second;
              stk::mesh::EntityId        owning_elementId    = subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
              stk::mesh::EntityRank      owning_element_rank = subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
              NodeIdsOnSubDimEntityType& nodeIds_onSE        = subDimCellData.get<SDC_DATA_GLOBAL_NODE_IDS>();

              VERIFY_OP_ON(owning_elementId , !=, 0, "owning_elementId=0");
              if (owning_elementId == 0) continue;
              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_element_rank, owning_elementId);
              if (!m_eMesh.is_valid(owning_element))
                {
                  std::cout << "owning_element is invalid, msg= " << msg << " owning_elementId= " << owning_elementId << " owning_element_rank= " << owning_element_rank
                            << " subDimCellData.size=" << subDimEntity.size()
                            << std::endl;
                  throw std::logic_error("logic: checkDB hmmm #11.0");
                }
              //bool isGhost = m_eMesh.isGhostElement(*owning_element);
              if (0 && !m_eMesh.isGhostElement(owning_element))
                {
                  for (unsigned i=0; i < nodeIds_onSE.size(); i++)
                    {
                      stk::mesh::Entity node = nodeIds_onSE[i];
                      stk::mesh::EntityId nodeId = nodeIds_onSE.m_entity_id_vector[i];
                      VERIFY_OP_ON(node, !=, stk::mesh::Entity(), "checkDB #11.1");
                      VERIFY_OP_ON(nodeId, !=, 0, "checkDB #11.1.1");
                      VERIFY_OP_ON(m_eMesh.identifier(node), ==, nodeId, "checkDB #11.2");
                      stk::mesh::Entity node_0 = m_eMesh.get_bulk_data()->get_entity(0, nodeId);

                      VERIFY_OP_ON(node, ==, node_0, "checkDB #11.3");
                      VERIFY_OP_ON(m_eMesh.identifier(node_0), ==, nodeId, "checkDB #11.4");
                      const int owner_rank = m_eMesh.owner_rank(node);
                      VERIFY_OP_ON(owner_rank, ==, m_eMesh.owner_rank(owning_element), "checkDB #11.6");
                    }
                }
            }
          std::cout << "P[" << m_eMesh.get_rank() << "] NodeRegistry::checkDB end 1 msg= " << msg << std::endl;

        }

    }

    /// allocate the send/recv buffers for all-to-all communication
    bool NodeRegistry::allocateBuffers()
    {
      stk::CommAll& comm_all = *m_comm_all;
      unsigned proc_size = comm_all.parallel_size();
      unsigned proc_rank = comm_all.parallel_rank();

      bool local = true;
      unsigned num_msg_bounds = proc_size < 4 ? proc_size : proc_size/4 ;
      bool global = comm_all.allocate_buffers(num_msg_bounds , false, local );
      if ( not global )
        {
          std::cout << "P[" << proc_rank << "] : not global" << std::endl;
          return false;
        }
      return true;
    }

    void NodeRegistry::communicate()
    {
      stk::CommAll& comm_all = *m_comm_all;
      comm_all.communicate();

      stk::ParallelMachine pm = m_eMesh.parallel();
      int failed = 0;
      stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );

      unpack();
    }

    void NodeRegistry::
    unpack()
    {
      stk::CommAll& comm_all = *m_comm_all;

      int failed = 0;
      std::string msg;

      stk::ParallelMachine pm = m_eMesh.parallel();
      unsigned proc_size = m_eMesh.get_parallel_size();
      unsigned proc_rank = comm_all.parallel_rank();

      vector<stk::mesh::EntityProc> nodes_to_ghost;

      CommDataType buffer_entry;
      NodeIdsOnSubDimEntityType nodeIds_onSE;
      try
        {
          for(unsigned from_proc = 0; from_proc < proc_size; ++from_proc )
            {
              stk::CommBuffer & recv_buffer = comm_all.recv_buffer( from_proc );

              while ( recv_buffer.remaining() )
                {
                  // Rank of sub-dim cells needing new nodes, which sub-dim entity, one non-owning element identifier, nodeId_elementOwnderId.first
                  recv_buffer.unpack< CommDataType >( buffer_entry );
                  nodeIds_onSE.unpack(m_eMesh, recv_buffer);
                  createNodeAndConnect(buffer_entry, nodeIds_onSE, from_proc, nodes_to_ghost);
                }
            }
        }
      catch ( std::exception &x )
        {
          failed = 1;
          msg = std::string("unpack error: ")+x.what();
        }

      stk::all_reduce( pm, stk::ReduceSum<1>( &failed ) );
      if ( failed )
        {
          throw std::runtime_error( msg+" from unpack error, rank = "+toString(proc_rank) );
        }

      if (nodes_to_ghost.size())
        {
          stk::mesh::Ghosting & ghosting = m_eMesh.get_bulk_data()->create_ghosting( std::string("new_nodes") );

          vector<stk::mesh::EntityKey> receive;
          ghosting.receive_list( receive );
          m_eMesh.get_bulk_data()->change_ghosting( ghosting, nodes_to_ghost, receive);
        }

    }// unpack


    /// after registering all needed nodes, this method is used to request new nodes on this processor
    void NodeRegistry::createNewNodesInParallel()
    {
      static int num_times_called = 0;
      ++num_times_called;

      stk::mesh::Part* new_nodes_part = m_eMesh.get_non_const_part("refine_new_nodes_part");
      if (new_nodes_part)
        {
          // first remove any existing nodes
          std::vector<stk::mesh::Part*> remove_parts(1, new_nodes_part);
          std::vector<stk::mesh::Part*> add_parts;
          std::vector<stk::mesh::Entity> node_vec;

          stk::mesh::Selector removePartSelector(*new_nodes_part & m_eMesh.get_fem_meta_data()->locally_owned_part() );
          const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              if (removePartSelector(bucket))
                {
                  const unsigned num_entity_in_bucket = bucket.size();
                  for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                    {
                      stk::mesh::Entity node = bucket[ientity];
                      node_vec.push_back(node);
                    }
                }
            }
          for (unsigned ii=0; ii < node_vec.size(); ii++)
            {
              m_eMesh.get_bulk_data()->change_entity_parts( node_vec[ii], add_parts, remove_parts );
            }
        }

      unsigned num_nodes_needed = local_size();

      // assert( bulk data is in modifiable mode)
      // create new entities on this proc
      vector<stk::mesh::Entity> new_nodes;

      if (m_useCustomGhosting)
        {
          m_eMesh.createEntities( stk::mesh::MetaData::NODE_RANK, num_nodes_needed, new_nodes);
        }

      std::vector<stk::mesh::EntityId> ids(num_nodes_needed);

      if (new_nodes_part)
        {
          stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->locally_owned_part() );
          std::vector<stk::mesh::Part*> add_parts(1, new_nodes_part);
          std::vector<stk::mesh::Part*> remove_parts;
          for (unsigned ind = 0; ind < new_nodes.size(); ind++)
            {
              m_eMesh.get_bulk_data()->change_entity_parts( new_nodes[ind], add_parts, remove_parts );
            }
        }

      // set map values to new node id's
      unsigned inode=0;

      for (SubDimCellToDataMap::iterator cell_iter = m_cell_2_data_map.begin(); cell_iter != m_cell_2_data_map.end(); ++cell_iter)
        {
          SubDimCellData& data = (*cell_iter).second;
          NodeIdsOnSubDimEntityType& nodeIds_onSE = data.get<SDC_DATA_GLOBAL_NODE_IDS>();
          if (!nodeIds_onSE.size())
            continue;

          stk::mesh::EntityId owning_elementId = data.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();

          if (!owning_elementId)
            {
              throw std::logic_error("logic: hmmm #5.4.0");
            }

          stk::mesh::EntityRank erank = data.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();
          stk::mesh::Entity owning_element = get_entity_element(*m_eMesh.get_bulk_data(), erank, owning_elementId);

          if (!m_eMesh.is_valid(owning_element))
            {
              throw std::logic_error("logic: hmmm #5.4");
            }
          if (!m_eMesh.isGhostElement(owning_element))
            {
              if (nodeIds_onSE.m_entity_id_vector.size() != nodeIds_onSE.size())
                {
                  throw std::logic_error("NodeRegistry:: createNewNodesInParallel logic err #0.0");
                }

              for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
                {
                  VERIFY_OP(inode, < , num_nodes_needed, "UniformRefiner::doBreak() too many nodes");
                  if ( DEBUG_NR_UNREF)
                    {
                      std::cout << "tmp createNewNodesInParallel: old node id= " << (m_eMesh.is_valid(nodeIds_onSE[ii]) ? toString(m_eMesh.identifier(nodeIds_onSE[ii])) : std::string("null")) << std::endl;
                      std::cout << "tmp createNewNodesInParallel: new node=";
                      m_eMesh.print_entity(std::cout, new_nodes[inode]);
                    }

                  // if already exists from a previous iteration/call to doBreak, don't reset it and just use the old node
                  if (m_eMesh.is_valid(nodeIds_onSE[ii]))
                    {
                      if (DEBUG_NR_UNREF)
                        {
                          std::cout << "tmp createNewNodesInParallel: old node id is no-null, re-using it= " << (m_eMesh.is_valid(nodeIds_onSE[ii]) ? toString(m_eMesh.identifier(nodeIds_onSE[ii])) : std::string("null")) << std::endl;
                          std::cout << "tmp createNewNodesInParallel: new node=";
                          m_eMesh.print_entity(std::cout, new_nodes[inode]);
                        }
                    }
                  else
                    {
                      nodeIds_onSE[ii] = new_nodes[inode];
                      nodeIds_onSE.m_entity_id_vector[ii] = m_eMesh.identifier(new_nodes[inode]);
                    }

                  inode++;
                }
            }
        }
    }


    /// unpacks the incoming information in @param buffer_entry and adds that information to my local node registry
    /// (i.e. the map of sub-dimensional entity to global node id is updated)
    void NodeRegistry::
    createNodeAndConnect(CommDataType& buffer_entry, NodeIdsOnSubDimEntityType& nodeIds_onSE, unsigned from_proc, vector<stk::mesh::EntityProc>& nodes_to_ghost)
    {
      stk::mesh::EntityRank& needed_entity_rank                    = buffer_entry.get<CDT_NEEDED_ENTITY_RANK>();
      unsigned    iSubDimOrd                                       = buffer_entry.get<CDT_SUB_DIM_ENTITY_ORDINAL>();
      stk::mesh::EntityKey&  non_owning_elementKey                 = buffer_entry.get<CDT_NON_OWNING_ELEMENT_KEY>();
      stk::mesh::EntityId    non_owning_elementId                  = non_owning_elementKey.id();
      stk::mesh::EntityRank  non_owning_elementRank                = non_owning_elementKey.rank();

      // create a new relation here?  no, we are going to delete this element, so we just register that the new node is attached to
      stk::mesh::EntityRank erank = non_owning_elementRank;
      stk::mesh::Entity element = get_entity_element(*m_eMesh.get_bulk_data(), erank, non_owning_elementId);

      for (unsigned iid = 0; iid < nodeIds_onSE.size(); iid++)
        {
          stk::mesh::Entity node = nodeIds_onSE[iid];

          // has to be null, right?
          if (m_eMesh.is_valid(node))
            {
              throw std::logic_error("logic: node should be null in createNodeAndConnect");
            }
        }
      if (!m_eMesh.is_valid(element))
        {
          throw std::logic_error("logic: element shouldn't be null in createNodeAndConnect");
        }

      SubDimCell_SDCEntityType subDimEntity(m_eMesh);
      getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
      SubDimCellData& subDimCellData = getNewNodeAndOwningElement(subDimEntity);
      subDimCellData.get<SDC_DATA_GLOBAL_NODE_IDS>() = nodeIds_onSE;

#ifndef NDEBUG
      stk::mesh::EntityId owning_element_id = subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
      stk::mesh::EntityRank owning_element_rank = subDimCellData.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

      if (owning_element_rank == non_owning_elementRank)
        {
          if (owning_element_id == non_owning_elementId)
            {
              NR_PRINT(owning_element_id);
              NR_PRINT(non_owning_elementId);
              NR_PRINT(owning_element_rank);

              {
                std::cout << "P[" << m_eMesh.get_rank() << "] unpacking edge, from_proc= " << from_proc << " with new node= " << m_eMesh.identifier(nodeIds_onSE[0])
                          << " buffer_entry= " << buffer_entry
                          << " new node id= " << nodeIds_onSE.m_entity_id_vector[0]
                          << " owning_element_id= " << owning_element_id << " owning_element_rank= " << owning_element_rank
                          << " subdim.size()= " << subDimEntity.size()
                          << " edge= " << m_eMesh.identifier(subDimEntity[0])
                          << " " << m_eMesh.identifier(subDimEntity[1])
                          << std::endl;
              }
            }
          VERIFY_OP(owning_element_id, !=, non_owning_elementId, "createNodeAndConnect:: bad elem ids");
          VERIFY_OP(owning_element_id, !=, non_owning_elementId, "createNodeAndConnect:: bad elem ids");
          VERIFY_OP(owning_element_id, < , non_owning_elementId, "createNodeAndConnect:: bad elem ids 2");
        }
#endif

    }

    // remove any sub-dim entities from the map that have a node in deleted_nodes
    void NodeRegistry::cleanDeletedNodes(std::set<stk::mesh::Entity, stk::mesh::EntityLess>& deleted_nodes,
                                         std::set<stk::mesh::Entity, stk::mesh::EntityLess>& kept_nodes_orig_minus_kept_nodes,
                                         SubDimCellToDataMap& to_save,
                                         bool debug)
    {
      std::set<stk::mesh::Entity, stk::mesh::EntityLess> deleted_nodes_copy = deleted_nodes;

      if (DEBUG_NR_UNREF)
        std::cout << "tmp cleanDeletedNodes deleted_nodes size: " << deleted_nodes_copy.size() << std::endl;

      SubDimCellToDataMap::iterator iter;
      std::vector<SubDimCellToDataMap::iterator> to_delete;

      SubDimCellToDataMap& map = m_cell_2_data_map;
      if (DEBUG_NR_UNREF)
        std::cout << "tmp cleanDeletedNodes map size: " << map.size() << std::endl;

      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;

          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

          unsigned jj = 0;
          bool found = false;
          for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
            {
              if (deleted_nodes.find(nodeIds_onSE[ii]) != deleted_nodes.end())
                {
                  if (kept_nodes_orig_minus_kept_nodes.find(nodeIds_onSE[ii]) != kept_nodes_orig_minus_kept_nodes.end())
                    {
                      to_save[subDimEntity] = nodeId_elementOwnderId;
                    }
                  found = true;
                  jj = ii;
                  deleted_nodes_copy.erase(nodeIds_onSE[ii]);
                  break;
                }
            }
          if (found)
            {
              if (DEBUG_NR_UNREF)
                {
                  std::cout << "tmp cleanDeletedNodes:: removing node id= " << m_eMesh.identifier(nodeIds_onSE[jj]) << std::endl;
                  std::cout << "Node: ";
                  m_eMesh.print_entity(std::cout, nodeIds_onSE[jj]);
                }
              if (!debug)
                {
                  to_delete.push_back(iter);
                }
            }
        }

      //std::cout << "tmp cleanDeletedNodes to_delete.size()= " << to_delete.size() << " map.size()= " << map.size() << std::endl;
      for (unsigned itd=0; itd < to_delete.size(); itd++)
        {
          map.erase(to_delete[itd]);
        }

      if (DEBUG_NR_UNREF && deleted_nodes_copy.size())
        {
          std::cout << "tmp cleanDeletedNodes some deleted nodes not found, size()=: " << deleted_nodes_copy.size() << " nodes= " << std::endl;
          std::set<stk::mesh::Entity, stk::mesh::EntityLess>::iterator it;
          for (it = deleted_nodes_copy.begin(); it != deleted_nodes_copy.end(); ++it)
            {
              stk::mesh::Entity node = *it;
              std::cout << "Node: ";
              m_eMesh.print_entity(std::cout, node);
            }

        }
    }

    // further cleanup of the NodeRegistry db - some elements get deleted on some procs but the ghost elements
    //   are still in the db - the easiest way to detect this is as used here: during modification_begin(),
    //   cleanup of cached transactions are performed, then we find these by seeing if our element id's are
    //   no longer in the stk_mesh db.

    void NodeRegistry::clear_element_owner_data_phase_2(bool resize_nodeId_data, bool mod_begin_end)
    {
      if (mod_begin_end)
        m_eMesh.get_bulk_data()->modification_begin();

      SubDimCellToDataMap::iterator iter;

      SubDimCellToDataMap& map = m_cell_2_data_map;

      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;

          stk::mesh::EntityId owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          stk::mesh::EntityRank owning_elementRank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

          if (owning_elementId)
            {
              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);

              if (!m_eMesh.is_valid(owning_element))
                {
                  NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                  if (resize_nodeId_data) nodeIds_onSE.resize(0);
                  nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = stk::mesh::EntityKey(owning_elementRank, 0u);
                }
            }
          else
            {
              NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
              if (resize_nodeId_data) nodeIds_onSE.resize(0);
            }
        }
      if (mod_begin_end)
        m_eMesh.get_bulk_data()->modification_end();

    }

    // remove/zero any data that points to a deleted element
    // called by Refiner with "children_to_be_removed_with_ghosts"
    void NodeRegistry::clear_element_owner_data( std::set<stk::mesh::Entity, stk::mesh::EntityLess>& elems_to_be_deleted, bool resize_nodeId_data)
    {
      SubDimCellToDataMap::iterator iter;
      SubDimCellToDataMap& map = m_cell_2_data_map;

      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;

          stk::mesh::EntityId owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          stk::mesh::EntityRank owning_elementRank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

          if (owning_elementId)
            {
              stk::mesh::Entity owning_element = get_entity_element(*m_eMesh.get_bulk_data(), owning_elementRank, owning_elementId);
              stk::mesh::Entity owning_element_1 = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);

              if (owning_element != owning_element_1)
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] NR::clear_1 error # 1= " << owning_element << " " << owning_element_1 <<  std::endl;
                  throw std::logic_error("NodeRegistry:: clear_element_owner_data_phase_1 error # 1");
                }

              if (m_eMesh.is_valid(owning_element))
                {
                  bool in_deleted_list = elems_to_be_deleted.find(owning_element) != elems_to_be_deleted.end();

                  if (in_deleted_list)
                    {
                      NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
                      if (resize_nodeId_data) nodeIds_onSE.resize(0);
                      nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>() = stk::mesh::EntityKey(owning_elementRank, 0u);
                    }
                }
            }
          else
            {
              NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
              if (resize_nodeId_data) nodeIds_onSE.resize(0);
            }
        }
    }


    void NodeRegistry::dumpDB(std::string msg)
    {
      bool dnru = DEBUG_NR_UNREF;
      if (!dnru) return;
      SubDimCellToDataMap::iterator iter;
      SubDimCellToDataMap& map = m_cell_2_data_map;
      std::cout << msg << " tmp dumpDB map size: " << map.size() << std::endl;

      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;

          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

          for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
            {
              std::cout << "tmp ddb:: node id= " << m_eMesh.identifier(nodeIds_onSE[ii]) << std::endl;
              std::cout << "subDimEntity= ";
              for (unsigned k=0; k < subDimEntity.size(); k++)
                {
                  std::cout << " " << m_eMesh.identifier(subDimEntity[k]) << " ";
                }
              std::cout << "Node: ";
              m_eMesh.print_entity(std::cout, nodeIds_onSE[ii]);
            }
        }
    }

    // estimate of memory used by this object
    unsigned NodeRegistry::get_memory_usage()
    {
      SubDimCellToDataMap::iterator iter;
      SubDimCellToDataMap& map = m_cell_2_data_map;

      unsigned mem=0;

      for (iter = map.begin(); iter != map.end(); ++iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*iter).second;

          mem += sizeof(SDCEntityType)*subDimEntity.size();
          mem += sizeof(stk::mesh::EntityKey);
          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

          unsigned mem1 = (sizeof(NodeIdsOnSubDimEntityTypeQuantum)+
                           sizeof(stk::mesh::EntityId))*nodeIds_onSE.size() +sizeof(unsigned);

          mem += mem1;
        }
      return mem;
    }

  }
}
