#include <stk_adapt/NodeRegistry.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_percept/mesh/mod/smoother/SpacingFieldUtil.hpp>

#include <set>

#if 0 && NODE_REGISTRY_MAP_TYPE_TEUCHOS_HASHTABLE
namespace Teuchos
{
  //template <class T> int hashCode(const T& x);

  template <> int hashCode(const SDCell_HashTable_Key& x)
  {
    return (int)x.getHash();
  }

}
#endif

namespace stk {
  namespace adapt {

    bool s_compare_using_entity_impl = false;

    static int s_nsz_parent = 0;
    static bool s_element_is_ghost = false;
    static stk::mesh::Selector *s_oldPartSelector = 0;

    // FIXME
    //static double m_min_spacing_factor = 0.5;  // reproduce (mostly) old behavior (centroids will be slightly diff)
    static double m_min_spacing_factor = 0.05;

#if !NODE_REGISTRY_MAP_ACCESSORS_INLINED
    SubDimCellData& NodeRegistry::getFromMap(SubDimCell_EntityId& subDimEntity)
    {
      //ftest(subDimEntity);
      return m_cell_2_data_map[subDimEntity];
    }
    void NodeRegistry::putInMap(SubDimCell_EntityId& subDimEntity, SubDimCellData& data)
    {
      //m_cell_2_data_map.insert(std::subDimEntity, data);
      //SubDimCellData& dataInMap = m_cell_2_data_map[subDimEntity];
      //dataInMap = data;
      //ftest(subDimEntity);
      m_cell_2_data_map[subDimEntity] = data;
    }
#endif

      /// fill
      ///    @param subDimEntity with the stk::mesh::EntityId's of
      ///    the ordinal @param iSubDimOrd sub-dimensional entity of
      ///    @param element of rank
      ///    @param needed_entity_rank
      ///
    void NodeRegistry::
      noInline_getSubDimEntity(SubDimCell_SDSEntityType& subDimEntity, const stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        subDimEntity.clear();
        // in the case of elements, we don't share any nodes so we just make a map of element id to node
        if (needed_entity_rank == stk::mesh::MetaData::ELEMENT_RANK)
          {
            subDimEntity.insert( element );
            //!!subDimEntity.insert(element.identifier());
            return;
          }

        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);

        //CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::mesh::MetaData::NODE_RANK);

        const unsigned *  inodes = 0;
        unsigned nSubDimNodes = 0;
        static const unsigned edge_nodes_2[2] = {0,1};
        static const unsigned face_nodes_3[3] = {0,1,2};
        static const unsigned face_nodes_4[4] = {0,1,2,3};

        // special case for faces in 3D
        if (needed_entity_rank == m_eMesh.face_rank() && needed_entity_rank == element.entity_rank())
          {
            nSubDimNodes = cell_topo_data->vertex_count;

            // note, some cells have sides with both 3 and 4 nodes (pyramid, prism)
            if (nSubDimNodes ==3 )
              inodes = face_nodes_3;
            else
              inodes = face_nodes_4;

          }
        // special case for edges in 2D
        else if (needed_entity_rank == m_eMesh.edge_rank() && needed_entity_rank == element.entity_rank())
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
        static SubDimCell_SDSEntityType subDimEntity(m_eMesh);
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
                        << " n0 =  " << nodes[iv0].identifier() << " n1= " << nodes[iv1].identifier()
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
        double alp1[2]={0.0,0.0};

        double alpsum=0.0;
        for (unsigned ipts=0; ipts < nsz; ipts++)
          {
            alp[ipts] = nodeId_elementOwnderId.get<SDC_DATA_SPACING>()[ipts];
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
                    if (0 && element.identifier() == 6659)
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
      void NodeRegistry::makeCentroid(stk::mesh::FieldBase *field, unsigned *subDimSize_in)
      {
        EXCEPTWATCH;
        bool do_respect_spacing = m_eMesh.get_respect_spacing();
        // called from main code
        if (do_respect_spacing && !subDimSize_in)
          {
            // recurse to specialize to compute edges first
            unsigned subDimSize_2 = 2;
            makeCentroid(field, &subDimSize_2);
            // then signal compute all non-edges
            subDimSize_2 = 0;
            makeCentroid(field, &subDimSize_2);
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
              field_rank = fr.entity_rank();
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
            const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            if (do_respect_spacing && *subDimSize_in == 2 && subDimEntity.size() != 2)
              continue;

            SubDimCellData& nodeId_elementOwnderId = (*iter).second;

            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            unsigned owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
            //unsigned owning_elementRank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            unsigned char owning_elementSubDimOrd = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_ORDINAL>();
            VERIFY_OP_ON(owning_elementSubDimOrd, >, 0, "hmm 2");
            --owning_elementSubDimOrd ;

            Double2& nodeId_spacing = nodeId_elementOwnderId.get<SDC_DATA_SPACING>();

            static const SubDimCellData empty_SubDimCellData;

            bool is_empty = (nodeId_elementOwnderId == empty_SubDimCellData);

            if (s_allow_empty_sub_dims && is_empty)
              {
                return;
              }

            if (is_empty)
              {
                throw std::runtime_error("makeCentroid(field) empty cell found");
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

            if (!nodeIds_onSE[0].is_valid())
              {
                continue;
              }

            stk::mesh::Entity c_node = nodeIds_onSE[0];

            if (!c_node.is_valid())
              {
                throw std::runtime_error("makeCentroid(field): bad node found 0.0");
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
                  SDSEntityType elementId = *subDimEntity.begin();
                  //!!element_p = get_entity_element(*m_eMesh.get_bulk_data(), stk::mesh::MetaData::ELEMENT_RANK, elementId);
                  element_p = elementId;
                  if (!element_p.is_valid())
                    {
                      throw std::runtime_error("makeCentroid(field): bad elem found 2");
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
                        if (!node.is_valid())
                          {
                            throw std::runtime_error("makeCentroid(field): bad node found 1.0");
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

                if (do_respect_spacing && nsz == 4 &&  spatialDim == 3)
                  {
                    //exit(1);
                    //element_side_nodes( const stk::mesh::Entity elem , int local_side_id, stk::mesh::EntityRank side_entity_rank, std::vector<stk::mesh::Entity>& side_node_entities )
                    stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(stk::mesh::MetaData::ELEMENT_RANK, owning_elementId);
                    VERIFY_OP_ON(owning_element, !=, stk::mesh::Entity(), "hmmm");
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
                    for (SubDimCell_SDSEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids, ++ipts)
                      {
                        SDSEntityType nodeId = *ids;
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
                      if (!node.is_valid())
                        {
                          throw std::runtime_error("makeCentroid(field): bad node found 2.0");
                        }
                      //double *  coord = m_eMesh.field_data(field, *node, null_u);
                      double *  field_data = m_eMesh.field_data_inlined(field, node);

                      if (doPrint && field_data)
                        {
                          //const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
                          //CellTopology cell_topo(cell_topo_data);

                          std::cout << "tmp NodeRegistry::makeCentroid(field) npts= " << subDimEntity.size() << " ipts= " << ipts
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
                    std::cout << "tmp NodeRegistry::makeCentroid(field) c_coord= " << c_coord[0] << " " << c_coord[1] << " " << c_coord[2] << std::endl;


                }
            }
          }
      } // makeCentroid(stk::mesh::FieldBase *)


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

          const CellTopologyData *const topology = stk::percept::PerceptMesh::get_cell_topology(part);
          (void)topology;
          stk::mesh::Selector selector(part);

          add_parts[0] = &part;

          SubDimCellToDataMap::iterator iter;

          // Step 1. loop over all pseudo-face/edges in the NodeRegistry

          for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
            {
              const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
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

              for (SubDimCell_SDSEntityType::const_iterator ids = subDimEntity.begin(); ids != subDimEntity.end(); ++ids)
                {
                  SDSEntityType nodeId = *ids;
                  stk::mesh::Entity node = nodeId;
                  found = false;
                  percept::MyPairIterRelation beams (m_eMesh, node, m_eMesh.element_rank());

                  if (DEBUG_ADD_RBARS > 1 && !m_eMesh.get_rank())
                    {
                      for (unsigned ii=0; ii < beams.size(); ii++)
                        {
                          if (selector(beams[ii].entity().bucket()))
                            {
                              percept::MyPairIterRelation beam_nodes (m_eMesh, beams[ii].entity(), m_eMesh.node_rank());
                              VERIFY_OP_ON(beam_nodes.size(), ==, 2, "rbar issue");
                              std::cout << "node= " << node.identifier() << " beam_nodes[" << beams[ii].entity().identifier() << "]= { "
                                        << std::setw(20) << beam_nodes[0].entity().identifier()
                                        << std::setw(20) << beam_nodes[1].entity().identifier()
                                        << std::endl;
                            }
                        }
                    }

                  // Step 3. for all beams attached to this node that belong to the current rbar block...
                  for (unsigned ii=0; ii < beams.size(); ii++)
                    {
                      if (selector(beams[ii].entity().bucket()))
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
                                      VERIFY_OP_ON(common_node.identifier(), ==, beam_nodes[jj].entity().identifier(), "rbar issue2: please rerun and exclude rbar blocks from refinement using --block_name option ");
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

                      if (!c_node.is_valid())
                        {
                          continue;
                        }

                      // only try to add element if I am the owner
                      if (c_node.owner_rank() == m_eMesh.get_parallel_rank())
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

                  UniformRefinerPatternBase::interpolateElementFields(m_eMesh, new_elems_attached_rbars[i], newElement);

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

  }
}
