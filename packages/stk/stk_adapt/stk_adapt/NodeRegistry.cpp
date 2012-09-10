#include <stk_adapt/NodeRegistry.hpp>

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
      noInline_getSubDimEntity(SubDimCell_SDSEntityType& subDimEntity, const stk::mesh::Entity& element, stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd)
      {
        subDimEntity.clear();
        // in the case of elements, we don't share any nodes so we just make a map of element id to node
        if (needed_entity_rank == m_eMesh.element_rank())
          {
            subDimEntity.insert( const_cast<stk::mesh::Entity*>(&element) );
            //!!subDimEntity.insert(element.identifier());
            return;
          }

        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);

        //CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

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

      static double spacing_edge(unsigned iv0, unsigned iv1, unsigned nsz, unsigned nsp,  double lspc[8][3], double den_xyz[3], double *coord[8])
      {
        unsigned iv[2]={iv0,iv1};
        double alp[2]={0.0,0.0};
        double alp1[2]={0.0,0.0};
        for (unsigned ipts=0; ipts < nsz; ipts++)
          {
            unsigned jpts = iv[ipts];
            double len = 0.0;
            alp[ipts]=0.0;

            for (unsigned isp = 0; isp < nsp; isp++)
              {
                alp[ipts] += (coord[iv[1]][isp] - coord[iv[0]][isp])*lspc[jpts][isp];
                len += (coord[iv[1]][isp] - coord[iv[0]][isp])*(coord[iv[1]][isp] - coord[iv[0]][isp]);
              }                            
            alp[ipts] = std::fabs(alp[ipts])/std::sqrt(len);
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
        return 1.0-alp1[0];
        //return alp1[0];
      }

      static void normalize_spacing_0(unsigned nsz, unsigned nsp, double spc[8][3], double den_xyz[3])
      {
        //double m_min_spacing_factor = 0.5;
        double m_min_spacing_factor = 0.0;
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

      static void normalize_spacing(unsigned nsz, unsigned nsp, double spc[8][3], double den_xyz[3], double *coord[8])
      {
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
                double alp01 = spacing_edge(0, 1, 2, nsp,  lspc, den_xyz, coord);
                double alp32 = spacing_edge(3, 2, 2, nsp,  lspc, den_xyz, coord);
                double x = 0.5*(alp01+alp32);
                double alp12 = spacing_edge(1, 2, 2, nsp,  lspc, den_xyz, coord);
                double alp03 = spacing_edge(0, 3, 2, nsp,  lspc, den_xyz, coord);
                double y = 0.5*(alp12+alp03);
                if (isp == 0)
                  {
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
                double alp01 = spacing_edge(0, 1, 2, nsp,  lspc, den_xyz, coord);
                double alp32 = spacing_edge(3, 2, 2, nsp,  lspc, den_xyz, coord);
                double alp45 = spacing_edge(4, 5, 2, nsp,  lspc, den_xyz, coord);
                double alp76 = spacing_edge(7, 6, 2, nsp,  lspc, den_xyz, coord);
                double x = 0.25*(alp01+alp32+alp45+alp76);
                double alp12 = spacing_edge(1, 2, 2, nsp,  lspc, den_xyz, coord);
                double alp03 = spacing_edge(0, 3, 2, nsp,  lspc, den_xyz, coord);
                double alp56 = spacing_edge(5, 6, 2, nsp,  lspc, den_xyz, coord);
                double alp47 = spacing_edge(4, 7, 2, nsp,  lspc, den_xyz, coord);
                double y = 0.25*(alp12+alp03+alp56+alp47);
                double alp04 = spacing_edge(0, 4, 2, nsp,  lspc, den_xyz, coord);
                double alp15 = spacing_edge(1, 5, 2, nsp,  lspc, den_xyz, coord);
                double alp26 = spacing_edge(2, 6, 2, nsp,  lspc, den_xyz, coord);
                double alp37 = spacing_edge(3, 7, 2, nsp,  lspc, den_xyz, coord);
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
            else
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
              }
            double sum=0.0;
            for (unsigned ipts=0; ipts < nsz; ipts++)
              {
                sum += spc[ipts][isp];
              }
            for (unsigned ipts=0; ipts < nsz; ipts++)
              {
                spc[ipts][isp] /= sum;
                if (0 && nsz == 8 && isp == 1) 
                  std::cout << "spc[" << ipts << "]= " << spc[ipts][isp] << std::endl;
              }
          }
      }

      /// makes coordinates of this new node be the centroid of its sub entity - this version does it for all new nodes
      void NodeRegistry::makeCentroid(stk::mesh::FieldBase *field)
      {
        EXCEPTWATCH;
        //unsigned *null_u = 0;
        stk::mesh::FieldBase *spacing_field    = m_eMesh.get_field("ref_spacing_field");

        int spatialDim = m_eMesh.get_spatial_dim();
        int fieldDim = spatialDim;
        stk::mesh::EntityRank field_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
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
        if (field_rank != stk::mesh::fem::FEMMetaData::NODE_RANK)
          {
            if (field_rank == m_eMesh.element_rank())
              {
                
              }
            return;
          }

        SubDimCellToDataMap::iterator iter;

        for (iter = m_cell_2_data_map.begin(); iter != m_cell_2_data_map.end(); ++iter)
          {
            const SubDimCell_SDSEntityType& subDimEntity = (*iter).first;
            SubDimCellData& nodeId_elementOwnderId = (*iter).second;

            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            unsigned owning_elementId = stk::mesh::entity_id(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            //unsigned owning_elementRank = stk::mesh::entity_rank(nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>());
            unsigned char owning_elementSubDimOrd = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_ORDINAL>();
            VERIFY_OP_ON(owning_elementSubDimOrd, >, 0, "hmm 2");
            --owning_elementSubDimOrd ;

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

            stk::mesh::EntityRank needed_entity_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
            // SPECIAL CASE
            // SPECIAL CASE
            // SPECIAL CASE
            if (subDimEntity.size() == 1)
              {
                needed_entity_rank = m_eMesh.element_rank();
              }

            if (nodeIds_onSE[0] == 0)
              {
                continue;
              }

            stk::mesh::Entity * c_node = nodeIds_onSE[0];

            if (!c_node)
              {
                throw std::runtime_error("makeCentroid(field): bad node found 0.0");
              }

            std::vector<double> c_p(fieldDim,0);
            bool doPrint = false;
            std::vector<stk::mesh::Entity *> nodes(8,(stk::mesh::Entity *)0);
            unsigned nsz = 0;
            bool do_spacing=true;

            if (needed_entity_rank == m_eMesh.element_rank())
              {
                EXCEPTWATCH;
                stk::mesh::Entity *element_p = 0;
                {
                  SDSEntityType elementId = *subDimEntity.begin();
                  //!!element_p = get_entity_element(*m_eMesh.get_bulk_data(), m_eMesh.element_rank(), elementId);
                  element_p = elementId;
                  if (!element_p)
                    {
                      throw std::runtime_error("makeCentroid(field): bad elem found 2");
                    }
                }

                stk::mesh::Entity& element = *element_p;
                bool element_is_ghost = m_eMesh.isGhostElement(element);
                if (element_is_ghost)
                  {
                    //std::cout << "tmp found ghost" << std::endl;
                  }
                else
                  {
                    const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
                    unsigned npts = elem_nodes.size();
                    nsz = npts;
                    nodes.resize(nsz, (stk::mesh::Entity *)0);
                    c_p.resize(fieldDim,0);
                    //if (npts == 2) doPrint=true;
                    //double dnpts = elem_nodes.size();
                    for (unsigned ipts = 0; ipts < npts; ipts++)
                      {
                        stk::mesh::Entity * node = elem_nodes[ipts].entity();
                        if (!node)
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
                nodes.resize(nsz, (stk::mesh::Entity *)0);
                c_p.resize(fieldDim,0);

                if (do_spacing && nsz == 4 &&  spatialDim == 3)
                  {
                    //exit(1);
                    //element_side_nodes( const stk::mesh::Entity & elem , int local_side_id, stk::mesh::EntityRank side_entity_rank, std::vector<stk::mesh::Entity *>& side_node_entities )
                    stk::mesh::Entity *owning_element = m_eMesh.get_bulk_data()->get_entity(m_eMesh.element_rank(), owning_elementId);
                    VERIFY_OP_ON(owning_element, !=, 0, "hmmm");
                    std::vector<stk::mesh::Entity *> side_node_entities;
                    PerceptMesh::element_side_nodes(*owning_element, owning_elementSubDimOrd, m_eMesh.face_rank(), side_node_entities);
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
                        stk::mesh::Entity * node = nodeId;
                        nodes[ipts]=node;
                      }
                  }
              }

            {
              //if ( (spacing_field && (spacing_field != field) && subDimEntity.size() == 2))
              // FIXME for quadratic elements
              if (do_spacing && (nsz <= 8 && spacing_field && (spacing_field != field) ) )
                {
                  EXCEPTWATCH;
                  unsigned ipts=0;

                  double * coord[8] = {0,0,0,0,0,0,0,0};
                  double * field_data[8] = {0,0,0,0,0,0,0,0};
                  double * spacing[8] = {0,0,0,0,0,0,0,0};

                  for (ipts=0; ipts < nsz; ipts++)
                    {
                      coord[ipts] = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), *nodes[ipts]);
                      field_data[ipts] = m_eMesh.field_data_inlined(field, *nodes[ipts]);
                      spacing[ipts] = m_eMesh.field_data_inlined(spacing_field, *nodes[ipts]);
                    }

                  double spc[8][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
                  double den = 0.0;
                  double den_xyz[3] = {0,0,0};
                  if (nsz == 2) 
                    {
                      den = 0.0;
                      for (ipts=0; ipts < nsz; ipts++)
                        {
                          double len = 0.0;
                          for (int isp = 0; isp < spatialDim; isp++)
                            {
                              spc[ipts][0] += (coord[1][isp] - coord[0][isp])*spacing[ipts][isp];
                              len += (coord[1][isp] - coord[0][isp])*(coord[1][isp] - coord[0][isp]);
                            }                            
                          spc[ipts][0] = std::fabs(spc[ipts][0])/std::sqrt(len);
                        }
                      for (ipts=0; ipts < nsz; ipts++)
                        {
                          for (int isp = 1; isp < spatialDim; isp++)
                            {
                              spc[ipts][isp] = spc[ipts][0];
                            }                            
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
                  normalize_spacing(nsz, spatialDim, spc, den_xyz, coord);
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

                }
              else
                {
                  EXCEPTWATCH;
                  double dnpts = nsz;
                  unsigned ipts=0;
                  for (ipts=0; ipts < nsz; ipts++)
                    {
                      stk::mesh::Entity * node = nodes[ipts];
                      if (!node)
                        {
                          throw std::runtime_error("makeCentroid(field): bad node found 2.0");
                        }
                      //double *  coord = m_eMesh.field_data(field, *node, null_u);
                      double *  field_data = m_eMesh.field_data_inlined(field, *node);

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
              double *  c_coord = m_eMesh.field_data_inlined(field, *c_node);

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


  }
}
