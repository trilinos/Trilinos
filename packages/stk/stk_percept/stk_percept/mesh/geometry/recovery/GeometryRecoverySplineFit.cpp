#include <stk_percept/PerceptMesh.hpp>

#include <stk_percept/mesh/geometry/recovery/GeometryRecoverySplineFit.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <stk_percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
#include <stk_percept/mesh/geometry/stk_geom/LocalCubicSplineFit.hpp>
#endif

namespace stk {
  namespace percept {

#if defined( STK_PERCEPT_HAS_GEOMETRY )
    stk::mesh::Part& GeometryRecoverySplineFit::create_clone_part(const std::string& clone_name, const stk::mesh::EntityRank rank_to_clone, bool make_part_io_part )
    {
      stk::mesh::Part& clone = m_eMesh.get_fem_meta_data()->declare_part(clone_name, rank_to_clone);
      if (make_part_io_part && clone.attribute<Ioss::GroupingEntity>() == NULL) {
        stk::io::put_io_part_attribute(clone);
      }
      return clone;
    }

    stk::mesh::Part& GeometryRecoverySplineFit::clone_part_bulk(const stk::mesh::Part& part, const std::string& clone_name, const stk::mesh::EntityRank rank_to_clone, bool make_part_io_part )
    {
      stk::mesh::Part * clone_p = m_eMesh.get_fem_meta_data()->get_part(clone_name);
      if (!clone_p)
        throw std::runtime_error("AdaptMain::clone_part: no part named: "+clone_name);
      stk::mesh::Part& clone = *clone_p;

      stk::mesh::Selector this_part(part);
      std::vector<stk::mesh::Entity> entities;
      stk::mesh::PartVector add_parts(1,&clone), remove_parts;
      const std::vector<stk::mesh::Bucket*> & entity_buckets = m_eMesh.get_bulk_data()->buckets( rank_to_clone );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = entity_buckets.begin() ; k != entity_buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          if (this_part(bucket))
            {
              const unsigned num_entities_in_bucket = bucket.size();
              for (unsigned i_entity = 0; i_entity < num_entities_in_bucket; i_entity++)
                {
                  stk::mesh::Entity entity = bucket[i_entity];
                  entities.push_back(entity);
                }
            }
        }

      for (unsigned ii=0; ii < entities.size(); ii++)
        {
          m_eMesh.get_bulk_data()->change_entity_parts( entities[ii], add_parts, remove_parts );
        }
      return clone;
    }

    /// returns true if closed found
    bool GeometryRecoverySplineFit::get_sorted_curve_node_entities(stk::mesh::Part& part, std::vector<stk::mesh::Entity>& sorted_entities)
    {
      bool debug_print = false;
#define APRINTLN(a) do { if (debug_print) std::cout << #a << " = " << a << std::endl; } while(0)
#define APRINTLN2(a,b) do { if (debug_print) std::cout << #a << " = " << a << " " << #b << " = " << b << std::endl; } while(0)

      bool closed = false;
      if (debug_print) std::cout << "AdaptMain::get_sorted_curve_node_entities: processing part = " << part.name() << std::endl;
      VERIFY_OP_ON(part.primary_entity_rank(), ==, m_eMesh.edge_rank(), "bad part");

      typedef std::set<stk::mesh::Entity, stk::mesh::EntityLess> SetOfEntities;
      SetOfEntities node_set(*m_eMesh.get_bulk_data());

      //typedef std::list<stk::mesh::Entity, stk::mesh::EntityLess> ListOfEntities;
      //ListOfEntities node_list(*m_eMesh.get_bulk_data());
      typedef std::list<stk::mesh::Entity> ListOfEntities;
      ListOfEntities node_list;

      stk::mesh::Selector this_part(part);
      stk::mesh::Entity edge_first = stk::mesh::Entity();
      const std::vector<stk::mesh::Bucket*> & edge_buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.edge_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = edge_buckets.begin() ; k != edge_buckets.end() ; ++k )
        {
          if (this_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              edge_first = bucket[0];
              break;
            }
        }

      int node_set_nnodes=0;
      if (debug_print)
        {
          const std::vector<stk::mesh::Bucket*> & node_buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator k = node_buckets.begin() ; k != node_buckets.end() ; ++k )
            {
              if (this_part(**k))
                {
                  stk::mesh::Bucket & bucket = **k ;
                  const unsigned num_nodes_in_bucket = bucket.size();
                  node_set_nnodes += num_nodes_in_bucket;
                  for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                    {
                      stk::mesh::Entity node = bucket[iNode];
                      node_set.insert(node);
                    }
                }
            }
        }

      const MyPairIterRelation edge_nodes(*m_eMesh.get_bulk_data(), edge_first, m_eMesh.node_rank() );
      VERIFY_OP_ON(edge_nodes.size(), ==, 2, "bad edge");

      for (unsigned idir=0; idir < 2; idir++)
        {
          stk::mesh::Entity current_edge = edge_first;
          stk::mesh::Entity last_node = stk::mesh::Entity();
          bool continue_proc=true;
          while(continue_proc)
            {
              const MyPairIterRelation current_edge_nodes(*m_eMesh.get_bulk_data(), current_edge, m_eMesh.node_rank() );
              unsigned JDIR=0;
              if (current_edge == edge_first)
                JDIR = idir;
              else
                {
                  for (unsigned jdir=0; jdir < 2; jdir++)
                    {
                      if (current_edge_nodes[jdir].entity() != last_node)
                        {
                          JDIR = jdir;
                          break;
                        }
                    }
                }
              last_node = current_edge_nodes[JDIR].entity();

              bool is_not_in_list = std::find(node_list.begin(), node_list.end(), last_node) == node_list.end();
              VERIFY_OP_ON(is_not_in_list, ==, true, "bad node list");

              bool is_node_in_part = this_part(m_eMesh.bucket(last_node));
              VERIFY_OP_ON(is_node_in_part, ==, true, "bad node not in part");
              if (idir==0)
                node_list.push_back(last_node);
              else
                node_list.push_front(last_node);

              is_not_in_list = std::find(node_list.begin(), node_list.end(), last_node) == node_list.end();
              VERIFY_OP_ON(is_not_in_list, ==, false, "bad node list 2");

              const MyPairIterRelation last_node_edges(*m_eMesh.get_bulk_data(), last_node, m_eMesh.edge_rank() );
              bool found_new_edge = false;
              for (unsigned kdir=0; kdir < 2; kdir++)
                {
                  if (last_node_edges[kdir].entity() != current_edge && this_part(m_eMesh.bucket(last_node_edges[kdir].entity() )))
                    {
                      current_edge = last_node_edges[kdir].entity();
                      // check for closed
                      if (idir == 0 && current_edge == edge_first)
                        {
                          if (debug_print) std::cout << "AdaptMain::get_sorted_curve_node_entities: found closed condition..." << std::endl;
                          VERIFY_OP_ON(last_node, ==, edge_nodes[1].entity(), "integrity check failed");
                          node_list.push_front(edge_nodes[1].entity());
                          ++idir; //force loop exit
                        }
                      else
                        {
                          found_new_edge = true;
                        }
                      break;
                    }
                }
              if (!found_new_edge)
                continue_proc = false;

              //VERIFY_OP_ON(found, ==, true, "hmmm");
            }
        }
      APRINTLN(node_list.size());
      APRINTLN2(node_set_nnodes,node_set.size());

      if (node_list.front() == node_list.back())
        {
          // closed case
          closed = true;
          sorted_entities.resize(0);
          sorted_entities.assign(node_list.begin(), node_list.end());
        }
      else
        {
          sorted_entities.resize(0);
          sorted_entities.assign(node_list.begin(), node_list.end());

          if (debug_print) VERIFY_OP_ON(node_list.size(), ==, node_set.size(), "node set/list error for part= "+part.name());
        }
      return closed;
    }

    void GeometryRecoverySplineFit::fit_geometry_create_parts_meta()
    {
      //bool debug_print = false;

      const stk::mesh::PartVector & parts = m_eMesh.get_fem_meta_data()->get_parts();
      m_surface_parts.resize(0);
      m_topo_parts.resize(0);

      int n_topo = 10000;
      unsigned nparts = parts.size();
      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];
          if (stk::mesh::is_auto_declared_part(part) )
            //|| (part.name().find("oldElem") != std::string::npos))
            continue;

          const STK_Adapt_Auto_Part *side_auto_part = part.attribute<STK_Adapt_Auto_Part>();
          if (side_auto_part)
            continue;

          if (part.primary_entity_rank() != m_eMesh.edge_rank())
            continue;
          if (part.subsets().size() == 0)  // skip parts like surface_quad4_edge2_4
            continue;

          std::string clone_name = "tbc_curve_block_"+toString(n_topo++);
          stk::mesh::Part& clone = create_clone_part(clone_name, m_eMesh.node_rank());
          m_topo_parts.push_back(&clone);
          m_surface_parts.push_back(&part);
        }
    }


    void GeometryRecoverySplineFit::fit_geometry(const std::string& filename)
    {
      using namespace stk::geom;

      bool debug_print = false;
      SplineFit::s_debug_print = debug_print;
      ON::Begin();

      FILE* file = ON::OpenFile( filename.c_str(), "wb" );
      if (!file)
        throw std::runtime_error("couldn't open file: "+filename+" in fit_geometry");

      // allow additional parts to be added
      ONX_Model model;

      // set revision history information
      model.m_properties.m_RevisionHistory.NewRevision();

      // set application information
      model.m_properties.m_Application.m_application_name = "Percept";
      model.m_properties.m_Application.m_application_URL = "http://www.sandia.gov";
      model.m_properties.m_Application.m_application_details = "Fit cubic splines.";

      // some notes
      model.m_properties.m_Notes.m_notes = "notes";
      model.m_properties.m_Notes.m_bVisible = true;

      // file settings (units, tolerances, views, ...)
      model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::inches; //ON::meters; //ON::inches;
      model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
      model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
      model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%

      // layer table
      {
        // OPTIONAL - define some layers
        ON_Layer layer[1];

        layer[0].SetLayerName("Default");
        layer[0].SetVisible(true);
        layer[0].SetLocked(false);
        layer[0].SetColor( ON_Color(0,0,0) );

        model.m_layer_table.Append(layer[0]);
      }

      // object table
      {
        //const stk::mesh::PartVector & parts = m_eMesh.get_fem_meta_data()->get_parts();
        //stk::mesh::PartVector m_surface_parts;
        //stk::mesh::PartVector m_topo_parts;

        m_eMesh.get_bulk_data()->modification_begin();
        unsigned nparts = m_topo_parts.size();
        for (unsigned ipart=0; ipart < nparts; ipart++)
          {
            clone_part_bulk(*m_surface_parts[ipart], m_topo_parts[ipart]->name(), m_eMesh.node_rank());
          }
        m_eMesh.get_bulk_data()->modification_end();

        nparts = m_topo_parts.size();
        for (unsigned ipart=0; ipart < nparts; ipart++)
          {
            stk::mesh::Part& part = *m_topo_parts[ipart];

            std::vector<stk::mesh::Entity> sorted_entities;
            bool isClosed = get_sorted_curve_node_entities(*m_surface_parts[ipart], sorted_entities);

            LocalCubicSplineFit cf;
            if (isClosed) {
              cf.setIsPeriodic(true);  // default is to assume a smooth seam at the closed/repeated node
            }
            int n = sorted_entities.size();
            Vectors2D Q(n);
            for (int i=0; i < n; i++)
              {
                double *cdata = m_eMesh.field_data(m_eMesh.get_coordinates_field(), sorted_entities[i]);
                Q[i] = Vector2D(cdata[0], cdata[1]);
              }

            DPRINTLN(Q);
            ON_Curve *curve = cf.fit(Q);
            if (debug_print)
              {
                std::cout << "Part = " << part.name() << std::endl;
                cf.print();
              }

            if ( curve->IsValid() )
              {
                ONX_Model_Object& mo = model.m_object_table.AppendNew();
                mo.m_object = curve;
                mo.m_bDeleteObject = true;
                mo.m_attributes.m_layer_index = 0;
                mo.m_attributes.m_name = part.name().c_str();
              }
            else
              delete curve;
          }
      }

      // start section comments
      const char* sStartSectionComment = __FILE__ "PerceptMesh::fit_geometry" __DATE__;

      ON_BinaryFile archive( ON::write3dm, file );

      // Set uuid's, indices, etc.
      model.Polish();
      // writes model to archive
      int version = 5; // File can be read by Rhino 5
      ON_TextLog error_log;
      bool ok = model.Write(archive, version, sStartSectionComment, &error_log );
      VERIFY_OP_ON(ok,==,true,"failed write of 3dm file");

      ON::End();

    }
#endif

  }
}
