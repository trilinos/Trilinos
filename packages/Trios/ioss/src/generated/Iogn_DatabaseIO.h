/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Iogn_DatabaseIO_h
#define SIERRA_Iogn_DatabaseIO_h

#include <Ioss_CodeTypes.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Field.h>
#include <Ioss_Map.h>
#include <Ioss_DBUsage.h>

#include <string>
#include <assert.h>
#include <iostream>

namespace Ioss {
  class GroupingEntity;
  class Region;
  class EntityBlock;
  class NodeBlock;
  class FaceBlock;
  class ElementBlock;
  class NodeSet;
  class EdgeSet;
  class FaceSet;
  class CommSet;
}

namespace Iogn {
  class GeneratedMesh;

  class IOFactory : public Ioss::IOFactory
    {
    public:
      static const IOFactory* factory();
    private:
      IOFactory();
      Ioss::DatabaseIO* make_IO(const std::string& filename,
				Ioss::DatabaseUsage db_usage,
				MPI_Comm communicator) const;
    };

  class DatabaseIO : public Ioss::DatabaseIO
    {
    public:
      DatabaseIO(Ioss::Region *region, const std::string& filename,
		 Ioss::DatabaseUsage db_usage, MPI_Comm communicator);
      ~DatabaseIO();

      int node_global_to_local(int /* global */, bool /* must_exist */) const {return 0;}
      int element_global_to_local(int /* global */) const {return 0;}

      // Check capabilities of input/output database...
      bool supports_nodal_fields()    const {return false;}
      bool supports_edge_fields()     const {return false;}
      bool supports_face_fields()     const {return false;}
      bool supports_element_fields()  const {return false;}
      bool supports_nodelist_fields() const {return false;}

      void read_meta_data();

      bool begin(Ioss::State state);
      bool   end(Ioss::State state);

      bool begin_state(Ioss::Region *region, int state, double time);
      bool   end_state(Ioss::Region *region, int state, double time);

      const GeneratedMesh* get_generated_mesh() const
      { return m_generatedMesh; }

      const std::vector<std::string>& get_faceset_names() const
      { return m_faceset_names; }
    private:
      void get_nodeblocks();
      void get_elemblocks();
      void get_nodesets();
      void get_facesets();
      void get_commsets();

      const Ioss::MapContainer& get_node_map() const;
      const Ioss::MapContainer& get_element_map() const;

      int get_field_internal(const Ioss::Region* reg, const Ioss::Field& field,
			     void *data, size_t data_size) const;

      int get_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int get_field_internal(const Ioss::FaceBlock* eb, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int get_field_internal(const Ioss::EdgeBlock* eb, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int get_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
			     void *data, size_t data_size) const;

      int get_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int get_field_internal(const Ioss::EdgeSet* es, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int get_field_internal(const Ioss::FaceSet* fs, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int get_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
			     void *data, size_t data_size) const;

      int put_field_internal(const Ioss::Region* region, const Ioss::Field& field,
			     void *data, size_t data_size) const;

      int put_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int put_field_internal(const Ioss::FaceBlock* fb, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int put_field_internal(const Ioss::EdgeBlock* nb, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int put_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
			     void *data, size_t data_size) const;

      int put_field_internal(const Ioss::NodeSet* ns, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int put_field_internal(const Ioss::EdgeSet* es, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int put_field_internal(const Ioss::FaceSet* fs, const Ioss::Field& field,
			     void *data, size_t data_size) const;
      int put_field_internal(const Ioss::CommSet* cs, const Ioss::Field& field,
			     void *data, size_t data_size) const;

      // Private member functions
      DatabaseIO(const DatabaseIO& from); // do not implement
      DatabaseIO& operator=(const DatabaseIO& from); // do not implement

      GeneratedMesh *m_generatedMesh;
      std::vector<std::string> m_faceset_names;

      int spatialDimension;
      int nodeCount;
      int elementCount;

      int nodeBlockCount;
      int elementBlockCount;
      int nodesetCount;
      int sidesetCount;

      // MAPS -- Used to convert from local exodusII ids/names to Sierra
      // database global ids/names

      //---Node Map -- Maps internal (1..NUMNP) ids to global ids used on the
      //               sierra side.   global = nodeMap[local]
      // nodeMap[0] contains: -1 if sequential, 0 if ordering unknown, 1
      // if nonsequential
      mutable Ioss::MapContainer        nodeMap;
      mutable Ioss::MapContainer        reorderNodeMap;
      mutable Ioss::ReverseMapContainer reverseNodeMap;
      // (local==global)

      //---Element Map -- Maps internal (1..NUMEL) ids to global ids used on the
      //               sierra side.   global = elementMap[local]
      // elementMap[0] contains: -1 if sequential, 0 if ordering unknown,
      // 1 if nonsequential
      mutable Ioss::MapContainer        elementMap;
      mutable Ioss::MapContainer        reorderElementMap;
      mutable Ioss::ReverseMapContainer reverseElementMap;

      mutable bool sequentialNG2L; // true if reverse node map is sequential
      mutable bool sequentialEG2L; // true if reverse element map is
				   // sequential (local==global)
    };
}
#endif // SIERRA_Iogn_DatabaseIO_h
