/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Iohb_DatabaseIO_h
#define SIERRA_Iohb_DatabaseIO_h

#include <Ioss_CodeTypes.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Field.h>
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

namespace Iohb {
  class Layout;
}

namespace Iohb {

  enum Format{DEFAULT=0,SPYHIS=1};

  class IOFactory : public Ioss::IOFactory
    {
    public:
      static const IOFactory* factory();
    private:
      IOFactory();
      Ioss::DatabaseIO* make_IO(const std::string& filename,
				Ioss::DatabaseUsage db_usage,
				MPI_Comm communicator) const;
      
      void register_library_versions() const {} // Nothing to register
      //  static const IOFactory registerThis;
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
      bool supports_nodal_fields()   const {return false;}
      bool supports_edge_fields()    const {return false;}
      bool supports_face_fields()    const {return false;}
      bool supports_element_fields() const {return false;}
      bool supports_nodelist_fields() const {return false;}

      void read_meta_data() {}

      bool begin(Ioss::State state);
      bool   end(Ioss::State state);

      bool begin_state(Ioss::Region *region, int state, double time);
      bool   end_state(Ioss::Region *region, int state, double time);

    private:
      void initialize(const Ioss::Region *region) const;
      
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

      std::ostream *logStream;
      Layout *layout_;
      Layout *legend_;

      std::string tsFormat;
      int precision_;
      bool showLabels;
      bool showLegend;

      bool initialized_;
      enum Format fileFormat;
    };
}
#endif // SIERRA_Iohb_DatabaseIO_h
