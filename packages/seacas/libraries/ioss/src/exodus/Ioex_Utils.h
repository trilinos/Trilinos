#ifndef IOEX_UTILS_H
#define IOEX_UTILS_H
#include <Ioss_Utils.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_ElementTopology.h>
#include <Ioss_CoordinateFrame.h>

#include <exodusII.h>
#include <string>
#include <vector>
#include <set>

// Contains code that is common between the file-per-processor (Iofx)
// and parallel exodus (Iopx) and base exodus (Ioex) classes.

namespace Ioss {
  class GroupingEntity;
  typedef std::vector<CoordinateFrame> CoordinateFrameContainer;
}

typedef std::vector<std::string> NameList;

namespace Ioex {
  typedef std::set<std::pair<int64_t, int64_t> > EntityIdSet;
  typedef std::set<std::string> SideSetSet;
  typedef std::map<std::string, const std::string, std::less<const std::string> > SideSetMap;

  struct TopologyMapCompare {
      bool operator() (const std::pair<std::string, const Ioss::ElementTopology*> &lhs,
                       const std::pair<std::string, const Ioss::ElementTopology*> &rhs) const
      {
        assert(lhs.second != NULL);
        assert(rhs.second != NULL);
        return lhs.first < rhs.first || (!(rhs.first < lhs.first) && lhs.second->name() < rhs.second->name());
      }
  };

  typedef std::map<std::pair<std::string, const Ioss::ElementTopology*>, int, TopologyMapCompare > TopologyMap;
  typedef TopologyMap::value_type TopoMapPair;

  const char *Version();
  bool check_processor_info(int exodusFilePtr, int processor_count, int processor_id);

  void update_last_time_attribute(int exodusFilePtr, double value);
  bool read_last_time_attribute(int exodusFilePtr, double *value);

  bool type_match(const std::string& type, const char *substring);
  int64_t extract_id(const std::string &name_id);
  bool set_id(const Ioss::GroupingEntity *entity, ex_entity_type type, Ioex::EntityIdSet *idset);
  int64_t  get_id(const Ioss::GroupingEntity *entity, ex_entity_type type,
		  Ioex::EntityIdSet *idset);
  void decode_surface_name(Ioex::SideSetMap &fs_map, Ioex::SideSetSet &fs_set,
			   const std::string &name);
  void fix_bad_name(char* name);

  void exodus_error(int exoid, int lineno, int /* processor */);

  void check_non_null(void *ptr, const char *type, const std::string &name);

  void add_map_fields(int exoid, Ioss::ElementBlock *block, int64_t my_element_count,
                      size_t name_length);

  void add_coordinate_frames(int exoid, Ioss::Region *region);
  void write_coordinate_frames(int exoid, const Ioss::CoordinateFrameContainer &frames);
  
  template <typename T>
  bool check_block_order(const std::vector<T*> &blocks);

  bool find_displacement_field(Ioss::NameList &fields,
                               const Ioss::GroupingEntity *block,
                               int ndim, std::string *disp_name);

  char ** get_exodus_names(size_t count, int size);
  void delete_exodus_names(char **names, int count);

  void get_fields(int64_t entity_count, char** names, size_t num_names,
                  Ioss::Field::RoleType fld_role,
                  const char suffix_separator, int *local_truth,
                  std::vector<Ioss::Field> &fields);

  std::string get_entity_name(int exoid, ex_entity_type type, int64_t id,
			      const std::string &basename, int length,
			      bool &db_has_name);

  void filter_element_list(Ioss::Region *region,
                           Ioss::Int64Vector &elements, Ioss::Int64Vector &sides,
                           bool remove_omitted_elements);

  bool filter_node_list(Ioss::Int64Vector &nodes,
                        const std::vector<unsigned char> &node_connectivity_status);

  template <typename T>
  void filter_node_list(T* data, std::vector<T> &dbvals, const std::vector<int64_t> &active_node_index)
  {
    for (size_t i=0; i < active_node_index.size(); i++) {
      data[i] = dbvals[active_node_index[i]];
    }
  }

  void filter_element_list(Ioss::Region *region,
                           Ioss::Int64Vector &elements, Ioss::Int64Vector &sides,
                           bool remove_omitted_elements);

  void separate_surface_element_sides(Ioss::Int64Vector &element,
				      Ioss::Int64Vector &sides,
				      Ioss::Region *region,
				      Ioex::TopologyMap &topo_map,
				      Ioex::TopologyMap &side_map,
				      Ioss::SurfaceSplitType split_type);
}
#endif
