// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef IOSS_Ioss_Region_h
#define IOSS_Ioss_Region_h

#include <Ioss_CodeTypes.h>
#include <string>
#include <Ioss_GroupingEntity.h>
#include <Ioss_State.h>
// Needed for node_global_to_local inline function.
#include <Ioss_DatabaseIO.h>
#include <Ioss_EntityType.h>
#include <vector>
#include <map>
#include <iosfwd>

namespace Ioss {

  class NodeBlock;
  class ElementBlock;
  class SideSet;
  class SideBlock;
  class NodeSet;
  class CommSet;
  class CoordinateFrame;

  typedef std::vector<NodeBlock*>    NodeBlockContainer;
  typedef std::vector<EdgeBlock*>    EdgeBlockContainer;
  typedef std::vector<FaceBlock*>    FaceBlockContainer;
  typedef std::vector<ElementBlock*> ElementBlockContainer;

  typedef std::vector<NodeSet*>      NodeSetContainer;
  typedef std::vector<EdgeSet*>      EdgeSetContainer;
  typedef std::vector<FaceSet*>      FaceSetContainer;
  typedef std::vector<ElementSet*>   ElementSetContainer;

  typedef std::vector<SideSet*>      SideSetContainer;
  typedef std::vector<CommSet*>      CommSetContainer;
  typedef std::vector<double>        StateTimeContainer;

  typedef std::vector<CoordinateFrame> CoordinateFrameContainer;
  
  typedef std::map<std::string, std::string, std::less<std::string> > AliasMap;
  typedef AliasMap::value_type IOAliasValuePair;

  class Region : public GroupingEntity {
  public:

    explicit Region(DatabaseIO *iodatabase = NULL, const std::string& name="");

    ~Region();

    std::string type_string() const {return "Region";}
    EntityType type() const {return REGION;}

    void output_summary(std::ostream &strm, bool do_transient=true);
    
    bool supports_field_type(Ioss::EntityType fld_type) const;

    // Helper function...
    int64_t node_global_to_local(int64_t global, bool must_exist=true) const;

    bool begin_mode(State new_state);
    bool end_mode(State current_state);

    // Add a new state at this time, return state number
    virtual int  add_state(double time);

    // Get time corresponding to specified state
    virtual double get_state_time(int state=-1) const;
    virtual double begin_state(int state);
    virtual double end_state(int state);
    bool model_defined() const {return modelDefined;}
    bool transient_defined() const {return transientDefined;}

    /**
     * Return a pair consisting of the step (1-based) corresponding to
     * the maximum time on the database and the corresponding maximum
     * time value. Note that this may not necessarily be the last step
     * on the database if cycle and overlay are being used.
     */
    std::pair<int, double> get_max_time() const;

    /**
     * Return a pair consisting of the step (1-based) corresponding to
     * the minimum time on the database and the corresponding minimum
     * time value. Note that this may not necessarily be the first step
     * on the database if cycle and overlay are being used.
     */
    std::pair<int, double> get_min_time() const;

    // Functions for an output region...
    bool add(NodeBlock    *block);
    bool add(EdgeBlock    *block);
    bool add(FaceBlock    *block);
    bool add(ElementBlock *block);
    bool add(SideSet      *set);
    bool add(NodeSet      *set);
    bool add(EdgeSet      *set);
    bool add(FaceSet      *set);
    bool add(ElementSet   *set);
    bool add(CommSet      *set);
    bool add(const CoordinateFrame &set);

    const NodeBlockContainer&    get_node_blocks() const;
    const EdgeBlockContainer&    get_edge_blocks() const;
    const FaceBlockContainer&    get_face_blocks() const;
    const ElementBlockContainer& get_element_blocks() const;
    const SideSetContainer&      get_sidesets() const;
    const NodeSetContainer&      get_nodesets() const;
    const EdgeSetContainer&      get_edgesets() const;
    const FaceSetContainer&      get_facesets() const;
    const ElementSetContainer&   get_elementsets() const;
    const CommSetContainer&      get_commsets() const;
    const CoordinateFrameContainer&  get_coordinate_frames() const;

    // Retrieve the Grouping Entity with the specified name.
    // Returns NULL if the entity does not exist
    GroupingEntity* get_entity(const std::string& name, EntityType io_type) const;
    GroupingEntity* get_entity(const std::string& name) const;
    NodeBlock*      get_node_block(const std::string& name) const;
    EdgeBlock*      get_edge_block(const std::string& name) const;
    FaceBlock*      get_face_block(const std::string& name) const;
    ElementBlock*   get_element_block(const std::string& name) const;
    SideSet*        get_sideset(const std::string& name) const;
    SideBlock*      get_sideblock(const std::string& name) const;
    NodeSet*        get_nodeset(const std::string& name) const;
    EdgeSet*        get_edgeset(const std::string& name) const;
    FaceSet*        get_faceset(const std::string& name) const;
    ElementSet*     get_elementset(const std::string& name) const;
    CommSet*        get_commset(const std::string& name) const;

    const CoordinateFrame& get_coordinate_frame(int64_t id) const;

    // Add the name 'alias' as an alias for the databae entity with the
    // name 'db_name'. Returns true if alias added; false if problems
    // adding alias.
    bool add_alias(const std::string &db_name, const std::string &alias);
    bool add_alias(const GroupingEntity *ge);
    std::string get_alias(const std::string &alias) const;

    const AliasMap & get_alias_map() const;

    /// Get a map containing all aliases defined for the entity with basename 'name' 
      int get_aliases(const std::string &name,
		      std::vector<std::string> &aliases) const;
      
    // This routine transfers all relavant aliases from the 'this'
    // region and applies them to the 'to' file.
    void transfer_mesh_aliases(Region *to) const;

    // Ensure that the 'this' region has the same ids and names as the 'from' region.
    void synchronize_id_and_name(const Region *from, bool synchronize_attribute_field_names=false);

    // Returns true if the passed in name refers to a known Entity
    // defined on this region.  If true, then 'type' (if non-NULL) is
    // filled in with the type of the entity; if false, then type (if
    // non-NULL) is set to 'INVALID' This function is defined to
    // consolidate several distinct implementations of this code in
    // client code. Because of this, the 'type' used in the client
    // code is repeated here instead of something more generic.
    bool is_valid_io_entity(const std::string& name, unsigned int io_type,
			    std::string *type = NULL) const;

    // Retrieve the element block that contains the specified element
    // The 'local_id' is the local database id (1-based), not the global id.
    // returns NULL if no element block contains this element (local_id <= 0
    // or greater than number of elements in database)
    ElementBlock* get_element_block(size_t local_id) const;

    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    Property get_implicit_property(const std::string& name) const;

    const std::vector<std::string> &get_information_records() const;
    void add_information_records(const std::vector<std::string> &info);
    void add_information_record(const std::string &info);
    
    const std::vector<std::string> &get_qa_records() const;
    void add_qa_record(const std::string &code, const std::string &code_qa,
		       const std::string &date="", const std::string &time="");

  protected:
    int64_t internal_get_field_data(const Field& field,
				void *data, size_t data_size) const;

    int64_t internal_put_field_data(const Field& field,
				void *data, size_t data_size) const;

  private:
    void delete_database();

    AliasMap aliases_; //< Stores alias mappings

    // Containers for all grouping entities
    NodeBlockContainer    nodeBlocks;
    EdgeBlockContainer    edgeBlocks;
    FaceBlockContainer    faceBlocks;
    ElementBlockContainer elementBlocks;
    
    NodeSetContainer      nodeSets;
    EdgeSetContainer      edgeSets;
    FaceSetContainer      faceSets;
    ElementSetContainer   elementSets;
    
    SideSetContainer      sideSets;
    CommSetContainer      commSets;
    CoordinateFrameContainer coordinateFrames;
    
    mutable StateTimeContainer stateTimes;

    int currentState;
    mutable int stateCount;
    bool modelDefined;
    bool transientDefined;
  };
}
inline bool Ioss::Region::supports_field_type(Ioss::EntityType fld_type) const
{return get_database()->entity_field_support() & fld_type;}

inline int64_t Ioss::Region::node_global_to_local(int64_t global, bool must_exist) const
{ return get_database()->node_global_to_local(global,must_exist); }

inline const std::vector<std::string> &Ioss::Region::get_information_records() const
{ return get_database()->get_information_records(); }

inline void Ioss::Region::add_information_records(const std::vector<std::string> &info)
{ return get_database()->add_information_records(info); }

inline void Ioss::Region::add_information_record(const std::string &info)
{ return get_database()->add_information_record(info); }
    
inline void Ioss::Region::add_qa_record(const std::string &code, const std::string &code_qa,
					const std::string &date, const std::string &time)
{ return get_database()->add_qa_record(code, code_qa, date, time); }

inline const std::vector<std::string> &Ioss::Region::get_qa_records() const
{ return get_database()->get_qa_records(); }

#endif
