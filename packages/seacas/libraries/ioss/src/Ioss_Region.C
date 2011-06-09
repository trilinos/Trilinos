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

#include <Ioss_Region.h>

#include <assert.h>

#include <sstream>
#include <cstring>
#include <string>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <limits.h>

#include <Ioss_DatabaseIO.h>
#include <Ioss_NodeBlock.h>
#include <Ioss_ElementBlock.h>
#include <Ioss_NodeSet.h>
#include <Ioss_SideSet.h>
#include <Ioss_SideBlock.h>
#include <Ioss_CommSet.h>
#include <Ioss_Utils.h>

namespace {
  const std::string id_str()           { return std::string("id");}
  const std::string orig_elem_str()    { return std::string("original_element_type");}
  const std::string orig_block_order() { return std::string("original_block_order");}

  bool lessOffset(const Ioss::ElementBlock *b1, const Ioss::ElementBlock *b2) {
    assert(b1->property_exists(orig_block_order()));
    assert(b2->property_exists(orig_block_order()));
    return b1->get_property(orig_block_order()).get_int() <
	   b2->get_property(orig_block_order()).get_int();
  }

  std::string uppercase(const std::string &name);

  template <typename T>
  void uniqify(std::vector<T> &container) {
    std::sort(container.begin(), container.end());
    container.erase(std::unique(container.begin(), container.end()), container.end());
    std::vector<T>(container).swap(container);
  }

}

namespace Ioss {
  Region::Region(DatabaseIO *iodatabase, const std::string& my_name)
    : GroupingEntity(iodatabase, my_name), currentState(-1), stateCount(0)
  {
    assert(iodatabase != NULL);
    iodatabase->set_region(this);

    properties.add(Property("entity_count", 1));

    if (iodatabase->is_input()) {
      // Read metadata -- populates GroupingEntity lists
      Region::begin_mode(STATE_DEFINE_MODEL);
      iodatabase->read_meta_data();
      Region::end_mode(STATE_DEFINE_MODEL);
      Region::begin_mode(STATE_READONLY);
    }
    properties.add(Property(this,
			    "node_block_count",    Property::INTEGER));
    properties.add(Property(this,
			    "element_block_count", Property::INTEGER));
    properties.add(Property(this,
			    "side_set_count",      Property::INTEGER));
    properties.add(Property(this,
			    "node_set_count",      Property::INTEGER));
    properties.add(Property(this,
			    "comm_set_count",      Property::INTEGER));
    properties.add(Property(this,
			    "node_count",          Property::INTEGER));
    properties.add(Property(this,
			    "element_count",       Property::INTEGER));
    properties.add(Property(this,
			    "state_count",         Property::INTEGER));
    properties.add(Property(this,
			    "current_state",       Property::INTEGER));
    properties.add(Property(this,
			    "database_name",       Property::STRING));
  }

  Region::~Region()
  {
    // Region owns all sub-grouping entities it contains...
    try {
      {
	NodeBlockContainer::const_iterator i = nodeBlocks.begin();
	while (i != nodeBlocks.end()) {
	  delete (*i++);
	}
      }

      // Element Blocks...
      {
	ElementBlockContainer::const_iterator i = elementBlocks.begin();
	while (i != elementBlocks.end()) {
	  delete (*i++);
	}
      }

      // SideSets...
      {
	SideSetContainer::const_iterator i = sideSets.begin();
	while (i != sideSets.end()) {
	  delete (*i++);
	}
      }

      // Node Sets
      {
	NodeSetContainer::const_iterator i = nodeSets.begin();
	while (i != nodeSets.end()) {
	  delete (*i++);
	}
      }

      // Communication Sets
      {
	CommSetContainer::const_iterator i = commSets.begin();
	while (i != commSets.end()) {
	  delete (*i++);
	}
      }
      // Region owns the database pointer even though other entities use it.
      GroupingEntity::really_delete_database();
    } catch (...) {
    }
  }

  void Region::delete_database()
  {
    GroupingEntity::really_delete_database();
  }

  void Region::output_summary(std::ostream &strm, bool do_transient)
  {
    strm << "\n Database: " << get_database()->get_filename() << "\n";
    
    strm << "\n Number of coordinates per node       =" << std::setw(9)
	 << get_property("spatial_dimension").get_int() << "\n";
    strm << " Number of nodes                      =" << std::setw(9)
	 << get_property("node_count").get_int() << "\n";
    strm << " Number of elements                   =" << std::setw(9)
	 << get_property("element_count").get_int() << "\n";
    strm << " Number of element blocks             =" << std::setw(9)
	 << get_property("element_block_count").get_int() << "\n";
    strm << " Number of nodal point sets           =" << std::setw(9)
	 << get_property("node_set_count").get_int() << "\n";
    strm << " Number of element side sets          =" << std::setw(9)
	 << get_property("side_set_count").get_int() << "\n\n";

    if (do_transient && get_property("state_count").get_int() > 0) {
      strm << " Number of global variables           =" << std::setw(9)
	   << field_count() << "\n";
      {
	Ioss::NameList names;
	nodeBlocks[0]->field_describe(Ioss::Field::TRANSIENT, &names);
	strm << " Number of variables at each node     =" << std::setw(9)
	     << names.size() << "\n";
      }

      {
	const Ioss::ElementBlockContainer &blocks = get_element_blocks();
	Ioss::ElementBlockContainer::const_iterator i = blocks.begin();
	Ioss::NameList names;
	while (i != blocks.end()) {
	  (*i)->field_describe(Ioss::Field::TRANSIENT, &names);
	  i++;
	}
	uniqify(names);
	strm << " Number of variables at each element  =" << std::setw(9)
	     << names.size() << "\n";
      }

      {
	const Ioss::NodeSetContainer &blocks = get_nodesets();
	Ioss::NodeSetContainer::const_iterator i = blocks.begin();
	Ioss::NameList names;
	while (i != blocks.end()) {
	  (*i)->field_describe(Ioss::Field::TRANSIENT, &names);
	  i++;
	}
	uniqify(names);
	strm << " Number of variables at each nodeset  =" << std::setw(9)
	     << names.size() << "\n";
      }

      {
	Ioss::NameList names;
	const Ioss::SideSetContainer fss = get_sidesets();
	Ioss::SideSetContainer::const_iterator i = fss.begin();
	while (i != fss.end()) {
	  const Ioss::SideBlockContainer fbs = (*i)->get_side_blocks();
	  Ioss::SideBlockContainer::const_iterator j = fbs.begin();
	  while (j != fbs.end()) {
	    (*j)->field_describe(Ioss::Field::TRANSIENT, &names);
	    ++j;
	  }
	  i++;
	}

	uniqify(names);
	strm << " Number of variables at each sideset  =" << std::setw(9)
	     << names.size() << "\n";
      }

      strm << "\n Number of time steps on the database =" << std::setw(9)
	   << get_property("state_count").get_int() << "\n";
    }
  }

  bool Region::begin_mode(State new_state)
  {
    // All transitions must begin from the 'STATE_CLOSED' state or be to
    // the 'STATE_CLOSED' state (There are no nested begin/end pairs at
    // this time.)

    bool success = false;
    if (new_state == STATE_CLOSED) {
      success = set_state(new_state);
    } else {
      switch (get_state()) {
      case STATE_CLOSED:
	// Make sure we can go to the specified state.
	switch (new_state) {
	default:
	  success = set_state(new_state);
	}
	break;

	// For the invalid transitions; provide a more meaningful
	// message in certain cases...
      case STATE_READONLY:
	{
	  std::ostringstream errmsg;
	  errmsg << "Cannot change state of an input (readonly) database in "
		 << get_database()->get_filename();
	  IOSS_ERROR(errmsg);
	}

	break;
      default:
	{
	  std::ostringstream errmsg;
	  errmsg << "Invalid nesting of begin/end pairs in "
		 << get_database()->get_filename();
	  IOSS_ERROR(errmsg);
	}
      }
    }
    // Pass the 'begin state' message on to the database so it can do any
    // cleanup/data checking/manipulations it needs to do.
    if (success) {
      DatabaseIO *db = (DatabaseIO*)get_database();
      success = db->begin(new_state);
    }

    return success;
  }


  bool Region::end_mode(State current_state)
  {
    // Check that 'current_state' matches the current state of the
    // Region (that is, we are leaving the state we are in).

    if (get_state() != current_state) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Specified end state does not match currently open state\n"
	     << "       [" << get_database()->get_filename() << "]\n";
      IOSS_ERROR(errmsg);
    }

    if (current_state == STATE_DEFINE_MODEL) {
      // Sort the element blocks based on the idOffset field...
      if (!get_database()->is_input()) {
	std::sort(elementBlocks.begin(), elementBlocks.end(), lessOffset);

	// Now update the block offsets based on this new order...
	int offset = 0;
	ElementBlockContainer::iterator i = elementBlocks.begin();
	while (i != elementBlocks.end()) {
	  (*i)->set_offset(offset);
	  offset += (*i++)->get_property("entity_count").get_int();
	}
      }
    }

    // Pass the 'end state' message on to the database so it can do any
    // cleanup/data checking/manipulations it needs to do.
    DatabaseIO *db = (DatabaseIO*)get_database();
    bool success = db->end(current_state);


    begin_mode(STATE_CLOSED);

    return success;
  }

  int Region::add_state(double time)
  {
    // NOTE:  For restart input databases, it is possible that the time
    //        is not monotonically increasing...
    double out_time = time;
    if ( !get_database()->is_input() && !get_database()->usage() == WRITE_HEARTBEAT &&
	 stateTimes.size() >= 1 && time <= stateTimes[stateTimes.size()-1]) {
      // Check that time is increasing...
      std::ostringstream errmsg;
      errmsg << "Current time, " << time
	     << ", is not greater than previous time, " << stateTimes[stateTimes.size()-1]
	     << " in\n"
	     << get_database()->get_filename() << std::endl
	     << "Modifying time to ensure monotonic increase.";
      IOSS_WARNING << errmsg.str();
      out_time = stateTimes[stateTimes.size()-1] + stateTimes[stateTimes.size()-1]/1.0e5;
    }

    if (get_database()->is_input() ||
	get_database()->usage() == WRITE_RESULTS ||
	get_database()->usage() == WRITE_RESTART ) {
      stateTimes.push_back(out_time);
      assert((int)stateTimes.size() == stateCount+1);

    } else {

      // Keep only the last time in the vector... This is to avoid memory growth for output
      // databases that write lots of steps (heartbeat, history).  There is no need to keep
      // a list of times that have been written since they are just streamed out and never read
      // We do sometimes need the list of times written to restart or results files though...
      if (stateTimes.empty()) {
	stateTimes.push_back(out_time);
      } else {
	stateTimes[0] = out_time;
      }
    }
    return ++stateCount;;
  }

  double Region::get_state_time(int state) const
  {

    // Return time corresponding to specified state.
    // If state == -1, return time for currently active state
    double time = 0.0;
    if (state == -1) {
      if (get_database()->is_input() ||
	  get_database()->usage() == WRITE_RESULTS ||
	  get_database()->usage() == WRITE_RESTART ) {
	if (currentState == -1) {
	  std::ostringstream errmsg;
	  errmsg << "ERROR: No currently active state.\n"
		 << "       [" << get_database()->get_filename() << "]\n";
	  IOSS_ERROR(errmsg);
	} else {
	  assert((int)stateTimes.size() >= currentState);

	  time = stateTimes[currentState-1];
	}
      } else {
	assert(stateTimes.size() >= 1);
	time = stateTimes[0];
      }
    } else if (state <= 0 || state > stateCount) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Requested state (" << state << ") is invalid.\n"
	     << "       State must be between 0 and " << stateCount << ".\n"
	     << "       [" << get_database()->get_filename() << "]\n";
      IOSS_ERROR(errmsg);
    } else {
      if (get_database()->is_input() ||
	  get_database()->usage() == WRITE_RESULTS ||
	  get_database()->usage() == WRITE_RESTART ) {
	assert((int)stateTimes.size() >= state);
	time = stateTimes[state-1];
      } else {
	assert(stateTimes.size() >= 1);
	time = stateTimes[0];
      }
    }
    return time;
  }

  std::pair<int, double> Region::get_max_time() const
  {
    if (!get_database()->is_input() &&
	get_database()->usage() != WRITE_RESULTS &&
	get_database()->usage() != WRITE_RESTART ) {
      return std::make_pair(currentState, stateTimes[0]);
    } else {
      // Cleanout the stateTimes vector and reload with current data in
      // case the database is being read and written at the same time.
      // This is rare, but is a supported use case.
      stateCount = 0;
      std::vector<double>().swap(stateTimes);
      DatabaseIO *db = (DatabaseIO*)get_database();
      db->get_step_times();

      int step = -1;
      double max_time = -1.0;
      for (int i=0; i < (int)stateTimes.size(); i++) {
	if (stateTimes[i] > max_time) {
	  step = i;
	  max_time = stateTimes[i];
	}
      }
      return std::make_pair(step+1, max_time);
    }
  }

  std::pair<int, double> Region::get_min_time() const 
  {
    if (!get_database()->is_input() &&
	get_database()->usage() != WRITE_RESULTS &&
	get_database()->usage() != WRITE_RESTART ) {
      return std::make_pair(currentState, stateTimes[0]);
    } else {
      // Cleanout the stateTimes vector and reload with current data in
      // case the database is being read and written at the same time.
      // This is rare, but is a supported use case.
      stateCount = 0;
      std::vector<double>().swap(stateTimes);
      DatabaseIO *db = (DatabaseIO*)get_database();
      db->get_step_times();

      int step = 0;
      double min_time = stateTimes[0];
      for (int i=1; i < (int)stateTimes.size(); i++) {
	if (stateTimes[i] < min_time) {
	  step = i;
	  min_time = stateTimes[i];
	}
      }
      return std::make_pair(step+1, min_time);
    }
  }

  double Region::begin_state(int state)
  {
    double time = 0.0;
    if (state > stateCount) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Requested state does not exist.\n"
	     << "       [" << get_database()->get_filename() << "]\n";
      IOSS_ERROR(errmsg);
    } else if (currentState != -1) {
      std::ostringstream errmsg;
      errmsg << "ERROR: State " << currentState
	     << " was not ended. Can not begin new state.\n"
	     << "       [" << get_database()->get_filename() << "]\n";
      IOSS_ERROR(errmsg);
    } else {
      assert(state <= stateCount);
      if (get_database()->is_input() ||
	  get_database()->usage() == WRITE_RESULTS ||
	  get_database()->usage() == WRITE_RESTART ) {
	assert((int)stateTimes.size() >= state);
	time = stateTimes[state-1];
      } else {
	assert(stateTimes.size() >= 1);
	time = stateTimes[0];
      }
      currentState = state;
      DatabaseIO *db = (DatabaseIO*)get_database();
      db->begin_state(this, state, time);
    }
    return time;
  }

  double Region::end_state(int state)
  {
    if (state != currentState) {
      std::ostringstream errmsg;
      errmsg << "ERROR: The current database state (" << currentState
	     << ") does not match the ending state (" << state << ").\n"
	     << "       [" << get_database()->get_filename() << "]\n";
      IOSS_ERROR(errmsg);
    }
    DatabaseIO *db = (DatabaseIO*)get_database();
    double time = 0.0;
    if (get_database()->is_input() ||
	get_database()->usage() == WRITE_RESULTS ||
	get_database()->usage() == WRITE_RESTART ) {
      assert((int)stateTimes.size() >= state);
      time = stateTimes[state-1];
    } else {
      assert(stateTimes.size() >= 1);
      time = stateTimes[0];
    }
    db->end_state(this, state, time);
    currentState = -1;
    return time;
  }

  bool Region::add(NodeBlock    *node_block)
  {
    // Check that region is in correct state for adding entities
    if (get_state() == STATE_DEFINE_MODEL) {
      // Add name as alias to itself to simplify later uses...
      add_alias(node_block);

      nodeBlocks.push_back(node_block);
      return true;
    } else {
      return false;
    }
  }

  bool Region::add(ElementBlock *element_block)
  {
    // Check that region is in correct state for adding entities
    if (get_state() == STATE_DEFINE_MODEL) {
      // Add name as alias to itself to simplify later uses...
      add_alias(element_block);

      // An input database defines these in the order matching the order
      // on the "file".  For output, we need to order based on the
      // "original_block_order" property and calculate the offset at that
      // point.  This is done in "end".

      if (get_database()->is_input()) {
	size_t nblocks = elementBlocks.size();
	int offset = 0;
	if (nblocks > 0) {
	  offset = elementBlocks[nblocks-1]->get_offset() +
	    elementBlocks[nblocks-1]->get_property("entity_count").get_int();
	}
	element_block->set_offset(offset);
      } else {
	// Check whether the "original_block_order" property exists on
	// this element block. If it isn't there, then add it with a
	// large value. If this is an element block read from the input
	// mesh, then the value will be updated during the
	// 'synchronize_id_and_name' function; if it is a block created
	// by the application during execution, then this value will
	// persist.  Add the property with a very large number such that
	// it will later be sorted after all "original" blocks.  Note
	// that it doesn't matter if two of the "new" blocks have the
	// same value since there is no ordering of new blocks that must
	// be preserved. (Use int_MAX/2 just to avoid some paranoia
	// about strange issue that might arise from int_MAX)
	if (!element_block->property_exists(orig_block_order())) {
	  element_block->property_add(Property(orig_block_order(), INT_MAX/2));
	}
      }
      elementBlocks.push_back(element_block);
      return true;
    } else {
      return false;
    }
  }

  bool Region::add(SideSet      *sideset)
  {
    // Check that region is in correct state for adding entities
    if (get_state() == STATE_DEFINE_MODEL) {
      // Add name as alias to itself to simplify later uses...
      add_alias(sideset);
      sideSets.push_back(sideset);
      return true;
    } else {
      return false;
    }
  }

  bool Region::add(NodeSet      *nodeset)
  {
    // Check that region is in correct state for adding entities
    if (get_state() == STATE_DEFINE_MODEL) {
      // Add name as alias to itself to simplify later uses...
      add_alias(nodeset);
      nodeSets.push_back(nodeset);
      return true;
    } else {
      return false;
    }
  }

  bool Region::add(CommSet      *commset)
  {
    // Check that region is in correct state for adding entities
    if (get_state() == STATE_DEFINE_MODEL) {
      // Add name as alias to itself to simplify later uses...
      add_alias(commset);
      commSets.push_back(commset);
      return true;
    } else {
      return false;
    }
  }

  bool Region::supports_field_type(EntityType fld_type) const
  {
    switch (fld_type) {
    case NODEBLOCK:
      return supports_nodal_fields();
    case ELEMENTBLOCK:
      return supports_element_fields();
    case NODESET:
      return supports_nodelist_fields();

    case SIDESET: // fall through
    case SIDEBLOCK:
      return supports_side_fields();
    case COMMSET:
      return false;
    default:
      return false;
    }
  }

  const NodeBlockContainer&  Region::get_node_blocks() const
  { return nodeBlocks; }

  const ElementBlockContainer& Region::get_element_blocks() const
  { return elementBlocks; }

  const SideSetContainer&  Region::get_sidesets() const
  { return sideSets; }

  const NodeSetContainer&  Region::get_nodesets() const
  { return nodeSets; }

  const CommSetContainer&  Region::get_commsets() const
  { return commSets; }

  bool Region::add_alias(const GroupingEntity *ge)
  {
    // Seeif an entity with this name already exists...
    std::string db_name = ge->name();
    const GroupingEntity *old_ge = get_entity(db_name);
    if (old_ge != NULL && ge != old_ge) {
      if (!((old_ge->type() == SIDEBLOCK &&     ge->type() == SIDESET) ||
	    (    ge->type() == SIDEBLOCK && old_ge->type() == SIDESET))) {
	int old_id = -1;
	int new_id = -1;
	if (old_ge->property_exists(id_str())) {
	  old_id = old_ge->get_property(id_str()).get_int();
	}
	if (ge->property_exists(id_str())) {
	  new_id = ge->get_property(id_str()).get_int();
	}
      
	std::ostringstream errmsg;
	errmsg << "\n\nERROR: Duplicate names detected.\n       The name '" << db_name << "' was found for both "
	       << old_ge->type_string() << " " << old_id << " and "
	       << ge->type_string() << " " << new_id
	       << ".\n       Names must be unique over all types in a finite element model.\n\n";
	IOSS_ERROR(errmsg);
      }
    }
    return add_alias(db_name, db_name);
  }

  bool Region::add_alias(const std::string &db_name, const std::string &alias)
  {
    // For use with the USTRING type in Sierra, create an uppercase
    // version of all aliases...
    std::string uname = uppercase(alias);
    if (uname != alias)
      aliases_.insert(IOAliasValuePair(uname, db_name));

    std::pair<AliasMap::iterator, bool> result = aliases_.insert(IOAliasValuePair(alias, db_name));
    return result.second;
  }

  std::string Region::get_alias(const std::string &alias) const
  {
    std::string ci_alias = uppercase(alias);
    AliasMap::const_iterator I = aliases_.find(ci_alias);
    if (I == aliases_.end()) {
      return "";
    } else {
      return (*I).second;
    }
  }

  int Region::get_aliases(const std::string& my_name, std::vector<std::string> &aliases) const
  {
    AliasMap::const_iterator I  = aliases_.begin();
    AliasMap::const_iterator IE = aliases_.end();
  
    size_t size = aliases.size();
    while (I != IE) {
      std::string alias = (*I).first;
      std::string base  = (*I).second;
      if (base == my_name) {
	aliases.push_back(alias);
      }
      ++I;
    }
    return static_cast<int>(aliases.size() - size);
  }

  const AliasMap & Region::get_alias_map() const
  {
    return aliases_;
  }

  GroupingEntity* Region::get_entity(const std::string& my_name, EntityType io_type) const
  {
    if (io_type == NODEBLOCK){
      return get_node_block(my_name);
    } else if (io_type == ELEMENTBLOCK) {
      return get_element_block(my_name);
    } else if (io_type == SIDESET) {
      return get_sideset(my_name);
    } else if (io_type == NODESET) {
      return get_nodeset(my_name);
    } else if (io_type == COMMSET) {
      return get_commset(my_name);
    } else if (io_type == SIDEBLOCK) {
      return get_sideblock(my_name);
    }
    return NULL;
  }

  GroupingEntity* Region::get_entity(const std::string& my_name) const
  {
    GroupingEntity *entity = NULL;
    entity = get_node_block(my_name);
    if (entity != NULL) { return entity;}
    entity = get_element_block(my_name);
    if (entity != NULL) { return entity;}
    entity = get_sideset(my_name);
    if (entity != NULL) { return entity;}
    entity = get_nodeset(my_name);
    if (entity != NULL) { return entity;}
    entity = get_commset(my_name);
    if (entity != NULL) { return entity;}
    entity = get_sideblock(my_name);
    if (entity != NULL) { return entity;}

    return entity;
  }

  NodeBlock*    Region::get_node_block(const std::string& my_name) const
  {
    const std::string db_name = get_alias(my_name);
    NodeBlock *ge = NULL;
    NodeBlockContainer::const_iterator i = nodeBlocks.begin();
    while (i != nodeBlocks.end()) {
      if ((*i)->name() == db_name) {
	ge = *i;
	break;
      }
      ++i;
    }
    return ge;
  }

  ElementBlock* Region::get_element_block(const std::string& my_name) const
  {
    const std::string db_name = get_alias(my_name);
    ElementBlock *ge = NULL;
    ElementBlockContainer::const_iterator i = elementBlocks.begin();
    while (i != elementBlocks.end()) {
      if ((*i)->name() == db_name) {
	ge = *i;
	break;
      }
      ++i;
    }
    return ge;
  }

  SideSet* Region::get_sideset(const std::string& my_name) const
  {
    const std::string db_name = get_alias(my_name);
    SideSet *ge = NULL;
    SideSetContainer::const_iterator i = sideSets.begin();
    while (i != sideSets.end()) {
      if ((*i)->name() == db_name) {
	ge = *i;
	break;
      }
      ++i;
    }
    return ge;
  }

  SideBlock* Region::get_sideblock(const std::string& my_name) const
  {
    SideBlock *ge = NULL;
    SideSetContainer::const_iterator i = sideSets.begin();
    while (i != sideSets.end()) {
      ge = (*i)->get_side_block(my_name);
      if (ge != NULL)
	break;
      ++i;
    }
    return ge;
  }

  NodeSet* Region::get_nodeset(const std::string& my_name) const
  {
    const std::string db_name = get_alias(my_name);
    NodeSet *ge = NULL;
    NodeSetContainer::const_iterator i = nodeSets.begin();
    while (i != nodeSets.end()) {
      if ((*i)->name() == db_name) {
	ge = *i;
	break;
      }
      ++i;
    }
    return ge;
  }

  CommSet* Region::get_commset(const std::string& my_name) const
  {
    const std::string db_name = get_alias(my_name);
    CommSet *ge = NULL;
    CommSetContainer::const_iterator i = commSets.begin();
    while (i != commSets.end()) {
      if ((*i)->name() == db_name) {
	ge = *i;
	break;
      }
      ++i;
    }
    return ge;
  }

  bool Region::is_valid_io_entity(const std::string& my_name, unsigned int io_type,
				  std::string *my_type) const
  {
    // Search all entities defined on this region for the name 'my_name'.
    // If found, then set 'type' (if non-NULL) to the type of the entity
    // (the 'type' values are from client code that was developed prior
    // to this function, so they are somewhat exodusII specific...).
    if ((io_type & NODEBLOCK) && get_node_block(my_name) != NULL) {
      if (my_type != NULL) *my_type = "NODE_BLOCK";
      return true;
    } else if ((io_type & ELEMENTBLOCK) && get_element_block(my_name) != NULL) {
      if (my_type != NULL) *my_type = "ELEMENT_BLOCK";
      return true;
    } else if ((io_type & SIDESET) && get_sideset(my_name) != NULL) {
      if (my_type != NULL) *my_type = "SURFACE";
      return true;
    } else if ((io_type & NODESET) && get_nodeset(my_name) != NULL) {
      if (my_type != NULL) *my_type = "NODESET";
      return true;
    } else if ((io_type & COMMSET) && get_commset(my_name) != NULL) {
      if (my_type != NULL) *my_type = "COMMSET";
      return true;
    }
    if (my_type != NULL) *my_type = "INVALID";
    return false;
  }

  // Retrieve the element block that contains the specified element
  // The 'local_id' is the local database id (1-based), not the global id.
  // returns NULL if no element block contains this element (local_id <= 0
  // or greater than number of elements in database)
  ElementBlock* Region::get_element_block(int local_id) const
  {
    ElementBlockContainer::const_iterator i = elementBlocks.begin();
    while (i != elementBlocks.end()) {
      if ((*i)->contains(local_id))
	return *i;
      ++i;
    }
    // Should not reach this point...
    std::ostringstream errmsg;
    errmsg << "Internal Program Error...Invalid local_id " << local_id << " specified";
    IOSS_ERROR(errmsg);
    return NULL;
  }

  Property Region::get_implicit_property(const std::string& my_name) const
  {
    if (my_name == "node_block_count")
      return Property(my_name, (int)nodeBlocks.size());

    if (my_name == "element_block_count")
      return Property(my_name, (int)elementBlocks.size());

    if (my_name == "side_set_count")
      return Property(my_name, (int)sideSets.size());

    if (my_name == "node_set_count")
      return Property(my_name, (int)nodeSets.size());

    if (my_name == "comm_set_count")
      return Property(my_name, (int)commSets.size());

    if (my_name == "state_count") {
      return Property(my_name, stateCount);
    }

    if (my_name == "current_state") {
      return Property(my_name, currentState);
    }

    if (my_name == "element_count") {
      int count = 0;
      ElementBlockContainer::const_iterator i = elementBlocks.begin();
      while (i != elementBlocks.end()) {
	count += (*i++)->get_property("entity_count").get_int();
      }
      return Property(my_name, count);
    }

    if (my_name == "node_count") {
      int count = 0;
      NodeBlockContainer::const_iterator i = nodeBlocks.begin();
      while (i != nodeBlocks.end()) {
	count += (*i++)->get_property("entity_count").get_int();
      }
      return Property(my_name, count);
    }

    if (my_name == "database_name") {
      std::string filename = get_database()->get_filename();
      return Property(my_name, filename);
    }

    else
      return GroupingEntity::get_implicit_property(my_name);
  }

  int Region::internal_get_field_data(const Field& field,
				      void *data, size_t data_size) const
  {
    return get_database()->get_field(this, field, data, data_size);
  }

  int Region::internal_put_field_data(const Field& field,
				      void *data, size_t data_size) const
  {
    return get_database()->put_field(this, field, data, data_size);
  }

  void Region::transfer_mesh_aliases(Region *to) const
  {
    // This routine transfers all relavant aliases from the 'this'
    // region and applies them to the 'to' file.

    // Iterate through list, [ returns <alias, base_entity_name> ], if
    // 'base_entity_name' is defined on the restart file, add 'alias' as
    // an alias for it...
    AliasMap::const_iterator I  = aliases_.begin();
    AliasMap::const_iterator IE = aliases_.end();

    while (I != IE) {
      std::string alias = (*I).first;
      std::string base  = (*I).second;
      if (alias != base && to->get_entity(base) != NULL) {
	to->add_alias(base, alias);
      }
      ++I;
    }
  }

  void Region::synchronize_id_and_name(const Region *from, bool sync_attribute_field_names)
  {
    // There is very little connection between an input (mesh) database
    // and an output (results/restart) database.  Basically, the entity
    // names are the same between the two files.  This works fine in the
    // case that an input database has 'generated' entity names of the
    // form 'block_10' or 'surface_32' since then the output database
    // can de-generate or decode the name and infer that the block
    // should have an id of 10 and the surface an id of 32.
    //
    // However, if alias or other renaming happens, then the output
    // block may have a name of the form 'fireset' and the underlying
    // database cannot infer that the id of the block should be 10.
    // Instead, it sets the id to an arbitrary number (1,2,...).  This
    // is annoying in the case of the results file since there is no
    // correspondence between the mesh numbering and the results
    // numbering. In the case of the restart output file, it can be
    // disastrous since when the file is used to restart the analysis,
    // there is no match between the mesh blocks and those found on the
    // restart file and the restart fails.
    //
    // So... We need to somehow ensure that the restart (and results)
    // files have the same ids.  To do this, we do the following:
    //
    // 1. The mesh database will set the property 'id' on input.
    // 2. The results/restart files will have a 'name' based on either
    //    the true name or an alias name.  Whichever, that alias will
    //    appear on the mesh database also, so we can query the mesh
    //    database aliases to get the entity.
    // 3. Once we have the entity, we can query the 'id' property and
    //    add the same 'id' property to the results/restart database.
    // 4. Also set the 'name' property to the base 'name' on the output file.
    // 5. Note that a property may already exist and must be removed
    //    before the 'correct' value is set.

    AliasMap::const_iterator I  = aliases_.begin();
    AliasMap::const_iterator IE = aliases_.end();

    while (I != IE) {
      std::string alias = (*I).first;
      std::string base  = (*I).second;
      ++I;

      if (alias == base) {

	// Query the 'from' database to get the entity (if any) referred
	// to by the 'alias'
	GroupingEntity *ge = from->get_entity(base);

	if (ge != NULL) {
	  // Get the entity from this region... Must be non-NULL
	  GroupingEntity *this_ge = get_entity(base);
	  if (this_ge == NULL) {
	    std::ostringstream errmsg;
	    errmsg << "INTERNAL ERROR: Could not find entity '" << base << "' in synchronize_id_and_name() "
		   << "                [" << get_database()->get_filename() << "]\n";
	    IOSS_ERROR(errmsg);
	  }

	  // See if there is an 'id' property...
	  if (ge->property_exists(id_str())) {
	    int id = ge->get_property(id_str()).get_int();

	    if (this_ge->property_exists(id_str())) {
	      // Remove the old property...
	      this_ge->property_erase(id_str());
	    }
	    // Set the new property
	    this_ge->property_add(Property(id_str(), id));
	  } else {
	    // No id, make sure the base name matches in both databases...
	    // There is always a 'name' property on an entity
	    if (this_ge->name() != base) {
	      this_ge->set_name(base);
	    }
	  }

	  // See if there is a 'original_element_type' property...
	  if (ge->property_exists(orig_elem_str())) {
	    std::string oes = ge->get_property(orig_elem_str()).get_string();

	    // Set the new property (erase if already exists; original file trumps...)
	    if (this_ge->property_exists(orig_elem_str())) {
	      this_ge->property_erase(orig_elem_str());
	    }
	    this_ge->property_add(Property(orig_elem_str(), oes));
	  }

	  // Specific to element blocks. Transfer the "original_block_order" property.
	  if (ge->property_exists(orig_block_order())) {
	    int offset = ge->get_property(orig_block_order()).get_int();
	    if (this_ge->property_exists(orig_block_order())) {
	      this_ge->property_erase(orig_block_order());
	    }
	    this_ge->property_add(Property(orig_block_order(), offset));
	  }

	  if (sync_attribute_field_names) {
	    // If there are any attribute fields, then copy those over
	    // to the new entity in order to maintain the same order
	    // since some codes access attributes by implicit order and
	    // not name... (typically, element blocks only)
	    Ioss::NameList attr_fields;
	    ge->field_describe(Ioss::Field::ATTRIBUTE, &attr_fields);
	    Ioss::NameList::const_iterator IF;
	    for (IF = attr_fields.begin(); IF != attr_fields.end(); ++IF) {
	      std::string field_name = *IF;
	      const Ioss::Field &field = ge->get_fieldref(field_name);
	      if (this_ge->field_exists(field_name)) {
		// If the field is already defined on the entity, make
		// sure that the attribute index matches...
		int index = field.get_index();
		const Ioss::Field &this_field = this_ge->get_fieldref(field_name);
		this_field.set_index(index);
	      } else {
		// If the field does not already exist, add it to the output node block
		this_ge->field_add(field);
	      }
	    }
	  }
	}
      }
    }

    while (I != IE) {
      std::string alias = (*I).first;
      std::string base  = (*I).second;
      ++I;

      if (alias != base) {
	GroupingEntity *ge = get_entity(base);
	if (ge != NULL) {
	  add_alias(base, alias);
	}
      }
    }
  }
}

  namespace {
    std::string uppercase(const std::string &my_name)
    {
      std::string s(my_name);
      std::transform(s.begin(), s.end(), s.begin(), toupper);
      return s;
    }
  }
