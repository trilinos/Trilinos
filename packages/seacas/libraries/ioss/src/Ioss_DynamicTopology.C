// Copyright(C) 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_Assembly.h"
#include "Ioss_Blob.h"
#include "Ioss_CodeTypes.h"
#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_DynamicTopology.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_EntityBlock.h"
#include "Ioss_EntityType.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h"
#include "Ioss_FileInfo.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_StructuredBlock.h"

#include <climits>
#include <cstddef>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <string>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <assert.h>

#include "Ioss_ParallelUtils.h"

namespace Ioss {

void DynamicTopologyObserver::check_region() const
{
  if(nullptr == m_region) {
    std::ostringstream errmsg;
    fmt::print(errmsg,
               "ERROR: A region has not been registered with the "
               "Dynamic Topology Observer.\n\n");
    IOSS_ERROR(errmsg);
  }
}

void DynamicTopologyObserver::register_region(Region *region)
{
  if(nullptr != region && nullptr != m_region && region != m_region) {
    std::ostringstream errmsg;
    fmt::print(errmsg,
               "ERROR: Attempt to re-register different region on "
               "Dynamic Topology Observer.\n\n");
    IOSS_ERROR(errmsg);
  }

  m_region = region;
}

void DynamicTopologyObserver::register_notifier(DynamicTopologyNotifier *notifier)
{
  if(nullptr != notifier && nullptr != m_notifier && notifier != m_notifier) {
    std::ostringstream errmsg;
    fmt::print(errmsg,
               "ERROR: Attempt to re-register different notifier on "
               "Dynamic Topology Observer.\n\n");
    IOSS_ERROR(errmsg);
  }

  m_notifier = notifier;
}

void DynamicTopologyObserver::set_cumulative_topology_modification(unsigned int type)
{ m_cumulativeTopologyModification = type; }

unsigned int DynamicTopologyObserver::get_cumulative_topology_modification() const
{ return m_cumulativeTopologyModification; }

unsigned int DynamicTopologyObserver::get_topology_modification() const
{ return m_topologyModification; }

void DynamicTopologyObserver::set_topology_modification_nl(unsigned int type)
{
  m_topologyModification |= type;
  m_cumulativeTopologyModification |= type;
}

void DynamicTopologyObserver::set_topology_modification(unsigned int type)
{
  if(!(m_topologyModification & type)) {
    set_topology_modification_nl(type);

    if(nullptr != m_notifier) {
      for(auto observer : m_notifier->get_observers()) {
        observer->set_topology_modification_nl(type);
      }
    }
  }
}

void DynamicTopologyObserver::reset_topology_modification()
{
  m_topologyModification = TOPOLOGY_SAME;
}

void DynamicTopologyObserver::reset_topology_modification_all()
{
  if(m_topologyModification != TOPOLOGY_SAME) {
    reset_topology_modification();

    if(nullptr != m_notifier) {
      for(auto observer : m_notifier->get_observers()) {
        observer->reset_topology_modification();
      }
    }
  }
}

bool DynamicTopologyObserver::is_topology_modified() const
{ return m_topologyModification != TOPOLOGY_SAME; }

const ParallelUtils &DynamicTopologyObserver::util() const
{
  check_region();
  return m_region->get_database()->util();
}

void DynamicTopologyObserver::synchronize_topology_modified_flags()
{
  check_region();
  int num_processors = m_region->get_database()->parallel_size();
  // Synchronize the topology flags between all processors in case
  // it has not been set consistently.
  if (num_processors > 1) {
    static unsigned int buffer[2];
    buffer[0] = m_cumulativeTopologyModification;
    buffer[1] = m_topologyModification;

    util().attribute_reduction(2*sizeof(unsigned int), reinterpret_cast<char*>(buffer));

    m_cumulativeTopologyModification = buffer[0];
    m_topologyModification = buffer[1];
  }
}

int DynamicTopologyObserver::get_cumulative_topology_modification_field()
{
  check_region();
  const std::string variable_name = topology_modification_change_name();

  int ivalue = 0;

  if (m_region->field_exists(variable_name)) {
    Field topo_field = m_region->get_field(variable_name);
    if (topo_field.get_type() == Field::INTEGER) {
      m_region->get_field_data(variable_name, &ivalue, sizeof(int));
    } else {
      double value;
      m_region->get_field_data(variable_name, &value, sizeof(double));
      ivalue = (int)value;
    }
  }

  int num_processors = m_region->get_database()->parallel_size();
  // Synchronize the value between all processors in case
  // it has not been set consistently.
  if (num_processors > 1) {
    unsigned int buffer[1];
    buffer[0] = ivalue;

    util().attribute_reduction(sizeof(unsigned int), reinterpret_cast<char*>(buffer));

    ivalue = (int)buffer[0];
  }

  m_cumulativeTopologyModification = ivalue;

  return ivalue;
}

void DynamicTopologyObserver::define_model()
{

}

void DynamicTopologyObserver::write_model()
{

}

void DynamicTopologyObserver::define_transient()
{

}


DynamicTopologyBroker* DynamicTopologyBroker::broker()
{
  static DynamicTopologyBroker broker_;
  return &broker_;
}

void DynamicTopologyBroker::register_model(const std::string& model_name)
{
  auto iter = m_notifiers.find(model_name);
  if(iter != m_notifiers.end()) {
    return;
  }

  m_notifiers[model_name] = std::make_shared<DynamicTopologyNotifier>(model_name);
}

std::shared_ptr<DynamicTopologyNotifier> DynamicTopologyBroker::get_notifier(const std::string& model_name) const
{
  auto iter = m_notifiers.find(model_name);
  if(iter != m_notifiers.end()) {
    return iter->second;
  }

  return {};
}

std::vector<std::shared_ptr<DynamicTopologyObserver>> DynamicTopologyBroker::get_observers(const std::string& model_name) const
{
  std::vector<std::shared_ptr<DynamicTopologyObserver>> observers;

  auto notifier = get_notifier(model_name);

  if(notifier) {
    return notifier->get_observers();
  }

  return observers;
}

void DynamicTopologyBroker::remove_model(const std::string& model_name)
{
  auto iter = m_notifiers.find(model_name);
  if(iter != m_notifiers.end()) {
    m_notifiers.erase(iter);
  }
}

void DynamicTopologyBroker::clear_models()
{
  m_notifiers.clear();
}

void DynamicTopologyBroker::register_observer(const std::string& model_name,
                                              std::shared_ptr<DynamicTopologyObserver> observer)
{
  auto notifier = get_notifier(model_name);

  if(!notifier) {
    register_model(model_name);
    notifier = get_notifier(model_name);
  }

  notifier->register_observer(observer);
}

void DynamicTopologyBroker::register_observer(const std::string& model_name,
                                              std::shared_ptr<DynamicTopologyObserver> observer,
                                              Region& region)
{
  region.register_mesh_modification_observer(observer);
  register_observer(model_name, observer);
}

void DynamicTopologyBroker::reset_topology_modification(const std::string& model_name)
{
  auto notifier = get_notifier(model_name);

  if(!notifier) return;

  notifier->reset_topology_modification();
}

void DynamicTopologyBroker::set_topology_modification(const std::string& model_name, unsigned int type)
{
  auto notifier = get_notifier(model_name);

  if(!notifier) return;

  notifier->set_topology_modification(type);
}


struct DynamicTopologyObserverCompare {
  bool operator()(const std::shared_ptr<DynamicTopologyObserver> & lhs,
                  const std::shared_ptr<DynamicTopologyObserver> & rhs) const {
    assert(lhs && (lhs->get_region() != nullptr));
    assert(rhs && (rhs->get_region() != nullptr));
    return (lhs->get_region() < rhs->get_region());
  }
};

void DynamicTopologyNotifier::register_observer(std::shared_ptr<DynamicTopologyObserver> observer)
{
  observer->register_notifier(this);
  m_observers.push_back(observer);
  std::sort(m_observers.begin(), m_observers.end(), DynamicTopologyObserverCompare());
}

void DynamicTopologyNotifier::unregister_observer(std::shared_ptr<DynamicTopologyObserver> observer)
{
  auto iter = std::find(m_observers.begin(), m_observers.end(), observer);
  if (iter != m_observers.end()) {
    (*iter)->register_notifier(nullptr);
    m_observers.erase(iter);
  }
}

void DynamicTopologyNotifier::reset_topology_modification()
{
  for(std::shared_ptr<DynamicTopologyObserver>& observer : m_observers) {
    observer->reset_topology_modification();
  }
}

void DynamicTopologyNotifier::set_topology_modification(unsigned int type)
{
  for(std::shared_ptr<DynamicTopologyObserver>& observer : m_observers) {
    observer->set_topology_modification(type);
  }
}


DynamicTopologyFileControl::DynamicTopologyFileControl(Region *region, unsigned int fileCyclicCount,
                                                       IfDatabaseExistsBehavior &ifDatabaseExists,
                                                       unsigned int &dbChangeCount)
  : m_region(region)
  , m_fileCyclicCount(fileCyclicCount)
  , m_ifDatabaseExists(ifDatabaseExists)
  , m_dbChangeCount(dbChangeCount)
{
  if(nullptr == region) {
    std::ostringstream errmsg;
    fmt::print(errmsg, "ERROR: null region passed in as argument to DynamicTopologyFileControl");
    IOSS_ERROR(errmsg);
  }

  m_ioDB   = region->get_property("base_filename").get_string();
  m_dbType = region->get_property("database_type").get_string();
}

const ParallelUtils &DynamicTopologyFileControl::util() const
{
  return m_region->get_database()->util();
}

bool DynamicTopologyFileControl::file_exists(const std::string &filename,
                                             const std::string &db_type,
                                             Ioss::DatabaseUsage db_usage)
{
  bool exists = false;
  int par_size = m_region->get_database()->parallel_size();
  int par_rank = m_region->get_database()->parallel_rank();
  bool is_parallel = par_size > 1;
  std::string full_filename = filename;
  if (is_parallel && db_type == "exodusII" && db_usage != Ioss::WRITE_HISTORY) {
    full_filename = Ioss::Utils::decode_filename(filename, par_rank, par_size);
  }

  if (!is_parallel || par_rank == 0) {
    // Now, see if this file exists...
    // Don't want to do a system call on all processors since it can take minutes
    // on some of the larger machines, filesystems, and processor counts...
    Ioss::FileInfo file = Ioss::FileInfo(full_filename);
    exists = file.exists();
  }

  if (is_parallel) {
    int iexists = exists ? 1 : 0;
    util().broadcast(iexists, 0);
    exists = iexists == 1;
  }
  return exists;
}

std::string DynamicTopologyFileControl::get_unique_filename(Ioss::DatabaseUsage db_usage)
{
  std::string filename = m_ioDB;

  do {
    // Run this loop at least once for all files.  If this is an automatic
    // restart, then make sure that the generated file does not already exist,
    // so keep running the loop until we generate a filename that doesn't exist...
    std::ostringstream tmp_filename;
    tmp_filename << m_ioDB;

    // Don't append the "-s000X" the first time in case the base filename doesn't
    // exist -- we want write to the name specified by the user if at all possible and
    // once that exists, then start adding on the suffix...
    if (m_dbChangeCount > 1) {
      tmp_filename << "-s" << std::setw(4) << std::setfill('0') << m_dbChangeCount;
    }
    filename = tmp_filename.str();
    ++m_dbChangeCount;
  } while(file_exists(filename, m_dbType, db_usage));
  --m_dbChangeCount;
  return filename;
}

std::string DynamicTopologyFileControl::construct_database_filename(int& step, Ioss::DatabaseUsage db_usage)
{
  // Filename will be of the form -- ioDB-sxxxx where xxxx is step
  // number.  Assume maximum of 9999 steps (will do more, but won't have
  // good lineup of step numbers.
  // Check database for validity (filename and a type)
  if(m_ioDB.empty() || m_dbType.empty())
  {
    std::string error_message;
    if(m_dbType.empty())
      error_message += "The database TYPE has not been defined\n";

    if(m_ioDB.empty())
    {
      error_message += "The database FILENAME has not been defined\n";
    }
    std::ostringstream errmsg;
    fmt::print(errmsg, fmt::runtime(error_message));
    IOSS_ERROR(errmsg);
  }
  assert(!m_ioDB.empty());
  assert(!m_dbType.empty());
  std::string filename = m_ioDB;
  if(m_fileCyclicCount > 0)
  {
    // In this mode, we close the old file and open a new file
    // every time this is called. The file suffix cycles through
    // the first fileCyclicCount'th entries in A,B,C,D,E,F,...
    if(step == 0)
      step++;

    static std::string suffix = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::string tmp = "-" + suffix.substr((step - 1) % m_fileCyclicCount, 1);
    filename += tmp;
    m_properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_OVERWRITE));
  }
  else
  {
    if(m_region->model_is_written())
    {
      // After the initial open, we want to add suffix if the topology changes
      // during the run
      m_ifDatabaseExists = Ioss::DB_ADD_SUFFIX_OVERWRITE;
    }

    // Handle complications of DB_APPEND mode...
    // If in DB_APPEND mode, then we don't output metadata
    // information, so some knowledge is needed at this level if
    // we are appending.  If user specified APPEND, but the file
    // doesn't yet exist OR it does exist and we are not
    // restarting, change the mode to OVERWRITE.
    // 0. Must be restarting; either manual or automatic.
    std::shared_ptr<DynamicTopologyObserver> observer = m_region->get_mesh_modification_observer();

    if(m_ifDatabaseExists == Ioss::DB_APPEND)
    {
      if(!observer->is_restart_requested())
      {
        // Not restarting
        m_ifDatabaseExists = Ioss::DB_OVERWRITE;
      }
      else if(!file_exists(m_ioDB, m_dbType, db_usage))
      {
        m_ifDatabaseExists = Ioss::DB_OVERWRITE;
      }
    }
    if(step > 1 || (m_dbChangeCount > 1))
    {
      // Use the !is_input_event test since restart input files already have the
      // -s000x extension...
      if(m_ifDatabaseExists == Ioss::DB_APPEND)
      {
        std::ostringstream tmp_filename;
        tmp_filename << m_ioDB;
        filename = m_ioDB;
        if(m_dbChangeCount > 1)
        {
          tmp_filename << "-s" << std::setw(4) << std::setfill('0') << m_dbChangeCount;
        }
        size_t inc = 0;
        while(file_exists(tmp_filename.str(), m_dbType, db_usage))
        {
          filename = tmp_filename.str();
          tmp_filename.clear();
          tmp_filename.str("");
          tmp_filename << m_ioDB << "-s" << std::setw(4) << std::setfill('0') << m_dbChangeCount + (++inc);
        }
        if(inc > 0)
        {
          m_dbChangeCount += (inc - 1);
        }
        else
        {
          m_ifDatabaseExists = Ioss::DB_OVERWRITE;
        }
      }
      else if(m_ifDatabaseExists == Ioss::DB_ADD_SUFFIX)
      {
        filename = get_unique_filename(db_usage);
      }
      else if(m_ifDatabaseExists == Ioss::DB_ADD_SUFFIX_OVERWRITE)
      {
        if(m_dbChangeCount > 0)
        {
          std::ostringstream tmp_filename;
          tmp_filename << m_ioDB << "-s" << std::setw(4) << std::setfill('0') << ++m_dbChangeCount;
          filename = tmp_filename.str();
        }
        else
        {
          filename = m_ioDB;
        }
      }
      else
      {
        filename = m_ioDB;
      }
    }
    else if(m_ifDatabaseExists == Ioss::DB_ADD_SUFFIX)
    {
      filename = get_unique_filename(db_usage);
    }
    else
    {
      filename = m_ioDB;
    }

    m_properties.add(Ioss::Property("APPEND_OUTPUT", m_ifDatabaseExists));
    // A little complicated on deciding whether we are actually
    // overwriting the database.  The 'validate' routine for Results and
    // History will call create_database once the parser block is
    // ended. This routine will then create the database and the
    // ioRegion_. However, the database will not really be opened or
    // written to at this time.  If the code is auto-restarting, then it will
    // detect that the database exists and create a database with the
    // -s000x extension.
    // At this point, we need to skip the 'abort_if_exists' test if we
    // are in this routine from the 'validate' and we are restarting
    // since we won't really write to the file.  So, the cases where we
    // *don't* check are:
    // -- is_input_event(db_usage)
    // -- ifExists_ == DB_OVERWRITE || DB_ADD_SUFFIX_OVERWRITE || DB_APPEND
    // -- is_automatic_restart() && step == 0 (coming from validate)
    if(m_ifDatabaseExists != DB_OVERWRITE &&
       m_ifDatabaseExists != DB_APPEND &&
       m_ifDatabaseExists != DB_ADD_SUFFIX_OVERWRITE &&
       !(step == 0 && observer->is_automatic_restart()))
    {
      abort_if_exists(filename, m_dbType, db_usage);
    }
  }
  return filename;
}

bool DynamicTopologyFileControl::abort_if_exists(const std::string &filename,
                                                 const std::string &db_type,
                                                 Ioss::DatabaseUsage db_usage)
{
  // Check whether file with same name as database already exists.  If so,
  // print error message and stop...
  // At the current time, only check on processor 0 and assume if it doesn't exist
  // there, then it doesn't exist on other processors.  Or, if it doesn't exist on
  // processor 0, then it doesn't matter if it doesn't exist on other processors
  // since we don't have all pieces...

  bool exists = file_exists(filename, db_type, db_usage);
  if (exists) {
    std::ostringstream errmsg;
    fmt::print(errmsg,
               "ERROR: The database file named '{} exists"
               "and would be overwritten if the code continued.\n\n"
               "Input options specified that this file *not* be overwritten,\n"
               "\tso you must rename or remove this file and restart the code.\n",
               filename);
    IOSS_ERROR(errmsg);
  }
  return exists;
}

Ioss::DatabaseIO * DynamicTopologyFileControl::clone_output_database(int steps)
{
  auto current_db = m_region->get_database();

  if (current_db->is_input())
    return nullptr;

  const Ioss::PropertyManager& current_properties = current_db->get_property_manager();
  Ioss::NameList names;
  current_properties.describe(&names);

  // Iterate through properties and transfer to new output database...
  Ioss::NameList::const_iterator I;
  for (I = names.begin(); I != names.end(); ++I) {
    if (!current_properties.exists(*I))
      m_properties.add(current_properties.get(*I));
  }

  auto db_usage = current_db->usage();

  std::string filename = construct_database_filename(steps, db_usage);

  Ioss::DatabaseIO *db = Ioss::IOFactory::create(m_dbType, filename, db_usage,
                                                 current_db->util().communicator(),
                                                 m_properties);

  if (nullptr == db) {
    std::ostringstream errmsg;
    fmt::print(errmsg,
               "ERROR: unable to create output database named '{}'"
               " of type '{}'", filename, m_dbType);
    IOSS_ERROR(errmsg);
  }

  assert(db != nullptr);
  if(!db->ok(true)) {
    std::ostringstream errmsg;
    fmt::print(errmsg,
               "ERROR: unable to validate output database named '{}'"
               " of type '{}'", filename, m_dbType);
    IOSS_ERROR(errmsg);
  }

  db->set_field_separator(current_db->get_field_separator());
  db->set_surface_split_type(current_db->get_surface_split_type());
  db->set_maximum_symbol_length(current_db->maximum_symbol_length());
  db->set_int_byte_size_api(current_db->int_byte_size_data_size());

  return db;
}

template<typename T>
void update_database_for_grouping_entities(const T& container, Ioss::DatabaseIO *db)
{
  for(auto * entity : container) {
    Ioss::GroupingEntity* ge = dynamic_cast<Ioss::GroupingEntity*>(entity);
    assert(ge != nullptr);

    if(ge->type() == Ioss::SIDESET) {
      Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(ge);
      assert(sset != nullptr);

      sset->reset_database(db);
      const auto &sblocks = sset->get_side_blocks();
      for (const auto &sblock : sblocks) {
        sblock->reset_database(db);
      }
    } else {
      ge->reset_database(db);
    }
  }
}

bool DynamicTopologyFileControl::replace_output_database(Ioss::DatabaseIO *db)
{
  auto current_db = m_region->get_database();

  if (current_db->is_input())
    return false;

  current_db->finalize_database();
  current_db->closeDatabase();
  delete current_db;

  m_region->reset_database(db);
  db->set_region(m_region);

  update_database_for_grouping_entities(m_region->get_node_blocks(), db);
  update_database_for_grouping_entities(m_region->get_edge_blocks(), db);
  update_database_for_grouping_entities(m_region->get_face_blocks(), db);
  update_database_for_grouping_entities(m_region->get_element_blocks(), db);
  update_database_for_grouping_entities(m_region->get_sidesets(), db);
  update_database_for_grouping_entities(m_region->get_nodesets(), db);
  update_database_for_grouping_entities(m_region->get_edgesets(), db);
  update_database_for_grouping_entities(m_region->get_facesets(), db);
  update_database_for_grouping_entities(m_region->get_elementsets(), db);
  update_database_for_grouping_entities(m_region->get_commsets(), db);
  update_database_for_grouping_entities(m_region->get_structured_blocks(), db);
  update_database_for_grouping_entities(m_region->get_assemblies(), db);
  update_database_for_grouping_entities(m_region->get_blobs(), db);

  return true;
}

void DynamicTopologyFileControl::clone_and_replace_output_database(int steps)
{
  auto db = clone_output_database(steps);

  if(nullptr != db)
    replace_output_database(db);
}

void DynamicTopologyFileControl::add_output_database_group(int steps)
{
  auto current_db = m_region->get_database();

  std::ostringstream oss;
  oss << group_prefix();
  oss << m_dbChangeCount;

  current_db->release_memory();
  current_db->open_root_group();
  current_db->create_subgroup(oss.str());

  m_dbChangeCount++;
}

}




