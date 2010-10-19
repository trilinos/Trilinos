/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef Ioi_IOBroker_h
#define Ioi_IOBroker_h

#include <string>

#include <Ioi_Types.h>
#include <Ioi_Observer.h>

#include <Ioss_EntityType.h>

namespace Ioss {
  class Region;
}

namespace sierra {
  template<typename T> class Listener;

  namespace Frio {
    class IOSchedUser;
    class MeshReader;
    class ReadRestart;
    class MeshField;
    class IOBase;
  }
}
  
namespace Ioi {
  class IOBroker 
  {
  public:

    typedef std::vector<std::pair<sierra::Frio::IOSchedUser*,InputOutputType> > OutputList;
    typedef std::vector<sierra::Frio::MeshField*> FieldList;

    explicit IOBroker(const std::string&) {}
    virtual ~IOBroker() {}

    // Specify that this region is special and does not need
    // to read a restart file if the calculation is being restarted.
    virtual void restart_file_needs(Ioi::RestartNeed need) = 0;
    virtual Ioi::RestartNeed restart_file_needs() const = 0;

    //! Tell the IO databases that the model topology has been modified.
    virtual void topology_modified(unsigned type /* Fmwk::TopologyModified type */ ) = 0;
    virtual void set_cumulative_topology(unsigned type) = 0;

    virtual std::pair<double,int> restart(double time, bool is_automatic_restart)= 0;
    virtual double check_restart_time()= 0;
    virtual double adjust_time_step_for_output(const double& provisional_dt, const double& time)= 0;

    //! Begin a new database if at the correct time.
    virtual void initialize_output_request(double time, int step, double termination_time)= 0;

    //! Output transient data to results database.
    virtual void process_output_request(double time, int step, double termination_time)= 0;

    //! If an output database has the "sychronize" option set, force output on it.
    virtual void synchronize_output(InputOutputType type)= 0;

    //! Prepare interpolation fields (if any) and read in initial values.
    virtual void initialize_interpolation_request(double time)= 0;
    virtual void initialize_interpolation_request(Ioi::FieldRoleType role, double time)= 0;
    
    //! Provide interpolated data 
    virtual void process_interpolation_request()= 0;

    /*!
     * Force the specified restart or results IO system to write at the
     * next call to process_IO. If request_name is empty, force
     * all IO systems to write.Returns the number of databases written to.
     */
    virtual int force_output(InputOutputType type, const std::string& request_name=std::string(""),
		       const std::string& info=std::string(""))= 0;

    /*!
     * Immediately output the specified restart or results IO system.
     * If request_name is empty, force all IO systems of the specified type
     * to write. Returns the number of databases written to.
     */
    virtual int output_now(double time, int step, InputOutputType type,
		     const std::string& request_name=std::string(""))= 0;

    //! Read mesh, populate regions database, then delete mesh
    virtual void process_mesh(double restart_time)= 0;

    virtual std::string name_of_input_mesh() const= 0;
    
    virtual bool mesh_is_defined() const= 0;
    
    // Returns true if the entity name "entity_name" refers to an entity
    // of type "entity_type" on the input mesh.  The valid types for
    // entity_type are listed in Ioss_EntityType.h
    virtual bool is_mesh_entity(const std::string &entity_name, unsigned int entity_type) const= 0;

    virtual Ioss::Region const* io_mesh_region() const= 0;
    virtual Ioss::Region*       io_mesh_region()= 0;

    /*! Create a new named suffix field type with the provided suffices.
     *  field type size is given by suffices.size()
     *  Returns true if field type can be created; false if it already exists.
     */
    virtual bool create_named_suffix_field_type(const std::string &type_name, std::vector<std::string> &suffices) const=0;

    /*! Add a custom mapping between the field named 'field' and the Ioss::StorageType 'type'
     *  The StorageType must exist at the time the mapping is added.
     *  This can be used to override the default field to type mapping provided by the IO system.
     *  An example is a fuego species field can use the named suffix storage type to suffix the field
     *  components with the species name "Ynd_O2, Ynd_H2, Ynd_N2" instead of Ynd_1, Ynd_2, Ynd_3...
     *  Returns true if mapping successful.  False if type doesn't exist or if a mapping for this field
     *  already exists.
     */
    virtual bool add_field_type_mapping(const std::string &field, const std::string &type) const=0 ;

    /*! Query database capabilities...
     * If request_name is empty, queries all databases of the requested
     * type and returns true if *any* of them support the specified field type.
     */
    virtual bool supports_field_type(Ioss::EntityType fld_type, InputOutputType db_type,
			       const std::string& request_name=std::string("")) const= 0;
		   
    //! Methods for observer/publish/subscribe implementation:
    virtual void subscribe(sierra::Listener<Ioi::IOEvent>& listener)= 0;
    virtual void unsubscribe(sierra::Listener<Ioi::IOEvent>& listener)= 0;

    /*!
     * This routine provides a mechanism for an application to request
     * the writing of specific variables on the results databases.
     */
    virtual bool register_output_request(int io_type,
				   const std::string &request_name,
				   const std::string &field_name,
				   const std::string &fmwk_field_name,
				   Ioi::FieldType type,
				   const std::string &entity_name = "ALL")= 0;
  
    /*!
     * This routine provides a mechanism for an application to request
     * the reading and/or writing of specific variables on the mesh,
     * or restart databases.  An example is the reading of
     * element attributes. Note that this is not the mechanism to be
     * used for registering output results variables since this is
     * accessed at model input/definition time, use
     * register_output_request for that function.
     */
    virtual void register_field_request(const sierra::Frio::MeshField &mesh_field)= 0;
    
    /*!
     * Return a list of which databases will perform output
     * at the end of the current timestep. 
     */
    virtual const OutputList will_output() const= 0;
    virtual const OutputList will_output(InputOutputType type) const= 0;

    /*!
     * Return a reference to the current list of outputs for
     * this region.
     */
    virtual const OutputList &outputs() const= 0;

    /*!
     * Return a pointer to the IOBase portion of an output request of
     * the specified type and with the specified name.  If name is
     * empty, then return the first output request of the specified type
     */
    virtual const sierra::Frio::IOBase* output(InputOutputType type,
				 const std::string& request_name=std::string("")) const = 0;
    
    /*! Inform all output databases that 'time' may possibly be reset to
     * a value less than the current time.  The databases should do what
     * is required to handle this (typically, close the current file and
     * output subsequent output to a different file.
     */
    virtual void reset_deja_vu() = 0;

    /*! Delete all output databases (results, heartbeat, history, *
     *  restart). Typically this will be called from the IOBroker *
     *  destructor; however, it may also be called if an exception is
     *  thrown so we can try to flush all data to disk...
     */
    virtual void delete_all_output_requests() = 0;
    
    virtual sierra::Frio::MeshReader *get_mesh() const = 0;
    virtual void add_mesh(sierra::Frio::MeshReader *io_region_mesh)= 0;
    virtual void replace_mesh(Ioss::Region *new_mesh, const std::string &model_name,
			bool owns_io_region = true, bool delete_mesh = true)= 0;
    
    virtual void output_list_add(InputOutputType type, sierra::Frio::IOSchedUser *)= 0;
    virtual void restart_list_add(sierra::Frio::ReadRestart *)= 0;

    virtual FieldList::const_iterator additional_fields_begin() const= 0;
    virtual FieldList::const_iterator additional_fields_end() const= 0;

    virtual OutputList::const_iterator output_list_begin() const = 0;
    virtual OutputList::const_iterator output_list_end() const = 0;

  };
}

#endif // Ioi_IOBroker_h
