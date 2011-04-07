#include <stk_io/util/IO_Fixture.hpp>
#include <stk_io/util/UseCase_mesh.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

namespace stk {
namespace io {
namespace util {

IO_Fixture::IO_Fixture(stk::ParallelMachine comm)
  : m_comm(comm)
  , m_fem_meta_data(NULL)
  , m_bulk_data(NULL)
  , m_ioss_input_region(NULL)
  , m_ioss_output_region(NULL)
  , m_mesh_type()
  , m_mesh_data()
{}

IO_Fixture::~IO_Fixture() {}

void IO_Fixture::create_output_mesh(
                                    const std::string & base_exodus_filename,
                                    bool  add_transient,
                                    bool  add_all_fields
                                   )
{
  // TODO: Check that the meta-data has a coordinates field and IO parts defined.

  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(
      "exodusII",
      base_exodus_filename,
      Ioss::WRITE_RESULTS,
      bulk_data().parallel()
      );

  ThrowErrorMsgIf(dbo == NULL || !dbo->ok(),
      "ERROR: Could not open results database '" << base_exodus_filename <<
      "' of type 'exodusII'");

  // NOTE: 'm_ioss_output_region' owns 'dbo' pointer at this time
  m_ioss_output_region = Teuchos::rcp(new Ioss::Region(dbo, "results_output"));

  /* Given the newly created Ioss::Region 'm_ioss_output_region', define the
   * model corresponding to the stk::mesh 'bulk_data'.  If the
   * optional 'input_region' is passed as an argument, then
   * synchronize all names and ids found on 'input_region' to the
   * output region 'm_ioss_output_region'.  The routine will query all parts
   * in 'bulk_data' and if they are io_parts (define by the existance
   * of the IOPartAttribute attribute on the part), then a
   * corresponding Ioss entity will be defined.  This routine only
   * deals with the non-transient portion of the model; no transient
   * fields are defined at this point.
   */
  stk::io::define_output_db( *m_ioss_output_region, bulk_data(), m_ioss_input_region.get() );

  /* Given an Ioss::Region 'm_ioss_output_region' which has already had its
   * metadata defined via 'define_output_db()' call; transfer all bulk
   * data (node coordinates, element connectivity, ...) to the
   * output database that corresponds to this Ioss::Region. At
   * return, all non-transient portions of the output database will
   * have been output.
   */
  stk::io::write_output_db( *m_ioss_output_region, bulk_data());

  if (add_transient) {
    m_ioss_output_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    // Special processing for nodeblock (all nodes in model)...
    stk::io::ioss_add_fields(
                              meta_data().universal_part(),
                              meta_data().node_rank(),
                              m_ioss_output_region->get_node_blocks()[0],
                              Ioss::Field::TRANSIENT,
                              add_all_fields
                            );

    const stk::mesh::PartVector & all_parts = meta_data().get_parts();
    for ( stk::mesh::PartVector::const_iterator ip = all_parts.begin();
          ip != all_parts.end();
          ++ip
        )
    {
      stk::mesh::Part & part = **ip;

      // Check whether this part should be output to results database.
      if (stk::io::is_part_io_part(part)) {
        // Get Ioss::GroupingEntity corresponding to this part...
        Ioss::GroupingEntity *entity = m_ioss_output_region->get_entity(part.name());
        if (entity != NULL && entity->type() == Ioss::ELEMENTBLOCK) {
          stk::io::ioss_add_fields(
                                    part,
                                    part.primary_entity_rank(),
                                    entity,
                                    Ioss::Field::TRANSIENT,
                                    add_all_fields
                                  );
        }
      }
    }
    m_ioss_output_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }
}

void IO_Fixture::add_timestep_to_output_mesh( double time )
{
  ThrowErrorMsgIf( Teuchos::is_null(m_ioss_output_region),
                   "Please call create_output_mesh before add_timestep_to_output_mesh" );
  m_ioss_output_region->begin_mode(Ioss::STATE_TRANSIENT);
  int out_step = m_ioss_output_region->add_state(time);
  stk::io::util::process_output_request(*m_ioss_output_region, bulk_data(), out_step);
  m_ioss_output_region->end_mode(Ioss::STATE_TRANSIENT);
}

void IO_Fixture::set_meta_data( Teuchos::RCP<stk::mesh::fem::FEMMetaData> arg_meta_data )
{
  ThrowErrorMsgIf( !Teuchos::is_null(m_fem_meta_data),
                   "Meta data already initialized" );
  m_fem_meta_data = arg_meta_data;
}

void IO_Fixture::set_bulk_data( Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data )
{
  ThrowErrorMsgIf( !Teuchos::is_null(m_bulk_data),
                   "Bulk data already initialized" );
  m_bulk_data = arg_bulk_data;
}

void IO_Fixture::initialize_meta_data( const std::string & base_filename, const std::string & mesh_type)
{
  ThrowErrorMsgIf( !Teuchos::is_null(m_fem_meta_data),
                   "Meta data already initialized" );
  ThrowErrorMsgIf( !Teuchos::is_null(m_ioss_input_region),
                   "Input region was already initialized");

  m_fem_meta_data = Teuchos::rcp( new stk::mesh::fem::FEMMetaData());
  m_mesh_type = mesh_type;

  std::string no_working_dir = "";

  Ioss::Init::Initializer init_db;

  stk::io::util::create_input_mesh(
                                    m_mesh_type,
                                    base_filename,
                                    no_working_dir,
                                    m_comm,
                                    meta_data(),
                                    m_mesh_data
                                  );

  // TODO: Restore this once m_mesh_data is fixed
  //m_ioss_input_region = Teuchos::rcp( m_mesh_data.m_region );
}

void IO_Fixture::initialize_bulk_data()
{
  ThrowErrorMsgIf( !Teuchos::is_null(m_bulk_data),
                   "Bulk data already initialized" );

  // TODO: Probable better to check m_ioss_input_region once that's fixed
  ThrowErrorMsgIf( m_mesh_type == "",
    "Can only use this method if meta-data was initialized with initialize_meta_data");

  m_bulk_data = Teuchos::rcp( new stk::mesh::BulkData(stk::mesh::fem::FEMMetaData::get_meta_data(meta_data()), m_comm));

  stk::io::util::populate_bulk_data(
                                    bulk_data(),
                                    m_mesh_data,
                                    m_mesh_type
                                   );
}

void IO_Fixture::set_input_ioss_region( Teuchos::RCP<Ioss::Region> input_region )
{
  ThrowErrorMsgIf( !Teuchos::is_null(m_ioss_input_region),
                   "Input region was already initialized");

  m_ioss_input_region = input_region;
}

} // namespace util
} // namespace io
} // namespace stk
