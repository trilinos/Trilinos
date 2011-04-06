
#include <stk_io/util/IO_Fixture.hpp>
#include <stk_io/util/UseCase_mesh.hpp>


namespace stk {
namespace io {
namespace util {

IO_Fixture::IO_Fixture()
  : m_fem_meta_data(NULL)
  , m_bulk_data(NULL)
  , m_ioss_output_region(NULL)
{}

IO_Fixture::~IO_Fixture() {}


void IO_Fixture::create_output_mesh(
                                    const std::string & base_exodus_filename,
                                    bool  add_transient,
                                    bool  add_all_fields
                                   )
{

}

void IO_Fixture::process_output_mesh( double time )
{

}


void IO_Fixture::initialize_meta_data( Teuchos::RCP<stk::mesh::fem::FEMMetaData> meta_data )
{

}

void IO_Fixture::serial_initialize_meta_data( const std::string & base_exodus_filename )
{

}

void IO_Fixture::parallel_initialize_meta_data( const std::string & base_exodus_filename )
{

}

void IO_Fixture::initialize_bulk_data( Teuchos::RCP<stk::mesh::BulkData> bulk_data )
{

}


void IO_Fixture::serial_initialize_bulk_data( const std::string & base_exodus_filename )
{

}

void IO_Fixture::parallel_initialize_bulk_data( const std::string & base_exodus_filename )
{

}


}//namespace util
}//namespace io
}//namespace stk

