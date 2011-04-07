


#ifndef STK_IO_UTIL_IO_FIXTURE_HPP
#define STK_IO_UTIL_IO_FIXTURE_HPP


#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_io/IossBridge.hpp>


#include <Teuchos_RCP.hpp>

namespace stk {
namespace io {
namespace util {

class IO_Fixture {
  public:

    IO_Fixture();
    ~IO_Fixture();

    void create_output_mesh(
                            const std::string & base_exodus_filename,
                            bool  add_transient = true,
                            bool  add_all_fields = false
                           );

    void process_output_mesh( double time );


    void initialize_meta_data( Teuchos::RCP<stk::mesh::fem::FEMMetaData> meta_data );
    void serial_initialize_meta_data( const std::string & base_exodus_filename );
    void parallel_initialize_meta_data( const std::string & base_exodus_filename );

    void initialize_bulk_data( Teuchos::RCP<stk::mesh::BulkData> bulk_data );
    void serial_initialize_bulk_data( const std::string & base_exodus_filename );
    void parallel_initialize_bulk_data( const std::string & base_exodus_filename );

    stk::mesh::fem::FEMMetaData & meta_data()   { return *m_fem_meta_data; }
    stk::mesh::BulkData         & bulk_data()   { return *m_bulk_data; }
    Teuchos::RCP<Ioss::Region>    ioss_region() { return  m_ioss_output_region; }


  private:

    Teuchos::RCP<stk::mesh::fem::FEMMetaData>  m_fem_meta_data;
    Teuchos::RCP<stk::mesh::BulkData>          m_bulk_data;

    Teuchos::RCP<Ioss::Region>                 m_ioss_output_region;

    //disallow copy constructor and assignment operator
    IO_Fixture( const IO_Fixture & );
    IO_Fixture & operator = ( const IO_Fixture & );
};

}//namespace util
}//namespace io
}//namespace stk

#endif //STK_IO_UTIL_IO_FIXTURE_HPP


