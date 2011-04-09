


#ifndef STK_IO_UTIL_IO_FIXTURE_HPP
#define STK_IO_UTIL_IO_FIXTURE_HPP


#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/util/UseCase_mesh.hpp>

#include <stk_util/environment/ReportHandler.hpp>


#include <Teuchos_RCP.hpp>
#include <string>

namespace stk {
namespace io {
namespace util {

class IO_Fixture {
  public:

    typedef stk::mesh::Field< double, stk::mesh::Cartesian> coord_field_type;

    IO_Fixture(stk::ParallelMachine comm);
    ~IO_Fixture();

    void create_output_mesh(
                            const std::string & base_exodus_filename,
                            bool  add_transient = true,
                            bool  add_all_fields = false
                           );

    void process_output_mesh( double time );


    void initialize_meta_data( Teuchos::RCP<stk::mesh::fem::FEMMetaData> arg_meta_data );
    void initialize_bulk_data( Teuchos::RCP<stk::mesh::BulkData> arg_bulk_data );

    void parallel_initialize_meta_data( const std::string & base_filename, const std::string & type = "exodusii" );
    void parallel_initialize_bulk_data();

    stk::mesh::fem::FEMMetaData & meta_data() {
      ThrowRequire( !Teuchos::is_null(m_fem_meta_data)) ;
      return *m_fem_meta_data;
    }

    stk::mesh::BulkData & bulk_data() {
      ThrowRequire( !Teuchos::is_null(m_bulk_data)) ;
      return *m_bulk_data;
    }

    coord_field_type & get_coordinate_field() {
      coord_field_type * coord_field = meta_data().get_field<coord_field_type>("coordinates");
      ThrowRequire( coord_field != NULL);
      return * coord_field;
    }

    Teuchos::RCP<Ioss::Region> input_ioss_region()  { return m_ioss_input_region; }
    Teuchos::RCP<Ioss::Region> output_ioss_region() { return m_ioss_output_region; }


  private:
    stk::ParallelMachine                       m_comm;
    Teuchos::RCP<stk::mesh::fem::FEMMetaData>  m_fem_meta_data;
    Teuchos::RCP<stk::mesh::BulkData>          m_bulk_data;

    Teuchos::RCP<Ioss::Region>                 m_ioss_input_region;
    Teuchos::RCP<Ioss::Region>                 m_ioss_output_region;

    std::string                                m_mesh_type;
    stk::io::util::MeshData                    m_mesh_data;

    //disallow copy constructor and assignment operator
    IO_Fixture( const IO_Fixture & );
    IO_Fixture & operator = ( const IO_Fixture & );
};

}//namespace util
}//namespace io
}//namespace stk

#endif //STK_IO_UTIL_IO_FIXTURE_HPP


