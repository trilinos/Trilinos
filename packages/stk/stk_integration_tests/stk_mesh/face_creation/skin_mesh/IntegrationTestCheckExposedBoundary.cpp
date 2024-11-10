/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <Ionit_Initializer.h>     // for Initializer
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_io/DatabasePurpose.hpp"
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_unit_test_utils/FaceCreationTestUtils.hpp>

namespace
{

class SkinnedMeshWithModifiedSkinPart : public stk::unit_test_util::MeshTestFixture
{
protected:
  virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    stk::unit_test_util::read_from_serial_file_and_decompose("ARA.e", get_bulk(), "cyclic");
    stk::mesh::Part& skinnedPart = SideTestUtil::run_skin_mesh(get_bulk(), get_things_to_skin(get_bulk()));

    run_modification(skinnedPart);
    EXPECT_FALSE(stk::mesh::check_exposed_block_boundary_sides(get_bulk(), get_things_to_skin(get_bulk()), skinnedPart));
  }

  virtual void run_modification(stk::mesh::Part &skin) = 0;

  stk::mesh::EntityVector get_faces(stk::mesh::Selector selector)
  {
    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(selector, get_bulk().buckets(stk::topology::FACE_RANK), faces);
    return faces;
  }

  stk::mesh::Selector get_things_to_skin(const stk::mesh::BulkData& bulkData)
  {
    return bulkData.mesh_meta_data().universal_part();
  }
};

class SkinnedMeshWithExtraFace: public SkinnedMeshWithModifiedSkinPart
{
protected:
  virtual void run_modification(stk::mesh::Part &skin)
  {
    get_bulk().modification_begin();
    add_extra_face_to_skin(skin);
    get_bulk().modification_end();
  }

private:
  void add_extra_face_to_skin(stk::mesh::Part &skin)
  {
    if(get_bulk().parallel_rank() == 0) {
      stk::mesh::EntityVector notSkinFaces = get_faces(!skin);
      ASSERT_EQ(1u, notSkinFaces.size());
      get_bulk().change_entity_parts(notSkinFaces[0], stk::mesh::ConstPartVector{&skin});
    }
  }
};

class SkinnedMeshWithMissingFace: public SkinnedMeshWithModifiedSkinPart
{
protected:
  virtual void run_modification(stk::mesh::Part &skin)
  {
    get_bulk().modification_begin();
    remove_face_from_skin(skin);
    get_bulk().modification_end();
  }

private:
  void remove_face_from_skin(stk::mesh::Part &skin)
  {
    stk::mesh::EntityVector skinFaces = get_faces(skin);
    ASSERT_EQ(5u, skinFaces.size());
    get_bulk().change_entity_parts(skinFaces[0], stk::mesh::ConstPartVector{}, stk::mesh::ConstPartVector{&skin});
  }
};

TEST_F(SkinnedMeshWithExtraFace, skin_no_aura)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(SkinnedMeshWithMissingFace, skin_no_aura)
{
  run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}

}
