
#ifndef FACETESTINGUTILS_HPP_
#define FACETESTINGUTILS_HPP_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>

inline unsigned read_file_count_sides(std::string filename, bool create_faces)
{
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);

    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    std::vector<unsigned> countVec;
    stk::mesh::count_entities(mesh.mesh_meta_data().universal_part(), mesh, countVec);
    return countVec[2];
}

inline bool fully_connected_elements_to_faces(stk::mesh::BulkData& mesh)
{
    bool fully_connected = true;
    stk::mesh::BucketVector const & elem_buckets = mesh.buckets(stk::topology::ELEMENT_RANK);
    for (size_t bucket_count=0, bucket_end=elem_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *elem_buckets[bucket_count];
        const unsigned num_expected_faces = bucket.topology().num_faces();
        for (size_t elem_count=0, elem_end=bucket.size(); elem_count < elem_end; ++elem_count) {
            stk::mesh::Entity elem = bucket[elem_count];
            if (num_expected_faces != mesh.num_faces(elem)) {
                fully_connected = false;
                break;
            }
        }
    }
    return fully_connected;
}

inline unsigned read_file_fully_connected_stk(std::string filename, bool create_faces) {
    stk::mesh::MetaData meta;
    stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
    stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);

    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return fully_connected_elements_to_faces(mesh);
}

inline unsigned count_shared_faces_between_different_elements(stk::mesh::BulkData& mesh) {
    unsigned shared_face_count = 0;
    stk::mesh::BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);
    for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
        for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
            stk::mesh::Entity face = bucket[face_count];
            bool is_face_shared = false;
            stk::mesh::Entity const * elements = mesh.begin_elements(face);
            for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                for (unsigned other_elem_count = elem_count;
                        other_elem_count < mesh.num_elements(face); ++other_elem_count) {
                    if ((elem_count != other_elem_count) &&
                            (elements[elem_count] != elements[other_elem_count])) {
                        is_face_shared = true;
                        break;
                    }
                }
            }
            if (is_face_shared) {
                ++shared_face_count;
            }
        }
    }
    return shared_face_count;
}

inline unsigned read_file_shared_faces_different_elements_stk(std::string filename, bool create_faces) {
    stk::mesh::MetaData meta;
        stk::mesh::BulkData mesh(meta, MPI_COMM_WORLD);
        stk::unit_test_util::fill_mesh_using_stk_io(filename, mesh, MPI_COMM_WORLD);

    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return count_shared_faces_between_different_elements(mesh);
}

inline unsigned count_shared_faces_between_same_element(stk::mesh::BulkData& mesh) {
    unsigned shared_face_count = 0;
    stk::mesh::BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);
    for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
        for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
            stk::mesh::Entity face = bucket[face_count];
            bool is_face_shared = false;
            stk::mesh::Entity const * elements = mesh.begin_elements(face);
            for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                for (unsigned other_elem_count = elem_count;
                        other_elem_count < mesh.num_elements(face); ++other_elem_count) {
                    if ((elem_count != other_elem_count) &&
                            (elements[elem_count] == elements[other_elem_count])) {
                        is_face_shared = true;
                        break;
                    }
                }
            }
            if (is_face_shared) {
                ++shared_face_count;
            }
        }
    }
    return shared_face_count;
}

inline unsigned read_file_shared_faces_same_elements_stk(std::string filename, bool create_faces) {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return count_shared_faces_between_same_element(mesh);

}



#endif
