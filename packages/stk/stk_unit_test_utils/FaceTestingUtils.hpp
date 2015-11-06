
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

inline bool check_face_elem_connectivity(stk::mesh::BulkData& mesh, const std::vector<int>& countsIn, bool debug=false) {
    std::vector<int> counts(6,-1);
    std::copy(countsIn.begin(),countsIn.end(),counts.begin());
    stk::mesh::FieldBase const * coord = mesh.mesh_meta_data().coordinate_field();
    stk::mesh::BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);
    bool extra_face_not_accounted_for = false;
    for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
        for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
            stk::mesh::Entity face = bucket[face_count];
            bool all_x_equal_half = true;
            stk::mesh::Entity const * node = mesh.begin_nodes(face);
            for (unsigned node_count = 0; node_count < mesh.num_nodes(face); ++node_count) {
                double *xyz = static_cast<double *>(stk::mesh::field_data(*coord, node[node_count]));
                if (xyz[0] != 0.5) {
                    all_x_equal_half = false;
                    break;
                }
            }
            if (all_x_equal_half) {
                if (counts[0] == static_cast<int>(mesh.num_elements(face))) {
                    counts[0] = -1;
                }
                else if (counts[1] == static_cast<int>(mesh.num_elements(face))) {
                    counts[1] = -1;
                }
                else if (counts[2] == static_cast<int>(mesh.num_elements(face))) {
                    counts[2] = -1;
                }
                else if (counts[3] == static_cast<int>(mesh.num_elements(face))) {
                    counts[3] = -1;
                }
                else if (counts[4] == static_cast<int>(mesh.num_elements(face))) {
                    counts[4] = -1;
                }
                else if (counts[5] == static_cast<int>(mesh.num_elements(face))) {
                    counts[5] = -1;
                }
                else {
                    extra_face_not_accounted_for = true;
                }
                if (debug) {
                    std::cout << "num_elements:" << mesh.num_elements(face) << std::endl;
                    stk::mesh::Entity const * elements = mesh.begin_elements(face);
                    for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                        std::cout << "elem_count " << elem_count << " id " << mesh.entity_key(elements[elem_count])  << " for face_count " << face_count << " for bucket count " << bucket_count << std::endl;
                    }
                    std::cout.flush();
                }
            }
        }
    }
    if (!debug && counts[5] == -1 && counts[4] == -1 &&  counts[3] == -1 && counts[2] == -1 && counts[1] == -1 && counts[0] == -1 && !extra_face_not_accounted_for) {
        return true;
    }
    else {
        if (debug) {
            return false;
        }
        else {
            return check_face_elem_connectivity(mesh, counts, true);
        }
    }
}

inline bool read_file_check_face_elem_connectivity_stk(std::string filename, bool create_faces, const std::vector<int>& counts) {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return check_face_elem_connectivity(mesh, counts);

}


inline bool check_face_elem_connectivity(stk::mesh::BulkData& mesh, int count1=-1, int count2=-1, int count3=-1, int count4=-1, int count5=-1, int count6=-1, bool debug=false) {
    stk::mesh::FieldBase const * coord = mesh.mesh_meta_data().coordinate_field();
    stk::mesh::BucketVector const & face_buckets = mesh.buckets(stk::topology::FACE_RANK);
    bool extra_face_not_accounted_for = false;
    for (size_t bucket_count=0, bucket_end=face_buckets.size(); bucket_count < bucket_end; ++bucket_count) {
        stk::mesh::Bucket & bucket = *face_buckets[bucket_count];
        for (size_t face_count=0, face_end=bucket.size(); face_count < face_end; ++face_count) {
            stk::mesh::Entity face = bucket[face_count];
            bool all_x_equal_half = true;
            stk::mesh::Entity const * node = mesh.begin_nodes(face);
            for (unsigned node_count = 0; node_count < mesh.num_nodes(face); ++node_count) {
                double *xyz = static_cast<double *>(stk::mesh::field_data(*coord, node[node_count]));
                if (xyz[0] != 0.5) {
                    all_x_equal_half = false;
                    break;
                }
            }
            if (all_x_equal_half) {
                if (count1 == static_cast<int>(mesh.num_elements(face))) {
                    count1 = -1;
                }
                else if (count2 == static_cast<int>(mesh.num_elements(face))) {
                    count2 = -1;
                }
                else if (count3 == static_cast<int>(mesh.num_elements(face))) {
                    count3 = -1;
                }
                else if (count4 == static_cast<int>(mesh.num_elements(face))) {
                    count4 = -1;
                }
                else if (count5 == static_cast<int>(mesh.num_elements(face))) {
                    count5 = -1;
                }
                else if (count6 == static_cast<int>(mesh.num_elements(face))) {
                    count6 = -1;
                }
                else {
                    extra_face_not_accounted_for = true;
                }
                if (debug) {
                    std::cout << "num_elements:" << mesh.num_elements(face) << std::endl;
                    stk::mesh::Entity const * elements = mesh.begin_elements(face);
                    for (unsigned elem_count = 0; elem_count < mesh.num_elements(face); ++elem_count) {
                        std::cout << "elem_count " << elem_count << " id " << mesh.entity_key(elements[elem_count])  << " for face_count " << face_count << " for bucket count " << bucket_count << std::endl;
                    }
                    std::cout.flush();
                }
            }
        }
    }
    if (!debug && count6 == -1 && count5 == -1 &&  count4 == -1 && count3 == -1 && count2 == -1 && count1 == -1 && !extra_face_not_accounted_for) {
        return true;
    }
    else {
        if (debug) {
            return false;
        }
        else {
            return check_face_elem_connectivity(mesh, count1, count2, count3, count4, count5, count6, true);
        }
    }
}

inline bool read_file_check_face_elem_connectivity_stk(std::string filename, bool create_faces, int count1=-1, int count2=-1, int count3=-1, int count4=-1, int count5=-1, int count6=-1) {
    stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
    stkMeshIoBroker.set_sideset_face_creation_behavior(stk::io::StkMeshIoBroker::STK_IO_SIDESET_FACE_CREATION_CURRENT);
    stkMeshIoBroker.add_mesh_database(filename, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();
    stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();
    if (create_faces) {
        stk::mesh::create_faces(mesh);
    }
    return check_face_elem_connectivity(mesh, count1, count2, count3, count4, count5, count6);

}


#endif
