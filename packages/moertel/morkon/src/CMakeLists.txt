TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(INCLUDES
    mrk_api_classes.hpp
    mrk_compute_normals.hpp
    mrk_data_types.hpp
	mrk_interface_host_side_adapter.hpp
    mrk_interface_impl.hpp
	mrk_int_types.hpp
	mrk_manager_impl.hpp
    mrk_search_for_pallet_generating_faces.hpp
    mrk_compute_pallets_from_candidate_face_pairs.hpp
    mrk_planar_geom.hpp
   )

SET(CPPSOURCES
   )


#TRIBITS_ADD_LIBRARY(morkon HEADERS ${INCLUDES})

#SET_PROPERTY(TARGET morkon_exp_exp_lib APPEND PROPERTY LINKER_LANGUAGE CXX)

