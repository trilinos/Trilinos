// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_POINTCLOUDSEARCH_HPP_
#define _COMPADRE_POINTCLOUDSEARCH_HPP_

#include "Compadre_Typedefs.hpp"
#include "Compadre_NeighborLists.hpp"
#include "nanoflann.hpp"
#include <Kokkos_Core.hpp>
#include <memory>

namespace Compadre {

//class sort_indices
//{
//   private:
//     double* mparr;
//   public:
//     sort_indices(double* parr) : mparr(parr) {}
//     bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
//};

//! Custom RadiusResultSet for nanoflann that uses pre-allocated space for indices and radii
//! instead of using std::vec for std::pairs
template <typename _DistanceType, typename _IndexType = size_t>
class RadiusResultSet {

  public:

    typedef _DistanceType DistanceType;
    typedef _IndexType IndexType;

    const DistanceType radius;
    IndexType count;
    DistanceType* r_dist;
    IndexType* i_dist;
    const IndexType max_size;

    RadiusResultSet(
        DistanceType radius_,
        DistanceType* r_dist_, IndexType* i_dist_, const IndexType max_size_)
        : radius(radius_), count(0), r_dist(r_dist_), i_dist(i_dist_), max_size(max_size_) {
        init();
    }

    void init() {}

    void clear() { count = 0; }

    size_t size() const { return count; }

    bool full() const { return true; }

    bool addPoint(DistanceType dist, IndexType index) {

        if (dist < radius) {
            // would throw an exception here if count>=max_size, but this code is 
            // called inside of a parallel region so only an abort is possible, 
            // but this situation is recoverable 
            //
            // instead, we increase count, but there is nowhere to store neighbors
            // since we are working with pre-allocated space
            // this will be discovered after returning from the search by comparing
            // the count against the pre-allocate space dimensions
            if (count<max_size) {
                i_dist[count] = index;
                r_dist[count] = dist;
            }
            count++;
        }
        return true;

    }

    DistanceType worstDist() const { return radius; }

    std::pair<IndexType, DistanceType> worst_item() const {
        // just to verify this isn't ever called
        compadre_kernel_assert_release(false && "worst_item() should not be called.");
    }

    void sort() { 
        // puts closest neighbor as the first entry in the neighbor list
        // leaves the rest unsorted
 
        if (count > 0) {

            // alternate sort for every entry, not currently being used
            //int indices[count];
            //for (int i=0; i<count; ++i) {
            //    indices[i] = i;
            //}
            //std::sort(indices, indices+count, sort_indices(r_dist));
            //std::vector<IndexType> tmp_indices(count);
            //std::vector<DistanceType> tmp_r(count);
            //for (int i=0; i<count; ++i) {
            //    tmp_indices[i] = i_dist[indices[i]];
            //    tmp_r[i] = r_dist[indices[i]];
            //}
            //for (int i=0; i<count; ++i) {
            //    i_dist[i] = tmp_indices[i];
            //    r_dist[i] = tmp_r[i];
            //}
            IndexType loop_max = (count < max_size) ? count : max_size;

            DistanceType best_distance = std::numeric_limits<DistanceType>::max();
            IndexType best_distance_index = 0;
            int best_index = -1;
            for (IndexType i=0; i<loop_max; ++i) {
                if (r_dist[i] < best_distance) {
                    best_distance = r_dist[i];
                    best_distance_index = i_dist[i];
                    best_index = i;
                }
            }

            if (best_index != 0) {
                auto tmp_ind = i_dist[0];
                i_dist[0] = best_distance_index;
                i_dist[best_index] = tmp_ind;
            }
        }
    }
};


//!  PointCloudSearch generates neighbor lists and window sizes for each target site
/*!
*  Search methods can be run in dry-run mode, or not.
*
*  #### When in dry-run mode:
*                                                                                                                         
*    `neighbors_list` will be populated with number of neighbors found for each target site.
*
*    This allows a user to know memory allocation needed before storage of neighbor indices.
*                                                                                                                         
*  #### When not in dry-run mode:
*
*    `neighbors_list_offsets` will be populated with offsets for values (dependent on method) determined by neighbor_list.
*    If a 2D view for `neighbors_list` is used, then \f$ N(i,j+1) \f$ will store the \f$ j^{th} \f$ neighbor of \f$ i \f$,
*    and \f$ N(i,0) \f$ will store the number of neighbors for target \f$ i \f$.
*
*/
template <typename view_type>
class PointCloudSearch {

    public:

        typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloudSearch<view_type> >, 
                PointCloudSearch<view_type>, 1> tree_type_1d;
        typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloudSearch<view_type> >, 
                PointCloudSearch<view_type>, 2> tree_type_2d;
        typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloudSearch<view_type> >, 
                PointCloudSearch<view_type>, 3> tree_type_3d;

    protected:

        //! source site coordinates
        view_type _src_pts_view;
        local_index_type _dim;
        local_index_type _max_leaf;

        std::shared_ptr<tree_type_1d> _tree_1d;
        std::shared_ptr<tree_type_2d> _tree_2d;
        std::shared_ptr<tree_type_3d> _tree_3d;

    public:

        PointCloudSearch(view_type src_pts_view, const local_index_type dimension = -1,
                const local_index_type max_leaf = -1) 
                : _src_pts_view(src_pts_view), 
                  _dim((dimension < 0) ? src_pts_view.extent(1) : dimension),
                  _max_leaf((max_leaf < 0) ? 10 : max_leaf) {
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename view_type::memory_space>::accessible==1)
                    && "Views passed to PointCloudSearch at construction should be accessible from the host.");
        };
    
        ~PointCloudSearch() {};

        //! Returns a liberal estimated upper bound on number of neighbors to be returned by a neighbor search
        //! for a given choice of dimension, basis size, and epsilon_multiplier. Assumes quasiuniform distribution
        //! of points. This result can be used to size a preallocated neighbor_lists kokkos view.
        static inline int getEstimatedNumberNeighborsUpperBound(int unisolvency_size, const int dimension, const double epsilon_multiplier) {
            int multiplier = 1;
            if (dimension==1) multiplier = 2;
            return multiplier * 2.0 * unisolvency_size * std::pow(epsilon_multiplier, dimension) + 1; // +1 is for the number of neighbors entry needed in the first entry of each row
        }
    
        //! Bounding box query method required by Nanoflann.
        template <class BBOX> bool kdtree_get_bbox(BBOX& bb) const {return false;}

        //! Returns the number of source sites
        inline int kdtree_get_point_count() const {return _src_pts_view.extent(0);}

        //! Returns the coordinate value of a point
        inline double kdtree_get_pt(const int idx, int dim) const {return _src_pts_view(idx,dim);}

        //! Returns the distance between a point and a source site, given its index
        inline double kdtree_distance(const double* queryPt, const int idx, long long sz) const {

            double distance = 0;
            for (int i=0; i<_dim; ++i) {
                distance += (_src_pts_view(idx,i)-queryPt[i])*(_src_pts_view(idx,i)-queryPt[i]);
            }
            return std::sqrt(distance);

        }

        void generateKDTree() {
            if (_dim==1) {
                _tree_1d = std::make_shared<tree_type_1d>(1, *this, nanoflann::KDTreeSingleIndexAdaptorParams(_max_leaf));
                _tree_1d->buildIndex();
            } else if (_dim==2) {
                _tree_2d = std::make_shared<tree_type_2d>(2, *this, nanoflann::KDTreeSingleIndexAdaptorParams(_max_leaf));
                _tree_2d->buildIndex();
            } else if (_dim==3) {
                _tree_3d = std::make_shared<tree_type_3d>(3, *this, nanoflann::KDTreeSingleIndexAdaptorParams(_max_leaf));
                _tree_3d->buildIndex();
            }
        }

        /*! \brief Generates neighbor lists of 2D view by performing a radius search 
            where the radius to be searched is in the epsilons view.
            If uniform_radius is given, then this overrides the epsilons view radii sizes.
            Accepts 2D neighbor_lists without number_of_neighbors_list.
            \param is_dry_run               [in] - whether to do a dry-run (find neighbors, but don't store)
            \param trg_pts_view             [in] - target coordinates from which to seek neighbors
            \param neighbor_lists           [out] - 2D view of neighbor lists to be populated from search
            \param epsilons                 [in/out] - radius to search, overwritten if uniform_radius != 0
            \param uniform_radius           [in] - double != 0 determines whether to overwrite all epsilons for uniform search
            \param max_search_radius        [in] - largest valid search (useful only for MPI jobs if halo size exists)
        */
        template <typename trg_view_type, typename neighbor_lists_view_type, typename epsilons_view_type>
        size_t generate2DNeighborListsFromRadiusSearch(bool is_dry_run, trg_view_type trg_pts_view, 
                neighbor_lists_view_type neighbor_lists, epsilons_view_type epsilons, 
                const double uniform_radius = 0.0, double max_search_radius = 0.0) {

            // function does not populate epsilons, they must be prepopulated

            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename trg_view_type::memory_space>::accessible==1) &&
                    "Target coordinates view passed to generate2DNeighborListsFromRadiusSearch should be accessible from the host.");
            compadre_assert_release((((int)trg_pts_view.extent(1))>=_dim) &&
                    "Target coordinates view passed to generate2DNeighborListsFromRadiusSearch must have \
                    second dimension as large as _dim.");
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename neighbor_lists_view_type::memory_space>::accessible==1) &&
                    "Views passed to generate2DNeighborListsFromRadiusSearch should be accessible from the host.");
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename epsilons_view_type::memory_space>::accessible==1) &&
                    "Views passed to generate2DNeighborListsFromRadiusSearch should be accessible from the host.");

            // loop size
            const int num_target_sites = trg_pts_view.extent(0);

            if ((!_tree_1d && _dim==1) || (!_tree_2d && _dim==2) || (!_tree_3d && _dim==3)) {
                this->generateKDTree();
            }

            // check neighbor lists and epsilons view sizes
            compadre_assert_release((neighbor_lists.extent(0)==(size_t)num_target_sites 
                        && neighbor_lists.extent(1)>=1)
                        && "neighbor lists View does not have large enough dimensions");
            compadre_assert_release((neighbor_lists_view_type::rank==2) && "neighbor_lists must be a 2D Kokkos view.");

            compadre_assert_release((epsilons.extent(0)==(size_t)num_target_sites)
                        && "epsilons View does not have the correct dimension");

            typedef Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
                    scratch_double_view;

            typedef Kokkos::View<size_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
                    scratch_int_view;

            // determine scratch space size needed
            int team_scratch_size = 0;
            team_scratch_size += scratch_double_view::shmem_size(neighbor_lists.extent(1)); // distances
            team_scratch_size += scratch_int_view::shmem_size(neighbor_lists.extent(1)); // indices
            team_scratch_size += scratch_double_view::shmem_size(_dim); // target coordinate

            // maximum number of neighbors found over all target sites' neighborhoods
            size_t max_num_neighbors = 0;
            // part 2. do radius search using window size from knn search
            // each row of neighbor lists is a neighbor list for the target site corresponding to that row
            Kokkos::parallel_reduce("radius search", host_team_policy(num_target_sites, Kokkos::AUTO)
                    .set_scratch_size(0 /*shared memory level*/, Kokkos::PerTeam(team_scratch_size)), 
                    KOKKOS_LAMBDA(const host_member_type& teamMember, size_t& t_max_num_neighbors) {

                // make unmanaged scratch views
                scratch_double_view neighbor_distances(teamMember.team_scratch(0 /*shared memory*/), neighbor_lists.extent(1));
                scratch_int_view neighbor_indices(teamMember.team_scratch(0 /*shared memory*/), neighbor_lists.extent(1));
                scratch_double_view this_target_coord(teamMember.team_scratch(0 /*shared memory*/), _dim);

                size_t neighbors_found = 0;

                const int i = teamMember.league_rank();

                // set epsilons if radius is specified
                if (uniform_radius > 0) epsilons(i) = uniform_radius;

                // needs furthest neighbor's distance for next portion
                compadre_kernel_assert_release((epsilons(i)<=max_search_radius || max_search_radius==0) && "max_search_radius given (generally derived from the size of a halo region), and search radius needed would exceed this max_search_radius.");

                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, neighbor_lists.extent(1)), [&](const int j) { 
                    neighbor_indices(j) = 0;
                    neighbor_distances(j) = -1.0;
                });
                teamMember.team_barrier();

                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    // target_coords is LayoutLeft on device and its HostMirror, so giving a pointer to 
                    // this data would lead to a wrong result if the device is a GPU

                    for (int j=0; j<_dim; ++j) {
                        this_target_coord(j) = trg_pts_view(i,j);
                    }

                    nanoflann::SearchParams sp; // default parameters
                    Compadre::RadiusResultSet<double> rrs(epsilons(i)*epsilons(i), neighbor_distances.data(), neighbor_indices.data(), neighbor_lists.extent(1));
                    if (_dim==1) {
                        neighbors_found = _tree_1d->template radiusSearchCustomCallback<Compadre::RadiusResultSet<double> >(this_target_coord.data(), rrs, sp) ;
                    } else if (_dim==2) {
                        neighbors_found = _tree_2d->template radiusSearchCustomCallback<Compadre::RadiusResultSet<double> >(this_target_coord.data(), rrs, sp) ;
                    } else if (_dim==3) {
                        neighbors_found = _tree_3d->template radiusSearchCustomCallback<Compadre::RadiusResultSet<double> >(this_target_coord.data(), rrs, sp) ;
                    }

                    t_max_num_neighbors = (neighbors_found > t_max_num_neighbors) ? neighbors_found : t_max_num_neighbors;
            
                    // the number of neighbors is stored in column zero of the neighbor lists 2D array
                    neighbor_lists(i,0) = neighbors_found;

                    // epsilons already scaled and then set by search radius
                });
                teamMember.team_barrier();

                // loop_bound so that we don't write into memory we don't have allocated
                int loop_bound = (neighbors_found < neighbor_lists.extent(1)-1) ? neighbors_found : neighbor_lists.extent(1)-1;
                // loop over each neighbor index and fill with a value
                if (!is_dry_run) {
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, loop_bound), [&](const int j) {
                        // cast to an whatever data type the 2D array of neighbor lists is using
                        neighbor_lists(i,j+1) = static_cast<typename std::remove_pointer<typename std::remove_pointer<typename neighbor_lists_view_type::data_type>::type>::type>(neighbor_indices(j));
                    });
                    teamMember.team_barrier();
                }

            }, Kokkos::Max<size_t>(max_num_neighbors) );
            Kokkos::fence();

            // check if max_num_neighbors will fit onto pre-allocated space
            compadre_assert_release((neighbor_lists.extent(1) >= (max_num_neighbors+1) || is_dry_run) 
                    && "neighbor_lists does not contain enough columns for the maximum number of neighbors needing to be stored.");

            return max_num_neighbors;
        }

        /*! \brief Generates compressed row neighbor lists by performing a radius search 
            where the radius to be searched is in the epsilons view.
            If uniform_radius is given, then this overrides the epsilons view radii sizes.
            Accepts 1D neighbor_lists with 1D number_of_neighbors_list.
            \param is_dry_run               [in] - whether to do a dry-run (find neighbors, but don't store)
            \param trg_pts_view             [in] - target coordinates from which to seek neighbors
            \param neighbor_lists           [out] - 1D view of neighbor lists to be populated from search
            \param number_of_neighbors_list [in/out] - number of neighbors for each target site
            \param epsilons                 [in/out] - radius to search, overwritten if uniform_radius != 0
            \param uniform_radius           [in] - double != 0 determines whether to overwrite all epsilons for uniform search
            \param max_search_radius        [in] - largest valid search (useful only for MPI jobs if halo size exists)
        */
        template <typename trg_view_type, typename neighbor_lists_view_type, typename epsilons_view_type>
        size_t generateCRNeighborListsFromRadiusSearch(bool is_dry_run, trg_view_type trg_pts_view, 
                neighbor_lists_view_type neighbor_lists, neighbor_lists_view_type number_of_neighbors_list, 
                epsilons_view_type epsilons, const double uniform_radius = 0.0, double max_search_radius = 0.0) {

            // function does not populate epsilons, they must be prepopulated

            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename trg_view_type::memory_space>::accessible==1) &&
                    "Target coordinates view passed to generateCRNeighborListsFromRadiusSearch should be accessible from the host.");
            compadre_assert_release((((int)trg_pts_view.extent(1))>=_dim) &&
                    "Target coordinates view passed to generateCRNeighborListsFromRadiusSearch must have \
                    second dimension as large as _dim.");
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename neighbor_lists_view_type::memory_space>::accessible==1) &&
                    "Views passed to generateCRNeighborListsFromRadiusSearch should be accessible from the host.");
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename epsilons_view_type::memory_space>::accessible==1) &&
                    "Views passed to generateCRNeighborListsFromRadiusSearch should be accessible from the host.");

            // loop size
            const int num_target_sites = trg_pts_view.extent(0);

            if ((!_tree_1d && _dim==1) || (!_tree_2d && _dim==2) || (!_tree_3d && _dim==3)) {
                this->generateKDTree();
            }

            compadre_assert_release((number_of_neighbors_list.extent(0)==(size_t)num_target_sites)
                        && "number_of_neighbors_list or neighbor lists View does not have large enough dimensions");
            compadre_assert_release((neighbor_lists_view_type::rank==1) && "neighbor_lists must be a 1D Kokkos view.");

            typedef Kokkos::View<global_index_type*, typename neighbor_lists_view_type::array_layout,
                    typename neighbor_lists_view_type::memory_space, typename neighbor_lists_view_type::memory_traits> row_offsets_view_type;
            row_offsets_view_type row_offsets;
            int max_neighbor_list_row_storage_size = 1;
            if (!is_dry_run) {
                auto nla = CreateNeighborLists(neighbor_lists, number_of_neighbors_list);
                max_neighbor_list_row_storage_size = nla.getMaxNumNeighbors();
                Kokkos::resize(row_offsets, num_target_sites);
                Kokkos::fence();
                Kokkos::parallel_for(Kokkos::RangePolicy<host_execution_space>(0,num_target_sites), [&](const int i) {
                    row_offsets(i) = nla.getRowOffsetHost(i); 
                });
                Kokkos::fence();
            }

            compadre_assert_release((epsilons.extent(0)==(size_t)num_target_sites)
                        && "epsilons View does not have the correct dimension");

            typedef Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
                    scratch_double_view;

            typedef Kokkos::View<size_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
                    scratch_int_view;

            // determine scratch space size needed
            int team_scratch_size = 0;
            team_scratch_size += scratch_double_view::shmem_size(max_neighbor_list_row_storage_size); // distances
            team_scratch_size += scratch_int_view::shmem_size(max_neighbor_list_row_storage_size); // indices
            team_scratch_size += scratch_double_view::shmem_size(_dim); // target coordinate

            // part 2. do radius search using window size from knn search
            // each row of neighbor lists is a neighbor list for the target site corresponding to that row
            Kokkos::parallel_for("radius search", host_team_policy(num_target_sites, Kokkos::AUTO)
                    .set_scratch_size(0 /*shared memory level*/, Kokkos::PerTeam(team_scratch_size)), 
                    KOKKOS_LAMBDA(const host_member_type& teamMember) {

                // make unmanaged scratch views
                scratch_double_view neighbor_distances(teamMember.team_scratch(0 /*shared memory*/), max_neighbor_list_row_storage_size);
                scratch_int_view neighbor_indices(teamMember.team_scratch(0 /*shared memory*/), max_neighbor_list_row_storage_size);
                scratch_double_view this_target_coord(teamMember.team_scratch(0 /*shared memory*/), _dim);

                size_t neighbors_found = 0;

                const int i = teamMember.league_rank();

                // set epsilons if radius is specified
                if (uniform_radius > 0) epsilons(i) = uniform_radius;

                // needs furthest neighbor's distance for next portion
                compadre_kernel_assert_release((epsilons(i)<=max_search_radius || max_search_radius==0) && "max_search_radius given (generally derived from the size of a halo region), and search radius needed would exceed this max_search_radius.");

                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, max_neighbor_list_row_storage_size), [&](const int j) { 
                    neighbor_indices(j) = 0;
                    neighbor_distances(j) = -1.0;
                });
                teamMember.team_barrier();

                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    // target_coords is LayoutLeft on device and its HostMirror, so giving a pointer to 
                    // this data would lead to a wrong result if the device is a GPU

                    for (int j=0; j<_dim; ++j) {
                        this_target_coord(j) = trg_pts_view(i,j);
                    }

                    nanoflann::SearchParams sp; // default parameters
                    Compadre::RadiusResultSet<double> rrs(epsilons(i)*epsilons(i), neighbor_distances.data(), neighbor_indices.data(), max_neighbor_list_row_storage_size);
                    if (_dim==1) {
                        neighbors_found = _tree_1d->template radiusSearchCustomCallback<Compadre::RadiusResultSet<double> >(this_target_coord.data(), rrs, sp) ;
                    } else if (_dim==2) {
                        neighbors_found = _tree_2d->template radiusSearchCustomCallback<Compadre::RadiusResultSet<double> >(this_target_coord.data(), rrs, sp) ;
                    } else if (_dim==3) {
                        neighbors_found = _tree_3d->template radiusSearchCustomCallback<Compadre::RadiusResultSet<double> >(this_target_coord.data(), rrs, sp) ;
                    }
            
                    // we check that neighbors found doesn't differ from dry-run or we store neighbors_found
                    // no check that neighbors found stay the same if uniform_radius specified (!=0)
                    if (is_dry_run || uniform_radius!=0.0) {
                        number_of_neighbors_list(i) = neighbors_found;
                    } else {
                        compadre_kernel_assert_debug((neighbors_found==(size_t)number_of_neighbors_list(i)) 
                                && "Number of neighbors found changed since dry-run.");
                    }

                    // epsilons already scaled and then set by search radius
                });
                teamMember.team_barrier();

                // loop_bound so that we don't write into memory we don't have allocated
                int loop_bound = neighbors_found;
                // loop over each neighbor index and fill with a value
                if (!is_dry_run) {
                    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, loop_bound), [&](const int j) {
                        // cast to an whatever data type the 2D array of neighbor lists is using
                        neighbor_lists(row_offsets(i)+j) = static_cast<typename std::remove_pointer<typename std::remove_pointer<typename neighbor_lists_view_type::data_type>::type>::type>(neighbor_indices(j));
                    });
                    teamMember.team_barrier();
                }

            });
            Kokkos::fence();
            auto nla = CreateNeighborLists(number_of_neighbors_list);
            return nla.getTotalNeighborsOverAllListsHost();
        }


        /*! \brief Generates neighbor lists as 2D view by performing a k-nearest neighbor search
            Only accepts 2D neighbor_lists without number_of_neighbors_list.
            \param is_dry_run               [in] - whether to do a dry-run (find neighbors, but don't store)
            \param trg_pts_view             [in] - target coordinates from which to seek neighbors
            \param neighbor_lists           [out] - 2D view of neighbor lists to be populated from search
            \param epsilons                 [in/out] - radius to search, overwritten if uniform_radius != 0
            \param neighbors_needed         [in] - k neighbors needed as a minimum
            \param epsilon_multiplier       [in] - distance to kth neighbor multiplied by epsilon_multiplier for follow-on radius search
            \param max_search_radius        [in] - largest valid search (useful only for MPI jobs if halo size exists)
        */
        template <typename trg_view_type, typename neighbor_lists_view_type, typename epsilons_view_type>
        size_t generate2DNeighborListsFromKNNSearch(bool is_dry_run, trg_view_type trg_pts_view, 
                neighbor_lists_view_type neighbor_lists, epsilons_view_type epsilons, 
                const int neighbors_needed, const double epsilon_multiplier = 1.6, 
                double max_search_radius = 0.0) {

            // First, do a knn search (removes need for guessing initial search radius)

            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename trg_view_type::memory_space>::accessible==1) &&
                    "Target coordinates view passed to generate2DNeighborListsFromKNNSearch should be accessible from the host.");
            compadre_assert_release((((int)trg_pts_view.extent(1))>=_dim) &&
                    "Target coordinates view passed to generate2DNeighborListsFromRadiusSearch must have \
                    second dimension as large as _dim.");
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename neighbor_lists_view_type::memory_space>::accessible==1) &&
                    "Views passed to generate2DNeighborListsFromKNNSearch should be accessible from the host.");
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename epsilons_view_type::memory_space>::accessible==1) &&
                    "Views passed to generate2DNeighborListsFromKNNSearch should be accessible from the host.");

            // loop size
            const int num_target_sites = trg_pts_view.extent(0);

            if ((!_tree_1d && _dim==1) || (!_tree_2d && _dim==2) || (!_tree_3d && _dim==3)) {
                this->generateKDTree();
            }
            Kokkos::fence();

            compadre_assert_release((num_target_sites==0 || // sizes don't matter when there are no targets
                    (neighbor_lists.extent(0)==(size_t)num_target_sites 
                        && neighbor_lists.extent(1)>=(size_t)(neighbors_needed+1)))
                        && "neighbor lists View does not have large enough dimensions");
            compadre_assert_release((neighbor_lists_view_type::rank==2) && "neighbor_lists must be a 2D Kokkos view.");

            compadre_assert_release((epsilons.extent(0)==(size_t)num_target_sites)
                        && "epsilons View does not have the correct dimension");

            typedef Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
                    scratch_double_view;

            typedef Kokkos::View<size_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
                    scratch_int_view;

            // determine scratch space size needed
            int team_scratch_size = 0;
            team_scratch_size += scratch_double_view::shmem_size(neighbor_lists.extent(1)); // distances
            team_scratch_size += scratch_int_view::shmem_size(neighbor_lists.extent(1)); // indices
            team_scratch_size += scratch_double_view::shmem_size(_dim); // target coordinate

            // minimum number of neighbors found over all target sites' neighborhoods
            size_t min_num_neighbors = 0;
            //
            // part 1. do knn search for neighbors needed for unisolvency
            // each row of neighbor lists is a neighbor list for the target site corresponding to that row
            //
            // as long as neighbor_lists can hold the number of neighbors_needed, we don't need to check
            // that the maximum number of neighbors will fit into neighbor_lists
            //
            Kokkos::parallel_reduce("knn search", host_team_policy(num_target_sites, Kokkos::AUTO)
                    .set_scratch_size(0 /*shared memory level*/, Kokkos::PerTeam(team_scratch_size)), 
                    KOKKOS_LAMBDA(const host_member_type& teamMember, size_t& t_min_num_neighbors) {

                // make unmanaged scratch views
                scratch_double_view neighbor_distances(teamMember.team_scratch(0 /*shared memory*/), neighbor_lists.extent(1));
                scratch_int_view neighbor_indices(teamMember.team_scratch(0 /*shared memory*/), neighbor_lists.extent(1));
                scratch_double_view this_target_coord(teamMember.team_scratch(0 /*shared memory*/), _dim);

                size_t neighbors_found = 0;

                const int i = teamMember.league_rank();

                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, neighbor_lists.extent(1)), [=](const int j) {
                    neighbor_indices(j) = 0;
                    neighbor_distances(j) = -1.0;
                });
            
                teamMember.team_barrier();
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    // target_coords is LayoutLeft on device and its HostMirror, so giving a pointer to 
                    // this data would lead to a wrong result if the device is a GPU

                    for (int j=0; j<_dim; ++j) { 
                        this_target_coord(j) = trg_pts_view(i,j);
                    }

                    if (_dim==1) {
                        neighbors_found = _tree_1d->knnSearch(this_target_coord.data(), neighbors_needed, 
                                neighbor_indices.data(), neighbor_distances.data()) ;
                    } else if (_dim==2) {
                        neighbors_found = _tree_2d->knnSearch(this_target_coord.data(), neighbors_needed, 
                                neighbor_indices.data(), neighbor_distances.data()) ;
                    } else if (_dim==3) {
                        neighbors_found = _tree_3d->knnSearch(this_target_coord.data(), neighbors_needed, 
                                neighbor_indices.data(), neighbor_distances.data()) ;
                    }

                    // get minimum number of neighbors found over all target sites' neighborhoods
                    t_min_num_neighbors = (neighbors_found < t_min_num_neighbors) ? neighbors_found : t_min_num_neighbors;
            
                    // scale by epsilon_multiplier to window from location where the last neighbor was found
                    epsilons(i) = (neighbor_distances(neighbors_found-1) > 0) ?
                        std::sqrt(neighbor_distances(neighbors_found-1))*epsilon_multiplier : 1e-14*epsilon_multiplier;
                    // the only time the second case using 1e-14 is used is when either zero neighbors or exactly one 
                    // neighbor (neighbor is target site) is found.  when the follow on radius search is conducted, the one
                    // neighbor (target site) will not be found if left at 0, so any positive amount will do, however 1e-14 
                    // should is small enough to ensure that other neighbors are not found

                    // needs furthest neighbor's distance for next portion
                    compadre_kernel_assert_release((neighbors_found<neighbor_lists.extent(1) || is_dry_run) 
                            && "Neighbors found exceeded storage capacity in neighbor list.");
                    compadre_kernel_assert_release((epsilons(i)<=max_search_radius || max_search_radius==0 || is_dry_run) 
                            && "max_search_radius given (generally derived from the size of a halo region), \
                                and search radius needed would exceed this max_search_radius.");
                    // neighbor_distances stores squared distances from neighbor to target, as returned by nanoflann
                });
            }, Kokkos::Min<size_t>(min_num_neighbors) );
            Kokkos::fence();
            
            // if no target sites, then min_num_neighbors is set to neighbors_needed
            // which also avoids min_num_neighbors being improperly set by min reduction
            if (num_target_sites==0) min_num_neighbors = neighbors_needed;

            // Next, check that we found the neighbors_needed number that we require for unisolvency
            compadre_assert_release((num_target_sites==0 || (min_num_neighbors>=(size_t)neighbors_needed))
                    && "Neighbor search failed to find number of neighbors needed for unisolvency.");
            
            // call a radius search using values now stored in epsilons
            size_t max_num_neighbors = generate2DNeighborListsFromRadiusSearch(is_dry_run, trg_pts_view, neighbor_lists, 
                    epsilons, 0.0 /*don't set uniform radius*/, max_search_radius);

            return max_num_neighbors;
        }

        /*! \brief Generates compressed row neighbor lists by performing a k-nearest neighbor search
            Only accepts 1D neighbor_lists with 1D number_of_neighbors_list.
            \param is_dry_run               [in] - whether to do a dry-run (find neighbors, but don't store)
            \param trg_pts_view             [in] - target coordinates from which to seek neighbors
            \param neighbor_lists           [out] - 1D view of neighbor lists to be populated from search
            \param number_of_neighbors_list [in/out] - number of neighbors for each target site
            \param epsilons                 [in/out] - radius to search, overwritten if uniform_radius != 0
            \param neighbors_needed         [in] - k neighbors needed as a minimum
            \param epsilon_multiplier       [in] - distance to kth neighbor multiplied by epsilon_multiplier for follow-on radius search
            \param max_search_radius        [in] - largest valid search (useful only for MPI jobs if halo size exists)
        */
        template <typename trg_view_type, typename neighbor_lists_view_type, typename epsilons_view_type>
        size_t generateCRNeighborListsFromKNNSearch(bool is_dry_run, trg_view_type trg_pts_view, 
                neighbor_lists_view_type neighbor_lists, neighbor_lists_view_type number_of_neighbors_list,
                epsilons_view_type epsilons, const int neighbors_needed, const double epsilon_multiplier = 1.6, 
                double max_search_radius = 0.0) {

            // First, do a knn search (removes need for guessing initial search radius)

            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename trg_view_type::memory_space>::accessible==1) &&
                    "Target coordinates view passed to generateCRNeighborListsFromKNNSearch should be accessible from the host.");
            compadre_assert_release((((int)trg_pts_view.extent(1))>=_dim) &&
                    "Target coordinates view passed to generateCRNeighborListsFromRadiusSearch must have \
                    second dimension as large as _dim.");
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename neighbor_lists_view_type::memory_space>::accessible==1) &&
                    "Views passed to generateCRNeighborListsFromKNNSearch should be accessible from the host.");
            compadre_assert_release((Kokkos::SpaceAccessibility<host_execution_space, typename epsilons_view_type::memory_space>::accessible==1) &&
                    "Views passed to generateCRNeighborListsFromKNNSearch should be accessible from the host.");

            // loop size
            const int num_target_sites = trg_pts_view.extent(0);

            if ((!_tree_1d && _dim==1) || (!_tree_2d && _dim==2) || (!_tree_3d && _dim==3)) {
                this->generateKDTree();
            }
            Kokkos::fence();

            compadre_assert_release((number_of_neighbors_list.extent(0)==(size_t)num_target_sites ) 
                        && "number_of_neighbors_list or neighbor lists View does not have large enough dimensions");
            compadre_assert_release((neighbor_lists_view_type::rank==1) && "neighbor_lists must be a 1D Kokkos view.");

            // if dry-run, neighbors_needed, else max over previous dry-run
            int max_neighbor_list_row_storage_size = neighbors_needed;
            if (!is_dry_run) {
                auto nla = CreateNeighborLists(neighbor_lists, number_of_neighbors_list);
                max_neighbor_list_row_storage_size = nla.getMaxNumNeighbors();
            }

            compadre_assert_release((epsilons.extent(0)==(size_t)num_target_sites)
                        && "epsilons View does not have the correct dimension");

            typedef Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
                    scratch_double_view;

            typedef Kokkos::View<size_t*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > 
                    scratch_int_view;

            // determine scratch space size needed
            int team_scratch_size = 0;
            team_scratch_size += scratch_double_view::shmem_size(max_neighbor_list_row_storage_size); // distances
            team_scratch_size += scratch_int_view::shmem_size(max_neighbor_list_row_storage_size); // indices
            team_scratch_size += scratch_double_view::shmem_size(_dim); // target coordinate

            // minimum number of neighbors found over all target sites' neighborhoods
            size_t min_num_neighbors = 0;
            //
            // part 1. do knn search for neighbors needed for unisolvency
            // each row of neighbor lists is a neighbor list for the target site corresponding to that row
            //
            // as long as neighbor_lists can hold the number of neighbors_needed, we don't need to check
            // that the maximum number of neighbors will fit into neighbor_lists
            //
            Kokkos::parallel_reduce("knn search", host_team_policy(num_target_sites, Kokkos::AUTO)
                    .set_scratch_size(0 /*shared memory level*/, Kokkos::PerTeam(team_scratch_size)), 
                    KOKKOS_LAMBDA(const host_member_type& teamMember, size_t& t_min_num_neighbors) {

                // make unmanaged scratch views
                scratch_double_view neighbor_distances(teamMember.team_scratch(0 /*shared memory*/), max_neighbor_list_row_storage_size);
                scratch_int_view neighbor_indices(teamMember.team_scratch(0 /*shared memory*/), max_neighbor_list_row_storage_size);
                scratch_double_view this_target_coord(teamMember.team_scratch(0 /*shared memory*/), _dim);

                size_t neighbors_found = 0;

                const int i = teamMember.league_rank();

                Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, max_neighbor_list_row_storage_size), [=](const int j) {
                    neighbor_indices(j) = 0;
                    neighbor_distances(j) = -1.0;
                });
            
                teamMember.team_barrier();
                Kokkos::single(Kokkos::PerTeam(teamMember), [&] () {
                    // target_coords is LayoutLeft on device and its HostMirror, so giving a pointer to 
                    // this data would lead to a wrong result if the device is a GPU

                    for (int j=0; j<_dim; ++j) { 
                        this_target_coord(j) = trg_pts_view(i,j);
                    }

                    if (_dim==1) {
                        neighbors_found = _tree_1d->knnSearch(this_target_coord.data(), neighbors_needed, 
                                neighbor_indices.data(), neighbor_distances.data()) ;
                    } else if (_dim==2) {
                        neighbors_found = _tree_2d->knnSearch(this_target_coord.data(), neighbors_needed, 
                                neighbor_indices.data(), neighbor_distances.data()) ;
                    } else if (_dim==3) {
                        neighbors_found = _tree_3d->knnSearch(this_target_coord.data(), neighbors_needed, 
                                neighbor_indices.data(), neighbor_distances.data()) ;
                    }

                    // get minimum number of neighbors found over all target sites' neighborhoods
                    t_min_num_neighbors = (neighbors_found < t_min_num_neighbors) ? neighbors_found : t_min_num_neighbors;
            
                    // scale by epsilon_multiplier to window from location where the last neighbor was found
                    epsilons(i) = (neighbor_distances(neighbors_found-1) > 0) ?
                        std::sqrt(neighbor_distances(neighbors_found-1))*epsilon_multiplier : 1e-14*epsilon_multiplier;
                    // the only time the second case using 1e-14 is used is when either zero neighbors or exactly one 
                    // neighbor (neighbor is target site) is found.  when the follow on radius search is conducted, the one
                    // neighbor (target site) will not be found if left at 0, so any positive amount will do, however 1e-14 
                    // should is small enough to ensure that other neighbors are not found

                    compadre_kernel_assert_release((epsilons(i)<=max_search_radius || max_search_radius==0 || is_dry_run) 
                            && "max_search_radius given (generally derived from the size of a halo region), \
                                and search radius needed would exceed this max_search_radius.");
                    // neighbor_distances stores squared distances from neighbor to target, as returned by nanoflann
                });
            }, Kokkos::Min<size_t>(min_num_neighbors) );
            Kokkos::fence();
            
            // if no target sites, then min_num_neighbors is set to neighbors_needed
            // which also avoids min_num_neighbors being improperly set by min reduction
            if (num_target_sites==0) min_num_neighbors = neighbors_needed;

            // Next, check that we found the neighbors_needed number that we require for unisolvency
            compadre_assert_release((num_target_sites==0 || (min_num_neighbors>=(size_t)neighbors_needed))
                    && "Neighbor search failed to find number of neighbors needed for unisolvency.");
            
            // call a radius search using values now stored in epsilons
            generateCRNeighborListsFromRadiusSearch(is_dry_run, trg_pts_view, neighbor_lists, 
                    number_of_neighbors_list, epsilons, 0.0 /*don't set uniform radius*/, max_search_radius);

            auto nla = CreateNeighborLists(number_of_neighbors_list);
            return nla.getTotalNeighborsOverAllListsHost();
        }
}; // PointCloudSearch

//! CreatePointCloudSearch allows for the construction of an object of type PointCloudSearch with template deduction
template <typename view_type>
PointCloudSearch<view_type> CreatePointCloudSearch(view_type src_view, const local_index_type dimensions = -1, const local_index_type max_leaf = -1) { 
    return PointCloudSearch<view_type>(src_view, dimensions, max_leaf);
}

} // Compadre

#endif
