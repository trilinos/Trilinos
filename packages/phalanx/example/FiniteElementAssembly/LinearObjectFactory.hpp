#ifndef PHALANX_EXAMPLE_LINEAR_OBJECT_FACTORY_HPP
#define PHALANX_EXAMPLE_LINEAR_OBJECT_FACTORY_HPP

#include <utility>
#include "Mesh.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Kokkos_View.hpp"
#include "Kokkos_StaticCrsGraph.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

namespace phx_example {

  //! Functor for creating the node set 
  class LinearObjectFactory {

    // Phalanx device spaces
    using exec_t = PHX::exec_space;
    using mem_t = PHX::mem_space;
    using team_t =  Kokkos::TeamPolicy<exec_t>::member_type;
    
    const int num_nodes_;     //! Total number of nodes in mesh
    const int num_equations_; //! Number of equations per node
    const int num_dofs_;      //! Total number of DOFs in the problem

    // Global node ids <cell,node>
    Kokkos::View<int**,PHX::Device> gids_;
    
    //! Stores the global node dependencies for Jacobian creation (mesh nodes, not dofs).
    using SetType = Kokkos::UnorderedMap<std::pair<int,int>,void,PHX::Device>;
    SetType global_node_set_;

    //! Number of entries in each row of Jacobian, sizeof(num_dofs_)
    Kokkos::View<int*,PHX::Device> row_count_;

    //! Row offsets for building graph, sizeof(num_dofs_+1)
    Kokkos::View<size_t*,PHX::Device> row_offsets_;

    //! Single value for teh size of the matrix
    Kokkos::View<size_t,PHX::Device> matrix_size_;
    
    //! Global ids for each column
    Kokkos::View<int*, PHX::Device> col_ids_;

    //! Crs Graph for Jacobian
    Kokkos::StaticCrsGraph<int,Kokkos::LayoutLeft,PHX::Device> graph_;

    public:
    
    struct TagFillNodeSet{};
    struct TagComputeRowOffsets{};
    struct TagFillGraphEntries{};
    struct TagSortGraphEntries{};

    LinearObjectFactory(const int& num_nodes,
                        const int& num_equations,
                        const Kokkos::View<int**,PHX::Device>& global_node_ids) :
      num_nodes_(num_nodes),
      num_equations_(num_equations),
      num_dofs_(num_nodes * num_equations),
      gids_(global_node_ids),
      row_count_("row_count",num_dofs_),
      row_offsets_("row_offsets",num_dofs_+1),
      matrix_size_("matrix_size"),
      graph_()
    {
      // Use Unordered_Map to get node/node dependencies for graph
      {
        size_t set_capacity = (28ull * num_nodes_);
        unsigned failed_insert_count = 0;
        
        do {
          // Zero the row count to restart the fill
          Kokkos::deep_copy(row_count_, 0);
          global_node_set_ = Kokkos::UnorderedMap<std::pair<int,int>,void,PHX::Device>( ( set_capacity += failed_insert_count ) );
          
          // May be larger than requested:
          set_capacity = global_node_set_.capacity();
          
          Kokkos::parallel_reduce(Kokkos::RangePolicy<exec_t,TagFillNodeSet>(0,gids_.extent(0)),
                                  *this, failed_insert_count);
          exec_t::fence();
          
        } while ( failed_insert_count ); 
      }

      // Build the offsets using exclusive scan of row_count_ into row_offsets_.
      // This also handles the final total in the 'row_count + 1' position.
      Kokkos::parallel_scan(Kokkos::RangePolicy<exec_t,TagComputeRowOffsets>(0,num_dofs_),*this);      
      exec_t::fence();

      // Fill the graph values
      size_t num_matrix_entries = 0;
      Kokkos::deep_copy(num_matrix_entries,matrix_size_);
      Kokkos::deep_copy(row_count_,0); // reuse this view
      Kokkos::fence();
      graph_.row_map = row_offsets_;
      graph_.entries = Kokkos::StaticCrsGraph<int,Kokkos::LayoutLeft,PHX::Device>::entries_type("graph_entries",num_matrix_entries);
      Kokkos::parallel_for(Kokkos::RangePolicy<exec_t,TagFillGraphEntries>(0,global_node_set_.capacity()), *this );
      exec_t::fence();
      
      // Free memory for temporaries
      matrix_size_ = Kokkos::View<size_t,PHX::Device>();
      row_count_ = Kokkos::View<int*,PHX::Device>();
      row_offsets_ = Kokkos::View<size_t*,PHX::Device>();
      global_node_set_.clear();
      
      // Sort the graph values in each row
      Kokkos::parallel_for(Kokkos::RangePolicy<exec_t,TagSortGraphEntries>(0,num_dofs_), *this );      
    }
    
    Kokkos::View<double*,PHX::Device> createSolutionVector(const std::string& name) const
    { return Kokkos::View<double*,PHX::Device>(name,num_dofs_); }
    
    KokkosSparse::CrsMatrix<double,int,PHX::Device> createJacobianMatrix(const std::string& name) const
    { return KokkosSparse::CrsMatrix<double,int,PHX::Device>(name,graph_); }


    // ****************************************************
    // TagFillNodeSet
    // ****************************************************
    
    KOKKOS_INLINE_FUNCTION
    void operator()(TagFillNodeSet, int cell , unsigned & failed_count) const
    {
      const int num_nodes = static_cast<int>(gids_.extent(1));
      
      for (int row=0; row < num_nodes; ++row) {
        const int row_gid = gids_(cell,row);
        
        for (int col=0; col < num_nodes; ++col) {
          const int col_gid = gids_(cell,col);

          const auto result = global_node_set_.insert( std::make_pair(row_gid,col_gid) );

          // Successful only if the key did not already exist
          if (result.success()) {
            // Add all dofs for the node/node pairing
            for (int row_eq_offset=0; row_eq_offset < num_equations_; ++row_eq_offset) {
              Kokkos::atomic_fetch_add(&row_count_( row_gid * num_equations_ + row_eq_offset ), num_equations_);
            }
          }
          else if (result.failed()) {
            ++failed_count ;
          }
        }
      }    
    }
    
    KOKKOS_INLINE_FUNCTION
    void init(const TagFillNodeSet& , unsigned& update) const
    { update = 0 ; }
    
    KOKKOS_INLINE_FUNCTION
    void join(const TagFillNodeSet& , volatile unsigned& update, volatile const unsigned& input ) const
    { update += input ; }  
    
    // ****************************************************
    // TagComputeRowOffsets: Parallel Scan for row_offsets_
    // ****************************************************
    KOKKOS_INLINE_FUNCTION
    void operator()(const TagComputeRowOffsets,
                    const unsigned irow,
                    unsigned& update,
                    const bool final ) const
    {
      // exclusive scan
      if (final)
        row_offsets_(irow) = update ;
      
      update += row_count_(irow);
      
      if (final) {
        if ((irow + 1) == row_count_.extent(0)) {
          row_offsets_(irow + 1) = update;
          matrix_size_() = update; 
        }
      }
    }
    
    KOKKOS_INLINE_FUNCTION
    void init(const TagComputeRowOffsets&, unsigned & update) const
    { update = 0; }
    
    KOKKOS_INLINE_FUNCTION
    void join(const TagComputeRowOffsets&,
              volatile unsigned& update,
              volatile const unsigned & input) const
    { update += input ; }

    // ****************************************************
    // TagFillGraphEntries
    // ****************************************************

    KOKKOS_INLINE_FUNCTION
    void operator()(const TagFillGraphEntries&, const int iset) const
    {
      if ( global_node_set_.valid_at(iset) ) {
        // Add each entry to the graph entries.
        
        const std::pair<int,int> key = global_node_set_.key_at(iset) ;
        const int row_node = key.first ;
        const int col_node = key.second ;
        
        // loop over equations
        for (int row_eq=0; row_eq < num_equations_; ++row_eq) {
          const int row_gid = row_node * num_equations_ + row_eq;  
          for (int col_eq=0; col_eq < num_equations_; ++col_eq) {
            const int col_gid = col_node * num_equations_ + col_eq;
            const size_t matrix_offset = graph_.row_map(row_gid) + Kokkos::atomic_fetch_add(&row_count_(row_gid), 1);
            graph_.entries( matrix_offset ) = col_gid;
          }
        }
        
      }
    }

    // ****************************************************
    // TagSortGraphEntries
    // ****************************************************

    KOKKOS_INLINE_FUNCTION
    void operator()(const TagSortGraphEntries&, const int irow) const
    {
      const size_t row_beg = graph_.row_map( irow );
      const size_t row_end = graph_.row_map( irow + 1 );
      for (size_t i = row_beg + 1; i < row_end; ++i ) {
        const int col = graph_.entries(i);
        size_t j = i ;
        for ( ; row_beg < j && col < graph_.entries(j-1) ; --j ) {
          graph_.entries(j) = graph_.entries(j-1);
        }
        graph_.entries(j) = col ;
      }
    }
    
    
  };

}

#endif
