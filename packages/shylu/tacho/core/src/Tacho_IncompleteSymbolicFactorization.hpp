#ifndef __TACHO_INCOMPLETE_SYMBOLIC_FACTORIZATION_HPP__
#define __TACHO_INCOMPLETE_SYMBOLIC_FACTORIZATION_HPP__

/// \file symbolic_factor_helper.hpp
/// \brief The class compute a nonzero pattern with a given level of fills
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

  template<typename CrsMatBaseType>
  class IncompleteSymbolicFactorization {
  public:

    typedef typename CrsMatBaseType::space_type         space_type;

    typedef typename CrsMatBaseType::ordinal_type       ordinal_type;
    typedef typename CrsMatBaseType::size_type          size_type;

    typedef typename CrsMatBaseType::size_type_array    size_type_array;
    typedef typename CrsMatBaseType::ordinal_type_array ordinal_type_array;
    typedef typename CrsMatBaseType::value_type_array   value_type_array;

    typedef Kokkos::TeamPolicy<space_type,Kokkos::Schedule<Kokkos::Static> > team_policy_type;      
    typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;

    typedef typename space_type::scratch_memory_space shmem_space;
    typedef Kokkos::View<ordinal_type*,shmem_space,Kokkos::MemoryUnmanaged> shared_ordinal_type_array;

    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

  private:

    class Queue {
    private:
      shared_ordinal_type_array _q;
      ordinal_type _begin, _end;

    public:
      KOKKOS_INLINE_FUNCTION
      Queue(shared_ordinal_type_array q)
        : _q(q),_begin(0),_end(0) { }

      KOKKOS_INLINE_FUNCTION
      ~Queue() = default;

      KOKKOS_INLINE_FUNCTION
      ordinal_type size() const { return _end - _begin; }

      KOKKOS_INLINE_FUNCTION
      bool empty() const { return !size(); }

      KOKKOS_INLINE_FUNCTION
      void push(const ordinal_type val) { _q(_end++) = val; }

      KOKKOS_INLINE_FUNCTION
      ordinal_type pop() { return _q(_begin++); }

      KOKKOS_INLINE_FUNCTION
      ordinal_type end() { return _end; }

      KOKKOS_INLINE_FUNCTION
      void reset() { _begin = 0; _end = 0; }
    };

    class FunctorComputeNonZeroPatternInRow {
    private:
      const ordinal_type _level;
      const CrsMatBaseType _A;
      
      size_type_array _ap;
      ordinal_type_array _aj;

      const ordinal_type _rows_per_team;
      const ordinal_type _phase;

    public:

      FunctorComputeNonZeroPatternInRow(size_type_array ap,
                                        const ordinal_type level,
                                        const CrsMatBaseType A,
                                        const ordinal_type rows_per_team)
        : _level(level), _A(A), _ap(ap), _aj(), _rows_per_team(rows_per_team), _phase(0) {}
      
      FunctorComputeNonZeroPatternInRow(size_type_array ap,
                                        ordinal_type_array aj,
                                        const ordinal_type level,
                                        const CrsMatBaseType A,
                                        const ordinal_type rows_per_team)
        : _level(level), _A(A), _ap(ap), _aj(aj), _rows_per_team(rows_per_team), _phase(1) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const typename team_policy_type::member_type &member) const {
        const ordinal_type lrank = member.league_rank();
        const ordinal_type lsize = member.league_size();

        const ordinal_type m = _A.NumRows(), n = _A.NumCols();

        shared_ordinal_type_array queue   (member.team_shmem(), n);
        shared_ordinal_type_array distance(member.team_shmem(), n);
        shared_ordinal_type_array visited (member.team_shmem(), n);

        // clean up visited array
        for (ordinal_type i=0;i<n;++i)
          visited(i) = 0;

        Queue q(queue); // fixed size queue

        // ** this approach is not suitable for the symbolic factorization
        // const ordinal_type ibegin = lrank*_rows_per_team;
        // const ordinal_type itemp  = ibegin + _rows_per_team;
        // const ordinal_type iend   = itemp < m ? itemp : m;

        //for (ordinal_type i=ibegin;i<iend;++i) {

        // ** rows should be shuffled for better load balance
        for (ordinal_type i=lrank;i<m;i+=lsize) {
          size_type cnt = 0;

          // account for the diagonal
          switch (_phase) {
          case 0:
            cnt = 1;        // pre-count diagonal
            break;
          case 1:
            cnt = _ap(i);
            _aj(cnt++) = i; // put diagonal 
            break;
          }
          
          {
            // initialize work space
            q.push(i);
            distance(i) = 0;

            const ordinal_type id = (i+1);
            visited(i) = id;

            // breath first search for i
            while (!q.empty()) {
              const ordinal_type h = q.pop();

              // loop over j adjancy
              const ordinal_type jbegin = _A.RowPtrBegin(h), jend = _A.RowPtrEnd(h);
              for (ordinal_type j=jbegin;j<jend;++j) {
                const ordinal_type t = _A.Col(j);

                // skip diagonal and visited ones
                if (h != t && visited(t) != id) {
                  visited(t) = id;

                  if (t < i && (_level < 0 || distance(h) < _level)) {
                    q.push(t);
                    distance(t) = distance(h) + 1;
                  }
                  if (t > i) {
                    switch (_phase) {
                    case 0:
                      ++cnt;
                      break;
                    case 1:
                      _aj(cnt++) = t;
                      break;
                    }
                  }
                }
              }
            }

            // clear work space
            for (ordinal_type j=0;j<q.end();++j) 
              distance(queue(j)) = 0;
            q.reset();
          }
          switch (_phase) {
          case 0:
            _ap(i+1) = cnt;
            break;
          case 1:
            Util::sort( _aj, _ap(i), _ap(i+1) );
            break;
          }
        }
      }
    };

    class FunctorCountOffsetsInRow {
    public:
      typedef size_type value_type;

    private:
      size_type_array _off_in_rows;

    public:
      KOKKOS_INLINE_FUNCTION      
      FunctorCountOffsetsInRow(const size_type_array off_in_rows)
        : _off_in_rows(off_in_rows) { }
      
      KOKKOS_INLINE_FUNCTION
      ~FunctorCountOffsetsInRow() = default;
      
      KOKKOS_INLINE_FUNCTION
      void init(value_type &update) const {
        update = 0;
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type i, value_type &update, const bool final) const {
        update += _off_in_rows(i);
        if (final)
          _off_in_rows(i) = update;
      }
      
      KOKKOS_INLINE_FUNCTION
      void join(volatile       value_type &update,
                volatile const value_type &input) const {
        update += input;
      }
    };

  public:
    
    static void createNonZeroPattern(CrsMatBaseType &F,
                                     const ordinal_type level,
                                     const int uplo,
                                     const CrsMatBaseType A,
                                     const ordinal_type rows_per_team) {
      // static assert for space type 

      if (level) {

        // create crs matrix
        const ordinal_type m = A.NumRows(), n = A.NumCols();
        size_type_array ap = size_type_array("Symbolic::Temp::RowPtrArray", m+1);
        
        // later determined
        ordinal_type_array aj;
        value_type_array ax;
        size_type nnz  = 0;
        
        // league size is determined by workset size per team; use single team
        const size_type league_size = (m + rows_per_team - 1)/rows_per_team;
        const size_type team_size   = 1;
        
        team_policy_type team_policy(league_size, team_size);
        
        // scratch space for team
        const size_type scratch_size_per_team   = 3*shared_ordinal_type_array::shmem_size(n);
        const size_type scratch_size_per_thread = 0; 
        
        team_policy = team_policy.set_scratch_size( 1, // memory hierarchy level
                                                    Kokkos::PerTeam  (scratch_size_per_team),
                                                    Kokkos::PerThread(scratch_size_per_thread) );

        // parallel for count # of nonzeros 
        Kokkos::parallel_for(team_policy, FunctorComputeNonZeroPatternInRow(ap, level, A, rows_per_team));

        // parallel scan with range policy: too slow
        // range_policy_type range_policy(0, m+1);
        // Kokkos::parallel_scan(range_policy, FunctorCountOffsetsInRow(ap));

        for (size_type i=0;i<m;++i) 
          ap(i+1) += ap(i);

        nnz = ap(m);
        aj  = ordinal_type_array(std::string(F.Label()) + 
                                 std::string("::ColIndexArray"), nnz);
        ax  = value_type_array(std::string(F.Label()) + 
                               std::string("::Output::ValueArray"), nnz);
        
        // parallel for to fill column indices
        Kokkos::parallel_for(team_policy, FunctorComputeNonZeroPatternInRow(ap, aj, level, A, rows_per_team));

        size_type_array ap_begin = size_type_array(std::string(F.Label()) + 
                                                   std::string("::RowPtrArrayBegin"), m);
        size_type_array ap_end   = size_type_array(std::string(F.Label()) + 
                                                   std::string("::RowPtrArrayEnd"),   m);
        
        // create a matrix 
        Kokkos::deep_copy(ap_begin, Kokkos::subview(ap, range_type(0, m  )));
        Kokkos::deep_copy(ap_end,   Kokkos::subview(ap, range_type(1, m+1)));
        
        F = CrsMatBaseType(F.Label(), m, n, nnz, ap_begin, ap_end, aj, ax);
        F.setNumNonZeros();
        
        // add A
        CrsMatrixTools::add(F, A);
      } else {
        F.createConfTo(A);
        CrsMatrixTools::copy(F, uplo, 0, A);
      }
    }
    void createNonZeroPattern(const int uplo,
                              CrsMatBaseType &F) {
      createNonZeroPattern(-1, uplo, F);
    }

  };

}

#endif



