#pragma once
#ifndef __SYMBOLIC_FACTOR_HELPER_HPP__
#define __SYMBOLIC_FACTOR_HELPER_HPP__

/// \file symbolic_factor_helper.hpp
/// \brief The class compute a nonzero pattern with a given level of fills
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "util.hpp"

namespace Example {

  using namespace std;

  template<class CrsMatrixType>
  class SymbolicFactorHelper : public Disp {
  public:
    typedef typename CrsMatrixType::ordinal_type ordinal_type;
    typedef typename CrsMatrixType::size_type    size_type;
    typedef typename CrsMatrixType::space_type   space_type;

    typedef typename CrsMatrixType::ordinal_type_array ordinal_type_array;
    typedef typename CrsMatrixType::size_type_array    size_type_array;
    typedef typename CrsMatrixType::value_type_array   value_type_array;

  private:
    string _label;                   // name of this class

    // matrix index base
    CrsMatrixType _A;                // input matrix
    ordinal_type _m, _n;             // matrix dimension

    struct crs_graph {
      ordinal_type_array _cidx;      // col index array
      size_type _nnz;                // # of nonzeros
      size_type_array _rptr;         // row ptr array
    };
    typedef struct crs_graph crs_graph_type;
    crs_graph_type _in, _out;

    typedef Kokkos::View<ordinal_type**, Kokkos::LayoutLeft, space_type> league_specific_ordinal_type_array;
    typedef typename league_specific_ordinal_type_array::value_type* league_specific_ordinal_type_array_ptr;

    int _lsize;
    league_specific_ordinal_type_array _queue, _visited, _distance;

    void createInternalWorkSpace() {
      _queue    = league_specific_ordinal_type_array(_label+"::QueueArray",    _m, _lsize);
      _visited  = league_specific_ordinal_type_array(_label+"::VisitedArray",  _m, _lsize);
      _distance = league_specific_ordinal_type_array(_label+"::DistanceArray", _m, _lsize);
    }

    void freeInternalWorkSpace() {
      _queue    = league_specific_ordinal_type_array();
      _visited  = league_specific_ordinal_type_array();
      _distance = league_specific_ordinal_type_array();
    }

  public:

    void setLabel(string label) { _label = label; }
    string Label() const { return _label; }

    SymbolicFactorHelper(const CrsMatrixType &A,
                         const int lsize = (space_type::thread_pool_size(0)/
                                            space_type::thread_pool_size(2)))  {

      _label = "SymbolicFactorHelper::" + A.Label();

      // matrix index base and the number of rows
      _A = A;

      _m = _A.NumRows();
      _n = _A.NumCols();

      // allocate memory for input crs matrix
      _in._nnz   = _A.NumNonZeros();
      _in._rptr  = size_type_array(_label+"::Input::RowPtrArray", _m+1);
      _in._cidx  = ordinal_type_array(_label+"::Input::ColIndexArray", _in._nnz);

      // adjust graph structure; A is assumed to have a graph without its diagonal
      A.convertGraph(_in._nnz, _in._rptr, _in._cidx);

      // league size
      _lsize = lsize;
    }
    virtual~SymbolicFactorHelper() { }

    class Queue {
    private:
      league_specific_ordinal_type_array_ptr _q;
      ordinal_type _begin, _end;

    public:
      Queue(league_specific_ordinal_type_array_ptr q)
        : _q(q),_begin(0),_end(0) { }

      ordinal_type size() const { return _end - _begin; }
      bool empty() const { return !size(); }

      void push(const ordinal_type val) { _q[_end++] = val; }
      ordinal_type pop() { return _q[_begin++]; }
      ordinal_type end() { return _end; }
    };

    class FunctorComputeNonZeroPatternInRow {
    public:
      typedef Kokkos::TeamPolicy<space_type> policy_type;

    private:
      ordinal_type _level, _m;
      crs_graph_type _graph;

      league_specific_ordinal_type_array _queue;
      league_specific_ordinal_type_array _visited;
      league_specific_ordinal_type_array _distance;

      size_type_array _rptr;
      ordinal_type_array _cidx;

      ordinal_type _phase;

    public:
      FunctorComputeNonZeroPatternInRow(const ordinal_type level,
                                        const ordinal_type m,
                                        const crs_graph_type &graph,
                                        league_specific_ordinal_type_array &queue,
                                        league_specific_ordinal_type_array &visited,
                                        league_specific_ordinal_type_array &distance,
                                        size_type_array &rptr,
                                        ordinal_type_array &cidx)
        : _level(level), _m(m), _graph(graph),
          _queue(queue), _visited(visited), _distance(distance),
          _rptr(rptr), _cidx(cidx), _phase(0)
      { }

      void setPhaseCountNumNonZeros() { _phase = 0; }
      void setPhaseComputeColIndex()  { _phase = 1; }

      KOKKOS_INLINE_FUNCTION
      void operator()(const typename policy_type::member_type &member) const {
        const int lrank = member.league_rank();
        const int lsize = member.league_size();

        // shuffle rows to get better load balance;
        // for instance, if ND is applied, more fills are generated in the last seperator.
        for (ordinal_type i=lrank;i<_m;i+=lsize) {
          league_specific_ordinal_type_array_ptr queue    = &_queue(0, lrank);
          league_specific_ordinal_type_array_ptr distance = &_distance(0, lrank);
          league_specific_ordinal_type_array_ptr visited  = &_visited(0, lrank);

          size_type cnt = 0;

          // account for the diagonal
          switch (_phase) {
          case 0:
            cnt = 1;
            break;
          case 1:
            cnt = _rptr[i];
            _cidx[cnt++] = i;
            break;
          }

          {
            Queue q(queue); // fixed size queue

            // initialize work space
            q.push(i);
            distance[i] = 0;

            const ordinal_type id = (i+1);
            visited[i] = id;

            // breath first search for i
            while (!q.empty()) {
              const ordinal_type h = q.pop();
              // loop over j adjancy
              const ordinal_type jbegin = _graph._rptr[h], jend = _graph._rptr[h+1];
              for (ordinal_type j=jbegin;j<jend;++j) {
                const ordinal_type t = _graph._cidx[j];
                if (visited[t] != id) {
                  visited[t] = id;

                  if (t < i && (_level < 0 || distance[h] < _level)) {
                    q.push(t);
                    distance[t] = distance[h] + 1;
                  }
                  if (t > i) {
                    switch (_phase) {
                    case 0:
                      ++cnt;
                      break;
                    case 1:
                      _cidx[cnt++] = t;
                      break;
                    }
                  }
                }
              }
            }

            // clear work space
            for (ordinal_type j=0;j<q.end();++j) {
              const ordinal_type jj = queue[j];
              distance[jj] = 0;
              visited[jj] = 0;
            }
          }
          switch (_phase) {
          case 0:
            _rptr[i+1] = cnt;
            break;
          case 1:
            sort(&_cidx[_rptr[i]], &_cidx[_rptr[i+1]]);
            break;
          }
        }
      }
    };

    class FunctorCountOffsetsInRow {
    public:
      typedef Kokkos::RangePolicy<space_type> policy_type;
      typedef size_type value_type;

    private:
      size_type_array _off_in_rows;

    public:
      FunctorCountOffsetsInRow(size_type_array &off_in_rows)
        : _off_in_rows(off_in_rows)
      { }

      KOKKOS_INLINE_FUNCTION
      void init(value_type &update) const {
        update = 0;
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const typename policy_type::member_type &i, value_type &update, const bool final) const {
        update += _off_in_rows(i);
        if (final)
          _off_in_rows(i) = update;
      }

      KOKKOS_INLINE_FUNCTION
      void join(volatile value_type &update,
                volatile const value_type &input) const {
        update += input;
      }
    };

    int createNonZeroPattern(const ordinal_type level,
                             const int uplo,
                             CrsMatrixType &F) {
      _out._nnz  = 0;
      _out._rptr = size_type_array(_label+"::Output::RowPtrArray", _m+1);

      createInternalWorkSpace();
      {
        FunctorComputeNonZeroPatternInRow functor(level, _m, _in,
                                                  _queue,
                                                  _visited,
                                                  _distance,
                                                  _out._rptr,
                                                  _out._cidx);

        functor.setPhaseCountNumNonZeros();
        Kokkos::parallel_for(typename FunctorComputeNonZeroPatternInRow::policy_type(_lsize, 1), functor);
      }
      {
        FunctorCountOffsetsInRow functor(_out._rptr);
        Kokkos::parallel_scan(typename FunctorCountOffsetsInRow::policy_type(0, _m+1), functor);

        _out._nnz  = _out._rptr[_m];
        _out._cidx = ordinal_type_array(_label+"::Output::ColIndexArray", _out._nnz);
      }
      {
        FunctorComputeNonZeroPatternInRow functor(level, _m, _in,
                                                  _queue,
                                                  _visited,
                                                  _distance,
                                                  _out._rptr,
                                                  _out._cidx);

        functor.setPhaseComputeColIndex();
        Kokkos::parallel_for(typename FunctorComputeNonZeroPatternInRow::policy_type(_lsize, 1), functor);
      }
      freeInternalWorkSpace();

      {
        value_type_array ax = value_type_array(_label+"::Output::ValueArray", _out._nnz);
        F = CrsMatrixType(_label+"::Output::Matrix", _m, _n, _out._nnz, _out._rptr, _out._cidx, ax);
        F.add(_A);
      }

      return 0;
    }

    ostream& showMe(ostream &os) const {
      streamsize prec = os.precision();
      os.precision(15);
      os << scientific;

      const int w = 6;

      os << " -- Matrix Dimension -- " << endl
         << "    # of Rows  = " << _m << endl
         << "    # of Cols  = " << _n << endl;

      os << endl;

      os << " -- Input Graph Without Diagonals -- " << endl
         << "    # of NonZeros  = " << _in._nnz << endl ;

      os << " -- Input Graph :: RowPtr -- " << endl;
      for (ordinal_type i=0;i<_in._rptr.dimension_0();++i)
        os << setw(w) << i
           << setw(w) << _in._rptr[i]
           << endl;

      os << endl;

      os << " -- Output Graph With Diagonals-- " << endl
         << "    # of NonZeros  = " << _out._nnz << endl ;

      os << " -- Output Graph :: RowPtr -- " << endl;
      for (ordinal_type i=0;i<_out._rptr.dimension_0();++i)
        os << setw(w) << i
           << setw(w) << _out._rptr[i]
           << endl;

      os.unsetf(ios::scientific);
      os.precision(prec);

      return os;
    }

  };

}

#endif



