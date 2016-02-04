#pragma once
#ifndef __GRAPH_HELPER_CAMD_HPP__
#define __GRAPH_HELPER_CAMD_HPP__

/// \file graph_helper_camd.hpp
/// \brief Interface to camd (suire sparse) ordering
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "camd.h"
#include "util.hpp"

namespace Tacho {

  using namespace std;

  template<class CrsMatBaseType>
  class GraphHelper_CAMD : public Disp {
  public:
    typedef typename CrsMatBaseType::ordinal_type ordinal_type;
    typedef typename CrsMatBaseType::size_type    size_type;

    typedef typename CrsMatBaseType::ordinal_type_array ordinal_type_array;
    typedef typename CrsMatBaseType::size_type_array    size_type_array;

  private:
    string _label;

    ordinal_type _m;
    size_type _nnz;
    size_type_array _rptr;
    ordinal_type_array _cidx, _perm, _peri, _cnst;

    double _control[CAMD_CONTROL], _info[CAMD_INFO];

    bool _is_ordered;

  public:

    void setLabel(string label) { _label = label; }
    string Label() const { return _label; }

    size_type NumNonZeros() const { return _nnz; }
    ordinal_type NumRows() const { return _m; }

    ordinal_type_array PermVector()       const { return _perm; }
    ordinal_type_array InvPermVector()    const { return _peri; }

    ordinal_type_array ConstraintVector() const { return _cnst; }

    GraphHelper_CAMD() = default;

    GraphHelper_CAMD(const CrsMatBaseType A,
                     const ordinal_type nblk,
                     const ordinal_type_array range) {
      _label = "GraphHelper_CAMD::" + A.Label();
      
      _m     = A.NumRows();
      _nnz   = A.NumNonZeros();

      // input matrix A should be strictly lower or upper wihtout diagonals
      // create a upper graph structure without diagonals
      _rptr  = size_type_array(_label+"::RowPtrArray", _m+1);
      _cidx  = ordinal_type_array(_label+"::ColIndexArray", _nnz);
      _cnst  = ordinal_type_array(_label+"::ConstraintArray", _m+1);
      _perm  = ordinal_type_array(_label+"::PermutationArray", _m);
      _peri  = ordinal_type_array(_label+"::InvPermutationArray", _m);

      auto rptr = A.RowPtr();
      for (ordinal_type i=0;i<(_m+1);++i)
        _rptr[i] = rptr[i];
      
      auto cidx = A.ColPtr();
      for (size_type i=0;i<_nnz;++i)
        _cidx[i] = cidx[i];

      for (ordinal_type i=0;i<nblk;++i)
        for (ordinal_type j=range[i];j<range[i+1];++j)
          _cnst[j] = i;
    }
    GraphHelper_CAMD(const GraphHelper_CAMD &b) = default;

    int computeOrdering() {
      int r_val = 0;

      camd_defaults(_control);
      camd_control(_control);

      ordinal_type *rptr = reinterpret_cast<ordinal_type*>(_rptr.ptr_on_device());
      ordinal_type *cidx = reinterpret_cast<ordinal_type*>(_cidx.ptr_on_device());
      ordinal_type *cnst = reinterpret_cast<ordinal_type*>(_cnst.ptr_on_device());
      ordinal_type *perm = reinterpret_cast<ordinal_type*>(_perm.ptr_on_device());
      
      r_val = camd_order(_m, rptr, cidx, perm, _control, _info, cnst) ;

      for (ordinal_type i=0;i<_m;++i)
        _peri[perm[i]] = i;

      _is_ordered = true;
      
      return r_val;
    }

    ostream& showMe(ostream &os) const {
      streamsize prec = os.precision();
      os.precision(15);
      os << scientific;

      os << " -- CAMD input -- " << endl
         << "    # of Rows      = " << _m << endl
         << "    # of NonZeros  = " << _nnz << endl;

      if (_is_ordered)
        os << " -- Ordering -- " << endl
           << "  PERM     PERI     CNST" << endl;

      const int w = 6;
      for (ordinal_type i=0;i<_m;++i)
        os << setw(w) << _perm[i] << "   "
           << setw(w) << _peri[i] << "   "
           << setw(w) << _cnst[i] << endl;

      os.unsetf(ios::scientific);
      os.precision(prec);

      return os;
    }

  };

}

#endif
