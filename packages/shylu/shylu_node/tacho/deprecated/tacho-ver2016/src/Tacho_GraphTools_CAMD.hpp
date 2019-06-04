#ifndef __TACHO_GRAPH_TOOLS_CAMD_HPP__
#define __TACHO_GRAPH_TOOLS_CAMD_HPP__

/// \file Tacho_GraphTools_CAMD.hpp
/// \brief Interface to camd (suire sparse) ordering
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

//#ifdef HAVE_SHYLU_NODETACHO_AMESOS
#include "trilinos_camd.h"
#define TACHO_CHOLMOD(run) trilinos_ ## run
typedef UF_long SuiteSparse_long;
// #else
// #include "camd.h"
// #define TACHO_CHOLMOD(run) run
// #endif

namespace Tacho {

  class CAMD {
  public:
    template<typename OrdinalType>
    static void run( OrdinalType n, OrdinalType Pe[], OrdinalType Iw[], 
                     OrdinalType Len[], OrdinalType iwlen, OrdinalType pfree,
                     OrdinalType Nv[], OrdinalType Next[], OrdinalType Last[],
                     OrdinalType Head[], OrdinalType Elen[], OrdinalType Degree[],
                     OrdinalType W[],
                     double Control[],
                     double Info[],
                     const OrdinalType C[],
                     OrdinalType BucketSet[] ) {
      TACHO_TEST_FOR_ABORT(true, "GraphTools_CAMD:: CAMD does not support the ordinal type");
    }
  };
  
  template<> 
  void CAMD::run<int>( int n, int Pe[], int Iw[], 
                       int Len[], int iwlen, int pfree,
                       int Nv[], int Next[], int Last[],
                       int Head[], int Elen[], int Degree[],
                       int W[],
                       double Control[],
                       double Info[],
                       const int C[],
                       int BucketSet[] ) {
    TACHO_CHOLMOD(camd_2)( n, Pe, Iw, Len, iwlen, pfree, 
                           Nv, Next, Last, Head, Elen, Degree, W, Control, Info, C, BucketSet );
  }
  
  template<> 
  void CAMD::run<SuiteSparse_long>( SuiteSparse_long n, SuiteSparse_long Pe[], SuiteSparse_long Iw[], 
                                    SuiteSparse_long Len[], SuiteSparse_long iwlen, SuiteSparse_long pfree,
                                    SuiteSparse_long Nv[], SuiteSparse_long Next[], SuiteSparse_long Last[],
                                    SuiteSparse_long Head[], SuiteSparse_long Elen[], SuiteSparse_long Degree[],
                                    SuiteSparse_long W[],
                                    double Control[],
                                    double Info[],
                                    const SuiteSparse_long C[],
                                    SuiteSparse_long BucketSet[] ) {
    TACHO_CHOLMOD(camd_l2)( n, Pe, Iw, Len, iwlen, pfree, 
                            Nv, Next, Last, Head, Elen, Degree, W, Control, Info, C, BucketSet );
  }
  
  template<typename OrdinalType,
           typename SizeType = OrdinalType,
           typename SpaceType = void>
  class GraphTools_CAMD {    
  public:
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;
    typedef SpaceType   space_type;

    typedef Kokkos::View<ordinal_type*, space_type> ordinal_type_array;
    typedef Kokkos::View<size_type*,    space_type> size_type_array;

  private:
    ordinal_type _m;
    size_type _nnz;
    size_type_array _rptr;
    ordinal_type_array _cidx, _cnst;
    
    // CAMD output
    ordinal_type_array _pe, _nv, _el, _next, _perm, _peri; // perm = last, peri = next
    
    double _control[CAMD_CONTROL], _info[CAMD_INFO];

    bool _is_ordered;

  public:

    size_type NumNonZeros() const { return _nnz; }
    ordinal_type NumRows() const { return _m; }

    ordinal_type_array PermVector()       const { return _perm; }
    ordinal_type_array InvPermVector()    const { return _peri; }

    ordinal_type_array ConstraintVector() const { return _cnst; }

    // static assert is necessary to enforce to use host space only
    GraphTools_CAMD() = default;
    GraphTools_CAMD(const GraphTools_CAMD &b) = default;
    virtual~GraphTools_CAMD() = default;

    void setGraph(const ordinal_type m,
                  const size_type_array rptr,
                  const ordinal_type_array cidx,
                  const ordinal_type nblk,
                  const ordinal_type_array range) {
      // graph information
      _m     = m;
      _nnz   = rptr[m];

      /// CAMD graph spec
      /// - no diagonals, symmetric

      _rptr  = rptr; 
      _cidx  = cidx; 

      // constraints are induced from range
      _cnst  = ordinal_type_array("CAMD::ConstraintArray", _m+1);
      for (ordinal_type i=0;i<nblk;++i)
        for (ordinal_type j=range[i];j<range[i+1];++j)
          _cnst[j] = i;

      // permutation vector
      _pe    = ordinal_type_array("CAMD::EliminationArray", _m);
      _nv    = ordinal_type_array("CAMD::SupernodesArray", _m);
      _el    = ordinal_type_array("CAMD::DegreeArray", _m);
      _next  = ordinal_type_array("CAMD::InvPermSupernodesArray", _m);
      _perm  = ordinal_type_array("CAMD::PermutationArray", _m);
      _peri  = ordinal_type_array("CAMD::InvPermutationArray", _m);
    }

    void computeOrdering() {
      TACHO_CHOLMOD(camd_defaults)(_control);
      TACHO_CHOLMOD(camd_control)(_control);

      ordinal_type *rptr = reinterpret_cast<ordinal_type*>(_rptr.data());
      ordinal_type *cidx = reinterpret_cast<ordinal_type*>(_cidx.data());
      ordinal_type *cnst = reinterpret_cast<ordinal_type*>(_cnst.data());

      ordinal_type *next = reinterpret_cast<ordinal_type*>(_next.data());
      ordinal_type *perm = reinterpret_cast<ordinal_type*>(_perm.data());

      // length array
      ordinal_type_array lwork("CAMD::LWorkArray", _m);
      ordinal_type *lwork_ptr = reinterpret_cast<ordinal_type*>(lwork.data());
      for (ordinal_type i=0;i<_m;++i)
        lwork_ptr[i] = rptr[i+1] - rptr[i];

      // workspace
      const size_type swlen = _nnz + _nnz/5 + 5*(_m+1);;
      ordinal_type_array swork("CAMD::SWorkArray", swlen);
      ordinal_type *swork_ptr = reinterpret_cast<ordinal_type*>(swork.data());

      ordinal_type *pe_ptr = reinterpret_cast<ordinal_type*>(_pe.data()); // 1) Pe
      size_type pfree = 0;
      for (ordinal_type i=0;i<_m;++i) {
        pe_ptr[i] = pfree;
        pfree += lwork_ptr[i];
      }
      TACHO_TEST_FOR_ABORT( _nnz != pfree, ">> nnz in the graph does not match to nnz count (pfree)");

      ordinal_type *nv_ptr = reinterpret_cast<ordinal_type*>(_nv.data()); // 2) Nv
      ordinal_type *hd_ptr = swork_ptr; swork_ptr += (_m+1);   // 3) Head
      ordinal_type *el_ptr = reinterpret_cast<ordinal_type*>(_el.data()); // 4) Elen
      ordinal_type *dg_ptr = swork_ptr; swork_ptr += _m;       // 5) Degree
      ordinal_type *wk_ptr = swork_ptr; swork_ptr += (_m+1);   // 6) W
      ordinal_type *bk_ptr = swork_ptr; swork_ptr += _m;       // 7) BucketSet

      const size_type iwlen = swlen - (4*_m+2);
      ordinal_type *iw_ptr = swork_ptr; swork_ptr += iwlen;    // Iw
      for (ordinal_type i=0;i<pfree;++i)
        iw_ptr[i] = cidx[i];

      CAMD::run(_m, pe_ptr, iw_ptr, lwork_ptr, iwlen, pfree,
                // output
                nv_ptr, next, perm, hd_ptr, el_ptr, dg_ptr, wk_ptr, 
                _control, _info, cnst, bk_ptr);
      
      TACHO_TEST_FOR_ABORT(_info[CAMD_STATUS] != CAMD_OK, "CAMD fails");

      for (ordinal_type i=0;i<_m;++i)
        _peri[_perm[i]] = i;

      _is_ordered = true;
    }

    std::ostream& showMe(std::ostream &os) const {
      std::streamsize prec = os.precision();
      os.precision(8);
      os << std::scientific;

      os << " -- CAMD input -- " << std::endl
         << "    # of Rows      = " << _m << std::endl
         << "    # of NonZeros  = " << _nnz << std::endl;

      if (_is_ordered)
        os << " -- Ordering -- " << std::endl
           << "  CNST     PERM     PERI       PE       NV     NEXT     ELEN" << std::endl;

      const int w = 6;
      for (ordinal_type i=0;i<_m;++i)
        os << std::setw(w) << _cnst[i] << "   "
           << std::setw(w) << _perm[i] << "   "
           << std::setw(w) << _peri[i] << "   "
           << std::setw(w) << _pe[i] << "   "
           << std::setw(w) << _nv[i] << "   "
           << std::setw(w) << _next[i] << "   "
           << std::setw(w) << _el[i] << "   "
           << std::endl;

      os.unsetf(std::ios::scientific);
      os.precision(prec);

      return os;
    }

  };

}

#endif
