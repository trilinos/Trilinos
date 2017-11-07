#ifndef __TACHO_DENSE_MATRIX_TOOLS_HPP__
#define __TACHO_DENSE_MATRIX_TOOLS_HPP__

/// \file Tacho_DenseMatrixTools.hpp
/// \brief Generic utility function for dense matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)  

#include "Tacho_Util.hpp"

namespace Tacho { 

  namespace Impl {
    
    // Serial and Team interface for dealing with dense matrices
    // - callable from device functor
    class DenseMatrixTools {
    public:

      struct Serial {
        /// \brief elementwise value set
        template<typename DenseMatrixTypeA, 
                 typename ValueType>
        KOKKOS_INLINE_FUNCTION
        static void
        set(/**/  DenseMatrixTypeA &A,
            const ValueType val) {
          
          const auto iend = A.NumRows();
          const auto jend = A.NumCols();
          for (auto j=0;j<jend;++j)
            for (auto i=0;i<iend;++i)
              A.Value(i,j) = val;
        }

        /// \brief elementwise copy 
        template<typename DenseMatrixTypeA, 
                 typename DenseMatrixTypeB>
        KOKKOS_INLINE_FUNCTION
        static void
        copy(/**/  DenseMatrixTypeA &B,
             const DenseMatrixTypeB &A) {
          const auto iend = B.NumRows();
          const auto jend = B.NumCols();
          
          for (auto j=0;j<jend;++j)
            for (auto i=0;i<iend;++i)
              B.Value(i,j) = A.Value(i,j);
        }

        /// \brief axpy
        template<typename DenseMatrixTypeA, 
                 typename DenseMatrixTypeB,
                 typename ValueType>
        KOKKOS_INLINE_FUNCTION
        static void
        axpy(/**/  DenseMatrixTypeA &B,
             const DenseMatrixTypeB &A,
             const ValueType alpha = 1.0) {
          const auto iend = B.NumRows();
          const auto jend = B.NumCols();
          
          for (auto j=0;j<jend;++j)
            for (auto i=0;i<iend;++i)
              B.Value(i,j) += alpha*A.Value(i,j);
        }

        template<typename DenseMatrixHierType,
                 typename DenseMatrixFlatType,
                 typename OrdinalType>
        KOKKOS_INLINE_FUNCTION
        static void
        getHierMatrix(DenseMatrixHierType &hier,
                      const DenseMatrixFlatType &flat,
                      const OrdinalType mb,
                      const OrdinalType nb) {
          const OrdinalType hm = hier.NumRows(), hn = hier.NumCols();
          const OrdinalType fm = flat.NumRows(), fn = flat.NumCols();
          
          typedef typename DenseMatrixHierType::ordinal_type ordinal_type;
          
          for (auto j=0;j<hn;++j) {
            const OrdinalType offn = nb*j;
            const OrdinalType ntmp = offn + nb; 
            const OrdinalType n    = ntmp < fn ? nb : (fn - offn); 
            
            for (auto i=0;i<hm;++i) {
              const OrdinalType offm = mb*i;
              const OrdinalType mtmp = offm + mb; 
              const OrdinalType m    = mtmp < fm ? mb : (fm - offm); 
              hier.Value(i, j).setView(flat, offm, m,
                                       /**/  offn, n);
            }
          }
        }
      };

    };
  }

  class DenseMatrixTools {
  public:

    /// Elementwise copy
    /// ------------------------------------------------------------------

    /// \brief elementwise copy 
    template<typename DenseMatrixTypeA, 
             typename DenseMatrixTypeB>
    inline
    static void
    copy(DenseMatrixTypeA &A,
         const DenseMatrixTypeB &B) {
      static_assert( Kokkos::Impl::is_same<
                     typename DenseMatrixTypeA::space_type,
                     typename DenseMatrixTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename DenseMatrixTypeA::space_type space_type;
      typedef typename DenseMatrixTypeA::ordinal_type ordinal_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;
      
      space_type::execution_space::fence();      

      Kokkos::parallel_for( range_policy(0, B.NumCols()), 
                            KOKKOS_LAMBDA(const ordinal_type j) 
                            {
                              //#pragma unroll
                              for (auto i=0;i<B.NumRows();++i)
                                A.Value(i,j) = B.Value(i,j);
                            } );

      space_type::execution_space::fence();
    }

    /// \brief elementwise copy 
    template<typename DenseMatrixTypeA, 
             typename DenseMatrixTypeB,
             typename PermutationVectorType>
    KOKKOS_INLINE_FUNCTION
    static void
    applyRowPermutation(DenseMatrixTypeA &A,
                        const DenseMatrixTypeB &B,
                        const PermutationVectorType &p) {
      static_assert( Kokkos::Impl::is_same<
                     typename DenseMatrixTypeA::space_type,
                     typename DenseMatrixTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename DenseMatrixTypeA::space_type space_type;
      typedef typename DenseMatrixTypeA::ordinal_type ordinal_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;
      
      space_type::execution_space::fence();      
      
      Kokkos::parallel_for( range_policy(0, B.NumCols()), 
                            KOKKOS_LAMBDA(const ordinal_type j) 
                            {
                              //#pragma unroll
                              for (auto i=0;i<B.NumRows();++i)
                                A.Value(p(i),j) = B.Value(i,j);
                            } );

      space_type::execution_space::fence();
    }

    /// \brief elementwise copy of lower/upper triangular of matrix including diagonals
    template<typename DenseMatrixTypeA, 
             typename DenseMatrixTypeB>
    KOKKOS_INLINE_FUNCTION
    static void
    copy(DenseMatrixTypeA &A,
         const int uplo,
         const int offset,
         const DenseMatrixTypeB &B) {
      static_assert( Kokkos::Impl
                     ::is_same<
                     typename DenseMatrixTypeA::space_type,
                     typename DenseMatrixTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename DenseMatrixTypeA::space_type space_type;
      typedef typename DenseMatrixTypeA::ordinal_type ordinal_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;
      
      space_type::execution_space::fence();      

      switch (uplo) {
      case Uplo::Lower: {
        Kokkos::parallel_for( range_policy(0, B.NumCols()),
                              KOKKOS_LAMBDA(const ordinal_type j)
                              {
                                //#pragma unroll
                                for (ordinal_type i=(j+offset);i<B.NumRows();++i)
                                  A.Value(i, j) = B.Value(i, j);
                              } );
        break;
      }
      case Uplo::Upper: {
        Kokkos::parallel_for( range_policy(0, B.NumCols()),
                              KOKKOS_LAMBDA(const ordinal_type j)
                              {
                                //#pragma unroll
                                for (ordinal_type i=0;i<=(j-offset);++i)
                                  A.Value(i, j) = B.Value(i, j);
                              } );
        break;
      }
      }

      space_type::execution_space::fence();
    }

    /// Flat to Hier
    /// ------------------------------------------------------------------

    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (o)

    /// \brief compute dimension of hier matrix when the flat is divided by mb x nb
    template<typename OrdinalType>
    KOKKOS_INLINE_FUNCTION
    static void
    getDimensionOfHierMatrix(OrdinalType &h,
                             const OrdinalType f,
                             const OrdinalType b) {
      const auto bb = (b == 0 ? f : b);

      h = f/bb + (f%bb > 0);
    }

    /// \brief compute dimension of hier matrix when the flat is divided by mb x nb
    template<typename DenseMatrixFlatType,
             typename OrdinalType>
    KOKKOS_INLINE_FUNCTION
    static void
    getDimensionOfHierMatrix(OrdinalType &hm,
                             OrdinalType &hn,
                             const DenseMatrixFlatType &flat,
                             const OrdinalType mb,
                             const OrdinalType nb) {
      getDimensionOfHierMatrix(hm, flat.NumRows(), mb);
      getDimensionOfHierMatrix(hn, flat.NumCols(), nb);
    }

      
    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (o)
    /// - Blocking with fence (o)

    template<typename DenseMatrixHierType,
             typename DenseMatrixFlatType,
             typename OrdinalType>
    inline
    static void
    getHierMatrix(DenseMatrixHierType &hier,
                  const DenseMatrixFlatType &flat,
                  const OrdinalType mb,
                  const OrdinalType nb) {
      static_assert( Kokkos::Impl::is_same<
                     typename DenseMatrixHierType::space_type,
                     typename DenseMatrixFlatType::space_type
                     >::value,
                     "Space type of input matrices does not match" );

      const OrdinalType hm = hier.NumRows(), hn = hier.NumCols();
      const OrdinalType fm = flat.NumRows(), fn = flat.NumCols();
      
      typedef typename DenseMatrixHierType::space_type space_type;
      typedef typename DenseMatrixHierType::ordinal_type ordinal_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;

      space_type::fence();
      
      Kokkos::parallel_for( range_policy(0, hn),
                            KOKKOS_LAMBDA(const ordinal_type j) 
        {
        const OrdinalType offn = nb*j;
        const OrdinalType ntmp = offn + nb; 
        const OrdinalType n    = ntmp < fn ? nb : (fn - offn); 
        
        //#pragma unroll
        for (auto i=0;i<hm;++i) {
          const OrdinalType offm = mb*i;
          const OrdinalType mtmp = offm + mb; 
          const OrdinalType m    = mtmp < fm ? mb : (fm - offm); 
          hier.Value(i, j).setView(flat, offm, m,
                                   /**/  offn, n);
        }
      });
    
      space_type::fence();
    }


    /// \brief fill hier matrix 
    // template<typename DenseMatrixHierType,
    //          typename DenseMatrixFlatType,
    //          typename OrdinalType>
    // KOKKOS_INLINE_FUNCTION
    // static void
    // getHierMatrix(DenseMatrixHierType &hier,
    //               const DenseMatrixFlatType &flat,
    //               const OrdinalType mb,
    //               const OrdinalType nb) {
    //   static_assert( Kokkos::Impl::is_same<
    //                  typename DenseMatrixHierType::space_type,
    //                  typename DenseMatrixFlatType::space_type
    //                  >::value,
    //                  "Space type of input matrices does not match" );
      
    //   typedef typename DenseMatrixHierType::space_type space_type;
    //   typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;

    //   const OrdinalType hm = hier.NumRows(), hn = hier.NumCols();
    //   const OrdinalType fm = flat.NumRows(), fn = flat.NumCols();

    //   space_type::execution_space::fence();

    //   Kokkos::parallel_for( range_policy(0, hn),
    //                         [&](const ordinal_type j)
    //                         {
    //                           const OrdinalType offn = nb*j;
    //                           const OrdinalType ntmp = offn + nb; 
    //                           const OrdinalType n    = ntmp < fn ? nb : (fn - offn); 

    //                           //#pragma unroll
    //                           for (ordinal_type i=0;i<hm;++i) {
    //                             const OrdinalType offm = mb*i;
    //                             const OrdinalType mtmp = offm + mb; 
    //                             const OrdinalType m    = mtmp < fm ? mb : (fm - offm); 
    //                             hier.Value(i, j).setView(flat, offm, m,
    //                                                      /**/  offn, n);
    //                           }
    //                         } );

    //   space_type::execution_space::fence();
    // }

    template<typename DenseMatrixHierType,
             typename DenseMatrixFlatType,
             typename OrdinalTypeArray,
             typename OrdinalType>
    KOKKOS_INLINE_FUNCTION
    static void
    getHierMatrix(DenseMatrixHierType &hier,
                  const DenseMatrixFlatType &flat,
                  const OrdinalTypeArray range,
                  const OrdinalType nb) {
      static_assert( Kokkos::Impl::is_same<
                     typename DenseMatrixHierType::space_type,
                     typename DenseMatrixFlatType::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename DenseMatrixHierType::ordinal_type ordinal_type;
      typedef typename DenseMatrixHierType::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;

      const OrdinalType hm = hier.NumRows(), hn = hier.NumCols();
      const OrdinalType fn = flat.NumCols();

      space_type::execution_space::fence();

      Kokkos::parallel_for( range_policy(0, hn),
                            KOKKOS_LAMBDA(const ordinal_type j)
                            {
                              const OrdinalType offn = nb*j;
                              const OrdinalType ntmp = offn + nb; 
                              const OrdinalType n    = ntmp < fn ? nb : (fn - offn); 

                              //#pragma unroll
                              for (ordinal_type i=0;i<hm;++i) {
                                const OrdinalType offm = range(i);
                                const OrdinalType m    = range(i+1) - range(i);
                                hier.Value(i, j).setView(flat, offm, m,
                                                         /**/  offn, n);
                              }
                            } );

      space_type::execution_space::fence();
    }

    /// \brief create hier matrix 
    template<typename DenseMatrixHierType,
             typename DenseMatrixFlatType,
             typename OrdinalType>
    inline
    static void
    createHierMatrix(DenseMatrixHierType &hier,
                     const DenseMatrixFlatType &flat,
                     const OrdinalType mb,
                     const OrdinalType nb) {
      OrdinalType hm, hn;
      getDimensionOfHierMatrix(hm, hn,
                               flat, 
                               mb, nb);
      
      hier.create(hm, hn);
      
      getHierMatrix(hier, flat, mb , nb);
    }

    /// \brief create hier matrix 
    template<typename DenseMatrixHierType,
             typename DenseMatrixFlatType,
             typename OrdinalType,
             typename OrdinalTypeArray>
    inline
    static void
    createHierMatrix(DenseMatrixHierType &hier,
                     const DenseMatrixFlatType &flat,
                     const OrdinalType nblks,
                     const OrdinalTypeArray range,
                     const OrdinalType nb) {
      OrdinalType hn;
      getDimensionOfHierMatrix(hn, flat.NumCols(), nb);
      
      hier.create(nblks, hn);
      
      getHierMatrix(hier, flat, range , nb);
    }
    
  };

}

#endif



//     /// \brief elementwise copy of matrix b
//     template<typename DenseMatrixTypeA, 
//              typename DenseMatrixTypeB,
//              typename OrdinalArrayTypeIp,
//              typename OrdinalArrayTypeJp>
//     KOKKOS_INLINE_FUNCTION
//     static void
//     copy(DenseMatrixTypeA &A,
//          const DenseMatrixTypeB &B,
//          const OrdinalArrayTypeIp &ip,         
//          const OrdinalArrayTypeJp &jp) { 
      
//       static_assert( Kokkos::Impl::is_same<DenseMatrixTypeA::space_type,DenseMatrixTypeB::space_type>::value,
//                      "Space type of input matrices does not match" );
      
//       typedef DenseMatrixTypeA::space_type space_type;
//       typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;
      
//       const int idx = ip.is_null() * 10 + jp.is_null();

//       space_type::execution_space::fence();      

//       switch (idx) {
//       case 11: { // ip is null and jp is null: no permutation
//         Kokkos::parallel_for( range_policy(0, B.NumCols()), 
//                               KOKKOS_LAMBDA(const ordinal_type j) 
//                               {
// #pragma unroll
//                                 for (auto i=0;i<B.NumRows();++i)
//                                   A.Value(i,j) = B.Value(i,j);
//                               } );
//         break;
//       }
//       case 0: { // ip is not null and jp is not null: row/column permutation 
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               KOKKOS_LAMBDA(const ordinal_type j) 
//                               {
// #pragma unroll 
//                                 for (auto i=0;i<B._m;++i)
//                                   A.Value(i, j) = B.Value(ip(i), jp(j));
//                               } );
//         break;
//       }
//       case 10: { // ip is not null and jp is null: row permutation
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               [&](const ordinal_type j) 
//                               {
// #pragma unroll 
//                                 for (auto i=0;i<B._m;++i)
//                                   A.Value(i, j) = B.Value(ip(i), j);
//                               } );
//         break;
//       } 
//       case 1: { // ip is null and jp is not null: column permutation
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               [&](const ordinal_type j) 
//                               {
//                                 const ordinal_type jj = jp(j);
// #pragma unroll 
//                                 for (auto i=0;i<B._m;++i)
//                                   A.Value(i, j) = B.Value(i, jj);
//                               } );
//         break;
//       }
//       }

//       space_type::execution_space::fence();
//     }

//     /// \brief elementwise copy of lower/upper triangular of matrix b
//     template<typename VT,
//              typename OT,
//              typename ST>
//     KOKKOS_INLINE_FUNCTION
//     void
//     copy(const int uplo, 
//          const DenseMatrixBase<VT,OT,ST,space_type> &b) { 

//       typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Dynamic> > range_policy;

//       space_type::execution_space::fence();

//       switch (uplo) {
//       case Uplo::Lower: {
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               [&](const ordinal_type j) 
//                               { 
// #pragma unroll 
//                                 for (ordinal_type i=j;i<B._m;++i) 
//                                   A.Value(i, j) = B.Value(i, j);
//                               } );
//         break;
//       }
//       case Uplo::Upper: {
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               [&](const ordinal_type j) 
//                               { 
// #pragma unroll 
//                                 for (ordinal_type i=0;i<(j+1);++i) 
//                                   A.Value(i, j) = B.Value(i, j);
//                               } );
//         break;
//       }
//       }

//       space_type::execution_space::fence();
//     }
    /// ------------------------------------------------------------------
