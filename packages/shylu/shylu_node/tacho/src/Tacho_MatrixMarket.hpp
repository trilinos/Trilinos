#ifndef __TACHO_MATRIX_MARKET_HPP__
#define __TACHO_MATRIX_MARKET_HPP__

/// \file Tacho_MatrixMarket.hpp
/// \brief IO utilities for matrix market format
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"

namespace Tacho {
    
    template<typename ValueType>
    inline
    typename std::enable_if<std::is_same<ValueType,double>::value ||
                            std::is_same<ValueType,float>::value>::type
    impl_read_value_from_file(std::ifstream &file,
                              ordinal_type &row,
                              ordinal_type &col,
                              ValueType &val) {
      file >> row >> col >> val;
    }  
    template<typename ValueType>
    inline
    typename std::enable_if<std::is_same<ValueType,Kokkos::complex<double> >::value ||
                            std::is_same<ValueType,Kokkos::complex<float> >::value>::type
    impl_read_value_from_file(std::ifstream &file,
                              ordinal_type &row,
                              ordinal_type &col,
                              ValueType &val) {
      typename ArithTraits<ValueType>::mag_type r, i;
      file >> row >> col >> r >> i;
      val = ValueType(r, i);
    }

    template<typename ValueType>
    inline
    typename std::enable_if<std::is_same<ValueType,double>::value ||
                            std::is_same<ValueType,float>::value>::type
    impl_write_value_to_file(std::ofstream &file,
                             const ordinal_type row,
                             const ordinal_type col,
                             const ValueType val,
                             const ordinal_type w = 10) {
      file << std::setw(w) << row << "  "
           << std::setw(w) << col << "  "
           << std::setw(w) << val << std::endl;
    }  


    template<typename ValueType>
    inline
    typename std::enable_if<std::is_same<ValueType,Kokkos::complex<double> >::value ||
                            std::is_same<ValueType,Kokkos::complex<float> >::value>::type
    impl_write_value_to_file(std::ofstream &file,
                             const ordinal_type row,
                             const ordinal_type col,
                             const ValueType val,
                             const ordinal_type w = 10) {
      file << std::setw(w) << row << "  "
           << std::setw(w) << col << "  "
           << std::setw(w) << val.real() << "  "
           << std::setw(w) << val.imag() << std::endl;
    }  

    template<typename ValueType>
    struct MatrixMarket {

      /// \brief matrix market reader
      template<typename ExecSpaceType>
      static void
      read(const std::string &filename, 
           CrsMatrixBase<ValueType,ExecSpaceType> &A,
           const ordinal_type verbose = 0) {
        static_assert(Kokkos::Impl::MemorySpaceAccess< 
                      Kokkos::HostSpace, 
                      typename ExecSpaceType::memory_space >::assignable, 
                      "ExecSpaceType is not assignable from HostSpace" );
        
        typedef ArithTraits<ValueType> arith_traits;
        Kokkos::Impl::Timer timer;

        timer.reset();

        std::ifstream file;
        file.open(filename);
        // TACHO_TEST_FOR_EXCEPTION(!file.good(), std::runtime_error,
        //                          std::string("failed to open" + filename));

        // reading mm header
        ordinal_type m, n;
        size_type nnz;
        bool symmetry = false, hermitian = false; //, cmplx = false;
        {
          std::string header;
          std::getline(file, header);
          while (file.good()) {
            char c = file.peek();
            if (c == '%' || c == '\n') {
              file.ignore(256, '\n');
              continue;
            }
            break;
          }
          std::transform(header.begin(), header.end(), header.begin(), ::tolower);
          symmetry = (header.find("symmetric") != std::string::npos || 
                      header.find("hermitian") != std::string::npos);

          hermitian = (header.find("hermitian") != std::string::npos);

          file >> m >> n >> nnz;
        }

        // read data into coo format
        const ordinal_type mm_base = 1;

        typedef ValueType value_type;
        typedef Coo<value_type> ijv_type;
        std::vector<ijv_type> mm;
        {
          mm.reserve(nnz*(symmetry ? 2 : 1));
          for (size_type i=0;i<nnz;++i) {
            ordinal_type row, col;
            value_type val;

            impl_read_value_from_file(file, row, col, val);

            row -= mm_base;
            col -= mm_base;

            mm.push_back(ijv_type(row, col, val));
            if (symmetry && row != col) 
              mm.push_back(ijv_type(col, row, hermitian ? arith_traits::conj(val) : val));
          }
          std::sort(mm.begin(), mm.end(), std::less<ijv_type>());

          // update nnz
          nnz = mm.size();
        }

        // change mm to crs
        Kokkos::View<size_type*,   ExecSpaceType> ap("ap", m+1);
        Kokkos::View<ordinal_type*,ExecSpaceType> aj("aj", nnz);
        Kokkos::View<value_type*,  ExecSpaceType> ax("ax", nnz);
        {
          ordinal_type icnt = 0;
          size_type jcnt = 0;
          ijv_type prev = mm[0];

          ap[icnt++] = 0;
          aj[jcnt]   = prev.j;
          ax[jcnt++] = prev.val;

          for (auto it=(mm.begin()+1);it<(mm.end());++it) {
            const ijv_type aij = (*it);

            if (aij.i != prev.i)
              ap[icnt++] = jcnt;

            if (aij == prev) {
              aj[jcnt]  = aij.j;
              ax[jcnt] += aij.val;
            } else {
              aj[jcnt]   = aij.j;
              ax[jcnt++] = aij.val;
            }
            prev = aij;
          }
          ap[icnt++] = jcnt;
          nnz = jcnt;
        }

        // create crs matrix view
        A.clear();
        A.setExternalMatrix(m, n, nnz, ap, aj, ax);

        const double t = timer.seconds();
        if (verbose) {

          printf("Summary: MatrixMarket\n");
          printf("=====================\n");
          printf("  File:      %s\n", filename.c_str());
          printf("  Time\n");
          printf("             time for reading A:                              %10.6f s\n", t);
          printf("\n");
          printf("  Sparse Matrix (%s) \n", (symmetry ? "symmetric" : "non-symmetric"));
          printf("             number of rows:                                  %10d\n", m);
          printf("             number of cols:                                  %10d\n", n);
          printf("             number of nonzeros:                              %10d\n", ordinal_type(nnz));
          printf("\n");
        }
      }

      /// \brief matrix marker writer
      template<typename ExecSpaceType, 
               typename uplo=void>
      static void
      write(std::ofstream &file,
            const CrsMatrixBase<ValueType,ExecSpaceType> &A,
            const std::string comment = "%% Tacho::MatrixMarket::Export") {
        static_assert(Kokkos::Impl::MemorySpaceAccess< 
                      Kokkos::HostSpace, 
                      typename ExecSpaceType::memory_space
                      >::assignable, 
                      "ExecSpaceType is not assignable from HostSpace" );

        typedef ArithTraits<ValueType> arith_traits;
        typedef ValueType value_type;

        std::streamsize prec = file.precision();
        file.precision(16);
        file << std::scientific;

        {
          file << "%%MatrixMarket matrix coordinate "
               << (arith_traits::is_complex ? "complex " : "real ")
               << (is_valid_uplo_tag<uplo>::value ? "symmetric " : "general ")
               << std::endl;
          file << comment << std::endl;
        }
        // cnt nnz
        size_type nnz = 0;
        {
          for (ordinal_type i=0;i<A.NumRows();++i) {
            const size_type jbegin = A.RowPtrBegin(i), jend = A.RowPtrEnd(i);
            for (size_type j=jbegin;j<jend;++j) {
              const auto aj = A.Col(j);
              if (std::is_same<uplo,Uplo::Upper>::value && i <= aj) ++nnz;
              if (std::is_same<uplo,Uplo::Lower>::value && i >= aj) ++nnz;
              if (std::is_same<uplo,void>::value) ++nnz;
            }
          }
          file << A.NumRows() << " " << A.NumCols() << " " << nnz << std::endl;
        }

        const int w = 10;
        {
          for (ordinal_type i=0;i<A.NumRows();++i) {
            const size_type jbegin = A.RowPtrBegin(i), jend = A.RowPtrEnd(i);
            for (size_type j=jbegin;j<jend;++j) {
              const auto aj = A.Col(j);
              bool flag = false;
              if (std::is_same<uplo,Uplo::Upper>::value && i <= aj) flag = true;
              if (std::is_same<uplo,Uplo::Lower>::value && i >= aj) flag = true;
              if (std::is_same<uplo,void>::value) flag = true;
              if (flag) {
                value_type val = A.Value(j);
                impl_write_value_to_file(file, i+1, aj+1, val, w);
              }
            }
          }
        }

        file.unsetf(std::ios::scientific);
        file.precision(prec);
      }
    };

}

#endif
