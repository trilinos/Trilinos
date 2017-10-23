#ifndef __TACHOEXP_MATRIX_MARKET_HPP__
#define __TACHOEXP_MATRIX_MARKET_HPP__

/// \file TachoExp_MatrixMarket.hpp
/// \brief IO utilities for matrix market format
/// \author Kyungjoo Kim (kyukim@sandia.gov)  

#include "TachoExp_Util.hpp"
#include "TachoExp_CrsMatrixBase.hpp"

namespace Tacho { 
  namespace Experimental {
    
    template<typename ValueType>
    struct MatrixMarket {
 
      /// \brief matrix market reader
      static CrsMatrixBase<ValueType,Kokkos::DefaultHostExecutionSpace>
      read(const std::string &filename) {
        std::ifstream file;
        file.open(filename);
        // TACHO_TEST_FOR_EXCEPTION(!file.good(), std::runtime_error, 
        //                          std::string("failed to open" + filename));
      
        // reading mm header
        ordinal_type m, n;
        size_type nnz;
        bool symmetry;
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
          symmetry = (header.find("symmetric") != std::string::npos);

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
            file >> row >> col >> val;
          
            row -= mm_base;
            col -= mm_base;
          
            mm.push_back(ijv_type(row, col, val));
            if (symmetry && row != col)
              mm.push_back(ijv_type(col, row, val));
          }
          std::sort(mm.begin(), mm.end(), std::less<ijv_type>());
        }
      
        // change mm to crs
        typedef Kokkos::DefaultHostExecutionSpace host_space;
        Kokkos::View<size_type*,   host_space> ap("ap", m+1);
        Kokkos::View<ordinal_type*,host_space> aj("aj", nnz);
        Kokkos::View<value_type*,  host_space> ax("ax", nnz);
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
        CrsMatrixBase<value_type,host_space> A;
        A.setExternalMatrix(m, n, nnz, ap, aj, ax);
        
        return A;
      }

      /// \brief matrix marker writer
      template<typename uplo=void>
      static void
      write(std::ofstream &file, 
            const CrsMatrixBase<ValueType,Kokkos::DefaultHostExecutionSpace> &A,
            const std::string comment = "%% Tacho::MatrixMarket::Export") {
        typedef ValueType value_type;

        std::streamsize prec = file.precision();
        file.precision(16);
        file << std::scientific;
      
        {
          file << "%%MatrixMarket matrix coordinate "
               << (is_complex_type<value_type>::value ? "complex " : "real ")
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
                file << std::setw(w) << ( i+1) << "  "
                     << std::setw(w) << (aj+1) << "  "
                     << std::setw(w) <<    val << std::endl;
              }
            }
          }
        }      

        file.unsetf(std::ios::scientific);
        file.precision(prec);
      }
    };
  }
}

#endif
