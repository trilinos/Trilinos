#ifndef __TACHO_CRS_DATA_HPP__
#define __TACHO_CRS_DATA_HPP__

/// \file Tacho_CrsData.hpp
/// \brief IO utilities for crs data
/// \author Kyungjoo Kim (kyukim@sandia.gov)  

#include "Tacho_Util.hpp"

namespace Tacho { 

  class CrsData {
  public:
    /// \brief read data file (base is 0)
    template<typename CrsMatrixType>
    static void
    read(CrsMatrixType &A, std::ifstream &file) {
      //typedef typename CrsMatrixType::value_type      value_type;
      typedef typename CrsMatrixType::ordinal_type    ordinal_type;
      typedef typename CrsMatrixType::size_type       size_type;
      typedef typename CrsMatrixType::size_type_array size_type_array;

      //typedef typename CrsMatrixType::space_type space_type;
      //typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;
      
      ordinal_type m, n;
      size_type nnz;
      
      {
        std::string header;
        if (file.is_open()) {
          std::getline(file, header);
          while (file.good()) {
            char c = file.peek();
            if (c == '%' || c == '\n') {
              file.ignore(256, '\n');
              continue;
            }
            break;
          }
        } else {
          TACHO_TEST_FOR_EXCEPTION(true, std::runtime_error, MSG_INVALID_INPUT(file));
        }

        // read matrix specification
        file >> m >> n >> nnz;
        
        // create crs matrix 
        A.create(m, n, nnz);        

        // read column pointers
        size_type_array ap("CrsData::read::ap", m+1);
        for (size_type i=0;i<=m;++i)
          file >> ap(i);

        Kokkos::deep_copy(A.RowPtrBegin(), Kokkos::subview(ap, Kokkos::pair<ordinal_type,ordinal_type>(0, m  )));
        Kokkos::deep_copy(A.RowPtrEnd(),   Kokkos::subview(ap, Kokkos::pair<ordinal_type,ordinal_type>(1, m+1)));
        
        auto aj = A.Cols();
        auto ax = A.Values();

        for (size_type i=0;i<nnz;++i)
          file >> aj(i) >> ax(i);
      }
    }

    /// \brief write
    template<typename CrsMatrixType>
    static void
    write(std::ofstream &file, 
          const CrsMatrixType A,
          const std::string comment = "%% Tacho::CrsData::Export") {
      typedef typename CrsMatrixType::value_type      value_type;
      typedef typename CrsMatrixType::ordinal_type    ordinal_type;
      typedef typename CrsMatrixType::size_type       size_type;
      typedef typename CrsMatrixType::size_type_array size_type_array;

      std::streamsize prec = file.precision();
      file.precision(16);
      file << std::scientific;
      
      file << "%%Tacho Crs data "
           << (Util::isComplex<value_type>() ? "complex " : "real ")
           << std::endl;
      
      file << comment << std::endl;

      // count nnz and contiguous row ptr
      const ordinal_type m = A.NumRows(), n = A.NumCols();
      size_type_array ap("CraData::write::ap", m + 1);

      size_type nnz = 0;
      ap[0] = 0;
      for (ordinal_type i=0;i<m;++i) {
        const auto nnz_in_row = A.RowPtrEnd(i) - A.RowPtrBegin(i);
        nnz += nnz_in_row;
        ap[i+1] = nnz;
      }

      // print out matrix specification
      file << m << " " << n << " " << nnz << std::endl;

      // print out row ptr
      const int w = 10;
      for (ordinal_type i=0;i<=m;++i) 
        file << std::setw(w) << ap[i] << std::endl;

      for (ordinal_type i=0;i<m;++i) {
        const size_type jbegin = A.RowPtrBegin(i), jend = A.RowPtrEnd(i);
        for (size_type j=jbegin;j<jend;++j) {
          const auto aj = A.Col(j);
          const auto ax = A.Value(j);
            file << std::setw(w) << aj << "  "
                 << std::setw(w) << ax << std::endl;
        }
      }
      file.unsetf(std::ios::scientific);
      file.precision(prec);
    }
    
  };

}

#endif
