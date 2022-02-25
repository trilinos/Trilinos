/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_applyDirichletBoundaryCondition.hpp"

namespace { // (anonymous)
  
  //
  // UNIT TESTS
  //
  
  template<class OutputType, class InputType>
  struct ToValue {
    static KOKKOS_INLINE_FUNCTION OutputType
    toValue (const InputType& x) {
      return static_cast<OutputType> (x);
    }
  };

  template<class OutputRealType, class InputType>
  struct ToValue<Kokkos::complex<OutputRealType>, InputType> {
    using output_type = Kokkos::complex<OutputRealType>;
    static KOKKOS_INLINE_FUNCTION output_type
    toValue (const InputType& x) {
      return static_cast<output_type> (static_cast<OutputRealType> (x));
    }
  };

  template<class OutputType, class InputType>
  KOKKOS_INLINE_FUNCTION OutputType toValue (const InputType& x) {
    return ToValue<OutputType, InputType>::toValue (x);
  }

   TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ApplyDirichletRows, SC, LO, GO, NT )
  {
    using Tpetra::createContigMapWithNode;
    using Teuchos::RCP;
    using std::endl;
    using crs_matrix_type = Tpetra::CrsMatrix<SC,LO,GO,NT>;
    using map_type = Tpetra::Map<LO,GO,NT>;
    using vec_type = Tpetra::Vector<SC,LO,GO,NT>;
    using mag_type = typename vec_type::mag_type;
    using IST = typename vec_type::impl_scalar_type;    
    using GST = Tpetra::global_size_t;
    using STS = Teuchos::ScalarTraits<SC>;
    using KAT = Kokkos::ArithTraits<IST>;    
    using local_matrix_device_type = 
          typename crs_matrix_type::local_matrix_device_type;
    using local_graph_device_type = 
          typename local_matrix_device_type::staticcrsgraph_type;
    using device_type = typename crs_matrix_type::device_type;
    using execution_space = typename crs_matrix_type::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const auto comm = Tpetra::getDefaultComm ();
    const LO lclNumRows = 10;
    RCP<const map_type> rowMap =
      createContigMapWithNode<LO,GO,NT> (INVALID, lclNumRows,
                                         comm);
    RCP<const map_type> colMap = rowMap;
    RCP<const map_type> domMap = rowMap;
    RCP<const map_type> ranMap = rowMap;
    
    vec_type vec1 (domMap);
    vec_type vec2 (ranMap), vec2_compare(ranMap);

    vec1.putScalar (STS::one ());
    vec2.putScalar (STS::zero ());
    vec2_compare.putScalar (STS::zero ());

    using row_offsets_type =
      typename local_graph_device_type::row_map_type::non_const_type;
    using lcl_col_inds_type =
      typename local_graph_device_type::entries_type::non_const_type;

    // Device filled BCs
    Kokkos::View<typename crs_matrix_type::local_ordinal_type*, device_type> lclRowInds ("lclRowInds", 1);
    Kokkos::parallel_for
      ("Fill lclRowInds",
       range_type (0, 1),
       KOKKOS_LAMBDA (const LO lclRow) {
	  lclRowInds(lclRow) = lclNumRows-1;
      });
    
    // Host filled BCs
    Kokkos::View<typename crs_matrix_type::local_ordinal_type*, Kokkos::HostSpace> lclRowInds_h ("lclRowsInds_h",1);
    for (LO k = 0; k < 1; ++k) {
      lclRowInds_h[0] = lclNumRows-1;
    }

    // Input matrix: Dense all 2's.
    local_matrix_device_type I_A_lcl;
    {
      row_offsets_type I_rowOffsets ("rowOffsets", lclNumRows+1);    
      lcl_col_inds_type I_lclColInds ("lclColInds", lclNumRows*lclNumRows);
      Kokkos::View<IST*, device_type> I_values ("values", lclNumRows*lclNumRows);
      Kokkos::parallel_for
        ("Input Matrix",
         range_type (0, lclNumRows),
         KOKKOS_LAMBDA (const LO lclRow) {
          I_rowOffsets(lclRow+1) = (lclRow+1)*lclNumRows;
          LO base = lclRow*lclNumRows;
        for(LO j=0; j<lclNumRows; j++) {
          I_lclColInds(base + j) = j;
          I_values(base+j) = KAT::one () + KAT::one ();
        }
        });
      local_graph_device_type I_G_lcl (I_lclColInds, I_rowOffsets);
      I_A_lcl = local_matrix_device_type("A_lcl", colMap->getLocalNumElements(),I_values,I_G_lcl);
    }
    crs_matrix_type input_matrix (I_A_lcl, rowMap, colMap, domMap, ranMap);

    // Output matrix: Dense all 2's except for the last row on the rank, which get cleaned up    
    local_matrix_device_type O_A_lcl;
    {
      row_offsets_type O_rowOffsets ("rowOffsets", lclNumRows+1);    
      lcl_col_inds_type O_lclColInds ("lclColInds", lclNumRows*lclNumRows);
      Kokkos::View<IST*, device_type> O_values ("values", lclNumRows*lclNumRows);
      Kokkos::parallel_for
        ("Output Matrix",
         range_type (0, lclNumRows),
         KOKKOS_LAMBDA (const LO lclRow) {
          O_rowOffsets(lclRow+1) = (lclRow+1)*lclNumRows;
          LO base = lclRow*lclNumRows;

          for(LO j=0; j<lclNumRows; j++) {
            O_lclColInds(base + j) = j;
            if(lclRow == lclNumRows-1)  {
              if(j==lclNumRows-1) O_values(base+j) = KAT::one();                  
              else O_values(base+j) = KAT::zero();
            }
            else {
              O_values(base+j) = KAT::one () + KAT::one ();
            }
          }            
        });
      local_graph_device_type O_G_lcl (O_lclColInds, O_rowOffsets);
      O_A_lcl = local_matrix_device_type("A_lcl", colMap->getLocalNumElements(),O_values,O_G_lcl);
    }
    crs_matrix_type output_matrix (O_A_lcl, rowMap, colMap, domMap, ranMap);

    // Generate a comparison vector
    output_matrix.apply(vec1, vec2_compare);      

    /*** Test device-filled boundary list, no execution space ***/
    {
      crs_matrix_type test_matrix = input_matrix;

      Tpetra::applyDirichletBoundaryConditionToLocalMatrixRows (test_matrix, lclRowInds);
      test_matrix.apply (vec1, vec2);
      vec2.update (-STS::one (), vec2_compare, STS::one ());
      
      const mag_type resNorm1 = vec2.norm1 ();
      TEST_EQUALITY( resNorm1, Teuchos::ScalarTraits<mag_type>::zero () );      
    }      

    /*** Test device-filled boundary list, provided execution space ***/
    {
      crs_matrix_type test_matrix = input_matrix;
      Tpetra::applyDirichletBoundaryConditionToLocalMatrixRows (execution_space(), test_matrix, lclRowInds);
      test_matrix.apply (vec1, vec2);
      vec2.update (-STS::one (), vec2_compare, STS::one ());
      
      const mag_type resNorm1 = vec2.norm1 ();
      TEST_EQUALITY( resNorm1, Teuchos::ScalarTraits<mag_type>::zero () );      
    }    

   
    /*** Test host-filled boundary list ***/
    {
      crs_matrix_type test_matrix = input_matrix;
      Tpetra::applyDirichletBoundaryConditionToLocalMatrixRows (test_matrix, lclRowInds_h);
      test_matrix.apply (vec1, vec2);
      vec2.update (-STS::one (), vec2_compare, STS::one ());
      
      const mag_type resNorm1 = vec2.norm1 ();
      TEST_EQUALITY( resNorm1, Teuchos::ScalarTraits<mag_type>::zero () );      
    }    

  }


   TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ApplyDirichletRowsAndCols, SC, LO, GO, NT )
  {
    using Tpetra::createContigMapWithNode;
    using Teuchos::RCP;
    using std::endl;
    using crs_matrix_type = Tpetra::CrsMatrix<SC,LO,GO,NT>;
    using map_type = Tpetra::Map<LO,GO,NT>;
    using vec_type = Tpetra::Vector<SC,LO,GO,NT>;
    using mag_type = typename vec_type::mag_type;
    using IST = typename vec_type::impl_scalar_type;    
    using GST = Tpetra::global_size_t;
    using STS = Teuchos::ScalarTraits<SC>;
    using KAT = Kokkos::ArithTraits<IST>;    
    using local_matrix_device_type = 
          typename crs_matrix_type::local_matrix_device_type;
    using local_graph_device_type = 
          typename local_matrix_device_type::staticcrsgraph_type;
    using device_type = typename crs_matrix_type::device_type;
    using execution_space = typename crs_matrix_type::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const auto comm = Tpetra::getDefaultComm ();
    const LO lclNumRows = 10;
    RCP<const map_type> rowMap =
      createContigMapWithNode<LO,GO,NT> (INVALID, lclNumRows,
                                         comm);
    RCP<const map_type> colMap = rowMap;
    RCP<const map_type> domMap = rowMap;
    RCP<const map_type> ranMap = rowMap;
    
    vec_type vec1 (domMap);
    vec_type vec2 (ranMap), vec2_compare(ranMap);

    vec1.putScalar (STS::one ());
    vec2.putScalar (STS::zero ());
    vec2_compare.putScalar (STS::zero ());

    using row_offsets_type =
      typename local_graph_device_type::row_map_type::non_const_type;
    using lcl_col_inds_type =
      typename local_graph_device_type::entries_type::non_const_type;

    // Device filled BCs
    Kokkos::View<typename crs_matrix_type::local_ordinal_type*, device_type> lclRowInds ("lclRowInds", 1);
    Kokkos::parallel_for
      ("Fill lclRowInds",
       range_type (0, 1),
       KOKKOS_LAMBDA (const LO lclRow) {
	  lclRowInds(lclRow) = lclNumRows-1;
      });
    
    // Host filled BCs
    Kokkos::View<typename crs_matrix_type::local_ordinal_type*, Kokkos::HostSpace> lclRowInds_h ("lclRowsInds_h",1);
    for (LO k = 0; k < 1; ++k) {
      lclRowInds_h[0] = lclNumRows-1;
    }

    // Input matrix: Dense all 2's.
    local_matrix_device_type I_A_lcl;
    {
      row_offsets_type I_rowOffsets ("rowOffsets", lclNumRows+1);    
      lcl_col_inds_type I_lclColInds ("lclColInds", lclNumRows*lclNumRows);
      Kokkos::View<IST*, device_type> I_values ("values", lclNumRows*lclNumRows);
      Kokkos::parallel_for
        ("Input Matrix",
         range_type (0, lclNumRows),
         KOKKOS_LAMBDA (const LO lclRow) {
          I_rowOffsets(lclRow+1) = (lclRow+1)*lclNumRows;
          LO base = lclRow*lclNumRows;
        for(LO j=0; j<lclNumRows; j++) {
          I_lclColInds(base + j) = j;
          I_values(base+j) = KAT::one () + KAT::one ();
        }
        });
      local_graph_device_type I_G_lcl (I_lclColInds, I_rowOffsets);
      I_A_lcl = local_matrix_device_type("A_lcl", colMap->getLocalNumElements(),I_values,I_G_lcl);
    }
    crs_matrix_type input_matrix (I_A_lcl, rowMap, colMap, domMap, ranMap);

    // Output matrix: Dense all 2's except for the last row on the rank, which get cleaned up    
    local_matrix_device_type O_A_lcl;
    {
      row_offsets_type O_rowOffsets ("rowOffsets", lclNumRows+1);    
      lcl_col_inds_type O_lclColInds ("lclColInds", lclNumRows*lclNumRows);
      Kokkos::View<IST*, device_type> O_values ("values", lclNumRows*lclNumRows);
      Kokkos::parallel_for
        ("Output Matrix",
         range_type (0, lclNumRows),
         KOKKOS_LAMBDA (const LO lclRow) {
          O_rowOffsets(lclRow+1) = (lclRow+1)*lclNumRows;
          LO base = lclRow*lclNumRows;

          for(LO j=0; j<lclNumRows; j++) {
            O_lclColInds(base + j) = j;
            if(lclRow == lclNumRows-1)  {
              if(j==lclNumRows-1) O_values(base+j) = KAT::one();
              else O_values(base+j) = KAT::zero();
            }
            else {
              if(j==lclNumRows-1) O_values(base+j) = KAT::zero();
              else O_values(base+j) = KAT::one () + KAT::one ();
            }
          }            
        });
      local_graph_device_type O_G_lcl (O_lclColInds, O_rowOffsets);
      O_A_lcl = local_matrix_device_type("A_lcl", colMap->getLocalNumElements(),O_values,O_G_lcl);
    }
    crs_matrix_type output_matrix (O_A_lcl, rowMap, colMap, domMap, ranMap);

    // Generate a comparison vector
    output_matrix.apply(vec1, vec2_compare);

    /*** Test device-filled boundary list, no execution space ***/
    {
      crs_matrix_type test_matrix = input_matrix;

      Tpetra::applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns (test_matrix, lclRowInds);
      test_matrix.apply (vec1, vec2);
      vec2.update (-STS::one (), vec2_compare, STS::one ());
      
      const mag_type resNorm1 = vec2.norm1 ();
      TEST_EQUALITY( resNorm1, Teuchos::ScalarTraits<mag_type>::zero () );      
    }      

    /*** Test device-filled boundary list, provided execution space ***/
    {
      crs_matrix_type test_matrix = input_matrix;
      Tpetra::applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns (execution_space(), test_matrix, lclRowInds);
      test_matrix.apply (vec1, vec2);
      vec2.update (-STS::one (), vec2_compare, STS::one ());
      
      const mag_type resNorm1 = vec2.norm1 ();
      TEST_EQUALITY( resNorm1, Teuchos::ScalarTraits<mag_type>::zero () );      
    }    

   
    /*** Test host-filled boundary list ***/
    {
      crs_matrix_type test_matrix = input_matrix;
      Tpetra::applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns (test_matrix, lclRowInds_h);
      test_matrix.apply (vec1, vec2);
      vec2.update (-STS::one (), vec2_compare, STS::one ());
      
      const mag_type resNorm1 = vec2.norm1 ();
      TEST_EQUALITY( resNorm1, Teuchos::ScalarTraits<mag_type>::zero () );      
    }    

  }


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ApplyDirichletRows, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ApplyDirichletRowsAndCols, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
