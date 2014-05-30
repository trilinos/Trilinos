#ifndef __Panzer_Intrepid_ConstBasis_hpp__
#define __Panzer_Intrepid_ConstBasis_hpp__

#include "Intrepid_Basis.hpp"

namespace panzer {
  
template<class Scalar, class ArrayScalar> 
class Basis_Constant: public Intrepid::Basis<Scalar, ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:
  
  /** \brief  Constructor.
   */
  Basis_Constant(const shards::CellTopology & ct);  
  
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Triangle</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Triangle</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec .

      \param  outputValues      [out] - variable rank array with the basis values
      \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
      \param  operatorType      [in]  - the operator acting on the basis functions    
   */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const Intrepid::EOperator        operatorType) const;
  
  
  /**  \brief  FVD basis evaluation: invocation of this method throws an exception.
   */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const ArrayScalar &    cellVertices,
                 const Intrepid::EOperator        operatorType = Intrepid::OPERATOR_VALUE) const;

};

}// namespace panzer

#include "Panzer_Intrepid_ConstBasis_impl.hpp"

#endif
