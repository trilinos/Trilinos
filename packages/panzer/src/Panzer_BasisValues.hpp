#ifndef PANZER_BASIS_VALUES_DECL_HPP
#define PANZER_BASIS_VALUES_DECL_HPP

#include "Teuchos_RCP.hpp"

#include "Intrepid_Basis.hpp"
#include "Panzer_BasisIRLayout.hpp"

namespace panzer {

  template<typename Scalar,typename Array>
  struct BasisValues { 
    static const Array dummyArray;    
 
    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<panzer::BasisIRLayout>& basis);

    //! Sizes/allocates memory for arrays
    template <typename ArrayFactory>
    void setupArrays(const Teuchos::RCP<panzer::BasisIRLayout>& basis,
                     const ArrayFactory & af);

    void evaluateValues(const Array& cub_points,
			const Array& jac,
			const Array& jac_det,
			const Array& jac_inv,
			const Array& node_coordinates);

    void evaluateValues(const Array& cub_points,
			const Array& jac,
			const Array& jac_det,
			const Array& jac_inv,
			const Array& weighted_measure,
			const Array& node_coordinates);

    PureBasis::EElementSpace getElementSpace() const; 

    Array basis_ref;           // <BASIS,IP>
    Array basis;               // <Cell,BASIS,IP>
    Array grad_basis_ref;      // <BASIS,IP,Dim>
    Array grad_basis;          // <Cell,BASIS,IP,Dim>
    Array curl_basis_ref;      // <BASIS,IP,Dim> (dimension dependent)
    Array curl_basis;          // <Cell,BASIS,IP,Dim> (dimension dependent)
    Array weighted_basis;      // <Cell,BASIS,IP>
    Array weighted_grad_basis; // <Cell,BASIS,IP,Dim>
    Array weighted_curl_basis; // <Cell,BASIS,IP,Dim> (dimension dependent)

    /** Carterisan coordinates for basis coefficients

        NOTE: This quantity is not always available.  Certain bases
        may not have a corresponding coordiante value
    */
    Array basis_coordinates_ref;     // <Cell,BASIS>

    /** Carterisan coordinates for basis coefficients

        NOTE: This quantity is not always available.  Certain bases
        may not have a corresponding coordiante value
    */
    Array basis_coordinates;         // <Cell,BASIS,Dim>

    Teuchos::RCP<panzer::BasisIRLayout> basis_layout;
    
    Teuchos::RCP<Intrepid::Basis<double,Array> > intrepid_basis;
  };

  
  /** Implementation for intrepid field container factory. This
    * is intended to be used only with the BasisValues object.
    */
  template <typename Scalar>
  class IntrepidFieldContainerFactory {
  public:
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0) const;
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1) const;
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2) const;
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3) const;
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const;
  };

} // namespace panzer

#include "Panzer_BasisValues_impl.hpp"

#endif
