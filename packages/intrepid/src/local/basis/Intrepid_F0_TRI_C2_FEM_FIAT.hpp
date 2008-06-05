#ifndef INTREPID_F0_TRI_C2_FEM_FIAT_HPP
#define INTREPID_F0_TRI_C2_FEM_FIAT_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_Tabulate.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS_types.hpp"

namespace Intrepid {

/** \class Intrepid::Basis_F0_TRI_C2_FEM_FIAT
  \brief Implementation of FIAT FEM basis functions of degree 2 for 0-forms on TRI cells. 
  Reconstruction space type is COMPLETE. Definition of the DoF set 
  for this basis, its enumeration and the associated local DoF tags are as follows,

\verbatim
  =================================================================================================
  |        |                degree-of-freedom-tag              |                                  |
  | DoF Id |---------------------------------------------------|          DoF definition          |
  |        |  subc dim  |  subc id   | subc DoFId |subc num DoF|                                  |
  |========|============|============|============|============|==================================|
  |    0   |     0      |     0      |     0      |     1      |               L_0(u)=u(0.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    1   |     0      |     1      |     0      |     1      |               L_1(u)=u(1.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    2   |     0      |     2      |     0      |     1      |               L_2(u)=u(0.0,1.0)  |
  |========|============|============|============|============|==================================|
  |    3   |     1      |     0      |     0      |     1      |               L_3(u)=u(0.5,0.0)  |
  |========|============|============|============|============|==================================|
  |    4   |     1      |     1      |     0      |     1      |               L_4(u)=u(0.5,0.5)  |
  |========|============|============|============|============|==================================|
  |    5   |     1      |     2      |     0      |     1      |               L_5(u)=u(0.0,0.5)  |
  |========|============|============|============|============|==================================|


\endverbatim

  The DefaultBasisFactory will select this basis
  if the following parameters are specified:
  \verbatim
  |=======================|===================================|
  |  EField               |  FIELD_FORM_0                     |
  |-----------------------|-----------------------------------|
  |  ECell                |  CELL_TRI                   |
  |-----------------------|-----------------------------------|
  |  EReconstructionSpace |  RECONSTRUCTION_SPACE_COMPLETE    |
  |-----------------------|-----------------------------------|
  |  degree               |  2                          |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim

*/

template<class Scalar>
class Basis_F0_TRI_C2_FEM_FIAT: public Basis<Scalar> {
  private:
  
  /** \brief Dimension of the space spanned by the basis = number of degrees of freedom.
  */
  int numDof_;
  
  /**\brief Lookup table for the DoF's local enumeration (DoF Id) by its local DoF tag
  */
  Teuchos::Array<Teuchos::Array<Teuchos::Array<int> > > tagToEnum_;
  
  /**\brief Lookup table for the DoF's local DoF tag by its local enumeration (DoF Id)
  */
  Teuchos::Array<LocalDofTag> enumToTag_;
  
  /**\brief "true" if both lookup arrays have been set by initialize()
  */
  bool isSet_;

  /**\brief coefficients of nodal basis in terms of orthogonal polynomials */
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > vdm_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats0_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats1_;

/* \brief Static data array that provides the Vandermonde and derivative matrices.  They need to
     be static data members because of the templating w.r.t type */
    static Scalar *get_vdm_data() {
    static Scalar vdm_data[] = { 1.1102230246251565e-16 , 0.0000000000000000e+00 , -2.7755575615628914e-17 , 3.3333333333333331e-01 , 3.3333333333333337e-01 , 3.3333333333333337e-01 , -2.9999999999999999e-01 , 2.9999999999999999e-01 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 4.0000000000000002e-01 , -4.0000000000000002e-01 , -9.9999999999999992e-02 , -9.9999999999999992e-02 , 1.9999999999999998e-01 , -2.6666666666666666e-01 , 1.3333333333333333e-01 , 1.3333333333333333e-01 , 3.3333333333333331e-01 , 3.3333333333333331e-01 , 0.0000000000000000e+00 , -6.6666666666666663e-01 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 2.0000000000000001e-01 , -2.0000000000000001e-01 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 4.0000000000000002e-01 , -4.0000000000000002e-01 , 6.6666666666666680e-02 , 6.6666666666666680e-02 , 2.0000000000000001e-01 , 6.6666666666666666e-02 , -2.0000000000000001e-01 , -2.0000000000000001e-01 };
    return vdm_data;
  }

  static Scalar *get_dmats0_data() {
    static Scalar dmats0_data[] = { 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 1.9999999999999996e+00 , 0.0000000000000000e+00 , -4.4408920985006262e-16 , -8.8817841970012523e-16 , -2.6645352591003757e-15 , 4.4408920985006262e-16 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 3.6637359812630144e-16 , 6.0000000000000018e+00 , 2.2204460492503126e-16 , 8.8817841970012523e-16 , -1.3322676295501878e-15 , 2.2204460492503091e-16 , 1.3333333333333330e+00 , 0.0000000000000000e+00 , 3.3333333333333335e+00 , 0.0000000000000000e+00 , -1.7763568394002505e-15 , 8.8817841970012523e-16 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 };
    return dmats0_data;
  }

  static Scalar *get_dmats1_data() {
    static Scalar dmats1_data[] = { 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 9.9999999999999978e-01 , 0.0000000000000000e+00 , -2.2204460492503131e-16 , -4.4408920985006262e-16 , -1.3322676295501878e-15 , 2.2204460492503131e-16 , 2.9999999999999996e+00 , 0.0000000000000000e+00 , -8.8817841970012523e-16 , -1.7763568394002505e-15 , -3.5527136788005009e-15 , -4.4408920985006262e-16 , 6.6666666666666652e-01 , 3.0000000000000000e+00 , -3.3333333333333337e-01 , 0.0000000000000000e+00 , -4.4408920985006262e-16 , 3.3306690738754696e-16 , 6.6666666666666674e-01 , 5.0000000000000000e+00 , 1.6666666666666667e+00 , 8.8817841970012523e-16 , -8.8817841970012523e-16 , 4.4408920985006262e-16 , -1.3333333333333330e+00 , 0.0000000000000000e+00 , 6.6666666666666679e+00 , 1.7763568394002505e-15 , 8.8817841970012484e-16 , -1.7763568394002505e-15 };
    return dmats1_data;
  }


  public:

  /** \brief Constructor.
  */
  Basis_F0_TRI_C2_FEM_FIAT() : numDof_(6), isSet_(false) , vdm_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C2_FEM_FIAT::get_vdm_data(),6,6,6) ) ),dmats0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C2_FEM_FIAT::get_dmats0_data(),6,6,6) ) ),dmats1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C2_FEM_FIAT::get_dmats1_data(),6,6,6) ) ) { }
  
  /** \brief Initializes arrays needed for the lookup of the local enumeration (DoF Id) of a 
    degree-of-freedom by its local DoF tag and the reverse lookup of the DoF tag by DoF Id.
   */
  void initialize();
  
  
  /** \brief Returns FieldContainer with multi-indexed quantity with the values of an operator 
  applied to a set of FEM basis functions, evaluated at a set of
  points in a <strong>reference</strong> cell. The rank of the return
  FieldContainer argument depends on the number of points, number of
  basis, functions, their rank, and the rank of the operator applied to them; see 
  FieldContainer::resize(int,int,EField,EOperator,int) for summary of
  the admissible index combinations.
  In particular, the admissible range of the <var>operatorType</var> argument depends on the rank of 
  the basis functions and the space dimension.Because FEM
  reconstruction relies on COMPLETE or INCOMPLETE polynomials, i.e.,
  smooth local spaces, the admissible range of <var>operatorType</var>
  always includes VALUE, D0, D1,...,D10. If derivative order exceeds
  polynomial degree, output container is filled with 0. 
 
      \param outputValues   [out]         - FieldContainer with the computed values (see
                                            implementation for index ordering)
      \param inputPoints     [in]         - evaluation points on the reference cell  
      \param operatorType    [in]         - the operator being applied to the basis function
  */    

  void getValues(FieldContainer<Scalar>&                  outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const EOperator                        operatorType) const;
  
  
  /** \brief This method is intended for FVD reconstructions and should not be used here. Its 
    invocation will throw an exception. 
  */  
  void getValues(FieldContainer<Scalar>&                  outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const Cell<Scalar>&                    cell) const;

  
  int getNumLocalDof() const;
  
  
  int getLocalDofEnumeration(const LocalDofTag dofTag);

  
  LocalDofTag getLocalDofTag(int id);

  
  const Teuchos::Array<LocalDofTag> & getAllLocalDofTags();

  
  ECell getCellType() const;

  
  EBasis getBasisType() const;

  
  ECoordinates getCoordinateSystem() const;

  
  int getDegree() const;

};

} // namespace Intrepid

#include "Intrepid_F0_TRI_C2_FEM_FIATDef.hpp"

#endif

