#ifndef INTREPID_F0_TET_C1_FEM_FIAT_HPP
#define INTREPID_F0_TET_C1_FEM_FIAT_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_Tabulate.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS_types.hpp"

namespace Intrepid {

/** \class Intrepid::Basis_F0_TET_C1_FEM_FIAT
  \brief Implementation of FIAT FEM basis functions of degree 1 for 0-forms on TET cells. 
  Reconstruction space type is COMPLETE. Definition of the DoF set 
  for this basis, its enumeration and the associated local DoF tags are as follows,

\verbatim
  =================================================================================================
  |        |                degree-of-freedom-tag              |                                  |
  | DoF Id |---------------------------------------------------|          DoF definition          |
  |        |  subc dim  |  subc id   | subc DoFId |subc num DoF|                                  |
  |========|============|============|============|============|==================================|
  |    0   |     0      |     0      |     0      |     1      |           L_0(u)=u(0.0,0.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    1   |     0      |     1      |     0      |     1      |           L_1(u)=u(1.0,0.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    2   |     0      |     2      |     0      |     1      |           L_2(u)=u(0.0,1.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    3   |     0      |     3      |     0      |     1      |           L_3(u)=u(0.0,0.0,1.0)  |
  |========|============|============|============|============|==================================|


\endverbatim

  The DefaultBasisFactory will select this basis
  if the following parameters are specified:
  \verbatim
  |=======================|===================================|
  |  EField               |  FIELD_FORM_0                     |
  |-----------------------|-----------------------------------|
  |  ECell                |  CELL_TET                   |
  |-----------------------|-----------------------------------|
  |  EReconstructionSpace |  RECONSTRUCTION_SPACE_COMPLETE    |
  |-----------------------|-----------------------------------|
  |  degree               |  1                          |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim

*/

template<class Scalar>
class Basis_F0_TET_C1_FEM_FIAT: public Basis<Scalar> {
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
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats2_;

/* \brief Static data array that provides the Vandermonde and derivative matrices.  They need to
     be static data members because of the templating w.r.t type */
  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  static Scalar *get_vdm_data() {
    static Scalar vdm_data[] = { 2.500000000000000e-01 , 2.500000000000000e-01 , 2.500000000000000e-01 , 2.500000000000000e-01 , -5.000000000000000e-01 , 5.000000000000000e-01 , 0.000000000000000e+00 , 0.000000000000000e+00 , -1.666666666666667e-01 , -1.666666666666667e-01 , 3.333333333333333e-01 , 0.000000000000000e+00 , -8.333333333333334e-02 , -8.333333333333334e-02 , -8.333333333333333e-02 , 2.500000000000000e-01 };
    return vdm_data;
  }

  static Scalar *get_dmats0_data() {
    static Scalar dmats0_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 5.656854249492381e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 };
    return dmats0_data;
  }

  static Scalar *get_dmats1_data() {
    static Scalar dmats1_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 2.828427124746190e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 8.485281374238571e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 };
    return dmats1_data;
  }

  static Scalar *get_dmats2_data() {
    static Scalar dmats2_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 2.828427124746190e+00 , -8.881784197001252e-16 , 0.000000000000000e+00 , 0.000000000000000e+00 , 2.828427124746190e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , -4.440892098500626e-16 , 1.131370849898476e+01 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 };
    return dmats2_data;
  }
#endif


  public:

  /** \brief Constructor.
  */
  Basis_F0_TET_C1_FEM_FIAT() : numDof_(4), isSet_(false) , vdm_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TET_C1_FEM_FIAT::get_vdm_data(),4,4,4) ) ),dmats0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TET_C1_FEM_FIAT::get_dmats0_data(),4,4,4) ) ),dmats1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TET_C1_FEM_FIAT::get_dmats1_data(),4,4,4) ) ),dmats2_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TET_C1_FEM_FIAT::get_dmats2_data(),4,4,4) ) ) { }
  
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

#include "Intrepid_F0_TET_C1_FEM_FIATDef.hpp"

#endif

