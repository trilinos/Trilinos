#ifndef INTREPID_F1_TET_I1_FEM_FIAT_HPP
#define INTREPID_F1_TET_I1_FEM_FIAT_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_Tabulate.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS_types.hpp"

namespace Intrepid {

/** \class Intrepid::Basis_F1_TET_I1_FEM_FIAT
  \brief Implementation of FIAT FEM basis functions of degree 1 for 1-forms on TET cells. 
  Reconstruction space type is INCOMPLETE. Definition of the DoF set 
  for this basis, its enumeration and the associated local DoF tags are as follows,

\verbatim
  =================================================================================================
  |        |                degree-of-freedom-tag              |                                  |
  | DoF Id |---------------------------------------------------|          DoF definition          |
  |        |  subc dim  |  subc id   | subc DoFId |subc num DoF|                                  |
  |========|============|============|============|============|==================================|
  |    0   |     1      |     0      |     0      |     1      |       L_0(u)=(u.t)(0.5,0.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    1   |     1      |     1      |     0      |     1      |       L_1(u)=(u.t)(0.5,0.5,0.0)  |
  |========|============|============|============|============|==================================|
  |    2   |     1      |     2      |     0      |     1      |       L_2(u)=(u.t)(0.0,0.5,0.0)  |
  |========|============|============|============|============|==================================|
  |    3   |     1      |     3      |     0      |     1      |       L_3(u)=(u.t)(0.0,0.0,0.5)  |
  |========|============|============|============|============|==================================|
  |    4   |     1      |     4      |     0      |     1      |       L_4(u)=(u.t)(0.5,0.0,0.5)  |
  |========|============|============|============|============|==================================|
  |    5   |     1      |     5      |     0      |     1      |       L_5(u)=(u.t)(0.0,0.5,0.5)  |
  |========|============|============|============|============|==================================|


\endverbatim

  The DefaultBasisFactory will select this basis
  if the following parameters are specified:
  \verbatim
  |=======================|===================================|
  |  EField               |  FIELD_FORM_1                     |
  |-----------------------|-----------------------------------|
  |  ECell                |  CELL_TET                         |
  |-----------------------|-----------------------------------|
  |  EReconstructionSpace |  RECONSTRUCTION_SPACE_INCOMPLETE  |
  |-----------------------|-----------------------------------|
  |  degree               |  1                                |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim

*/

template<class Scalar>
class Basis_F1_TET_I1_FEM_FIAT: public Basis<Scalar> {
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
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > vdm0_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > vdm1_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > vdm2_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats0_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats1_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats2_;

/* \brief Static data array that provides the Vandermonde and derivative matrices.  They need to
     be static data members because of the templating w.r.t type */
  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  static Scalar *get_vdm0_data () {
    static Scalar vdm0_data [] = { 2.142857142857144e-01 , -3.571428571428571e-01 , -3.571428571428569e-01 , 4.285714285714287e-01 , -4.285714285714285e-01 , 1.375581506340984e-16 , -3.078079127378572e-16 , -1.982518532415071e-16 , -1.982518532415071e-16 , 1.095560594963500e-16 , -1.095560594963501e-16 , 3.697785493223493e-32 , -4.999999999999998e-01 , -4.999999999999997e-01 , -4.999999999999997e-01 , -9.714451465470120e-17 , -2.775557561562891e-17 , 5.551115123125783e-17 , -2.857142857142857e-01 , 1.428571428571429e-01 , 1.428571428571428e-01 , 4.285714285714285e-01 , -4.285714285714285e-01 , 6.418476861114186e-17 };
    return vdm0_data ;
  }

  static Scalar *get_vdm1_data () {
    static Scalar vdm1_data [] = { 1.071428571428571e-01 , 1.547619047619048e-01 , -5.119047619047618e-01 , 3.809523809523809e-01 , -4.761904761904757e-02 , -3.333333333333332e-01 , 2.500000000000002e-01 , 2.500000000000002e-01 , 2.500000000000002e-01 , 7.632783294297951e-17 , 1.387778780781446e-17 , -5.551115123125783e-17 , -2.499999999999999e-01 , -2.499999999999999e-01 , -2.499999999999998e-01 , -4.857225732735060e-17 , -1.387778780781446e-17 , 2.775557561562891e-17 , -1.428571428571428e-01 , -9.523809523809523e-02 , 2.380952380952381e-01 , 3.809523809523809e-01 , -4.761904761904758e-02 , -3.333333333333333e-01 };
    return vdm1_data ;
  }

  static Scalar *get_vdm2_data () {
    static Scalar vdm2_data [] = { 1.071428571428570e-01 , -6.746031746031750e-02 , -2.896825396825395e-01 , 6.031746031746034e-01 , 1.746031746031745e-01 , 2.222222222222222e-01 , 1.785714285714286e-01 , 3.571428571428563e-02 , 3.571428571428562e-02 , -1.428571428571429e-01 , 1.428571428571429e-01 , 2.255140518769849e-17 , -2.499999999999998e-01 , -1.388888888888888e-01 , -3.611111111111109e-01 , -1.111111111111112e-01 , -1.111111111111111e-01 , 2.222222222222222e-01 , -1.428571428571428e-01 , 1.587301587301589e-02 , 1.269841269841269e-01 , 2.698412698412697e-01 , -1.587301587301587e-01 , -1.111111111111110e-01 };
    return vdm2_data ;
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
  Basis_F1_TET_I1_FEM_FIAT() : numDof_(6), isSet_(false) , vdm0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TET_I1_FEM_FIAT::get_vdm0_data(),6,6,4) ) ),vdm1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TET_I1_FEM_FIAT::get_vdm1_data(),6,6,4) ) ),vdm2_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TET_I1_FEM_FIAT::get_vdm2_data(),6,6,4) ) ),dmats0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TET_I1_FEM_FIAT::get_dmats0_data(),4,4,4) ) ),dmats1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TET_I1_FEM_FIAT::get_dmats1_data(),4,4,4) ) ),dmats2_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TET_I1_FEM_FIAT::get_dmats2_data(),4,4,4) ) ) { }
  
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

#include "Intrepid_F1_TET_I1_FEM_FIATDef.hpp"

#endif

