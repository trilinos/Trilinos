#ifndef INTREPID_F1_TRI_I2_FEM_FIAT_HPP
#define INTREPID_F1_TRI_I2_FEM_FIAT_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_Tabulate.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS_types.hpp"

namespace Intrepid {

/** \class Intrepid::Basis_F1_TRI_I2_FEM_FIAT
  \brief Implementation of FIAT FEM basis functions of degree 2 for 1-forms on TRI cells. 
  Reconstruction space type is INCOMPLETE. Definition of the DoF set 
  for this basis, its enumeration and the associated local DoF tags are as follows,

\verbatim
  =================================================================================================
  |        |                degree-of-freedom-tag              |                                  |
  | DoF Id |---------------------------------------------------|          DoF definition          |
  |        |  subc dim  |  subc id   | subc DoFId |subc num DoF|                                  |
  |========|============|============|============|============|==================================|
  |    0   |     1      |     0      |     0      |     2      |L_0(u)=(u.t)(0.333333333333,0.0)  |
  |========|============|============|============|============|==================================|
  |    1   |     1      |     0      |     1      |     2      |L_1(u)=(u.t)(0.666666666667,0.0)  |
  |========|============|============|============|============|==================================|
  |    2   |     1      |     1      |     0      |     2      |L_2(u)=(u.t)(0.666666666667,0.333333333333)  |
  |========|============|============|============|============|==================================|
  |    3   |     1      |     1      |     1      |     2      |L_3(u)=(u.t)(0.333333333333,0.666666666667)  |
  |========|============|============|============|============|==================================|
  |    4   |     1      |     2      |     0      |     2      |L_4(u)=(u.t)(0.0,0.666666666667)  |
  |========|============|============|============|============|==================================|
  |    5   |     1      |     2      |     1      |     2      |L_5(u)=(u.t)(0.0,0.333333333333)  |
  |========|============|============|============|============|==================================|
  |    6   |     2      |     0      |     0      |     2      |           L_6(u)=IntegralMoment  |
  |========|============|============|============|============|==================================|
  |    7   |     2      |     0      |     1      |     2      |           L_7(u)=IntegralMoment  |
  |========|============|============|============|============|==================================|


\endverbatim

  The DefaultBasisFactory will select this basis
  if the following parameters are specified:
  \verbatim
  |=======================|===================================|
  |  EField               |  FIELD_FORM_1                     |
  |-----------------------|-----------------------------------|
  |  ECell                |  CELL_TRI                         |
  |-----------------------|-----------------------------------|
  |  EReconstructionSpace |  RECONSTRUCTION_SPACE_INCOMPLETE  |
  |-----------------------|-----------------------------------|
  |  degree               |  2                                |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim

*/

template<class Scalar>
class Basis_F1_TRI_I2_FEM_FIAT: public Basis<Scalar> {
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
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats0_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats1_;

/* \brief Static data array that provides the Vandermonde and derivative matrices.  They need to
     be static data members because of the templating w.r.t type */
  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  static Scalar *get_vdm0_data () {
    static Scalar vdm0_data [] = { 1.130958391173942e-16 , 3.630722025110833e-17 , -1.155352940054866e-16 , -3.268880137955637e-17 , 1.965371606305337e-16 , 7.325140927855189e-17 , 2.000000000000001e+00 , 4.117799851099946e-16 , -9.193548387096774e-01 , 9.193548387096773e-01 , -5.806451612903225e-01 , -1.920744616688060e-16 , -1.123263096848136e-16 , 5.806451612903227e-01 , -1.161290322580647e+00 , 2.322580645161291e+00 , -2.157894736842104e-01 , -2.157894736842105e-01 , 2.842105263157896e-01 , -5.684210526315793e-01 , -5.684210526315788e-01 , 2.842105263157894e-01 , 2.947368421052635e-01 , -1.025113154085044e-16 , -1.973227633212326e-15 , 1.992530516116116e-15 , 1.786674847982189e-15 , 7.149586169314260e-17 , 1.098436094932587e-16 , -1.543665843278697e-15 , 3.635556958506501e-15 , -6.322531556630950e-15 , 5.806451612903227e-01 , -5.806451612903227e-01 , -5.806451612903225e-01 , -1.752070710736575e-16 , -2.949029909160572e-16 , 5.806451612903227e-01 , -1.161290322580647e+00 , 2.322580645161291e+00 , 2.842105263157887e-01 , 2.842105263157901e-01 , 2.842105263157903e-01 , -5.684210526315792e-01 , -5.684210526315788e-01 , 2.842105263157889e-01 , -1.705263157894737e+00 , -2.220446049250313e-15 };
    return vdm0_data ;
  }

  static Scalar *get_vdm1_data () {
    static Scalar vdm1_data [] = { -4.987329993433320e-17 , 2.005774019098183e-16 , 3.599551212651875e-17 , -2.585483368199026e-17 , 1.045552059102291e-16 , -4.542807102714264e-17 , 1.487525380650112e-16 , 2.000000000000001e+00 , -5.149405772495756e-01 , 4.044142614601020e-01 , 4.044142614601018e-01 , 1.105263157894738e-01 , 1.105263157894737e-01 , 9.850594227504245e-01 , 1.250933786078100e+00 , 1.161290322580645e+00 , 4.533106960950760e-02 , -2.611205432937182e-01 , -2.611205432937179e-01 , 2.157894736842105e-01 , -7.842105263157896e-01 , 5.453310696095072e-01 , 3.409168081494067e-01 , -3.870967741935500e-01 , -5.806451612903222e-01 , 5.806451612903221e-01 , 5.806451612903220e-01 , 1.361757928641794e-16 , 2.437286483747414e-16 , -5.806451612903221e-01 , 1.161290322580645e+00 , -2.322580645161289e+00 , 1.955857385398981e-01 , -3.850594227504243e-01 , -3.850594227504242e-01 , 1.894736842105262e-01 , 1.894736842105262e-01 , 1.955857385398979e-01 , -1.222410865874335e-02 , 1.161290322580644e+00 , -3.056027164685615e-03 , 2.872665534804749e-01 , 2.872665534804748e-01 , -2.842105263157894e-01 , -2.842105263157893e-01 , -3.056027164685407e-03 , -5.623089983022083e-01 , -5.806451612903203e-01 };
    return vdm1_data ;
  }

  static Scalar *get_dmats0_data() {
    static Scalar dmats0_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 2.000000000000000e+00 , 0.000000000000000e+00 , -4.440892098500626e-16 , -8.881784197001252e-16 , -2.664535259100376e-15 , 4.440892098500626e-16 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 3.663735981263014e-16 , 6.000000000000002e+00 , 2.220446049250313e-16 , 8.881784197001252e-16 , -1.332267629550188e-15 , 2.220446049250309e-16 , 1.333333333333333e+00 , 0.000000000000000e+00 , 3.333333333333333e+00 , 0.000000000000000e+00 , -1.776356839400250e-15 , 8.881784197001252e-16 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 };
    return dmats0_data;
  }

  static Scalar *get_dmats1_data() {
    static Scalar dmats1_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 9.999999999999998e-01 , 0.000000000000000e+00 , -2.220446049250313e-16 , -4.440892098500626e-16 , -1.332267629550188e-15 , 2.220446049250313e-16 , 3.000000000000000e+00 , 0.000000000000000e+00 , -8.881784197001252e-16 , -1.776356839400250e-15 , -3.552713678800501e-15 , -4.440892098500626e-16 , 6.666666666666665e-01 , 3.000000000000000e+00 , -3.333333333333334e-01 , 0.000000000000000e+00 , -4.440892098500626e-16 , 3.330669073875470e-16 , 6.666666666666667e-01 , 5.000000000000000e+00 , 1.666666666666667e+00 , 8.881784197001252e-16 , -8.881784197001252e-16 , 4.440892098500626e-16 , -1.333333333333333e+00 , 0.000000000000000e+00 , 6.666666666666668e+00 , 1.776356839400250e-15 , 8.881784197001248e-16 , -1.776356839400250e-15 };
    return dmats1_data;
  }
#endif


  public:

  /** \brief Constructor.
  */
  Basis_F1_TRI_I2_FEM_FIAT() : numDof_(8), isSet_(false) , vdm0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TRI_I2_FEM_FIAT::get_vdm0_data(),8,8,6) ) ),vdm1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TRI_I2_FEM_FIAT::get_vdm1_data(),8,8,6) ) ),dmats0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TRI_I2_FEM_FIAT::get_dmats0_data(),6,6,6) ) ),dmats1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F1_TRI_I2_FEM_FIAT::get_dmats1_data(),6,6,6) ) ) { }
  
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

#include "Intrepid_F1_TRI_I2_FEM_FIATDef.hpp"

#endif

