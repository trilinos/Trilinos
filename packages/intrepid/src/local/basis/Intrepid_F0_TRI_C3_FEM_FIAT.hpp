#ifndef INTREPID_F0_TRI_C3_FEM_FIAT_HPP
#define INTREPID_F0_TRI_C3_FEM_FIAT_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_Tabulate.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS_types.hpp"

namespace Intrepid {

/** \class Intrepid::Basis_F0_TRI_C3_FEM_FIAT
  \brief Implementation of FIAT FEM basis functions of degree 3 for 0-forms on TRI cells. 
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
  |    3   |     1      |     0      |     0      |     2      |    L_3(u)=u(0.333333333333,0.0)  |
  |========|============|============|============|============|==================================|
  |    4   |     1      |     0      |     1      |     2      |    L_4(u)=u(0.666666666667,0.0)  |
  |========|============|============|============|============|==================================|
  |    5   |     1      |     1      |     0      |     2      |L_5(u)=u(0.666666666667,0.333333333333)  |
  |========|============|============|============|============|==================================|
  |    6   |     1      |     1      |     1      |     2      |L_6(u)=u(0.333333333333,0.666666666667)  |
  |========|============|============|============|============|==================================|
  |    7   |     1      |     2      |     0      |     2      |    L_7(u)=u(0.0,0.666666666667)  |
  |========|============|============|============|============|==================================|
  |    8   |     1      |     2      |     1      |     2      |    L_8(u)=u(0.0,0.333333333333)  |
  |========|============|============|============|============|==================================|
  |    9   |     2      |     0      |     0      |     1      |L_9(u)=u(0.333333333333,0.333333333333)  |
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
  |  degree               |  3                          |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim

*/

template<class Scalar>
class Basis_F0_TRI_C3_FEM_FIAT: public Basis<Scalar> {
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
  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  static Scalar *get_vdm_data() {
    static Scalar vdm_data[] = { 3.333333333333331e-02 , 3.333333333333335e-02 , 3.333333333333335e-02 , 7.499999999999996e-02 , 7.499999999999996e-02 , 7.500000000000001e-02 , 7.500000000000004e-02 , 7.500000000000001e-02 , 7.499999999999990e-02 , 4.500000000000000e-01 , -4.999999999999993e-02 , 5.000000000000010e-02 , 1.427429603089487e-17 , -4.500000000000001e-01 , 4.500000000000000e-01 , 4.500000000000000e-01 , 0.000000000000000e+00 , 0.000000000000000e+00 , -4.499999999999999e-01 , 3.568574007723717e-17 , -1.666666666666668e-02 , -1.666666666666666e-02 , 3.333333333333333e-02 , -1.499999999999999e-01 , -1.500000000000001e-01 , -1.500000000000001e-01 , 3.000000000000000e-01 , 3.000000000000000e-01 , -1.500000000000000e-01 , 0.000000000000000e+00 , 2.142857142857143e-01 , 2.142857142857142e-01 , -8.489813632400870e-34 , -2.142857142857142e-01 , -2.142857142857144e-01 , 3.214285714285714e-01 , -7.646944302265103e-18 , 7.646944302265107e-18 , 3.214285714285714e-01 , -6.428571428571428e-01 , 1.285714285714286e-01 , -1.285714285714286e-01 , 4.996003610813204e-17 , 1.928571428571429e-01 , -1.928571428571429e-01 , 1.285714285714285e-01 , 3.214285714285714e-01 , -3.214285714285715e-01 , -1.285714285714285e-01 , 0.000000000000000e+00 , 4.285714285714286e-02 , 4.285714285714284e-02 , 1.285714285714286e-01 , 1.178571428571428e-01 , 1.178571428571429e-01 , -9.642857142857147e-02 , 3.214285714285713e-02 , 3.214285714285713e-02 , -9.642857142857143e-02 , -3.214285714285713e-01 , -2.250000000000001e-01 , 2.249999999999999e-01 , 2.674291294206276e-33 , 6.750000000000000e-01 , -6.749999999999999e-01 , -2.408787455213509e-17 , 2.408787455213509e-17 , -2.408787455213510e-17 , 2.408787455213509e-17 , -2.674291294206276e-33 , -1.607142857142857e-01 , -1.607142857142858e-01 , 6.367360224300658e-34 , 1.607142857142858e-01 , 1.607142857142857e-01 , 3.214285714285714e-01 , 5.735208226698832e-18 , -5.735208226698834e-18 , 3.214285714285714e-01 , -6.428571428571428e-01 , -9.642857142857141e-02 , 9.642857142857143e-02 , 3.568574007723717e-17 , -3.214285714285708e-02 , 3.214285714285710e-02 , -3.214285714285714e-01 , 3.214285714285714e-01 , -3.214285714285715e-01 , 3.214285714285714e-01 , -3.568574007723717e-17 , -3.214285714285713e-02 , -3.214285714285714e-02 , 1.285714285714286e-01 , -3.214285714285714e-02 , -3.214285714285712e-02 , 1.285714285714286e-01 , -1.928571428571429e-01 , -1.928571428571428e-01 , 1.285714285714286e-01 , 1.285714285714285e-01 };
    return vdm_data;
  }

  static Scalar *get_dmats0_data() {
    static Scalar dmats0_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 2.000000000000003e+00 , 2.601093943407511e-15 , 1.332267629550188e-15 , -3.357858209172413e-16 , -8.881784197001252e-16 , 4.440892098500626e-16 , 3.332182674458116e-31 , 9.044919008782904e-16 , 3.806478941571979e-16 , 4.440892098500626e-16 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 2.803313137178519e-16 , 6.000000000000002e+00 , -2.287059430727823e-15 , 2.808127879887354e-15 , -1.110223024625157e-16 , 8.691460249922652e-16 , 9.827852817271120e-16 , -6.647007207688176e-15 , 1.332267629550189e-15 , -1.795389234108110e-15 , 1.333333333333335e+00 , -4.440892098500591e-16 , 3.333333333333334e+00 , 1.160069772669550e-16 , -8.881784197001252e-16 , -8.881784197001252e-16 , -1.110223024625148e-16 , -6.960418636017305e-16 , -2.664535259100372e-15 , -2.664535259100376e-15 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 9.999999999999998e-01 , 4.988073446351596e-15 , -7.999999999999989e-01 , 9.999999999999998e+00 , -3.497202527569243e-15 , 2.000000000000002e-01 , 3.354459567260294e-15 , 2.439771740645651e-15 , 7.375052949295683e-15 , 8.326672684688674e-16 , -1.206442353426004e-15 , 2.400000000000002e+00 , -2.205643075588645e-15 , 2.246360693830211e-15 , 8.399999999999999e+00 , 8.183929724379720e-16 , 4.785457744357509e-16 , -8.179823032069250e-15 , 2.789186772031433e-30 , -2.512276101437497e-15 , 1.000000000000000e+00 , -2.775557561562887e-15 , 3.199999999999998e+00 , 1.885113380588022e-16 , 1.110223024625157e-16 , 4.199999999999999e+00 , -4.440892098500619e-16 , -1.797201843127905e-15 , -4.218847493575590e-15 , -2.664535259100376e-15 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 };
    return dmats0_data;
  }

  static Scalar *get_dmats1_data() {
    static Scalar dmats1_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 1.000000000000001e+00 , 1.300546971703756e-15 , 6.661338147750939e-16 , -1.678929104586207e-16 , -4.440892098500626e-16 , 2.220446049250313e-16 , 1.666091337229058e-31 , 4.522459504391452e-16 , 1.903239470785990e-16 , 2.220446049250313e-16 , 3.000000000000004e+00 , 3.901640915111266e-15 , 1.776356839400250e-15 , -6.702121850696354e-16 , -1.332267629550188e-15 , 8.881784197001252e-16 , -1.776356839400250e-15 , 2.447022176724795e-17 , -7.612957883143911e-16 , -8.881784197001252e-16 , 6.666666666666679e-01 , 3.000000000000002e+00 , -3.333333333333336e-01 , 7.252418620300098e-16 , 1.110223024625157e-16 , 5.551115123125783e-16 , -3.967857788365693e-16 , -2.803284815162586e-15 , 9.516197353929919e-16 , -5.551115123125783e-17 , 6.666666666666674e-01 , 5.000000000000003e+00 , 1.666666666666664e+00 , 2.509132357668788e-15 , -6.661338147750939e-16 , 0.000000000000000e+00 , -6.919068492753150e-17 , -6.035223341491033e-15 , -8.881784197001221e-16 , -2.664535259100376e-15 , -1.333333333333336e+00 , -6.534455516365203e-15 , 6.666666666666663e+00 , 1.902786318531034e-15 , 1.332267629550188e-15 , -2.664535259100376e-15 , 5.598900026800757e-31 , -2.312889109259917e-15 , -4.314009467114889e-15 , -4.440892098500626e-15 , 4.999999999999999e-01 , 2.400000000000003e+00 , -4.000000000000002e-01 , 5.000000000000001e+00 , -6.000000000000018e-01 , 1.000000000000004e-01 , 3.174464659392941e-15 , 7.793793954884532e-16 , 4.464682591885451e-15 , -5.551115123125783e-17 , 4.999999999999992e-01 , 1.200000000000005e+00 , 7.999999999999996e-01 , 6.999999999999998e+00 , 4.199999999999996e+00 , -6.999999999999998e-01 , 4.230524660745312e-15 , 2.843672009884914e-16 , 5.265629202507887e-15 , -1.554312234475219e-15 , 4.999999999999987e-01 , -3.600000000000002e+00 , 1.600000000000000e+00 , 1.113904886696612e-16 , 8.399999999999999e+00 , 2.100000000000000e+00 , -5.042395072913611e-16 , -2.241092798866429e-15 , -1.776356839400246e-15 , -8.881784197001252e-16 , 2.500000000000001e+00 , 9.516197353930126e-17 , -2.000000000000000e+00 , -1.380981496957210e-15 , 4.440892098500626e-16 , 1.050000000000000e+01 , 2.664535259100376e-15 , -5.958952152579898e-16 , -3.679596310186229e-15 , 8.881784197001252e-16 };
    return dmats1_data;
  }
#endif


  public:

  /** \brief Constructor.
  */
  Basis_F0_TRI_C3_FEM_FIAT() : numDof_(10), isSet_(false) , vdm_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C3_FEM_FIAT::get_vdm_data(),10,10,10) ) ),dmats0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C3_FEM_FIAT::get_dmats0_data(),10,10,10) ) ),dmats1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C3_FEM_FIAT::get_dmats1_data(),10,10,10) ) ) { }
  
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

#include "Intrepid_F0_TRI_C3_FEM_FIATDef.hpp"

#endif

