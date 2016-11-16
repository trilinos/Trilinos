#ifndef Tempus_RKButcherTableauBuilder_decl_hpp
#define Tempus_RKButcherTableauBuilder_decl_hpp

#include "Tempus_RKButcherTableau.hpp"
#include "Teuchos_ObjectBuilder.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"


namespace Tempus {

/** \brief Runge-Kutta Builder class.
 *  This factory creates RKButcherTableau objects given the description
 *  string from the RKButcherTableau object.
 *
 *  This was taken and modified from Rythmos' RKButcherTableauBuilder class.
 */
template<class Scalar>
class RKButcherTableauBuilder :
  virtual public Teuchos::ParameterListAcceptor
{
  public:
    RKButcherTableauBuilder();
    virtual ~RKButcherTableauBuilder() {}

    void setRKButcherTableauFactory(
      const Teuchos::RCP<const Teuchos::AbstractFactory<RKButcherTableau<Scalar> > >
        &rkbtFactory,
      const std::string &rkbtFactoryName
      );

    Teuchos::RCP<RKButcherTableau<Scalar> > create(
      const std::string &rkbt_name = ""
      ) const;

    /** \name Overridden from Teuchos::ParameterListAcceptor */
    //@{
      /** \brief . */
      void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);
      /** \brief . */
      Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
      /** \brief . */
      Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
      /** \brief. */
      Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
      /** \brief. */
      Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    //@}
  private:
    Teuchos::ObjectBuilder<RKButcherTableau<Scalar> > builder_;

    void initializeDefaults_();
};

// Nonmember constructor
template<class Scalar>
Teuchos::RCP<RKButcherTableauBuilder<Scalar> > rKButcherTableauBuilder();

// Nonmember helper function
template<class Scalar>
Teuchos::RCP<RKButcherTableau<Scalar> >
createRKBT(const std::string& rkbt_name,
           Teuchos::RCP<Teuchos::ParameterList> pl);

} // namespace Tempus


#endif // Tempus_RKButcherTableauBuilder_decl_hpp
