/*
 * MueLu_DummyFactory.hpp
 *
 *  Created on: 11.09.2011
 *      Author: tobias
 */

#ifndef MUELU_DUMMYFACTORY_HPP_
#define MUELU_DUMMYFACTORY_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu
{
    class Level;

    /*!
      @class DummyFactory class.
      @brief Dummy Factory that is used for data stored in level class for that no generating factory is available/necessary.
    */
    class DummyFactory : public SingleLevelFactoryBase {
    public:
        //@{ Constructors/Destructors.

        //! Constructor.
        DummyFactory() {}

        //! Destructor.
        virtual ~DummyFactory() {}
        //@}

        //! Input
        //@{

        void DeclareInput(Level &currentLevel) const
        { std::cout << "Declare input of DummyFactory " << std::endl;     }

        //@}

        //@{
        //! @name Build methods.

        //! Build. The DummyFactory has no support for a Build routine. throw error
        bool Build(Level & currentLevel) const
        {
            TEST_FOR_EXCEPTION(1, MueLu::Exceptions::RuntimeError, "DummyFactory::Build(): You cannot call the Build method for DummyFactory!");
            return false;
        }

        //!
        bool NewBuild(Level & requestedLevel) const {
          return Build(requestedLevel);
        }

        //!
        void callDeclareInput(Level & requestedLevel) const {
          DeclareInput(requestedLevel); }

    }; // end DummyFactory

} // end namespace MueLu

#endif /* MUELU_DUMMYFACTORY_HPP_ */
