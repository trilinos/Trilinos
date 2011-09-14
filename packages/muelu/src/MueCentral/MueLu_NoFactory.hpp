/*
 * MueLu_NoFactory.hpp
 *
 *  Created on: Sep 13, 2011
 *      Author: wiesner
 */

#ifndef MUELU_NOFACTORY_HPP_
#define MUELU_NOFACTORY_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu
{
    class Level;

    /*!
      @class NoFactory class.
      @brief NoFactory that is used for data stored in level class for that no generating factory is available/necessary.
    */
    class NoFactory : public SingleLevelFactoryBase {
    private:
        static bool instanceFlag_; // flag for singleton instance
        static NoFactory* UserDefined; // static NoFactory instance for user defined "factories"

        //! Constructor.
        NoFactory() {}

    public:
        //@{ Destructor.


        //! Destructor.
        virtual ~NoFactory() { NoFactory::instanceFlag_ = false; }
        //@}

        //! Input
        //@{

        void DeclareInput(Level &currentLevel) const
        { std::cout << "Declare input of DummyFactory " << std::endl;     }

        //@}

        //@{
        //! @name Build methods.

        //! Build. The NoFactory has no support for a Build routine. throw error
        bool Build(Level & currentLevel) const
        {
            TEST_FOR_EXCEPTION(1, MueLu::Exceptions::RuntimeError, "NoFactory::Build(): You cannot call the Build method for NoFactory!");
            return false;
        }

        //!
        bool NewBuild(Level & requestedLevel) const {
          return Build(requestedLevel);
        }

        //!
        void callDeclareInput(Level & requestedLevel) const {
          DeclareInput(requestedLevel); }

        static const MueLu::NoFactory* get()
        {
        	if(!instanceFlag_)
        	{
        		NoFactory::UserDefined = new NoFactory();
        		NoFactory::instanceFlag_ = true;
        		return NoFactory::UserDefined;
        	}
        	else
        		return MueLu::NoFactory::UserDefined;
        }

        typedef const NoFactory* (*ptrGetInstance)();

    }; // end NoFactory






} // end namespace MueLu

#endif /* MUELU_NOFACTORY_HPP_ */
