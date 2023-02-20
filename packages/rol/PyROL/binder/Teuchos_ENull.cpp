#include <PyROL_Teuchos_Custom.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Dependency.hpp>
#include <Teuchos_DependencySheet.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_Exceptions.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListModifier.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPNode.hpp>
#include <Teuchos_StringIndexedOrderedValueObjectContainer.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_XMLObjectImplem.hpp>
#include <Teuchos_any.hpp>
#include <Teuchos_basic_oblackholestream.hpp>
#include <Teuchos_dyn_cast.hpp>
#include <Teuchos_toString.hpp>
#include <cwchar>
#include <deque>
#include <functional>
#include <ios>
#include <iterator>
#include <locale>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
#include <sstream> // __str__
#include <stdexcept>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>

#include <functional>
#include "PyROL_Smart_Holder.hpp"
#include <string>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

// Teuchos::ExceptionBase file:Teuchos_Exceptions.hpp line:57
struct PyCallBack_Teuchos_ExceptionBase : public Teuchos::ExceptionBase {
	using Teuchos::ExceptionBase::ExceptionBase;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::ExceptionBase *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateOwningRCPError file:Teuchos_Exceptions.hpp line:69
struct PyCallBack_Teuchos_DuplicateOwningRCPError : public Teuchos::DuplicateOwningRCPError {
	using Teuchos::DuplicateOwningRCPError::DuplicateOwningRCPError;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateOwningRCPError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::NullReferenceError file:Teuchos_Exceptions.hpp line:77
struct PyCallBack_Teuchos_NullReferenceError : public Teuchos::NullReferenceError {
	using Teuchos::NullReferenceError::NullReferenceError;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::NullReferenceError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::NonconstAccessError file:Teuchos_Exceptions.hpp line:85
struct PyCallBack_Teuchos_NonconstAccessError : public Teuchos::NonconstAccessError {
	using Teuchos::NonconstAccessError::NonconstAccessError;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::NonconstAccessError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::RangeError file:Teuchos_Exceptions.hpp line:93
struct PyCallBack_Teuchos_RangeError : public Teuchos::RangeError {
	using Teuchos::RangeError::RangeError;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RangeError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DanglingReferenceError file:Teuchos_Exceptions.hpp line:101
struct PyCallBack_Teuchos_DanglingReferenceError : public Teuchos::DanglingReferenceError {
	using Teuchos::DanglingReferenceError::DanglingReferenceError;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DanglingReferenceError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::IncompatibleIteratorsError file:Teuchos_Exceptions.hpp line:109
struct PyCallBack_Teuchos_IncompatibleIteratorsError : public Teuchos::IncompatibleIteratorsError {
	using Teuchos::IncompatibleIteratorsError::IncompatibleIteratorsError;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::IncompatibleIteratorsError *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateParameterSublist file:Teuchos_Exceptions.hpp line:118
struct PyCallBack_Teuchos_DuplicateParameterSublist : public Teuchos::DuplicateParameterSublist {
	using Teuchos::DuplicateParameterSublist::DuplicateParameterSublist;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateParameterSublist *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateParameterEntryException file:Teuchos_Exceptions.hpp line:132
struct PyCallBack_Teuchos_DuplicateParameterEntryException : public Teuchos::DuplicateParameterEntryException {
	using Teuchos::DuplicateParameterEntryException::DuplicateParameterEntryException;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateParameterEntryException *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateParameterEntryIDException file:Teuchos_Exceptions.hpp line:145
struct PyCallBack_Teuchos_DuplicateParameterEntryIDException : public Teuchos::DuplicateParameterEntryIDException {
	using Teuchos::DuplicateParameterEntryIDException::DuplicateParameterEntryIDException;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateParameterEntryIDException *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::DuplicateValidatorIDException file:Teuchos_Exceptions.hpp line:158
struct PyCallBack_Teuchos_DuplicateValidatorIDException : public Teuchos::DuplicateValidatorIDException {
	using Teuchos::DuplicateValidatorIDException::DuplicateValidatorIDException;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::DuplicateValidatorIDException *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return logic_error::what();
	}
};

// Teuchos::RCPNode file:Teuchos_RCPNode.hpp line:153
struct PyCallBack_Teuchos_RCPNode : public Teuchos::RCPNode {
	using Teuchos::RCPNode::RCPNode;

	bool is_valid_ptr() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RCPNode *>(this), "is_valid_ptr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"RCPNode::is_valid_ptr\"");
	}
	void delete_obj() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RCPNode *>(this), "delete_obj");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"RCPNode::delete_obj\"");
	}
	void throw_invalid_obj_exception(const std::string & a0, const void * a1, const class Teuchos::RCPNode * a2, const void * a3) const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RCPNode *>(this), "throw_invalid_obj_exception");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"RCPNode::throw_invalid_obj_exception\"");
	}
	const std::string get_base_obj_type_name() const override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::RCPNode *>(this), "get_base_obj_type_name");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const std::string>::value) {
				static pybind11::detail::override_caster_t<const std::string> caster;
				return pybind11::detail::cast_ref<const std::string>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const std::string>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"RCPNode::get_base_obj_type_name\"");
	}
};

// Teuchos::m_bad_cast file:Teuchos_dyn_cast.hpp line:60
struct PyCallBack_Teuchos_m_bad_cast : public Teuchos::m_bad_cast {
	using Teuchos::m_bad_cast::m_bad_cast;

	const char * what() const throw() override {
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Teuchos::m_bad_cast *>(this), "what");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<const char *>::value) {
				static pybind11::detail::override_caster_t<const char *> caster;
				return pybind11::detail::cast_ref<const char *>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<const char *>(std::move(o));
		}
		return m_bad_cast::what();
	}
};

void bind_Teuchos_ENull(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	// Teuchos::ENull file:Teuchos_ENull.hpp line:54
	pybind11::enum_<Teuchos::ENull>(M("Teuchos"), "ENull", pybind11::arithmetic(), "Used to initialize a RCP object to NULL using an\n implicit conversion!\n\n \n\n ", pybind11::module_local())
		.value("null", Teuchos::null)
		.export_values();

;

	{ // Teuchos::ExceptionBase file:Teuchos_Exceptions.hpp line:57
		pybind11::class_<Teuchos::ExceptionBase, std::shared_ptr<Teuchos::ExceptionBase>, PyCallBack_Teuchos_ExceptionBase, std::logic_error> cl(M("Teuchos"), "ExceptionBase", "Base exception class for Teuchos\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_ExceptionBase const &o){ return new PyCallBack_Teuchos_ExceptionBase(o); } ) );
		cl.def( pybind11::init( [](Teuchos::ExceptionBase const &o){ return new Teuchos::ExceptionBase(o); } ) );
		cl.def("assign", (class Teuchos::ExceptionBase & (Teuchos::ExceptionBase::*)(const class Teuchos::ExceptionBase &)) &Teuchos::ExceptionBase::operator=, "C++: Teuchos::ExceptionBase::operator=(const class Teuchos::ExceptionBase &) --> class Teuchos::ExceptionBase &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateOwningRCPError file:Teuchos_Exceptions.hpp line:69
		pybind11::class_<Teuchos::DuplicateOwningRCPError, std::shared_ptr<Teuchos::DuplicateOwningRCPError>, PyCallBack_Teuchos_DuplicateOwningRCPError, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateOwningRCPError", "Thrown if a duplicate owning RCP is creatd the the same object.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateOwningRCPError & (Teuchos::DuplicateOwningRCPError::*)(const class Teuchos::DuplicateOwningRCPError &)) &Teuchos::DuplicateOwningRCPError::operator=, "C++: Teuchos::DuplicateOwningRCPError::operator=(const class Teuchos::DuplicateOwningRCPError &) --> class Teuchos::DuplicateOwningRCPError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::NullReferenceError file:Teuchos_Exceptions.hpp line:77
		pybind11::class_<Teuchos::NullReferenceError, std::shared_ptr<Teuchos::NullReferenceError>, PyCallBack_Teuchos_NullReferenceError, Teuchos::ExceptionBase> cl(M("Teuchos"), "NullReferenceError", "Null reference error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_NullReferenceError const &o){ return new PyCallBack_Teuchos_NullReferenceError(o); } ) );
		cl.def( pybind11::init( [](Teuchos::NullReferenceError const &o){ return new Teuchos::NullReferenceError(o); } ) );
		cl.def("assign", (class Teuchos::NullReferenceError & (Teuchos::NullReferenceError::*)(const class Teuchos::NullReferenceError &)) &Teuchos::NullReferenceError::operator=, "C++: Teuchos::NullReferenceError::operator=(const class Teuchos::NullReferenceError &) --> class Teuchos::NullReferenceError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::NonconstAccessError file:Teuchos_Exceptions.hpp line:85
		pybind11::class_<Teuchos::NonconstAccessError, std::shared_ptr<Teuchos::NonconstAccessError>, PyCallBack_Teuchos_NonconstAccessError, Teuchos::ExceptionBase> cl(M("Teuchos"), "NonconstAccessError", "Null reference error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::NonconstAccessError & (Teuchos::NonconstAccessError::*)(const class Teuchos::NonconstAccessError &)) &Teuchos::NonconstAccessError::operator=, "C++: Teuchos::NonconstAccessError::operator=(const class Teuchos::NonconstAccessError &) --> class Teuchos::NonconstAccessError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::RangeError file:Teuchos_Exceptions.hpp line:93
		pybind11::class_<Teuchos::RangeError, std::shared_ptr<Teuchos::RangeError>, PyCallBack_Teuchos_RangeError, Teuchos::ExceptionBase> cl(M("Teuchos"), "RangeError", "Range error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_RangeError const &o){ return new PyCallBack_Teuchos_RangeError(o); } ) );
		cl.def( pybind11::init( [](Teuchos::RangeError const &o){ return new Teuchos::RangeError(o); } ) );
		cl.def("assign", (class Teuchos::RangeError & (Teuchos::RangeError::*)(const class Teuchos::RangeError &)) &Teuchos::RangeError::operator=, "C++: Teuchos::RangeError::operator=(const class Teuchos::RangeError &) --> class Teuchos::RangeError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DanglingReferenceError file:Teuchos_Exceptions.hpp line:101
		pybind11::class_<Teuchos::DanglingReferenceError, std::shared_ptr<Teuchos::DanglingReferenceError>, PyCallBack_Teuchos_DanglingReferenceError, Teuchos::ExceptionBase> cl(M("Teuchos"), "DanglingReferenceError", "Dangling reference error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_DanglingReferenceError const &o){ return new PyCallBack_Teuchos_DanglingReferenceError(o); } ) );
		cl.def( pybind11::init( [](Teuchos::DanglingReferenceError const &o){ return new Teuchos::DanglingReferenceError(o); } ) );
		cl.def("assign", (class Teuchos::DanglingReferenceError & (Teuchos::DanglingReferenceError::*)(const class Teuchos::DanglingReferenceError &)) &Teuchos::DanglingReferenceError::operator=, "C++: Teuchos::DanglingReferenceError::operator=(const class Teuchos::DanglingReferenceError &) --> class Teuchos::DanglingReferenceError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::IncompatibleIteratorsError file:Teuchos_Exceptions.hpp line:109
		pybind11::class_<Teuchos::IncompatibleIteratorsError, std::shared_ptr<Teuchos::IncompatibleIteratorsError>, PyCallBack_Teuchos_IncompatibleIteratorsError, Teuchos::ExceptionBase> cl(M("Teuchos"), "IncompatibleIteratorsError", "Incompatiable iterators error exception class.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::IncompatibleIteratorsError & (Teuchos::IncompatibleIteratorsError::*)(const class Teuchos::IncompatibleIteratorsError &)) &Teuchos::IncompatibleIteratorsError::operator=, "C++: Teuchos::IncompatibleIteratorsError::operator=(const class Teuchos::IncompatibleIteratorsError &) --> class Teuchos::IncompatibleIteratorsError &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateParameterSublist file:Teuchos_Exceptions.hpp line:118
		pybind11::class_<Teuchos::DuplicateParameterSublist, std::shared_ptr<Teuchos::DuplicateParameterSublist>, PyCallBack_Teuchos_DuplicateParameterSublist, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateParameterSublist", "Optionally thrown when a sublist is set twice by either\n updateParametersFromXmlFile(), updateParametersFromXmlFileAndUpdate() or\n updateParametersFromXmlString()\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateParameterSublist & (Teuchos::DuplicateParameterSublist::*)(const class Teuchos::DuplicateParameterSublist &)) &Teuchos::DuplicateParameterSublist::operator=, "C++: Teuchos::DuplicateParameterSublist::operator=(const class Teuchos::DuplicateParameterSublist &) --> class Teuchos::DuplicateParameterSublist &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateParameterEntryException file:Teuchos_Exceptions.hpp line:132
		pybind11::class_<Teuchos::DuplicateParameterEntryException, std::shared_ptr<Teuchos::DuplicateParameterEntryException>, PyCallBack_Teuchos_DuplicateParameterEntryException, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateParameterEntryException", "Thrown when a Parameter Entry that is already being tracked\n is attempted to be inserted again into the masterParameterEntryMap\n and masterIDMap\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateParameterEntryException & (Teuchos::DuplicateParameterEntryException::*)(const class Teuchos::DuplicateParameterEntryException &)) &Teuchos::DuplicateParameterEntryException::operator=, "C++: Teuchos::DuplicateParameterEntryException::operator=(const class Teuchos::DuplicateParameterEntryException &) --> class Teuchos::DuplicateParameterEntryException &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateParameterEntryIDException file:Teuchos_Exceptions.hpp line:145
		pybind11::class_<Teuchos::DuplicateParameterEntryIDException, std::shared_ptr<Teuchos::DuplicateParameterEntryIDException>, PyCallBack_Teuchos_DuplicateParameterEntryIDException, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateParameterEntryIDException", "Thrown when a Parameter Entry ID that is already being used\n is attempted to be reused again.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateParameterEntryIDException & (Teuchos::DuplicateParameterEntryIDException::*)(const class Teuchos::DuplicateParameterEntryIDException &)) &Teuchos::DuplicateParameterEntryIDException::operator=, "C++: Teuchos::DuplicateParameterEntryIDException::operator=(const class Teuchos::DuplicateParameterEntryIDException &) --> class Teuchos::DuplicateParameterEntryIDException &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::DuplicateValidatorIDException file:Teuchos_Exceptions.hpp line:158
		pybind11::class_<Teuchos::DuplicateValidatorIDException, std::shared_ptr<Teuchos::DuplicateValidatorIDException>, PyCallBack_Teuchos_DuplicateValidatorIDException, Teuchos::ExceptionBase> cl(M("Teuchos"), "DuplicateValidatorIDException", "Thrown when a ParameterEntryValidatorID that\n is already being used is attempted to be reused again.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def("assign", (class Teuchos::DuplicateValidatorIDException & (Teuchos::DuplicateValidatorIDException::*)(const class Teuchos::DuplicateValidatorIDException &)) &Teuchos::DuplicateValidatorIDException::operator=, "C++: Teuchos::DuplicateValidatorIDException::operator=(const class Teuchos::DuplicateValidatorIDException &) --> class Teuchos::DuplicateValidatorIDException &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:60
		pybind11::class_<Teuchos::ToStringTraits<int>, std::shared_ptr<Teuchos::ToStringTraits<int>>> cl(M("Teuchos"), "ToStringTraits_int_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<int>(); } ) );
		cl.def_static("toString", (std::string (*)(const int &)) &Teuchos::ToStringTraits<int>::toString, "C++: Teuchos::ToStringTraits<int>::toString(const int &) --> std::string", pybind11::arg("t"));
	}
	// Teuchos::toString(const double &) file:Teuchos_toString.hpp line:82
	M("Teuchos").def("toString", (std::string (*)(const double &)) &Teuchos::toString<double>, "C++: Teuchos::toString(const double &) --> std::string", pybind11::arg("t"));

	// Teuchos::toString(const int &) file:Teuchos_toString.hpp line:82
	M("Teuchos").def("toString", (std::string (*)(const int &)) &Teuchos::toString<int>, "C++: Teuchos::toString(const int &) --> std::string", pybind11::arg("t"));

	// Teuchos::toString(const bool &) file:Teuchos_toString.hpp line:82
	M("Teuchos").def("toString", (std::string (*)(const bool &)) &Teuchos::toString<bool>, "C++: Teuchos::toString(const bool &) --> std::string", pybind11::arg("t"));

	// Teuchos::toString(const std::string &) file:Teuchos_toString.hpp line:82
	M("Teuchos").def("toString", (std::string (*)(const std::string &)) &Teuchos::toString<std::string>, "C++: Teuchos::toString(const std::string &) --> std::string", pybind11::arg("t"));

	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:90
		pybind11::class_<Teuchos::ToStringTraits<bool>, std::shared_ptr<Teuchos::ToStringTraits<bool>>> cl(M("Teuchos"), "ToStringTraits_bool_t", "Specialization for bool. ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<bool>(); } ) );
		cl.def_static("toString", (std::string (*)(const bool &)) &Teuchos::ToStringTraits<bool>::toString, "C++: Teuchos::ToStringTraits<bool>::toString(const bool &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:103
		pybind11::class_<Teuchos::ToStringTraits<std::string>, std::shared_ptr<Teuchos::ToStringTraits<std::string>>> cl(M("Teuchos"), "ToStringTraits_std_string_t", "Specialization for std::string. ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<std::string>(); } ) );
		cl.def_static("toString", (std::string (*)(const std::string &)) &Teuchos::ToStringTraits<std::string >::toString, "C++: Teuchos::ToStringTraits<std::string >::toString(const std::string &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:113
		pybind11::class_<Teuchos::ToStringTraits<double>, std::shared_ptr<Teuchos::ToStringTraits<double>>> cl(M("Teuchos"), "ToStringTraits_double_t", "Specialization for double. ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<double>(); } ) );
		cl.def_static("toString", (std::string (*)(const double &)) &Teuchos::ToStringTraits<double>::toString, "C++: Teuchos::ToStringTraits<double>::toString(const double &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::ToStringTraits file:Teuchos_toString.hpp line:149
		pybind11::class_<Teuchos::ToStringTraits<float>, std::shared_ptr<Teuchos::ToStringTraits<float>>> cl(M("Teuchos"), "ToStringTraits_float_t", "Specialization for float. ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<float>(); } ) );
		cl.def_static("toString", (std::string (*)(const float &)) &Teuchos::ToStringTraits<float>::toString, "C++: Teuchos::ToStringTraits<float>::toString(const float &) --> std::string", pybind11::arg("t"));
	}
	// Teuchos::EPrePostDestruction file:Teuchos_RCPNode.hpp line:79
	pybind11::enum_<Teuchos::EPrePostDestruction>(M("Teuchos"), "EPrePostDestruction", pybind11::arithmetic(), "Used to specify a pre or post destruction of extra data\n\n \n\n ", pybind11::module_local())
		.value("PRE_DESTROY", Teuchos::PRE_DESTROY)
		.value("POST_DESTROY", Teuchos::POST_DESTROY)
		.export_values();

;

	// Teuchos::ERCPStrength file:Teuchos_RCPNode.hpp line:85
	pybind11::enum_<Teuchos::ERCPStrength>(M("Teuchos"), "ERCPStrength", pybind11::arithmetic(), "Used to specify if the pointer is weak or strong.\n\n \n\n ", pybind11::module_local())
		.value("RCP_STRONG", Teuchos::RCP_STRONG)
		.value("RCP_WEAK", Teuchos::RCP_WEAK)
		.export_values();

;

	// Teuchos::ERCPNodeLookup file:Teuchos_RCPNode.hpp line:91
	pybind11::enum_<Teuchos::ERCPNodeLookup>(M("Teuchos"), "ERCPNodeLookup", pybind11::arithmetic(), "Used to determine if RCPNode lookup is performed or not.\n\n \n\n ", pybind11::module_local())
		.value("RCP_ENABLE_NODE_LOOKUP", Teuchos::RCP_ENABLE_NODE_LOOKUP)
		.value("RCP_DISABLE_NODE_LOOKUP", Teuchos::RCP_DISABLE_NODE_LOOKUP)
		.export_values();

;

	// Teuchos::debugAssertStrength(enum Teuchos::ERCPStrength) file:Teuchos_RCPNode.hpp line:94
	M("Teuchos").def("debugAssertStrength", (void (*)(enum Teuchos::ERCPStrength)) &Teuchos::debugAssertStrength, ". \n\nC++: Teuchos::debugAssertStrength(enum Teuchos::ERCPStrength) --> void", pybind11::arg("strength"));

	{ // Teuchos::ToStringTraits file:Teuchos_RCPNode.hpp line:119
		pybind11::class_<Teuchos::ToStringTraits<Teuchos::ERCPStrength>, std::shared_ptr<Teuchos::ToStringTraits<Teuchos::ERCPStrength>>> cl(M("Teuchos"), "ToStringTraits_Teuchos_ERCPStrength_t", "Traits class specialization for toString(...) function for\n converting from ERCPStrength to std::string.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ToStringTraits<Teuchos::ERCPStrength>(); } ) );
		cl.def_static("toString", (std::string (*)(const enum Teuchos::ERCPStrength &)) &Teuchos::ToStringTraits<Teuchos::ERCPStrength>::toString, "C++: Teuchos::ToStringTraits<Teuchos::ERCPStrength>::toString(const enum Teuchos::ERCPStrength &) --> std::string", pybind11::arg("t"));
	}
	{ // Teuchos::RCPNode file:Teuchos_RCPNode.hpp line:153
		pybind11::class_<Teuchos::RCPNode, std::shared_ptr<Teuchos::RCPNode>, PyCallBack_Teuchos_RCPNode> cl(M("Teuchos"), "RCPNode", "Node class to keep track of address and the reference count for a\n reference-counted utility class and delete the object.\n\n This is not a general user-level class.  This is used in the implementation\n of all of the reference-counting utility classes.\n\n NOTE: The reference counts all start a 0 so the client (i.e. RCPNodeHandle)\n must increment them from 0 after creation.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init<bool>(), pybind11::arg("has_ownership_in") );

		cl.def("attemptIncrementStrongCountFromNonZeroValue", (bool (Teuchos::RCPNode::*)()) &Teuchos::RCPNode::attemptIncrementStrongCountFromNonZeroValue, "attemptIncrementStrongCountFromNonZeroValue() supports weak\n to strong conversion but this is forward looking code.\n\nC++: Teuchos::RCPNode::attemptIncrementStrongCountFromNonZeroValue() --> bool");
		cl.def("strong_count", (int (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::strong_count, ". \n\nC++: Teuchos::RCPNode::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::weak_count, ". \n\nC++: Teuchos::RCPNode::weak_count() const --> int");
		cl.def("incr_count", (void (Teuchos::RCPNode::*)(const enum Teuchos::ERCPStrength)) &Teuchos::RCPNode::incr_count, ". \n\nC++: Teuchos::RCPNode::incr_count(const enum Teuchos::ERCPStrength) --> void", pybind11::arg("strength"));
		cl.def("deincr_count", (int (Teuchos::RCPNode::*)(const enum Teuchos::ERCPStrength)) &Teuchos::RCPNode::deincr_count, ". \n\nC++: Teuchos::RCPNode::deincr_count(const enum Teuchos::ERCPStrength) --> int", pybind11::arg("strength"));
		cl.def("has_ownership", (void (Teuchos::RCPNode::*)(bool)) &Teuchos::RCPNode::has_ownership, ". \n\nC++: Teuchos::RCPNode::has_ownership(bool) --> void", pybind11::arg("has_ownership_in"));
		cl.def("has_ownership", (bool (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::has_ownership, ". \n\nC++: Teuchos::RCPNode::has_ownership() const --> bool");
		cl.def("set_extra_data", (void (Teuchos::RCPNode::*)(const class Teuchos::any &, const std::string &, enum Teuchos::EPrePostDestruction, bool)) &Teuchos::RCPNode::set_extra_data, ". \n\nC++: Teuchos::RCPNode::set_extra_data(const class Teuchos::any &, const std::string &, enum Teuchos::EPrePostDestruction, bool) --> void", pybind11::arg("extra_data"), pybind11::arg("name"), pybind11::arg("destroy_when"), pybind11::arg("force_unique"));
		cl.def("get_extra_data", (class Teuchos::any & (Teuchos::RCPNode::*)(const std::string &, const std::string &)) &Teuchos::RCPNode::get_extra_data, ". \n\nC++: Teuchos::RCPNode::get_extra_data(const std::string &, const std::string &) --> class Teuchos::any &", pybind11::return_value_policy::automatic, pybind11::arg("type_name"), pybind11::arg("name"));
		cl.def("get_optional_extra_data", (class Teuchos::any * (Teuchos::RCPNode::*)(const std::string &, const std::string &)) &Teuchos::RCPNode::get_optional_extra_data, ". \n\nC++: Teuchos::RCPNode::get_optional_extra_data(const std::string &, const std::string &) --> class Teuchos::any *", pybind11::return_value_policy::automatic, pybind11::arg("type_name"), pybind11::arg("name"));
		cl.def("is_valid_ptr", (bool (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::is_valid_ptr, ". \n\nC++: Teuchos::RCPNode::is_valid_ptr() const --> bool");
		cl.def("delete_obj", (void (Teuchos::RCPNode::*)()) &Teuchos::RCPNode::delete_obj, ". \n\nC++: Teuchos::RCPNode::delete_obj() --> void");
		cl.def("throw_invalid_obj_exception", (void (Teuchos::RCPNode::*)(const std::string &, const void *, const class Teuchos::RCPNode *, const void *) const) &Teuchos::RCPNode::throw_invalid_obj_exception, ". \n\nC++: Teuchos::RCPNode::throw_invalid_obj_exception(const std::string &, const void *, const class Teuchos::RCPNode *, const void *) const --> void", pybind11::arg("rcp_type_name"), pybind11::arg("rcp_ptr"), pybind11::arg("rcp_node_ptr"), pybind11::arg("rcp_obj_ptr"));
		cl.def("get_base_obj_type_name", (const std::string (Teuchos::RCPNode::*)() const) &Teuchos::RCPNode::get_base_obj_type_name, ". \n\nC++: Teuchos::RCPNode::get_base_obj_type_name() const --> const std::string");
	}
	// Teuchos::throw_null_ptr_error(const std::string &) file:Teuchos_RCPNode.hpp line:334
	M("Teuchos").def("throw_null_ptr_error", (void (*)(const std::string &)) &Teuchos::throw_null_ptr_error, "Throw that a pointer passed into an RCP object is null.\n\n \n\n \n\nC++: Teuchos::throw_null_ptr_error(const std::string &) --> void", pybind11::arg("type_name"));

	{ // Teuchos::RCPNodeTracer file:Teuchos_RCPNode.hpp line:368
		pybind11::class_<Teuchos::RCPNodeTracer, std::shared_ptr<Teuchos::RCPNodeTracer>> cl(M("Teuchos"), "RCPNodeTracer", "Debug-mode RCPNode tracing class.\n\n This is a static class that is used to trace all RCP nodes that are created\n and destroyed and to look-up RCPNodes given an an object's address.  This\n database is used for several different types of debug-mode runtime checking\n including a) the detection of cicular references, b) detecting the creation\n of duplicate owning RCPNode objects for the same reference-counted object,\n and c) to create weak RCP objects for existing RCPNode objects.\n\n This is primarily an internal implementation class but there are a few\n functions (maked as such below) that can be called by general users to turn\n on and off node tracing and to print the active RCPNode objects at any\n time.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCPNodeTracer(); } ) );
		cl.def_static("isTracingActiveRCPNodes", (bool (*)()) &Teuchos::RCPNodeTracer::isTracingActiveRCPNodes, "Return if we are tracing active nodes or not.\n\n NOTE: This will always return false when TEUCHOS_DEBUG is\n not defined.\n\nC++: Teuchos::RCPNodeTracer::isTracingActiveRCPNodes() --> bool");
		cl.def_static("numActiveRCPNodes", (int (*)()) &Teuchos::RCPNodeTracer::numActiveRCPNodes, "Print the number of active RCPNode objects currently being\n tracked.\n\nC++: Teuchos::RCPNodeTracer::numActiveRCPNodes() --> int");
		cl.def_static("getRCPNodeStatistics", (struct Teuchos::RCPNodeTracer::RCPNodeStatistics (*)()) &Teuchos::RCPNodeTracer::getRCPNodeStatistics, "Return the statistics on RCPNode allocations. \n\nC++: Teuchos::RCPNodeTracer::getRCPNodeStatistics() --> struct Teuchos::RCPNodeTracer::RCPNodeStatistics");
		cl.def_static("setPrintRCPNodeStatisticsOnExit", (void (*)(bool)) &Teuchos::RCPNodeTracer::setPrintRCPNodeStatisticsOnExit, "Set if RCPNode usage statistics will be printed when the program\n ends or not.\n\nC++: Teuchos::RCPNodeTracer::setPrintRCPNodeStatisticsOnExit(bool) --> void", pybind11::arg("printRCPNodeStatisticsOnExit"));
		cl.def_static("getPrintRCPNodeStatisticsOnExit", (bool (*)()) &Teuchos::RCPNodeTracer::getPrintRCPNodeStatisticsOnExit, "Return if RCPNode usage statistics will be printed when the\n program ends or not.\n\nC++: Teuchos::RCPNodeTracer::getPrintRCPNodeStatisticsOnExit() --> bool");
		cl.def_static("setPrintActiveRcpNodesOnExit", (void (*)(bool)) &Teuchos::RCPNodeTracer::setPrintActiveRcpNodesOnExit, "Set if printActiveRCPNodes() is called on exit from the\n program.\n\nC++: Teuchos::RCPNodeTracer::setPrintActiveRcpNodesOnExit(bool) --> void", pybind11::arg("printActiveRcpNodesOnExit"));
		cl.def_static("getPrintActiveRcpNodesOnExit", (bool (*)()) &Teuchos::RCPNodeTracer::getPrintActiveRcpNodesOnExit, "Return if printActiveRCPNodes() is called on exit from the\n program.\n\nC++: Teuchos::RCPNodeTracer::getPrintActiveRcpNodesOnExit() --> bool");
		cl.def_static("addNewRCPNode", (void (*)(class Teuchos::RCPNode *, const std::string &)) &Teuchos::RCPNodeTracer::addNewRCPNode, "Add new RCPNode to the global list.\n\n Only gets called when RCPNode tracing has been activated.\n\nC++: Teuchos::RCPNodeTracer::addNewRCPNode(class Teuchos::RCPNode *, const std::string &) --> void", pybind11::arg("rcp_node"), pybind11::arg("info"));
		cl.def_static("removeRCPNode", (void (*)(class Teuchos::RCPNode *)) &Teuchos::RCPNodeTracer::removeRCPNode, "Remove an RCPNode from global list.\n\n Always gets called in a debug build (TEUCHOS_DEBUG defined) when\n node tracing is enabled.\n\nC++: Teuchos::RCPNodeTracer::removeRCPNode(class Teuchos::RCPNode *) --> void", pybind11::arg("rcp_node"));
		cl.def_static("getExistingRCPNodeGivenLookupKey", (class Teuchos::RCPNode * (*)(const void *)) &Teuchos::RCPNodeTracer::getExistingRCPNodeGivenLookupKey, "Return a raw pointer to an existing owning RCPNode given its\n lookup key.\n\n \n returnVal != 0 if an owning RCPNode exists, 0\n otherwsise.\n\nC++: Teuchos::RCPNodeTracer::getExistingRCPNodeGivenLookupKey(const void *) --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic, pybind11::arg("lookupKey"));
		cl.def_static("getActiveRCPNodeHeaderString", (std::string (*)()) &Teuchos::RCPNodeTracer::getActiveRCPNodeHeaderString, "Header string used in printActiveRCPNodes(). \n\nC++: Teuchos::RCPNodeTracer::getActiveRCPNodeHeaderString() --> std::string");
		cl.def_static("getCommonDebugNotesString", (std::string (*)()) &Teuchos::RCPNodeTracer::getCommonDebugNotesString, "Common error message string on how to debug RCPNode problems. \n\nC++: Teuchos::RCPNodeTracer::getCommonDebugNotesString() --> std::string");

		{ // Teuchos::RCPNodeTracer::RCPNodeStatistics file:Teuchos_RCPNode.hpp line:375
			auto & enclosing_class = cl;
			pybind11::class_<Teuchos::RCPNodeTracer::RCPNodeStatistics, std::shared_ptr<Teuchos::RCPNodeTracer::RCPNodeStatistics>> cl(enclosing_class, "RCPNodeStatistics", "RCP statistics struct. ", pybind11::module_local());
			cl.def( pybind11::init( [](){ return new Teuchos::RCPNodeTracer::RCPNodeStatistics(); } ) );
			cl.def_readwrite("maxNumRCPNodes", &Teuchos::RCPNodeTracer::RCPNodeStatistics::maxNumRCPNodes);
			cl.def_readwrite("totalNumRCPNodeAllocations", &Teuchos::RCPNodeTracer::RCPNodeStatistics::totalNumRCPNodeAllocations);
			cl.def_readwrite("totalNumRCPNodeDeletions", &Teuchos::RCPNodeTracer::RCPNodeStatistics::totalNumRCPNodeDeletions);
		}

	}
	{ // Teuchos::ActiveRCPNodesSetup file:Teuchos_RCPNode.hpp line:703
		pybind11::class_<Teuchos::ActiveRCPNodesSetup, std::shared_ptr<Teuchos::ActiveRCPNodesSetup>> cl(M("Teuchos"), "ActiveRCPNodesSetup", "Sets up node tracing and prints remaining RCPNodes on destruction.\n\n This class is used by automataic code that sets up support for RCPNode\n tracing and for printing of remaining nodes on destruction.\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::ActiveRCPNodesSetup(); } ) );
		cl.def( pybind11::init( [](Teuchos::ActiveRCPNodesSetup const &o){ return new Teuchos::ActiveRCPNodesSetup(o); } ) );
		cl.def("foo", (void (Teuchos::ActiveRCPNodesSetup::*)()) &Teuchos::ActiveRCPNodesSetup::foo, ". \n\nC++: Teuchos::ActiveRCPNodesSetup::foo() --> void");
	}
	{ // Teuchos::RCPNodeHandle file:Teuchos_RCPNode.hpp line:748
		pybind11::class_<Teuchos::RCPNodeHandle, std::shared_ptr<Teuchos::RCPNodeHandle>> cl(M("Teuchos"), "RCPNodeHandle", "Handle class that manages the RCPNode's reference counting.\n\n \n This class is not intended for Teuchos users.  It\n   is an implementation detail of Teuchos' reference-counting\n   \"smart\" pointer (RCP) and array (ArrayRCP) classes.\n\n NOTE: I (Ross Bartlett) am not generally a big fan of handle classes and\n greatly prefer smart pointers.  However, this is one case where a handle\n class makes sense.  First, I want special behavior in some functions when\n the wrapped RCPNode pointer is null.  Second, I can't use one of the\n smart-pointer classes because this class is used to implement all of those\n smart-pointer classes!\n\n \n\n ", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCPNodeHandle(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class Teuchos::RCPNode * a0){ return new Teuchos::RCPNodeHandle(a0); } ), "doc" , pybind11::arg("node"));
		cl.def( pybind11::init( [](class Teuchos::RCPNode * a0, enum Teuchos::ERCPStrength const & a1){ return new Teuchos::RCPNodeHandle(a0, a1); } ), "doc" , pybind11::arg("node"), pybind11::arg("strength_in"));
		cl.def( pybind11::init<class Teuchos::RCPNode *, enum Teuchos::ERCPStrength, bool>(), pybind11::arg("node"), pybind11::arg("strength_in"), pybind11::arg("newNode") );

		cl.def( pybind11::init( [](Teuchos::RCPNodeHandle const &o){ return new Teuchos::RCPNodeHandle(o); } ) );
		cl.def("swap", (void (Teuchos::RCPNodeHandle::*)(class Teuchos::RCPNodeHandle &)) &Teuchos::RCPNodeHandle::swap, "Swap the contents of  with \n\nC++: Teuchos::RCPNodeHandle::swap(class Teuchos::RCPNodeHandle &) --> void", pybind11::arg("node_ref"));
		cl.def("assign", (class Teuchos::RCPNodeHandle & (Teuchos::RCPNodeHandle::*)(enum Teuchos::ENull)) &Teuchos::RCPNodeHandle::operator=, "Null assignment.\n\n This method satisfies the strong exception guarantee: It either\n returns successfully, or throws an exception without modifying\n any user-visible state.\n\nC++: Teuchos::RCPNodeHandle::operator=(enum Teuchos::ENull) --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("assign", (class Teuchos::RCPNodeHandle & (Teuchos::RCPNodeHandle::*)(const class Teuchos::RCPNodeHandle &)) &Teuchos::RCPNodeHandle::operator=, "Copy assignment operator.\n\n This method satisfies the strong exception guarantee: It either\n returns successfully, or throws an exception without modifying\n any user-visible state.\n\nC++: Teuchos::RCPNodeHandle::operator=(const class Teuchos::RCPNodeHandle &) --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic, pybind11::arg("node_ref"));
		cl.def("create_strong_lock", (class Teuchos::RCPNodeHandle (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::create_strong_lock, "Return a strong handle if possible using thread safe atomics\n\nC++: Teuchos::RCPNodeHandle::create_strong_lock() const --> class Teuchos::RCPNodeHandle");
		cl.def("create_weak", (class Teuchos::RCPNodeHandle (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::create_weak, "Return a weak handle.\n\nC++: Teuchos::RCPNodeHandle::create_weak() const --> class Teuchos::RCPNodeHandle");
		cl.def("create_strong", (class Teuchos::RCPNodeHandle (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::create_strong, "Return a strong handle.\n\nC++: Teuchos::RCPNodeHandle::create_strong() const --> class Teuchos::RCPNodeHandle");
		cl.def("node_ptr", (class Teuchos::RCPNode * (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::node_ptr, "Return a pointer to the underlying RCPNode.\n\nC++: Teuchos::RCPNodeHandle::node_ptr() const --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic);
		cl.def("is_node_null", (bool (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::is_node_null, "Whether the underlying RCPNode is NULL.\n\nC++: Teuchos::RCPNodeHandle::is_node_null() const --> bool");
		cl.def("is_valid_ptr", (bool (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::is_valid_ptr, "Whether the underlying pointer is valid.\n\n \n NULL is a valid pointer; this method returns true in that case.\n\nC++: Teuchos::RCPNodeHandle::is_valid_ptr() const --> bool");
		cl.def("same_node", (bool (Teuchos::RCPNodeHandle::*)(const class Teuchos::RCPNodeHandle &) const) &Teuchos::RCPNodeHandle::same_node, "Whether the RCPNode for which  is a handle is the\n   same RCPNode as this object's RCPNode.\n\nC++: Teuchos::RCPNodeHandle::same_node(const class Teuchos::RCPNodeHandle &) const --> bool", pybind11::arg("node2"));
		cl.def("strong_count", (int (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::strong_count, "The strong count for this RCPNode, or 0 if the node is NULL.\n\nC++: Teuchos::RCPNodeHandle::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::weak_count, "The weak count for this RCPNode, or 0 if the node is NULL.\n\nC++: Teuchos::RCPNodeHandle::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::total_count, "The sum of the weak and string counts.\n\nC++: Teuchos::RCPNodeHandle::total_count() const --> int");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::strength, "The strength of this handle.\n\nC++: Teuchos::RCPNodeHandle::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("has_ownership", (void (Teuchos::RCPNodeHandle::*)(bool)) &Teuchos::RCPNodeHandle::has_ownership, ". \n\nC++: Teuchos::RCPNodeHandle::has_ownership(bool) --> void", pybind11::arg("has_ownership_in"));
		cl.def("has_ownership", (bool (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::has_ownership, ". \n\nC++: Teuchos::RCPNodeHandle::has_ownership() const --> bool");
		cl.def("set_extra_data", (void (Teuchos::RCPNodeHandle::*)(const class Teuchos::any &, const std::string &, enum Teuchos::EPrePostDestruction, bool)) &Teuchos::RCPNodeHandle::set_extra_data, ". \n\nC++: Teuchos::RCPNodeHandle::set_extra_data(const class Teuchos::any &, const std::string &, enum Teuchos::EPrePostDestruction, bool) --> void", pybind11::arg("extra_data"), pybind11::arg("name"), pybind11::arg("destroy_when"), pybind11::arg("force_unique"));
		cl.def("get_extra_data", (class Teuchos::any & (Teuchos::RCPNodeHandle::*)(const std::string &, const std::string &)) &Teuchos::RCPNodeHandle::get_extra_data, ". \n\nC++: Teuchos::RCPNodeHandle::get_extra_data(const std::string &, const std::string &) --> class Teuchos::any &", pybind11::return_value_policy::automatic, pybind11::arg("type_name"), pybind11::arg("name"));
		cl.def("get_optional_extra_data", (class Teuchos::any * (Teuchos::RCPNodeHandle::*)(const std::string &, const std::string &)) &Teuchos::RCPNodeHandle::get_optional_extra_data, ". \n\nC++: Teuchos::RCPNodeHandle::get_optional_extra_data(const std::string &, const std::string &) --> class Teuchos::any *", pybind11::return_value_policy::automatic, pybind11::arg("type_name"), pybind11::arg("name"));
		cl.def("debug_assert_not_null", (void (Teuchos::RCPNodeHandle::*)() const) &Teuchos::RCPNodeHandle::debug_assert_not_null, ". \n\nC++: Teuchos::RCPNodeHandle::debug_assert_not_null() const --> void");

		cl.def("__str__", [](Teuchos::RCPNodeHandle const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCPNodeThrowDeleter file:Teuchos_RCPNode.hpp line:1110
		pybind11::class_<Teuchos::RCPNodeThrowDeleter, std::shared_ptr<Teuchos::RCPNodeThrowDeleter>> cl(M("Teuchos"), "RCPNodeThrowDeleter", "Deletes a (non-owning) RCPNode but not it's underlying object in\n case of a throw.\n\n This class is used in contexts where RCPNodeTracer::addNewRCPNode(...)\n might thrown an exception for a duplicate node being added.  The assumption\n is that there must already be an owning (or non-owning) RCP object that\n will delete the underlying object and therefore this class should *not*\n call delete_obj()!", pybind11::module_local());
		cl.def( pybind11::init<class Teuchos::RCPNode *>(), pybind11::arg("node") );

		cl.def("get", (class Teuchos::RCPNode * (Teuchos::RCPNodeThrowDeleter::*)() const) &Teuchos::RCPNodeThrowDeleter::get, ". \n\nC++: Teuchos::RCPNodeThrowDeleter::get() const --> class Teuchos::RCPNode *", pybind11::return_value_policy::automatic);
		cl.def("release", (void (Teuchos::RCPNodeThrowDeleter::*)()) &Teuchos::RCPNodeThrowDeleter::release, "Releaes the RCPNode pointer before the destructor is called. \n\nC++: Teuchos::RCPNodeThrowDeleter::release() --> void");
	}
	{ // Teuchos::Ptr file:Teuchos_PtrDecl.hpp line:104
		pybind11::class_<Teuchos::Ptr<const Teuchos::ParameterEntry>, std::shared_ptr<Teuchos::Ptr<const Teuchos::ParameterEntry>>> cl(M("Teuchos"), "Ptr_const_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::Ptr<const Teuchos::ParameterEntry>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_in") );

		cl.def( pybind11::init<const class Teuchos::ParameterEntry *>(), pybind11::arg("ptr_in") );

		cl.def( pybind11::init( [](Teuchos::Ptr<const Teuchos::ParameterEntry> const &o){ return new Teuchos::Ptr<const Teuchos::ParameterEntry>(o); } ) );
		cl.def("assign", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> & (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &)) &Teuchos::Ptr<const Teuchos::ParameterEntry>::operator=, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::operator=(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) --> class Teuchos::Ptr<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic, pybind11::arg("ptr"));
		cl.def("arrow", (const class Teuchos::ParameterEntry * (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::operator->, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::operator->() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (const class Teuchos::ParameterEntry & (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::operator*, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::operator*() const --> const class Teuchos::ParameterEntry &", pybind11::return_value_policy::automatic);
		cl.def("get", (const class Teuchos::ParameterEntry * (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::get, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::get() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (const class Teuchos::ParameterEntry * (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::getRawPtr, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::getRawPtr() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("is_null", (bool (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::is_null, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::is_null() const --> bool");
		cl.def("assert_not_null", (const class Teuchos::Ptr<const class Teuchos::ParameterEntry> & (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::assert_not_null, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::assert_not_null() const --> const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("ptr", (const class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::ptr, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::ptr() const --> const class Teuchos::Ptr<const class Teuchos::ParameterEntry>");
		cl.def("getConst", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::Ptr<const Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<const Teuchos::ParameterEntry>::getConst, "C++: Teuchos::Ptr<const Teuchos::ParameterEntry>::getConst() const --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>");

		cl.def("__str__", [](Teuchos::Ptr<const Teuchos::ParameterEntry> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::Ptr file:Teuchos_PtrDecl.hpp line:104
		pybind11::class_<Teuchos::Ptr<Teuchos::ParameterEntry>, std::shared_ptr<Teuchos::Ptr<Teuchos::ParameterEntry>>> cl(M("Teuchos"), "Ptr_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::Ptr<Teuchos::ParameterEntry>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_in") );

		cl.def( pybind11::init<class Teuchos::ParameterEntry *>(), pybind11::arg("ptr_in") );

		cl.def( pybind11::init( [](Teuchos::Ptr<Teuchos::ParameterEntry> const &o){ return new Teuchos::Ptr<Teuchos::ParameterEntry>(o); } ) );
		cl.def("assign", (class Teuchos::Ptr<class Teuchos::ParameterEntry> & (Teuchos::Ptr<Teuchos::ParameterEntry>::*)(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &)) &Teuchos::Ptr<Teuchos::ParameterEntry>::operator=, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::operator=(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) --> class Teuchos::Ptr<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic, pybind11::arg("ptr"));
		cl.def("arrow", (class Teuchos::ParameterEntry * (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::operator->, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::operator->() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::ParameterEntry & (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::operator*, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::operator*() const --> class Teuchos::ParameterEntry &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::ParameterEntry * (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::get, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::get() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::ParameterEntry * (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::getRawPtr, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::getRawPtr() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("is_null", (bool (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::is_null, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::is_null() const --> bool");
		cl.def("assert_not_null", (const class Teuchos::Ptr<class Teuchos::ParameterEntry> & (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::assert_not_null, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::assert_not_null() const --> const class Teuchos::Ptr<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("ptr", (const class Teuchos::Ptr<class Teuchos::ParameterEntry> (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::ptr, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::ptr() const --> const class Teuchos::Ptr<class Teuchos::ParameterEntry>");
		cl.def("getConst", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::Ptr<Teuchos::ParameterEntry>::*)() const) &Teuchos::Ptr<Teuchos::ParameterEntry>::getConst, "C++: Teuchos::Ptr<Teuchos::ParameterEntry>::getConst() const --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>");

		cl.def("__str__", [](Teuchos::Ptr<Teuchos::ParameterEntry> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::Ptr file:Teuchos_PtrDecl.hpp line:104
		pybind11::class_<Teuchos::Ptr<Teuchos::ParameterList>, std::shared_ptr<Teuchos::Ptr<Teuchos::ParameterList>>> cl(M("Teuchos"), "Ptr_Teuchos_ParameterList_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::Ptr<Teuchos::ParameterList>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_in") );

		cl.def( pybind11::init<class Teuchos::ParameterList *>(), pybind11::arg("ptr_in") );

		cl.def( pybind11::init( [](Teuchos::Ptr<Teuchos::ParameterList> const &o){ return new Teuchos::Ptr<Teuchos::ParameterList>(o); } ) );
		cl.def("assign", (class Teuchos::Ptr<class Teuchos::ParameterList> & (Teuchos::Ptr<Teuchos::ParameterList>::*)(const class Teuchos::Ptr<class Teuchos::ParameterList> &)) &Teuchos::Ptr<Teuchos::ParameterList>::operator=, "C++: Teuchos::Ptr<Teuchos::ParameterList>::operator=(const class Teuchos::Ptr<class Teuchos::ParameterList> &) --> class Teuchos::Ptr<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic, pybind11::arg("ptr"));
		cl.def("arrow", (class Teuchos::ParameterList * (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::operator->, "C++: Teuchos::Ptr<Teuchos::ParameterList>::operator->() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::ParameterList & (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::operator*, "C++: Teuchos::Ptr<Teuchos::ParameterList>::operator*() const --> class Teuchos::ParameterList &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::ParameterList * (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::get, "C++: Teuchos::Ptr<Teuchos::ParameterList>::get() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::ParameterList * (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::getRawPtr, "C++: Teuchos::Ptr<Teuchos::ParameterList>::getRawPtr() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("is_null", (bool (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::is_null, "C++: Teuchos::Ptr<Teuchos::ParameterList>::is_null() const --> bool");
		cl.def("assert_not_null", (const class Teuchos::Ptr<class Teuchos::ParameterList> & (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::assert_not_null, "C++: Teuchos::Ptr<Teuchos::ParameterList>::assert_not_null() const --> const class Teuchos::Ptr<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("ptr", (const class Teuchos::Ptr<class Teuchos::ParameterList> (Teuchos::Ptr<Teuchos::ParameterList>::*)() const) &Teuchos::Ptr<Teuchos::ParameterList>::ptr, "C++: Teuchos::Ptr<Teuchos::ParameterList>::ptr() const --> const class Teuchos::Ptr<class Teuchos::ParameterList>");

		cl.def("__str__", [](Teuchos::Ptr<Teuchos::ParameterList> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	// Teuchos::ERCPWeakNoDealloc file:Teuchos_RCPDecl.hpp line:75
	pybind11::enum_<Teuchos::ERCPWeakNoDealloc>(M("Teuchos"), "ERCPWeakNoDealloc", pybind11::arithmetic(), "", pybind11::module_local())
		.value("RCP_WEAK_NO_DEALLOC", Teuchos::RCP_WEAK_NO_DEALLOC)
		.export_values();

;

	// Teuchos::ERCPUndefinedWeakNoDealloc file:Teuchos_RCPDecl.hpp line:76
	pybind11::enum_<Teuchos::ERCPUndefinedWeakNoDealloc>(M("Teuchos"), "ERCPUndefinedWeakNoDealloc", pybind11::arithmetic(), "", pybind11::module_local())
		.value("RCP_UNDEFINED_WEAK_NO_DEALLOC", Teuchos::RCP_UNDEFINED_WEAK_NO_DEALLOC)
		.export_values();

;

	// Teuchos::ERCPUndefinedWithDealloc file:Teuchos_RCPDecl.hpp line:77
	pybind11::enum_<Teuchos::ERCPUndefinedWithDealloc>(M("Teuchos"), "ERCPUndefinedWithDealloc", pybind11::arithmetic(), "", pybind11::module_local())
		.value("RCP_UNDEFINED_WITH_DEALLOC", Teuchos::RCP_UNDEFINED_WITH_DEALLOC)
		.export_values();

;

	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<std::basic_ostringstream<char>>, std::shared_ptr<Teuchos::RCP<std::basic_ostringstream<char>>>> cl(M("Teuchos"), "RCP_std_basic_ostringstream_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<std::basic_ostringstream<char>>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](class std::basic_ostringstream<char> * a0){ return new Teuchos::RCP<std::basic_ostringstream<char>>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class std::basic_ostringstream<char> *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<std::basic_ostringstream<char>> const &o){ return new Teuchos::RCP<std::basic_ostringstream<char>>(o); } ) );
		cl.def( pybind11::init<class std::basic_ostringstream<char> *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class std::basic_ostringstream<char> *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class std::basic_ostringstream<char> *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class std::basic_ostringstream<char> > & (Teuchos::RCP<std::basic_ostringstream<char>>::*)(const class Teuchos::RCP<class std::basic_ostringstream<char> > &)) &Teuchos::RCP<std::basic_ostringstream<char> >::operator=, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::operator=(const class Teuchos::RCP<class std::basic_ostringstream<char> > &) --> class Teuchos::RCP<class std::basic_ostringstream<char> > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class std::basic_ostringstream<char> > & (Teuchos::RCP<std::basic_ostringstream<char>>::*)(enum Teuchos::ENull)) &Teuchos::RCP<std::basic_ostringstream<char> >::operator=, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class std::basic_ostringstream<char> > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<std::basic_ostringstream<char>>::*)(class Teuchos::RCP<class std::basic_ostringstream<char> > &)) &Teuchos::RCP<std::basic_ostringstream<char> >::swap, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::swap(class Teuchos::RCP<class std::basic_ostringstream<char> > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::is_null, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::is_null() const --> bool");
		cl.def("arrow", (class std::basic_ostringstream<char> * (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::operator->, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::operator->() const --> class std::basic_ostringstream<char> *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class std::basic_ostringstream<char> & (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::operator*, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::operator*() const --> class std::basic_ostringstream<char> &", pybind11::return_value_policy::automatic);
		cl.def("get", (class std::basic_ostringstream<char> * (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::get, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::get() const --> class std::basic_ostringstream<char> *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class std::basic_ostringstream<char> * (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::getRawPtr, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::getRawPtr() const --> class std::basic_ostringstream<char> *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::strength, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::is_valid_ptr, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::strong_count, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::weak_count, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::total_count, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<std::basic_ostringstream<char>>::*)()) &Teuchos::RCP<std::basic_ostringstream<char> >::set_has_ownership, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::has_ownership, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class std::basic_ostringstream<char> > (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::create_weak, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::create_weak() const --> class Teuchos::RCP<class std::basic_ostringstream<char> >");
		cl.def("create_strong", (class Teuchos::RCP<class std::basic_ostringstream<char> > (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::create_strong, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::create_strong() const --> class Teuchos::RCP<class std::basic_ostringstream<char> >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class std::basic_ostringstream<char> > & (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::assert_not_null, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::assert_not_null() const --> const class Teuchos::RCP<class std::basic_ostringstream<char> > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class std::basic_ostringstream<char> > & (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::assert_valid_ptr, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::assert_valid_ptr() const --> const class Teuchos::RCP<class std::basic_ostringstream<char> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class std::basic_ostringstream<char> > & (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::debug_assert_not_null, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::debug_assert_not_null() const --> const class Teuchos::RCP<class std::basic_ostringstream<char> > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class std::basic_ostringstream<char> > & (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::debug_assert_valid_ptr, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class std::basic_ostringstream<char> > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<std::basic_ostringstream<char>>::*)()) &Teuchos::RCP<std::basic_ostringstream<char> >::reset, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::reset() --> void");
		cl.def("access_private_ptr", (class std::basic_ostringstream<char> * (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::access_private_ptr, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::access_private_ptr() const --> class std::basic_ostringstream<char> *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<std::basic_ostringstream<char>>::*)()) &Teuchos::RCP<std::basic_ostringstream<char> >::nonconst_access_private_node, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<std::basic_ostringstream<char>>::*)() const) &Teuchos::RCP<std::basic_ostringstream<char> >::access_private_node, "C++: Teuchos::RCP<std::basic_ostringstream<char> >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<std::basic_ostringstream<char>> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>, std::shared_ptr<Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>>> cl(M("Teuchos"), "RCP_Teuchos_basic_FancyOStream_char_std_char_traits_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > * a0){ return new Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >> const &o){ return new Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>(o); } ) );
		cl.def( pybind11::init<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &)) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::operator=, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::operator=(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &) --> class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)(enum Teuchos::ENull)) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::operator=, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)(class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &)) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::swap, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::swap(class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::is_null, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::is_null() const --> bool");
		cl.def("arrow", (class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > * (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::operator->, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::operator->() const --> class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::operator*, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::operator*() const --> class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > * (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::get, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::get() const --> class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > * (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::getRawPtr, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::getRawPtr() const --> class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::strength, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::is_valid_ptr, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::strong_count, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::weak_count, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::total_count, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)()) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::set_has_ownership, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::has_ownership, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::create_weak, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::create_weak() const --> class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > >");
		cl.def("create_strong", (class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::create_strong, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::create_strong() const --> class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::assert_not_null, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::assert_not_null() const --> const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::debug_assert_not_null, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::debug_assert_not_null() const --> const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::debug_assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)()) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::reset, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::reset() --> void");
		cl.def("access_private_ptr", (class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > * (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::access_private_ptr, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::access_private_ptr() const --> class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)()) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::nonconst_access_private_node, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::access_private_node, "C++: Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> > >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<Teuchos::basic_FancyOStream<char, std::char_traits<char> >> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>, std::shared_ptr<Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>>> cl(M("Teuchos"), "RCP_Teuchos_basic_oblackholestream_char_std_char_traits_char_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > * a0){ return new Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >> const &o){ return new Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>(o); } ) );
		cl.def( pybind11::init<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)(const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &)) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::operator=, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::operator=(const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &) --> class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)(enum Teuchos::ENull)) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::operator=, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)(class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &)) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::swap, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::swap(class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::is_null, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::is_null() const --> bool");
		cl.def("arrow", (class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > * (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::operator->, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::operator->() const --> class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::operator*, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::operator*() const --> class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > * (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::get, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::get() const --> class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > * (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::getRawPtr, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::getRawPtr() const --> class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::strength, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::is_valid_ptr, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::strong_count, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::weak_count, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::total_count, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)()) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::set_has_ownership, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::has_ownership, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::create_weak, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::create_weak() const --> class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > >");
		cl.def("create_strong", (class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::create_strong, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::create_strong() const --> class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > >");
		cl.def("assert_not_null", (const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::assert_not_null, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::assert_not_null() const --> const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::debug_assert_not_null, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::debug_assert_not_null() const --> const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::debug_assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)()) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::reset, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::reset() --> void");
		cl.def("access_private_ptr", (class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > * (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::access_private_ptr, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::access_private_ptr() const --> class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)()) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::nonconst_access_private_node, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>::*)() const) &Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::access_private_node, "C++: Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> > >::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<Teuchos::basic_oblackholestream<char, std::char_traits<char> >> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<Teuchos::XMLObjectImplem>, std::shared_ptr<Teuchos::RCP<Teuchos::XMLObjectImplem>>> cl(M("Teuchos"), "RCP_Teuchos_XMLObjectImplem_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<Teuchos::XMLObjectImplem>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](class Teuchos::XMLObjectImplem * a0){ return new Teuchos::RCP<Teuchos::XMLObjectImplem>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class Teuchos::XMLObjectImplem *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<Teuchos::XMLObjectImplem> const &o){ return new Teuchos::RCP<Teuchos::XMLObjectImplem>(o); } ) );
		cl.def( pybind11::init<class Teuchos::XMLObjectImplem *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::XMLObjectImplem *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::XMLObjectImplem *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class Teuchos::XMLObjectImplem> & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)(const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &)) &Teuchos::RCP<Teuchos::XMLObjectImplem>::operator=, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::operator=(const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &) --> class Teuchos::RCP<class Teuchos::XMLObjectImplem> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class Teuchos::XMLObjectImplem> & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)(enum Teuchos::ENull)) &Teuchos::RCP<Teuchos::XMLObjectImplem>::operator=, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class Teuchos::XMLObjectImplem> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)(class Teuchos::RCP<class Teuchos::XMLObjectImplem> &)) &Teuchos::RCP<Teuchos::XMLObjectImplem>::swap, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::swap(class Teuchos::RCP<class Teuchos::XMLObjectImplem> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::is_null, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::is_null() const --> bool");
		cl.def("arrow", (class Teuchos::XMLObjectImplem * (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::operator->, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::operator->() const --> class Teuchos::XMLObjectImplem *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::XMLObjectImplem & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::operator*, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::operator*() const --> class Teuchos::XMLObjectImplem &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::XMLObjectImplem * (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::get, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::get() const --> class Teuchos::XMLObjectImplem *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::XMLObjectImplem * (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::getRawPtr, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::getRawPtr() const --> class Teuchos::XMLObjectImplem *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::strength, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::is_valid_ptr, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::strong_count, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::weak_count, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::total_count, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)()) &Teuchos::RCP<Teuchos::XMLObjectImplem>::set_has_ownership, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::has_ownership, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class Teuchos::XMLObjectImplem> (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::create_weak, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::create_weak() const --> class Teuchos::RCP<class Teuchos::XMLObjectImplem>");
		cl.def("create_strong", (class Teuchos::RCP<class Teuchos::XMLObjectImplem> (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::create_strong, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::create_strong() const --> class Teuchos::RCP<class Teuchos::XMLObjectImplem>");
		cl.def("assert_not_null", (const class Teuchos::RCP<class Teuchos::XMLObjectImplem> & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::assert_not_null, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::assert_not_null() const --> const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::XMLObjectImplem> & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class Teuchos::XMLObjectImplem> & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::debug_assert_not_null, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::debug_assert_not_null() const --> const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::XMLObjectImplem> & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::debug_assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)()) &Teuchos::RCP<Teuchos::XMLObjectImplem>::reset, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::reset() --> void");
		cl.def("access_private_ptr", (class Teuchos::XMLObjectImplem * (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::access_private_ptr, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::access_private_ptr() const --> class Teuchos::XMLObjectImplem *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)()) &Teuchos::RCP<Teuchos::XMLObjectImplem>::nonconst_access_private_node, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::XMLObjectImplem>::*)() const) &Teuchos::RCP<Teuchos::XMLObjectImplem>::access_private_node, "C++: Teuchos::RCP<Teuchos::XMLObjectImplem>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<Teuchos::XMLObjectImplem> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<const Teuchos::ParameterEntryValidator>, std::shared_ptr<Teuchos::RCP<const Teuchos::ParameterEntryValidator>>> cl(M("Teuchos"), "RCP_const_Teuchos_ParameterEntryValidator_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<const Teuchos::ParameterEntryValidator>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](const class Teuchos::ParameterEntryValidator * a0){ return new Teuchos::RCP<const Teuchos::ParameterEntryValidator>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<const class Teuchos::ParameterEntryValidator *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<const Teuchos::ParameterEntryValidator> const &o){ return new Teuchos::RCP<const Teuchos::ParameterEntryValidator>(o); } ) );
		cl.def( pybind11::init<const class Teuchos::ParameterEntryValidator *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::ParameterEntryValidator *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::ParameterEntryValidator *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &)) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::operator=, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::operator=(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &) --> class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)(enum Teuchos::ENull)) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::operator=, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)(class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &)) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::swap, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::swap(class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::is_null, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::is_null() const --> bool");
		cl.def("arrow", (const class Teuchos::ParameterEntryValidator * (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::operator->, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::operator->() const --> const class Teuchos::ParameterEntryValidator *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (const class Teuchos::ParameterEntryValidator & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::operator*, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::operator*() const --> const class Teuchos::ParameterEntryValidator &", pybind11::return_value_policy::automatic);
		cl.def("get", (const class Teuchos::ParameterEntryValidator * (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::get, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::get() const --> const class Teuchos::ParameterEntryValidator *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (const class Teuchos::ParameterEntryValidator * (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::getRawPtr, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::getRawPtr() const --> const class Teuchos::ParameterEntryValidator *", pybind11::return_value_policy::automatic);
		cl.def("getConst", (class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::getConst, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::getConst() const --> class Teuchos::RCP<const class Teuchos::ParameterEntryValidator>");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::strength, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::is_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::strong_count, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::weak_count, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::total_count, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)()) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::set_has_ownership, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::has_ownership, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::create_weak, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::create_weak() const --> class Teuchos::RCP<const class Teuchos::ParameterEntryValidator>");
		cl.def("create_strong", (class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::create_strong, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::create_strong() const --> class Teuchos::RCP<const class Teuchos::ParameterEntryValidator>");
		cl.def("assert_not_null", (const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::assert_not_null, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::debug_assert_not_null, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::debug_assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::debug_assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)()) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::reset, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::reset() --> void");
		cl.def("access_private_ptr", (const class Teuchos::ParameterEntryValidator * (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::access_private_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::access_private_ptr() const --> const class Teuchos::ParameterEntryValidator *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)()) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::nonconst_access_private_node, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::ParameterEntryValidator>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntryValidator>::access_private_node, "C++: Teuchos::RCP<const Teuchos::ParameterEntryValidator>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<const Teuchos::ParameterEntryValidator> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<const Teuchos::ParameterEntry>, std::shared_ptr<Teuchos::RCP<const Teuchos::ParameterEntry>>> cl(M("Teuchos"), "RCP_const_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<const Teuchos::ParameterEntry>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](const class Teuchos::ParameterEntry * a0){ return new Teuchos::RCP<const Teuchos::ParameterEntry>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<const class Teuchos::ParameterEntry *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<const Teuchos::ParameterEntry> const &o){ return new Teuchos::RCP<const Teuchos::ParameterEntry>(o); } ) );
		cl.def( pybind11::init<const class Teuchos::ParameterEntry *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::ParameterEntry *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::ParameterEntry *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<const class Teuchos::ParameterEntry> & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)(const class Teuchos::RCP<const class Teuchos::ParameterEntry> &)) &Teuchos::RCP<const Teuchos::ParameterEntry>::operator=, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::operator=(const class Teuchos::RCP<const class Teuchos::ParameterEntry> &) --> class Teuchos::RCP<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<const class Teuchos::ParameterEntry> & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)(enum Teuchos::ENull)) &Teuchos::RCP<const Teuchos::ParameterEntry>::operator=, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<const Teuchos::ParameterEntry>::*)(class Teuchos::RCP<const class Teuchos::ParameterEntry> &)) &Teuchos::RCP<const Teuchos::ParameterEntry>::swap, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::swap(class Teuchos::RCP<const class Teuchos::ParameterEntry> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::is_null, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::is_null() const --> bool");
		cl.def("arrow", (const class Teuchos::ParameterEntry * (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::operator->, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::operator->() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (const class Teuchos::ParameterEntry & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::operator*, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::operator*() const --> const class Teuchos::ParameterEntry &", pybind11::return_value_policy::automatic);
		cl.def("get", (const class Teuchos::ParameterEntry * (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::get, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::get() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (const class Teuchos::ParameterEntry * (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::getRawPtr, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::getRawPtr() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("ptr", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::ptr() const --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>");
		cl.def("__call__", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::operator(), "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::operator()() const --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>");
		cl.def("getConst", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::getConst, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::getConst() const --> class Teuchos::RCP<const class Teuchos::ParameterEntry>");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::strength, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::is_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::strong_count, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::weak_count, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::total_count, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<const Teuchos::ParameterEntry>::*)()) &Teuchos::RCP<const Teuchos::ParameterEntry>::set_has_ownership, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::has_ownership, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::has_ownership() const --> bool");
		cl.def("release", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (Teuchos::RCP<const Teuchos::ParameterEntry>::*)()) &Teuchos::RCP<const Teuchos::ParameterEntry>::release, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::release() --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>");
		cl.def("create_weak", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::create_weak, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::create_weak() const --> class Teuchos::RCP<const class Teuchos::ParameterEntry>");
		cl.def("create_strong", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::create_strong, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::create_strong() const --> class Teuchos::RCP<const class Teuchos::ParameterEntry>");
		cl.def("assert_not_null", (const class Teuchos::RCP<const class Teuchos::ParameterEntry> & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::assert_not_null, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::ParameterEntry> & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<const class Teuchos::ParameterEntry> & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::debug_assert_not_null, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::debug_assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::ParameterEntry> & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::debug_assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<const Teuchos::ParameterEntry>::*)()) &Teuchos::RCP<const Teuchos::ParameterEntry>::reset, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::reset() --> void");
		cl.def("access_private_ptr", (const class Teuchos::ParameterEntry * (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::access_private_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::access_private_ptr() const --> const class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)()) &Teuchos::RCP<const Teuchos::ParameterEntry>::nonconst_access_private_node, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<const Teuchos::ParameterEntry>::access_private_node, "C++: Teuchos::RCP<const Teuchos::ParameterEntry>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<const Teuchos::ParameterEntry> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<const Teuchos::ParameterListModifier>, std::shared_ptr<Teuchos::RCP<const Teuchos::ParameterListModifier>>> cl(M("Teuchos"), "RCP_const_Teuchos_ParameterListModifier_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<const Teuchos::ParameterListModifier>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](const class Teuchos::ParameterListModifier * a0){ return new Teuchos::RCP<const Teuchos::ParameterListModifier>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<const class Teuchos::ParameterListModifier *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<const Teuchos::ParameterListModifier> const &o){ return new Teuchos::RCP<const Teuchos::ParameterListModifier>(o); } ) );
		cl.def( pybind11::init<const class Teuchos::ParameterListModifier *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::ParameterListModifier *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::ParameterListModifier *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<const class Teuchos::ParameterListModifier> & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)(const class Teuchos::RCP<const class Teuchos::ParameterListModifier> &)) &Teuchos::RCP<const Teuchos::ParameterListModifier>::operator=, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::operator=(const class Teuchos::RCP<const class Teuchos::ParameterListModifier> &) --> class Teuchos::RCP<const class Teuchos::ParameterListModifier> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<const class Teuchos::ParameterListModifier> & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)(enum Teuchos::ENull)) &Teuchos::RCP<const Teuchos::ParameterListModifier>::operator=, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<const class Teuchos::ParameterListModifier> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)(class Teuchos::RCP<const class Teuchos::ParameterListModifier> &)) &Teuchos::RCP<const Teuchos::ParameterListModifier>::swap, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::swap(class Teuchos::RCP<const class Teuchos::ParameterListModifier> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::is_null, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::is_null() const --> bool");
		cl.def("arrow", (const class Teuchos::ParameterListModifier * (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::operator->, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::operator->() const --> const class Teuchos::ParameterListModifier *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (const class Teuchos::ParameterListModifier & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::operator*, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::operator*() const --> const class Teuchos::ParameterListModifier &", pybind11::return_value_policy::automatic);
		cl.def("get", (const class Teuchos::ParameterListModifier * (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::get, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::get() const --> const class Teuchos::ParameterListModifier *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (const class Teuchos::ParameterListModifier * (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::getRawPtr, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::getRawPtr() const --> const class Teuchos::ParameterListModifier *", pybind11::return_value_policy::automatic);
		cl.def("getConst", (class Teuchos::RCP<const class Teuchos::ParameterListModifier> (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::getConst, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::getConst() const --> class Teuchos::RCP<const class Teuchos::ParameterListModifier>");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::strength, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::is_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::strong_count, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::weak_count, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::total_count, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)()) &Teuchos::RCP<const Teuchos::ParameterListModifier>::set_has_ownership, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::has_ownership, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<const class Teuchos::ParameterListModifier> (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::create_weak, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::create_weak() const --> class Teuchos::RCP<const class Teuchos::ParameterListModifier>");
		cl.def("create_strong", (class Teuchos::RCP<const class Teuchos::ParameterListModifier> (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::create_strong, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::create_strong() const --> class Teuchos::RCP<const class Teuchos::ParameterListModifier>");
		cl.def("assert_not_null", (const class Teuchos::RCP<const class Teuchos::ParameterListModifier> & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::assert_not_null, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::ParameterListModifier> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::ParameterListModifier> & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::ParameterListModifier> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<const class Teuchos::ParameterListModifier> & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::debug_assert_not_null, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::debug_assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::ParameterListModifier> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::ParameterListModifier> & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::debug_assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::ParameterListModifier> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)()) &Teuchos::RCP<const Teuchos::ParameterListModifier>::reset, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::reset() --> void");
		cl.def("access_private_ptr", (const class Teuchos::ParameterListModifier * (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::access_private_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::access_private_ptr() const --> const class Teuchos::ParameterListModifier *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)()) &Teuchos::RCP<const Teuchos::ParameterListModifier>::nonconst_access_private_node, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::ParameterListModifier>::*)() const) &Teuchos::RCP<const Teuchos::ParameterListModifier>::access_private_node, "C++: Teuchos::RCP<const Teuchos::ParameterListModifier>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<const Teuchos::ParameterListModifier> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<Teuchos::ParameterEntry>, std::shared_ptr<Teuchos::RCP<Teuchos::ParameterEntry>>> cl(M("Teuchos"), "RCP_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<Teuchos::ParameterEntry>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](class Teuchos::ParameterEntry * a0){ return new Teuchos::RCP<Teuchos::ParameterEntry>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class Teuchos::ParameterEntry *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<Teuchos::ParameterEntry> const &o){ return new Teuchos::RCP<Teuchos::ParameterEntry>(o); } ) );
		cl.def( pybind11::init<class Teuchos::ParameterEntry *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::ParameterEntry *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::ParameterEntry *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class Teuchos::ParameterEntry> & (Teuchos::RCP<Teuchos::ParameterEntry>::*)(const class Teuchos::RCP<class Teuchos::ParameterEntry> &)) &Teuchos::RCP<Teuchos::ParameterEntry>::operator=, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::operator=(const class Teuchos::RCP<class Teuchos::ParameterEntry> &) --> class Teuchos::RCP<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class Teuchos::ParameterEntry> & (Teuchos::RCP<Teuchos::ParameterEntry>::*)(enum Teuchos::ENull)) &Teuchos::RCP<Teuchos::ParameterEntry>::operator=, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<Teuchos::ParameterEntry>::*)(class Teuchos::RCP<class Teuchos::ParameterEntry> &)) &Teuchos::RCP<Teuchos::ParameterEntry>::swap, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::swap(class Teuchos::RCP<class Teuchos::ParameterEntry> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::is_null, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::is_null() const --> bool");
		cl.def("arrow", (class Teuchos::ParameterEntry * (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::operator->, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::operator->() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::ParameterEntry & (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::operator*, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::operator*() const --> class Teuchos::ParameterEntry &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::ParameterEntry * (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::get, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::get() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::ParameterEntry * (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::getRawPtr, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::getRawPtr() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("ptr", (class Teuchos::Ptr<class Teuchos::ParameterEntry> (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::ptr, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::ptr() const --> class Teuchos::Ptr<class Teuchos::ParameterEntry>");
		cl.def("__call__", (class Teuchos::Ptr<class Teuchos::ParameterEntry> (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::operator(), "C++: Teuchos::RCP<Teuchos::ParameterEntry>::operator()() const --> class Teuchos::Ptr<class Teuchos::ParameterEntry>");
		cl.def("getConst", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::getConst, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::getConst() const --> class Teuchos::RCP<const class Teuchos::ParameterEntry>");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::strength, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::is_valid_ptr, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::strong_count, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::weak_count, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::total_count, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<Teuchos::ParameterEntry>::*)()) &Teuchos::RCP<Teuchos::ParameterEntry>::set_has_ownership, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::has_ownership, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::has_ownership() const --> bool");
		cl.def("release", (class Teuchos::Ptr<class Teuchos::ParameterEntry> (Teuchos::RCP<Teuchos::ParameterEntry>::*)()) &Teuchos::RCP<Teuchos::ParameterEntry>::release, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::release() --> class Teuchos::Ptr<class Teuchos::ParameterEntry>");
		cl.def("create_weak", (class Teuchos::RCP<class Teuchos::ParameterEntry> (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::create_weak, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::create_weak() const --> class Teuchos::RCP<class Teuchos::ParameterEntry>");
		cl.def("create_strong", (class Teuchos::RCP<class Teuchos::ParameterEntry> (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::create_strong, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::create_strong() const --> class Teuchos::RCP<class Teuchos::ParameterEntry>");
		cl.def("assert_not_null", (const class Teuchos::RCP<class Teuchos::ParameterEntry> & (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::assert_not_null, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::assert_not_null() const --> const class Teuchos::RCP<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::ParameterEntry> & (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class Teuchos::ParameterEntry> & (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::debug_assert_not_null, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::debug_assert_not_null() const --> const class Teuchos::RCP<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::ParameterEntry> & (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::debug_assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::ParameterEntry> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<Teuchos::ParameterEntry>::*)()) &Teuchos::RCP<Teuchos::ParameterEntry>::reset, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::reset() --> void");
		cl.def("access_private_ptr", (class Teuchos::ParameterEntry * (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::access_private_ptr, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::access_private_ptr() const --> class Teuchos::ParameterEntry *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::ParameterEntry>::*)()) &Teuchos::RCP<Teuchos::ParameterEntry>::nonconst_access_private_node, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::ParameterEntry>::*)() const) &Teuchos::RCP<Teuchos::ParameterEntry>::access_private_node, "C++: Teuchos::RCP<Teuchos::ParameterEntry>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<Teuchos::ParameterEntry> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<Teuchos::ParameterList>, std::shared_ptr<Teuchos::RCP<Teuchos::ParameterList>>> cl(M("Teuchos"), "RCP_Teuchos_ParameterList_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<Teuchos::ParameterList>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](class Teuchos::ParameterList * a0){ return new Teuchos::RCP<Teuchos::ParameterList>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class Teuchos::ParameterList *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership_in") );

		cl.def( pybind11::init( [](Teuchos::RCP<Teuchos::ParameterList> const &o){ return new Teuchos::RCP<Teuchos::ParameterList>(o); } ) );
		cl.def( pybind11::init<class Teuchos::ParameterList *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::ParameterList *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::ParameterList *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class Teuchos::ParameterList> & (Teuchos::RCP<Teuchos::ParameterList>::*)(const class Teuchos::RCP<class Teuchos::ParameterList> &)) &Teuchos::RCP<Teuchos::ParameterList>::operator=, "C++: Teuchos::RCP<Teuchos::ParameterList>::operator=(const class Teuchos::RCP<class Teuchos::ParameterList> &) --> class Teuchos::RCP<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class Teuchos::ParameterList> & (Teuchos::RCP<Teuchos::ParameterList>::*)(enum Teuchos::ENull)) &Teuchos::RCP<Teuchos::ParameterList>::operator=, "C++: Teuchos::RCP<Teuchos::ParameterList>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<Teuchos::ParameterList>::*)(class Teuchos::RCP<class Teuchos::ParameterList> &)) &Teuchos::RCP<Teuchos::ParameterList>::swap, "C++: Teuchos::RCP<Teuchos::ParameterList>::swap(class Teuchos::RCP<class Teuchos::ParameterList> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::is_null, "C++: Teuchos::RCP<Teuchos::ParameterList>::is_null() const --> bool");
		cl.def("arrow", (class Teuchos::ParameterList * (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::operator->, "C++: Teuchos::RCP<Teuchos::ParameterList>::operator->() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::ParameterList & (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::operator*, "C++: Teuchos::RCP<Teuchos::ParameterList>::operator*() const --> class Teuchos::ParameterList &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::ParameterList * (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::get, "C++: Teuchos::RCP<Teuchos::ParameterList>::get() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::ParameterList * (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::getRawPtr, "C++: Teuchos::RCP<Teuchos::ParameterList>::getRawPtr() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("ptr", (class Teuchos::Ptr<class Teuchos::ParameterList> (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::ptr, "C++: Teuchos::RCP<Teuchos::ParameterList>::ptr() const --> class Teuchos::Ptr<class Teuchos::ParameterList>");
		cl.def("__call__", (class Teuchos::Ptr<class Teuchos::ParameterList> (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::operator(), "C++: Teuchos::RCP<Teuchos::ParameterList>::operator()() const --> class Teuchos::Ptr<class Teuchos::ParameterList>");
		cl.def("getConst", (class Teuchos::RCP<const class Teuchos::ParameterList> (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::getConst, "C++: Teuchos::RCP<Teuchos::ParameterList>::getConst() const --> class Teuchos::RCP<const class Teuchos::ParameterList>");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::strength, "C++: Teuchos::RCP<Teuchos::ParameterList>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::is_valid_ptr, "C++: Teuchos::RCP<Teuchos::ParameterList>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::strong_count, "C++: Teuchos::RCP<Teuchos::ParameterList>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::weak_count, "C++: Teuchos::RCP<Teuchos::ParameterList>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::total_count, "C++: Teuchos::RCP<Teuchos::ParameterList>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<Teuchos::ParameterList>::*)()) &Teuchos::RCP<Teuchos::ParameterList>::set_has_ownership, "C++: Teuchos::RCP<Teuchos::ParameterList>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::has_ownership, "C++: Teuchos::RCP<Teuchos::ParameterList>::has_ownership() const --> bool");
		cl.def("release", (class Teuchos::Ptr<class Teuchos::ParameterList> (Teuchos::RCP<Teuchos::ParameterList>::*)()) &Teuchos::RCP<Teuchos::ParameterList>::release, "C++: Teuchos::RCP<Teuchos::ParameterList>::release() --> class Teuchos::Ptr<class Teuchos::ParameterList>");
		cl.def("create_weak", (class Teuchos::RCP<class Teuchos::ParameterList> (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::create_weak, "C++: Teuchos::RCP<Teuchos::ParameterList>::create_weak() const --> class Teuchos::RCP<class Teuchos::ParameterList>");
		cl.def("create_strong", (class Teuchos::RCP<class Teuchos::ParameterList> (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::create_strong, "C++: Teuchos::RCP<Teuchos::ParameterList>::create_strong() const --> class Teuchos::RCP<class Teuchos::ParameterList>");
		cl.def("assert_not_null", (const class Teuchos::RCP<class Teuchos::ParameterList> & (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::assert_not_null, "C++: Teuchos::RCP<Teuchos::ParameterList>::assert_not_null() const --> const class Teuchos::RCP<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::ParameterList> & (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::ParameterList>::assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class Teuchos::ParameterList> & (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::debug_assert_not_null, "C++: Teuchos::RCP<Teuchos::ParameterList>::debug_assert_not_null() const --> const class Teuchos::RCP<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::ParameterList> & (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::debug_assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::ParameterList>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<Teuchos::ParameterList>::*)()) &Teuchos::RCP<Teuchos::ParameterList>::reset, "C++: Teuchos::RCP<Teuchos::ParameterList>::reset() --> void");
		cl.def("access_private_ptr", (class Teuchos::ParameterList * (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::access_private_ptr, "C++: Teuchos::RCP<Teuchos::ParameterList>::access_private_ptr() const --> class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::ParameterList>::*)()) &Teuchos::RCP<Teuchos::ParameterList>::nonconst_access_private_node, "C++: Teuchos::RCP<Teuchos::ParameterList>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::ParameterList>::*)() const) &Teuchos::RCP<Teuchos::ParameterList>::access_private_node, "C++: Teuchos::RCP<Teuchos::ParameterList>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<Teuchos::ParameterList> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<const Teuchos::ParameterList>, std::shared_ptr<Teuchos::RCP<const Teuchos::ParameterList>>> cl(M("Teuchos"), "RCP_const_Teuchos_ParameterList_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<const Teuchos::ParameterList>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](const class Teuchos::ParameterList * a0){ return new Teuchos::RCP<const Teuchos::ParameterList>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<const class Teuchos::ParameterList *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<const Teuchos::ParameterList> const &o){ return new Teuchos::RCP<const Teuchos::ParameterList>(o); } ) );
		cl.def( pybind11::init<const class Teuchos::ParameterList *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::ParameterList *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::ParameterList *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<const class Teuchos::ParameterList> & (Teuchos::RCP<const Teuchos::ParameterList>::*)(const class Teuchos::RCP<const class Teuchos::ParameterList> &)) &Teuchos::RCP<const Teuchos::ParameterList>::operator=, "C++: Teuchos::RCP<const Teuchos::ParameterList>::operator=(const class Teuchos::RCP<const class Teuchos::ParameterList> &) --> class Teuchos::RCP<const class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<const class Teuchos::ParameterList> & (Teuchos::RCP<const Teuchos::ParameterList>::*)(enum Teuchos::ENull)) &Teuchos::RCP<const Teuchos::ParameterList>::operator=, "C++: Teuchos::RCP<const Teuchos::ParameterList>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<const class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<const Teuchos::ParameterList>::*)(class Teuchos::RCP<const class Teuchos::ParameterList> &)) &Teuchos::RCP<const Teuchos::ParameterList>::swap, "C++: Teuchos::RCP<const Teuchos::ParameterList>::swap(class Teuchos::RCP<const class Teuchos::ParameterList> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::is_null, "C++: Teuchos::RCP<const Teuchos::ParameterList>::is_null() const --> bool");
		cl.def("arrow", (const class Teuchos::ParameterList * (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::operator->, "C++: Teuchos::RCP<const Teuchos::ParameterList>::operator->() const --> const class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (const class Teuchos::ParameterList & (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::operator*, "C++: Teuchos::RCP<const Teuchos::ParameterList>::operator*() const --> const class Teuchos::ParameterList &", pybind11::return_value_policy::automatic);
		cl.def("get", (const class Teuchos::ParameterList * (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::get, "C++: Teuchos::RCP<const Teuchos::ParameterList>::get() const --> const class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (const class Teuchos::ParameterList * (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::getRawPtr, "C++: Teuchos::RCP<const Teuchos::ParameterList>::getRawPtr() const --> const class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("getConst", (class Teuchos::RCP<const class Teuchos::ParameterList> (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::getConst, "C++: Teuchos::RCP<const Teuchos::ParameterList>::getConst() const --> class Teuchos::RCP<const class Teuchos::ParameterList>");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::strength, "C++: Teuchos::RCP<const Teuchos::ParameterList>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::is_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterList>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::strong_count, "C++: Teuchos::RCP<const Teuchos::ParameterList>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::weak_count, "C++: Teuchos::RCP<const Teuchos::ParameterList>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::total_count, "C++: Teuchos::RCP<const Teuchos::ParameterList>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<const Teuchos::ParameterList>::*)()) &Teuchos::RCP<const Teuchos::ParameterList>::set_has_ownership, "C++: Teuchos::RCP<const Teuchos::ParameterList>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::has_ownership, "C++: Teuchos::RCP<const Teuchos::ParameterList>::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<const class Teuchos::ParameterList> (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::create_weak, "C++: Teuchos::RCP<const Teuchos::ParameterList>::create_weak() const --> class Teuchos::RCP<const class Teuchos::ParameterList>");
		cl.def("create_strong", (class Teuchos::RCP<const class Teuchos::ParameterList> (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::create_strong, "C++: Teuchos::RCP<const Teuchos::ParameterList>::create_strong() const --> class Teuchos::RCP<const class Teuchos::ParameterList>");
		cl.def("assert_not_null", (const class Teuchos::RCP<const class Teuchos::ParameterList> & (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::assert_not_null, "C++: Teuchos::RCP<const Teuchos::ParameterList>::assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::ParameterList> & (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterList>::assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<const class Teuchos::ParameterList> & (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::debug_assert_not_null, "C++: Teuchos::RCP<const Teuchos::ParameterList>::debug_assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::ParameterList> & (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::debug_assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterList>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::ParameterList> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<const Teuchos::ParameterList>::*)()) &Teuchos::RCP<const Teuchos::ParameterList>::reset, "C++: Teuchos::RCP<const Teuchos::ParameterList>::reset() --> void");
		cl.def("access_private_ptr", (const class Teuchos::ParameterList * (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::access_private_ptr, "C++: Teuchos::RCP<const Teuchos::ParameterList>::access_private_ptr() const --> const class Teuchos::ParameterList *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::ParameterList>::*)()) &Teuchos::RCP<const Teuchos::ParameterList>::nonconst_access_private_node, "C++: Teuchos::RCP<const Teuchos::ParameterList>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::ParameterList>::*)() const) &Teuchos::RCP<const Teuchos::ParameterList>::access_private_node, "C++: Teuchos::RCP<const Teuchos::ParameterList>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<const Teuchos::ParameterList> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<Teuchos::Dependency>, std::shared_ptr<Teuchos::RCP<Teuchos::Dependency>>> cl(M("Teuchos"), "RCP_Teuchos_Dependency_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<Teuchos::Dependency>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("null_arg") );

		cl.def( pybind11::init( [](class Teuchos::Dependency * a0){ return new Teuchos::RCP<Teuchos::Dependency>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<class Teuchos::Dependency *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<Teuchos::Dependency> const &o){ return new Teuchos::RCP<Teuchos::Dependency>(o); } ) );
		cl.def( pybind11::init<class Teuchos::Dependency *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::Dependency *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<class Teuchos::Dependency *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<class Teuchos::Dependency> & (Teuchos::RCP<Teuchos::Dependency>::*)(const class Teuchos::RCP<class Teuchos::Dependency> &)) &Teuchos::RCP<Teuchos::Dependency>::operator=, "C++: Teuchos::RCP<Teuchos::Dependency>::operator=(const class Teuchos::RCP<class Teuchos::Dependency> &) --> class Teuchos::RCP<class Teuchos::Dependency> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<class Teuchos::Dependency> & (Teuchos::RCP<Teuchos::Dependency>::*)(enum Teuchos::ENull)) &Teuchos::RCP<Teuchos::Dependency>::operator=, "C++: Teuchos::RCP<Teuchos::Dependency>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<class Teuchos::Dependency> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<Teuchos::Dependency>::*)(class Teuchos::RCP<class Teuchos::Dependency> &)) &Teuchos::RCP<Teuchos::Dependency>::swap, "C++: Teuchos::RCP<Teuchos::Dependency>::swap(class Teuchos::RCP<class Teuchos::Dependency> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::is_null, "C++: Teuchos::RCP<Teuchos::Dependency>::is_null() const --> bool");
		cl.def("arrow", (class Teuchos::Dependency * (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::operator->, "C++: Teuchos::RCP<Teuchos::Dependency>::operator->() const --> class Teuchos::Dependency *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (class Teuchos::Dependency & (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::operator*, "C++: Teuchos::RCP<Teuchos::Dependency>::operator*() const --> class Teuchos::Dependency &", pybind11::return_value_policy::automatic);
		cl.def("get", (class Teuchos::Dependency * (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::get, "C++: Teuchos::RCP<Teuchos::Dependency>::get() const --> class Teuchos::Dependency *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (class Teuchos::Dependency * (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::getRawPtr, "C++: Teuchos::RCP<Teuchos::Dependency>::getRawPtr() const --> class Teuchos::Dependency *", pybind11::return_value_policy::automatic);
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::strength, "C++: Teuchos::RCP<Teuchos::Dependency>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::is_valid_ptr, "C++: Teuchos::RCP<Teuchos::Dependency>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::strong_count, "C++: Teuchos::RCP<Teuchos::Dependency>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::weak_count, "C++: Teuchos::RCP<Teuchos::Dependency>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::total_count, "C++: Teuchos::RCP<Teuchos::Dependency>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<Teuchos::Dependency>::*)()) &Teuchos::RCP<Teuchos::Dependency>::set_has_ownership, "C++: Teuchos::RCP<Teuchos::Dependency>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::has_ownership, "C++: Teuchos::RCP<Teuchos::Dependency>::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<class Teuchos::Dependency> (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::create_weak, "C++: Teuchos::RCP<Teuchos::Dependency>::create_weak() const --> class Teuchos::RCP<class Teuchos::Dependency>");
		cl.def("create_strong", (class Teuchos::RCP<class Teuchos::Dependency> (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::create_strong, "C++: Teuchos::RCP<Teuchos::Dependency>::create_strong() const --> class Teuchos::RCP<class Teuchos::Dependency>");
		cl.def("assert_not_null", (const class Teuchos::RCP<class Teuchos::Dependency> & (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::assert_not_null, "C++: Teuchos::RCP<Teuchos::Dependency>::assert_not_null() const --> const class Teuchos::RCP<class Teuchos::Dependency> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::Dependency> & (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::Dependency>::assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::Dependency> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<class Teuchos::Dependency> & (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::debug_assert_not_null, "C++: Teuchos::RCP<Teuchos::Dependency>::debug_assert_not_null() const --> const class Teuchos::RCP<class Teuchos::Dependency> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<class Teuchos::Dependency> & (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::debug_assert_valid_ptr, "C++: Teuchos::RCP<Teuchos::Dependency>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<class Teuchos::Dependency> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<Teuchos::Dependency>::*)()) &Teuchos::RCP<Teuchos::Dependency>::reset, "C++: Teuchos::RCP<Teuchos::Dependency>::reset() --> void");
		cl.def("access_private_ptr", (class Teuchos::Dependency * (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::access_private_ptr, "C++: Teuchos::RCP<Teuchos::Dependency>::access_private_ptr() const --> class Teuchos::Dependency *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::Dependency>::*)()) &Teuchos::RCP<Teuchos::Dependency>::nonconst_access_private_node, "C++: Teuchos::RCP<Teuchos::Dependency>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<Teuchos::Dependency>::*)() const) &Teuchos::RCP<Teuchos::Dependency>::access_private_node, "C++: Teuchos::RCP<Teuchos::Dependency>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<Teuchos::Dependency> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCP file:Teuchos_RCPDecl.hpp line:429
		pybind11::class_<Teuchos::RCP<const Teuchos::DependencySheet>, std::shared_ptr<Teuchos::RCP<const Teuchos::DependencySheet>>> cl(M("Teuchos"), "RCP_const_Teuchos_DependencySheet_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::RCP<const Teuchos::DependencySheet>(); } ), "doc" );
		cl.def( pybind11::init<enum Teuchos::ENull>(), pybind11::arg("") );

		cl.def( pybind11::init( [](const class Teuchos::DependencySheet * a0){ return new Teuchos::RCP<const Teuchos::DependencySheet>(a0); } ), "doc" , pybind11::arg("p"));
		cl.def( pybind11::init<const class Teuchos::DependencySheet *, bool>(), pybind11::arg("p"), pybind11::arg("has_ownership") );

		cl.def( pybind11::init( [](Teuchos::RCP<const Teuchos::DependencySheet> const &o){ return new Teuchos::RCP<const Teuchos::DependencySheet>(o); } ) );
		cl.def( pybind11::init<const class Teuchos::DependencySheet *, enum Teuchos::ERCPWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::DependencySheet *, enum Teuchos::ERCPUndefinedWeakNoDealloc>(), pybind11::arg("p"), pybind11::arg("") );

		cl.def( pybind11::init<const class Teuchos::DependencySheet *, const class Teuchos::RCPNodeHandle &>(), pybind11::arg("p"), pybind11::arg("node") );

		cl.def("assign", (class Teuchos::RCP<const class Teuchos::DependencySheet> & (Teuchos::RCP<const Teuchos::DependencySheet>::*)(const class Teuchos::RCP<const class Teuchos::DependencySheet> &)) &Teuchos::RCP<const Teuchos::DependencySheet>::operator=, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::operator=(const class Teuchos::RCP<const class Teuchos::DependencySheet> &) --> class Teuchos::RCP<const class Teuchos::DependencySheet> &", pybind11::return_value_policy::automatic, pybind11::arg("r_ptr"));
		cl.def("assign", (class Teuchos::RCP<const class Teuchos::DependencySheet> & (Teuchos::RCP<const Teuchos::DependencySheet>::*)(enum Teuchos::ENull)) &Teuchos::RCP<const Teuchos::DependencySheet>::operator=, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::operator=(enum Teuchos::ENull) --> class Teuchos::RCP<const class Teuchos::DependencySheet> &", pybind11::return_value_policy::automatic, pybind11::arg(""));
		cl.def("swap", (void (Teuchos::RCP<const Teuchos::DependencySheet>::*)(class Teuchos::RCP<const class Teuchos::DependencySheet> &)) &Teuchos::RCP<const Teuchos::DependencySheet>::swap, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::swap(class Teuchos::RCP<const class Teuchos::DependencySheet> &) --> void", pybind11::arg("r_ptr"));
		cl.def("is_null", (bool (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::is_null, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::is_null() const --> bool");
		cl.def("arrow", (const class Teuchos::DependencySheet * (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::operator->, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::operator->() const --> const class Teuchos::DependencySheet *", pybind11::return_value_policy::automatic);
		cl.def("dereference", (const class Teuchos::DependencySheet & (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::operator*, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::operator*() const --> const class Teuchos::DependencySheet &", pybind11::return_value_policy::automatic);
		cl.def("get", (const class Teuchos::DependencySheet * (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::get, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::get() const --> const class Teuchos::DependencySheet *", pybind11::return_value_policy::automatic);
		cl.def("getRawPtr", (const class Teuchos::DependencySheet * (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::getRawPtr, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::getRawPtr() const --> const class Teuchos::DependencySheet *", pybind11::return_value_policy::automatic);
		cl.def("getConst", (class Teuchos::RCP<const class Teuchos::DependencySheet> (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::getConst, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::getConst() const --> class Teuchos::RCP<const class Teuchos::DependencySheet>");
		cl.def("strength", (enum Teuchos::ERCPStrength (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::strength, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::strength() const --> enum Teuchos::ERCPStrength");
		cl.def("is_valid_ptr", (bool (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::is_valid_ptr, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::is_valid_ptr() const --> bool");
		cl.def("strong_count", (int (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::strong_count, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::strong_count() const --> int");
		cl.def("weak_count", (int (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::weak_count, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::weak_count() const --> int");
		cl.def("total_count", (int (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::total_count, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::total_count() const --> int");
		cl.def("set_has_ownership", (void (Teuchos::RCP<const Teuchos::DependencySheet>::*)()) &Teuchos::RCP<const Teuchos::DependencySheet>::set_has_ownership, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::set_has_ownership() --> void");
		cl.def("has_ownership", (bool (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::has_ownership, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::has_ownership() const --> bool");
		cl.def("create_weak", (class Teuchos::RCP<const class Teuchos::DependencySheet> (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::create_weak, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::create_weak() const --> class Teuchos::RCP<const class Teuchos::DependencySheet>");
		cl.def("create_strong", (class Teuchos::RCP<const class Teuchos::DependencySheet> (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::create_strong, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::create_strong() const --> class Teuchos::RCP<const class Teuchos::DependencySheet>");
		cl.def("assert_not_null", (const class Teuchos::RCP<const class Teuchos::DependencySheet> & (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::assert_not_null, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::DependencySheet> &", pybind11::return_value_policy::automatic);
		cl.def("assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::DependencySheet> & (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::DependencySheet> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_not_null", (const class Teuchos::RCP<const class Teuchos::DependencySheet> & (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::debug_assert_not_null, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::debug_assert_not_null() const --> const class Teuchos::RCP<const class Teuchos::DependencySheet> &", pybind11::return_value_policy::automatic);
		cl.def("debug_assert_valid_ptr", (const class Teuchos::RCP<const class Teuchos::DependencySheet> & (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::debug_assert_valid_ptr, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::debug_assert_valid_ptr() const --> const class Teuchos::RCP<const class Teuchos::DependencySheet> &", pybind11::return_value_policy::automatic);
		cl.def("reset", (void (Teuchos::RCP<const Teuchos::DependencySheet>::*)()) &Teuchos::RCP<const Teuchos::DependencySheet>::reset, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::reset() --> void");
		cl.def("access_private_ptr", (const class Teuchos::DependencySheet * (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::access_private_ptr, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::access_private_ptr() const --> const class Teuchos::DependencySheet *", pybind11::return_value_policy::automatic);
		cl.def("nonconst_access_private_node", (class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::DependencySheet>::*)()) &Teuchos::RCP<const Teuchos::DependencySheet>::nonconst_access_private_node, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::nonconst_access_private_node() --> class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);
		cl.def("access_private_node", (const class Teuchos::RCPNodeHandle & (Teuchos::RCP<const Teuchos::DependencySheet>::*)() const) &Teuchos::RCP<const Teuchos::DependencySheet>::access_private_node, "C++: Teuchos::RCP<const Teuchos::DependencySheet>::access_private_node() const --> const class Teuchos::RCPNodeHandle &", pybind11::return_value_policy::automatic);

		cl.def("__str__", [](Teuchos::RCP<const Teuchos::DependencySheet> const &o) -> std::string { std::ostringstream s; Teuchos::operator<<(s, o); return s.str(); } );
	}
	{ // Teuchos::RCPComp file:Teuchos_RCPDecl.hpp line:955
		pybind11::class_<Teuchos::RCPComp, std::shared_ptr<Teuchos::RCPComp>> cl(M("Teuchos"), "RCPComp", "Struct for comparing two RCPs. Simply compares\n the raw pointers contained within the RCPs", pybind11::module_local());
		cl.def( pybind11::init( [](Teuchos::RCPComp const &o){ return new Teuchos::RCPComp(o); } ) );
		cl.def( pybind11::init( [](){ return new Teuchos::RCPComp(); } ) );
	}
	{ // Teuchos::RCPConstComp file:Teuchos_RCPDecl.hpp line:965
		pybind11::class_<Teuchos::RCPConstComp, std::shared_ptr<Teuchos::RCPConstComp>> cl(M("Teuchos"), "RCPConstComp", "Struct for comparing two RCPs. Simply compares\n the raw pointers contained within the RCPs", pybind11::module_local());
		cl.def( pybind11::init( [](Teuchos::RCPConstComp const &o){ return new Teuchos::RCPConstComp(o); } ) );
		cl.def( pybind11::init( [](){ return new Teuchos::RCPConstComp(); } ) );
		cl.def("__call__", (bool (Teuchos::RCPConstComp::*)(const class Teuchos::RCP<const class Teuchos::ParameterEntry>, const class Teuchos::RCP<const class Teuchos::ParameterEntry>) const) &Teuchos::RCPConstComp::operator()<Teuchos::ParameterEntry,Teuchos::ParameterEntry>, "C++: Teuchos::RCPConstComp::operator()(const class Teuchos::RCP<const class Teuchos::ParameterEntry>, const class Teuchos::RCP<const class Teuchos::ParameterEntry>) const --> bool", pybind11::arg("p1"), pybind11::arg("p2"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<Teuchos::ParameterEntry>, std::shared_ptr<Teuchos::DeallocNull<Teuchos::ParameterEntry>>> cl(M("Teuchos"), "DeallocNull_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<Teuchos::ParameterEntry>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<Teuchos::ParameterEntry> const &o){ return new Teuchos::DeallocNull<Teuchos::ParameterEntry>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<Teuchos::ParameterEntry>::*)(class Teuchos::ParameterEntry *)) &Teuchos::DeallocNull<Teuchos::ParameterEntry>::free, "C++: Teuchos::DeallocNull<Teuchos::ParameterEntry>::free(class Teuchos::ParameterEntry *) --> void", pybind11::arg("ptr"));
	}
	{ // Teuchos::DeallocNull file:Teuchos_RCPDecl.hpp line:996
		pybind11::class_<Teuchos::DeallocNull<const Teuchos::ParameterEntry>, std::shared_ptr<Teuchos::DeallocNull<const Teuchos::ParameterEntry>>> cl(M("Teuchos"), "DeallocNull_const_Teuchos_ParameterEntry_t", "", pybind11::module_local());
		cl.def( pybind11::init( [](){ return new Teuchos::DeallocNull<const Teuchos::ParameterEntry>(); } ) );
		cl.def( pybind11::init( [](Teuchos::DeallocNull<const Teuchos::ParameterEntry> const &o){ return new Teuchos::DeallocNull<const Teuchos::ParameterEntry>(o); } ) );
		cl.def("free", (void (Teuchos::DeallocNull<const Teuchos::ParameterEntry>::*)(const class Teuchos::ParameterEntry *)) &Teuchos::DeallocNull<const Teuchos::ParameterEntry>::free, "C++: Teuchos::DeallocNull<const Teuchos::ParameterEntry>::free(const class Teuchos::ParameterEntry *) --> void", pybind11::arg("ptr"));
	}
	// Teuchos::rcp(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > * a0) -> Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > (*)(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool)) &Teuchos::rcp<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>, "C++: Teuchos::rcp(class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > *, bool) --> class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > * a0) -> Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > > (*)(class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, bool)) &Teuchos::rcp<Teuchos::basic_oblackholestream<char, std::char_traits<char> >>, "C++: Teuchos::rcp(class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > *, bool) --> class Teuchos::RCP<class Teuchos::basic_oblackholestream<char, struct std::char_traits<char> > >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class std::basic_ostringstream<char> *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class std::basic_ostringstream<char> * a0) -> Teuchos::RCP<class std::basic_ostringstream<char> > { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class std::basic_ostringstream<char> > (*)(class std::basic_ostringstream<char> *, bool)) &Teuchos::rcp<std::basic_ostringstream<char>>, "C++: Teuchos::rcp(class std::basic_ostringstream<char> *, bool) --> class Teuchos::RCP<class std::basic_ostringstream<char> >", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcp(class Teuchos::ParameterList *, bool) file:Teuchos_RCP.hpp line:622
	M("Teuchos").def("rcp", [](class Teuchos::ParameterList * a0) -> Teuchos::RCP<class Teuchos::ParameterList> { return Teuchos::rcp(a0); }, "", pybind11::arg("p"));
	M("Teuchos").def("rcp", (class Teuchos::RCP<class Teuchos::ParameterList> (*)(class Teuchos::ParameterList *, bool)) &Teuchos::rcp<Teuchos::ParameterList>, "C++: Teuchos::rcp(class Teuchos::ParameterList *, bool) --> class Teuchos::RCP<class Teuchos::ParameterList>", pybind11::arg("p"), pybind11::arg("owns_mem"));

	// Teuchos::rcpFromRef(class Teuchos::ParameterEntry &) file:Teuchos_RCPDecl.hpp line:1297
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<class Teuchos::ParameterEntry> (*)(class Teuchos::ParameterEntry &)) &Teuchos::rcpFromRef<Teuchos::ParameterEntry>, "C++: Teuchos::rcpFromRef(class Teuchos::ParameterEntry &) --> class Teuchos::RCP<class Teuchos::ParameterEntry>", pybind11::arg("r"));

	// Teuchos::rcpFromRef(const class Teuchos::ParameterEntry &) file:Teuchos_RCPDecl.hpp line:1297
	M("Teuchos").def("rcpFromRef", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (*)(const class Teuchos::ParameterEntry &)) &Teuchos::rcpFromRef<const Teuchos::ParameterEntry>, "C++: Teuchos::rcpFromRef(const class Teuchos::ParameterEntry &) --> class Teuchos::RCP<const class Teuchos::ParameterEntry>", pybind11::arg("r"));

	// Teuchos::rcpWithEmbeddedObjPostDestroy(class Teuchos::ParameterList *, const class Teuchos::RCP<class Teuchos::ParameterList> &, bool) file:Teuchos_RCP.hpp line:676
	M("Teuchos").def("rcpWithEmbeddedObjPostDestroy", [](class Teuchos::ParameterList * a0, const class Teuchos::RCP<class Teuchos::ParameterList> & a1) -> Teuchos::RCP<class Teuchos::ParameterList> { return Teuchos::rcpWithEmbeddedObjPostDestroy(a0, a1); }, "", pybind11::arg("p"), pybind11::arg("embedded"));
	M("Teuchos").def("rcpWithEmbeddedObjPostDestroy", (class Teuchos::RCP<class Teuchos::ParameterList> (*)(class Teuchos::ParameterList *, const class Teuchos::RCP<class Teuchos::ParameterList> &, bool)) &Teuchos::rcpWithEmbeddedObjPostDestroy<Teuchos::ParameterList,Teuchos::RCP<Teuchos::ParameterList>>, "C++: Teuchos::rcpWithEmbeddedObjPostDestroy(class Teuchos::ParameterList *, const class Teuchos::RCP<class Teuchos::ParameterList> &, bool) --> class Teuchos::RCP<class Teuchos::ParameterList>", pybind11::arg("p"), pybind11::arg("embedded"), pybind11::arg("owns_mem"));

	// Teuchos::rcpWithEmbeddedObjPostDestroy(const class Teuchos::ParameterList *, const class Teuchos::RCP<const class Teuchos::ParameterList> &, bool) file:Teuchos_RCP.hpp line:676
	M("Teuchos").def("rcpWithEmbeddedObjPostDestroy", [](const class Teuchos::ParameterList * a0, const class Teuchos::RCP<const class Teuchos::ParameterList> & a1) -> Teuchos::RCP<const class Teuchos::ParameterList> { return Teuchos::rcpWithEmbeddedObjPostDestroy(a0, a1); }, "", pybind11::arg("p"), pybind11::arg("embedded"));
	M("Teuchos").def("rcpWithEmbeddedObjPostDestroy", (class Teuchos::RCP<const class Teuchos::ParameterList> (*)(const class Teuchos::ParameterList *, const class Teuchos::RCP<const class Teuchos::ParameterList> &, bool)) &Teuchos::rcpWithEmbeddedObjPostDestroy<const Teuchos::ParameterList,Teuchos::RCP<const Teuchos::ParameterList>>, "C++: Teuchos::rcpWithEmbeddedObjPostDestroy(const class Teuchos::ParameterList *, const class Teuchos::RCP<const class Teuchos::ParameterList> &, bool) --> class Teuchos::RCP<const class Teuchos::ParameterList>", pybind11::arg("p"), pybind11::arg("embedded"), pybind11::arg("owns_mem"));

	// Teuchos::is_null(const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &) file:Teuchos_RCP.hpp line:715
	M("Teuchos").def("is_null", (bool (*)(const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &)) &Teuchos::is_null<Teuchos::XMLObjectImplem>, "C++: Teuchos::is_null(const class Teuchos::RCP<class Teuchos::XMLObjectImplem> &) --> bool", pybind11::arg("p"));

	// Teuchos::nonnull(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &) file:Teuchos_RCP.hpp line:723
	M("Teuchos").def("nonnull", (bool (*)(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &)) &Teuchos::nonnull<Teuchos::basic_FancyOStream<char, std::char_traits<char> >>, "C++: Teuchos::nonnull(const class Teuchos::RCP<class Teuchos::basic_FancyOStream<char, struct std::char_traits<char> > > &) --> bool", pybind11::arg("p"));

	// Teuchos::nonnull(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &) file:Teuchos_RCP.hpp line:723
	M("Teuchos").def("nonnull", (bool (*)(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &)) &Teuchos::nonnull<const Teuchos::ParameterEntryValidator>, "C++: Teuchos::nonnull(const class Teuchos::RCP<const class Teuchos::ParameterEntryValidator> &) --> bool", pybind11::arg("p"));

	{ // Teuchos::m_bad_cast file:Teuchos_dyn_cast.hpp line:60
		pybind11::class_<Teuchos::m_bad_cast, std::shared_ptr<Teuchos::m_bad_cast>, PyCallBack_Teuchos_m_bad_cast, std::bad_cast> cl(M("Teuchos"), "m_bad_cast", "Exception class for bad cast.\n\n\nWe create this class so that we may throw a bad_cast when appropriate and\nstill use the TEUCHOS_TEST_FOR_EXCEPTION macro.  We recommend users try to catch a\nbad_cast.", pybind11::module_local());
		cl.def( pybind11::init<const std::string &>(), pybind11::arg("what_arg") );

		cl.def( pybind11::init( [](PyCallBack_Teuchos_m_bad_cast const &o){ return new PyCallBack_Teuchos_m_bad_cast(o); } ) );
		cl.def( pybind11::init( [](Teuchos::m_bad_cast const &o){ return new Teuchos::m_bad_cast(o); } ) );
		cl.def("what", (const char * (Teuchos::m_bad_cast::*)() const) &Teuchos::m_bad_cast::what, "C++: Teuchos::m_bad_cast::what() const --> const char *", pybind11::return_value_policy::automatic);
		cl.def("assign", (class Teuchos::m_bad_cast & (Teuchos::m_bad_cast::*)(const class Teuchos::m_bad_cast &)) &Teuchos::m_bad_cast::operator=, "C++: Teuchos::m_bad_cast::operator=(const class Teuchos::m_bad_cast &) --> class Teuchos::m_bad_cast &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// Teuchos::dyn_cast_throw_exception(const std::string &, const std::string &, const std::string &) file:Teuchos_dyn_cast.hpp line:70
	M("Teuchos").def("dyn_cast_throw_exception", (void (*)(const std::string &, const std::string &, const std::string &)) &Teuchos::dyn_cast_throw_exception, "C++: Teuchos::dyn_cast_throw_exception(const std::string &, const std::string &, const std::string &) --> void", pybind11::arg("T_from"), pybind11::arg("T_from_concr"), pybind11::arg("T_to"));

	// Teuchos::ptrFromRef(class Teuchos::ParameterEntry &) file:Teuchos_PtrDecl.hpp line:298
	M("Teuchos").def("ptrFromRef", (class Teuchos::Ptr<class Teuchos::ParameterEntry> (*)(class Teuchos::ParameterEntry &)) &Teuchos::ptrFromRef<Teuchos::ParameterEntry>, "C++: Teuchos::ptrFromRef(class Teuchos::ParameterEntry &) --> class Teuchos::Ptr<class Teuchos::ParameterEntry>", pybind11::arg("arg"));

	// Teuchos::ptrFromRef(const class Teuchos::ParameterEntry &) file:Teuchos_PtrDecl.hpp line:298
	M("Teuchos").def("ptrFromRef", (class Teuchos::Ptr<const class Teuchos::ParameterEntry> (*)(const class Teuchos::ParameterEntry &)) &Teuchos::ptrFromRef<const Teuchos::ParameterEntry>, "C++: Teuchos::ptrFromRef(const class Teuchos::ParameterEntry &) --> class Teuchos::Ptr<const class Teuchos::ParameterEntry>", pybind11::arg("arg"));

	// Teuchos::rcpFromPtr(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) file:Teuchos_PtrDecl.hpp line:309
	M("Teuchos").def("rcpFromPtr", (class Teuchos::RCP<class Teuchos::ParameterEntry> (*)(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &)) &Teuchos::rcpFromPtr<Teuchos::ParameterEntry>, "C++: Teuchos::rcpFromPtr(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) --> class Teuchos::RCP<class Teuchos::ParameterEntry>", pybind11::arg("ptr"));

	// Teuchos::rcpFromPtr(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) file:Teuchos_PtrDecl.hpp line:309
	M("Teuchos").def("rcpFromPtr", (class Teuchos::RCP<const class Teuchos::ParameterEntry> (*)(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &)) &Teuchos::rcpFromPtr<const Teuchos::ParameterEntry>, "C++: Teuchos::rcpFromPtr(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) --> class Teuchos::RCP<const class Teuchos::ParameterEntry>", pybind11::arg("ptr"));

	// Teuchos::is_null(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) file:Teuchos_PtrDecl.hpp line:355
	M("Teuchos").def("is_null", (bool (*)(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &)) &Teuchos::is_null<Teuchos::ParameterEntry>, "C++: Teuchos::is_null(const class Teuchos::Ptr<class Teuchos::ParameterEntry> &) --> bool", pybind11::arg("p"));

	// Teuchos::is_null(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) file:Teuchos_PtrDecl.hpp line:355
	M("Teuchos").def("is_null", (bool (*)(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &)) &Teuchos::is_null<const Teuchos::ParameterEntry>, "C++: Teuchos::is_null(const class Teuchos::Ptr<const class Teuchos::ParameterEntry> &) --> bool", pybind11::arg("p"));

}
