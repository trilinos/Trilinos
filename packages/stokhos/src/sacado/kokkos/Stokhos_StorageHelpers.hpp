// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_STORAGE_HELPERS_HPP
#define STOKHOS_STORAGE_HELPERS_HPP

#define STOKHOS_STORAGE_HELPER_STRINGNAME_DYNAMIC(__storagename__)     \
  namespace Sacado                                                     \
  {                                                                    \
    template <typename ordinal_t, typename value_t, typename device_t> \
    struct StringName<Stokhos::__storagename__<ordinal_t,              \
                                               value_t,                \
                                               device_t>>              \
    {                                                                  \
      static std::string eval()                                        \
      {                                                                \
        std::stringstream ss;                                          \
        ss << "Stokhos::" #__storagename__ "<"                         \
           << StringName<ordinal_t>::eval() << ","                     \
           << StringName<value_t>::eval() << ","                       \
           << StringName<device_t>::eval() << ">";                     \
        return ss.str();                                               \
      }                                                                \
    };                                                                 \
  }

#define STOKHOS_STORAGE_HELPER_STRINGNAME_STATIC(__storagename__)                   \
    namespace Sacado                                                                \
    {                                                                               \
        template <typename ordinal_t, typename value_t, int Num, typename device_t> \
        struct StringName<Stokhos::__storagename__<ordinal_t,                       \
                                                   value_t,                         \
                                                   Num,                             \
                                                   device_t>>                       \
        {                                                                           \
            static std::string eval()                                               \
            {                                                                       \
                std::stringstream ss;                                               \
                ss << "Stokhos::" #__storagename__ "<"                              \
                   << StringName<ordinal_t>::eval() << ","                          \
                   << StringName<value_t>::eval() << ","                            \
                   << Num << ","                                                    \
                   << StringName<device_t>::eval() << ">";                          \
                return ss.str();                                                    \
            }                                                                       \
        };                                                                          \
    }

#endif // STOKHOS_STORAGE_HELPERS_HPP
