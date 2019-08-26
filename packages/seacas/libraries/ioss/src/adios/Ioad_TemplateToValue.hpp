// Copyright(C) 1999-2010 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef IOSS_Ioad_TemplateToValue_hpp
#define IOSS_Ioad_TemplateToValue_hpp

#include "Ioss_Field.h" // for Field, etc

namespace Ioad {

template <typename T>
constexpr Ioss::Field::BasicType template_to_basic_type() noexcept
{
    return Ioss::Field::BasicType::INVALID;
}

template <>
constexpr Ioss::Field::BasicType template_to_basic_type<double>() noexcept
{
    return Ioss::Field::BasicType::DOUBLE;
}

template <>
constexpr Ioss::Field::BasicType template_to_basic_type<int32_t>() noexcept
{
    return Ioss::Field::BasicType::INT32;
}

template <>
constexpr Ioss::Field::BasicType template_to_basic_type<int64_t>() noexcept
{
    return Ioss::Field::BasicType::INT64;
}

template <>
constexpr Ioss::Field::BasicType template_to_basic_type<Complex>() noexcept
{
    return Ioss::Field::BasicType::COMPLEX;
}

template <>
constexpr Ioss::Field::BasicType template_to_basic_type<std::string>() noexcept
{
    return Ioss::Field::BasicType::STRING;
}

template <>
constexpr Ioss::Field::BasicType template_to_basic_type<char>() noexcept
{
    return Ioss::Field::BasicType::CHARACTER;
}

template <>
constexpr char const *get_entity_type<Ioss::SideBlock>() noexcept
{
    return "SideBlock";
}

template <>
constexpr char const *get_entity_type<Ioss::SideSet>() noexcept
{
    return "SideSet";
}

template <>
constexpr char const *get_entity_type<Ioss::NodeBlock>() noexcept
{
    return "NodeBlock";
}

template <>
constexpr char const *get_entity_type<Ioss::EdgeBlock>() noexcept
{
    return "EdgeBlock";
}

template <>
constexpr char const *get_entity_type<Ioss::FaceBlock>() noexcept
{
    return "FaceBlock";
}

template <>
constexpr char const *get_entity_type<Ioss::ElementBlock>() noexcept
{
    return "ElementBlock";
}

template <>
constexpr char const *get_entity_type<Ioss::NodeSet>() noexcept
{
    return "NodeSet";
}

template <>
constexpr char const *get_entity_type<Ioss::EdgeSet>() noexcept
{
    return "EdgeSet";
}

template <>
constexpr char const *get_entity_type<Ioss::FaceSet>() noexcept
{
    return "FaceSet";
}

template <>
constexpr char const *get_entity_type<Ioss::ElementSet>() noexcept
{
    return "ElementSet";
}

template <>
constexpr char const *get_entity_type<Ioss::CommSet>() noexcept
{
    return "CommSet";
}

template <class T>
inline std::string GetType() noexcept
{
    return "compound";
}
template <>
inline std::string GetType<void>() noexcept
{
    return "unknown";
}

template <>
inline std::string GetType<std::string>() noexcept
{
    return "string";
}

template <>
inline std::string GetType<char>() noexcept
{
    return "char";
}
template <>
inline std::string GetType<signed char>() noexcept
{
    return "signed char";
}
template <>
inline std::string GetType<unsigned char>() noexcept
{
    return "unsigned char";
}
template <>
inline std::string GetType<short>() noexcept
{
    return "short";
}
template <>
inline std::string GetType<unsigned short>() noexcept
{
    return "unsigned short";
}
template <>
inline std::string GetType<int>() noexcept
{
    return "int";
}
template <>
inline std::string GetType<unsigned int>() noexcept
{
    return "unsigned int";
}
template <>
inline std::string GetType<long int>() noexcept
{
    return "long int";
}
template <>
inline std::string GetType<unsigned long int>() noexcept
{
    return "unsigned long int";
}
template <>
inline std::string GetType<long long int>() noexcept
{
    return "long long int";
}
template <>
inline std::string GetType<unsigned long long int>() noexcept
{
    return "unsigned long long int";
}
template <>
inline std::string GetType<float>() noexcept
{
    return "float";
}
template <>
inline std::string GetType<double>() noexcept
{
    return "double";
}
template <>
inline std::string GetType<long double>() noexcept
{
    return "long double";
}
template <>
inline std::string GetType<std::complex<float>>() noexcept
{
    return "float complex";
}
template <>
inline std::string GetType<std::complex<double>>() noexcept
{
    return "double complex";
}

} // namespace Ioad

#endif
