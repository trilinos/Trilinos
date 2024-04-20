// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
// 

#include "stk_util/util/Marshal.hpp"

#include <stdint.h>                     // for uint32_t
#include <algorithm>                    // for fill
#include <cstring>                      // for strlen
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for runtime_error
#include <utility>                      // for make_pair, pair


namespace stk {

namespace {

typedef std::map<uint32_t, std::string> TypeNameMap;

bool            s_crcInitialized = false;               ///< CRC table has been initialized
uint32_t        s_crcTable[256];                        ///< CRC lookup table
TypeNameMap     s_typeNameMap;                          ///< CRC'd type name map

uint32_t
crc32_reflect(
  uint32_t      seed,
  const char    c)
{
  uint32_t value = 0;

  // Swap bit 0 for bit 7, bit 1 For bit 6, etc....
  for(int i = 1; i < (c + 1); i++) {
    if (seed & 1) 
      value |= (1 << (c - i));
    seed >>= 1;
  }

  return value;
}


void
crc32_initialize()
{
  s_crcInitialized = true;
  
  //0x04C11DB7 is the official polynomial used by PKZip, WinZip and Ethernet.
  uint32_t polynomial = 0x04C11DB7;

  std::fill(s_crcTable, s_crcTable + 256, 0);
  
  // 256 values representing ASCII character codes.
  for (int i = 0; i <= 0xFF; i++) {
    s_crcTable[i] = crc32_reflect(i, 8) << 24;

    for (int j = 0; j < 8; j++)
      s_crcTable[i] = (s_crcTable[i] << 1) ^ ((s_crcTable[i] & (1u << 31)) ? polynomial : 0);
    
    s_crcTable[i] = crc32_reflect(s_crcTable[i], 32);
  }
}


void
crc32_part(
  uint32_t &            crc,
  const unsigned char * s,
  unsigned              l)
{
  while (l--)
    crc = (crc >> 8)^s_crcTable[(crc & 0xFF)^*s++];
}


uint32_t
crc32(
  const unsigned char * s,
  unsigned              l)
{
  uint32_t crc = 0xffffffff;
  crc32_part(crc, s, l);
  return crc ^ 0xffffffff;
}


void
crc32_write(
  std::stringstream &           os,
  const std::type_info &        typeinfo)
{
  if (!s_crcInitialized)
    crc32_initialize();
  
  const char *name = typeinfo.name();
  unsigned length = std::strlen(name);
  
  uint32_t      crc = crc32(reinterpret_cast<const unsigned char *>(name), length);

  TypeNameMap::iterator it = s_typeNameMap.find(crc);
  if (it == s_typeNameMap.end())
    s_typeNameMap.emplace(crc, std::string(name));
                         
  os.write(reinterpret_cast<char*>(&crc), sizeof(uint32_t));
}


void
crc32_check(
  std::stringstream &           is,
  const std::type_info &        typeinfo)
{
  if (!s_crcInitialized)
    crc32_initialize();
  
  const char *name = typeinfo.name();
  unsigned length = std::strlen(name);
  
  uint32_t      crc_check = crc32(reinterpret_cast<const unsigned char *>(name), length);
  uint32_t      crc = 0;
  
  {
    TypeNameMap::iterator it = s_typeNameMap.find(crc_check);
    if (it == s_typeNameMap.end())
      s_typeNameMap.emplace(crc_check, std::string(name));
  }
                         
  is.read(reinterpret_cast<char *>(&crc), sizeof(uint32_t));

  if (crc_check != crc) {
    std::ostringstream ss;
    ss << "Marshaller encountered type ";
    TypeNameMap::const_iterator it = s_typeNameMap.find(crc);
    if (it == s_typeNameMap.end())
      ss << "code " << std::hex << crc;
    else
      ss << (*it).second;
    
    ss << " when expecting type ";
    it = s_typeNameMap.find(crc_check);
    if (it == s_typeNameMap.end())
      ss << "code " << std::hex << crc_check;
    else
      ss << (*it).second;
   
    throw std::runtime_error(ss.str());
  }
}

} // namespace <empty>

Marshal::Marshal(
  unsigned              type_check)
  : stream(std::ios_base::out),
    m_typeCheck(TYPE_CHECK_NONE)
{
  (*this) << type_check;
  m_typeCheck = type_check;
}


Marshal::Marshal(
  const std::string &   s)
  : stream(s, std::ios_base::in),
    m_typeCheck(TYPE_CHECK_NONE)
{
  (*this) >> m_typeCheck;
}


std::string
Marshal::str() const
{
  return stream.str();
}


size_t
Marshal::size() const
{
  return stream.str().size();
}


Marshal::operator void * () const
{
  return static_cast<void *>(const_cast<std::stringstream*>(&stream));
}


void
Marshal::write(
  const char *          address, 
  size_t                byte_count)
{
  stream.write(address, byte_count);
}

  
void
Marshal::read(
  char *                address,
  size_t                byte_count)
{
  stream.read(address, byte_count);
}


template<>
Marshal &operator<<(Marshal &mout, const std::type_info &t) {
  crc32_write(mout.stream, t);
  return mout;
}

template<>
Marshal &operator>>(Marshal &min, const std::type_info &t) {
  crc32_check(min.stream, t);
  return min;
}

template<>
Marshal &operator<<(Marshal &mout, const signed char &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned char &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const char &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const short &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned short &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const int &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned int &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const long &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned long &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const long long &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const unsigned long long &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const float &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const double &t) {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(t);
  return write(mout, t);
}

template<>
Marshal &operator<<(Marshal &mout, const std::string &s)  {
  if (mout.m_typeCheck & Marshal::TYPE_CHECK_POD)
    mout << typeid(s);

  size_t ul = s.size();
  mout << ul;
  mout.stream.write(s.data(), s.size());
  return mout;
}


template<>
Marshal &operator>>(Marshal &min, signed char &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned char &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, char &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, short &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned short &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, int &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned int &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, long &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned long &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, long long &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, unsigned long long &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, float &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, double &t) {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(t);
  return read(min, t);
}

template<>
Marshal &operator>>(Marshal &min, std::string &s)  {
  if (min.m_typeCheck & Marshal::TYPE_CHECK_POD)
    min >> typeid(s);

  size_t size = 0;
  min >> size;
  std::vector<char> c(size);

  min.stream.read(c.data(), size);
  s.assign(c.data(), size);

  return min;
}

} // namespace stk

