// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
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

#ifndef PACKAGES_STK_STK_TOOLS_STK_TOOLS_BLOCK_EXTRACTOR_PARSECSV_HPP_
#define PACKAGES_STK_STK_TOOLS_STK_TOOLS_BLOCK_EXTRACTOR_PARSECSV_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

namespace stk {
namespace tools {

std::string strip_string(const std::string &token)
{
    std::string tmp(token);
    tmp.erase(0, tmp.find_first_not_of(" "));
    tmp.erase(tmp.find_last_not_of(" ")+1);
    return tmp;
}

std::vector<std::string> get_csv(const std::string &input)
{
    std::vector<std::string> separated;
    std::istringstream iss(input);
    std::string token;
    while(std::getline(iss, token, ','))
        separated.push_back(strip_string(token));
    return separated;
}

std::vector<std::string> get_block_names_given_ids(const std::vector<std::string> &ids)
{
    std::vector<std::string> names(ids.size());
    for(size_t i=0; i<ids.size(); i++)
    {
        names[i] = "block_" + ids[i];
    }
    return names;
}

}}

#endif /* PACKAGES_STK_STK_TOOLS_STK_TOOLS_BLOCK_EXTRACTOR_PARSECSV_HPP_ */
