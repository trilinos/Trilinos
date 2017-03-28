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
