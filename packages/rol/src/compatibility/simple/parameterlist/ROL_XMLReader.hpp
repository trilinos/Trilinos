// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_XMLREADER_HPP
#define ROL_XMLREADER_HPP

#include <pugixml.hpp>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <type_traits>
#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

namespace ROL {

class XMLReader {
private:
    enum class ArrayType { DOUBLE, INT, STRING, UNKNOWN };

    // Helper function to detect if string looks like an array
    static bool isArrayString(const std::string& value) {
        size_t start = value.find('{');
        size_t end = value.find('}');
        return (start != std::string::npos && end != std::string::npos && end > start);
    }

    // Helper function to detect array element type
    static ArrayType detectArrayType(const std::string& value) {
        std::string cleaned = value;
        size_t start = cleaned.find('{');
        size_t end = cleaned.find('}');
        if (start != std::string::npos && end != std::string::npos && end > start) {
            cleaned = cleaned.substr(start + 1, end - start - 1);
        }

        // Check first non-whitespace element to determine type
        std::istringstream iss(cleaned);
        std::string token;

        if (std::getline(iss, token, ',')) {
            // Trim whitespace
            size_t first = token.find_first_not_of(" \t");
            if (first != std::string::npos) {
                size_t last = token.find_last_not_of(" \t");
                token = token.substr(first, last - first + 1);

                // Try to parse as double, then int
                try {
                    std::stod(token);
                    // Check if it contains decimal point
                    if (token.find('.') != std::string::npos) {
                        return ArrayType::DOUBLE;
                    } else {
                        return ArrayType::INT;
                    }
                } catch (...) {
                    return ArrayType::STRING;
                }
            }
        }

        return ArrayType::UNKNOWN;
    }

public:
    static ROL::Ptr<ParameterList> readParametersFromFile(const std::string& filename) {
        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load_file(filename.c_str());

        if (!result) {
            throw std::runtime_error("XMLReader: Failed to load XML file: " + filename +
                                   " Error: " + result.description());
        }

        ROL::Ptr<ParameterList> params = ROL::makePtr<ParameterList>();

        // Find the root ParameterList element
        pugi::xml_node root = doc.child("ParameterList");
        if (!root) {
            throw std::runtime_error("XMLReader: No ParameterList element found in XML file");
        }

        parseParameterList(root, *params);
        return params;
    }

private:
    static void parseParameterList(const pugi::xml_node& node, ParameterList& params) {
        for (pugi::xml_node child : node.children()) {
            std::string name = child.name();

            if (name == "Parameter") {
                parseParameter(child, params);
            } else if (name == "ParameterList") {
                parseSublist(child, params);
            }
        }
    }

    static void parseParameter(const pugi::xml_node& node, ParameterList& params) {
        std::string name = node.attribute("name").as_string();
        std::string type = node.attribute("type").as_string();
        std::string value = node.attribute("value").as_string();

        if (name.empty()) {
            throw std::runtime_error("XMLReader: Parameter element missing 'name' attribute");
        }

        if (type.empty()) {
            throw std::runtime_error("XMLReader: Parameter element missing 'type' attribute");
        }

        // Parse value based on type
        if (type == "double") {
            params.set(name, std::stod(value));
        } else if (type == "int") {
            params.set(name, std::stoi(value));
        } else if (type == "bool") {
            bool bval = (value == "true" || value == "1");
            params.set(name, bval);
        } else if (type == "string") {
            // Check if this is an array-like string and try to parse it as an array
            if (isArrayString(value)) {
                if (detectArrayType(value) == ArrayType::DOUBLE) {
                    try {
                        std::vector<double> arr = fromStringToArray<double>(value);
                        params.set(name, arr);
                        return;
                    } catch (...) {
                        // Fall back to string if parsing fails
                    }
                } else if (detectArrayType(value) == ArrayType::INT) {
                    try {
                        std::vector<int> arr = fromStringToArray<int>(value);
                        params.set(name, arr);
                        return;
                    } catch (...) {
                        // Fall back to string if parsing fails
                    }
                }
            }
            params.set(name, value);
        } else {
            // Default to string for unknown types, but check for arrays
            if (isArrayString(value)) {
                if (detectArrayType(value) == ArrayType::DOUBLE) {
                    try {
                        std::vector<double> arr = fromStringToArray<double>(value);
                        params.set(name, arr);
                        return;
                    } catch (...) {
                        // Fall back to string if parsing fails
                    }
                } else if (detectArrayType(value) == ArrayType::INT) {
                    try {
                        std::vector<int> arr = fromStringToArray<int>(value);
                        params.set(name, arr);
                        return;
                    } catch (...) {
                        // Fall back to string if parsing fails
                    }
                }
            }
            params.set(name, value);
        }
    }

    static void parseSublist(const pugi::xml_node& node, ParameterList& params) {
        std::string name = node.attribute("name").as_string();
        if (name.empty()) {
            throw std::runtime_error("XMLReader: ParameterList element missing 'name' attribute");
        }

        ParameterList& sublist = params.sublist(name);
        parseParameterList(node, sublist);
    }
};

class XMLWriter {
public:
    static void writeParametersToFile(const ParameterList& params, const std::string& filename) {
        pugi::xml_document doc;

        // Create root ParameterList element
        pugi::xml_node root = doc.append_child("ParameterList");
        root.append_attribute("name") = "Main";

        // Write parameters to XML
        writeParameterList(params, root);

        // Save to file
        if (!doc.save_file(filename.c_str(), "  ")) {
            throw std::runtime_error("XMLWriter: Failed to write XML file: " + filename);
        }
    }

private:
    static void writeParameterList(const ParameterList& params, pugi::xml_node& node) {
        // Write all parameters using the iterator interface
        for (auto it = params.begin(); it != params.end(); ++it) {
            std::string key = params.name(it);
            writeParameter(params, key, node);
        }

        // Note: Sublists would need additional implementation if supported
        // For now, this handles the basic parameter types
    }

    static void writeParameter(const ParameterList& params, const std::string& key, pugi::xml_node& node) {
        pugi::xml_node param_node = node.append_child("Parameter");
        param_node.append_attribute("name") = key.c_str();

        // Check type and write accordingly
        if (params.isType<bool>(key)) {
            param_node.append_attribute("type") = "bool";
            bool value = params.get<bool>(key);
            param_node.append_attribute("value") = value ? "true" : "false";
        }
        else if (params.isType<int>(key)) {
            param_node.append_attribute("type") = "int";
            int value = params.get<int>(key);
            param_node.append_attribute("value") = std::to_string(value).c_str();
        }
        else if (params.isType<double>(key)) {
            param_node.append_attribute("type") = "double";
            double value = params.get<double>(key);
            param_node.append_attribute("value") = std::to_string(value).c_str();
        }
        else if (params.isType<std::string>(key)) {
            param_node.append_attribute("type") = "string";
            std::string value = params.get<std::string>(key);
            param_node.append_attribute("value") = value.c_str();
        }
        else if (params.isType<std::vector<double>>(key)) {
            param_node.append_attribute("type") = "string";
            std::vector<double> value = params.get<std::vector<double>>(key);
            param_node.append_attribute("value") = vectorToString(value).c_str();
        }
        else if (params.isType<std::vector<int>>(key)) {
            param_node.append_attribute("type") = "string";
            std::vector<int> value = params.get<std::vector<int>>(key);
            param_node.append_attribute("value") = vectorToString(value).c_str();
        }
        else if (params.isType<std::vector<std::string>>(key)) {
            param_node.append_attribute("type") = "string";
            std::vector<std::string> value = params.get<std::vector<std::string>>(key);
            param_node.append_attribute("value") = vectorToString(value).c_str();
        }
        else {
            // Default case - try as string
            param_node.append_attribute("type") = "string";
            param_node.append_attribute("value") = "unknown_type";
        }
    }

    // Helper function to convert vector to string representation
    template<typename T>
    static std::string vectorToString(const std::vector<T>& vec) {
        if (vec.empty()) return "{ }";

        std::ostringstream oss;
        oss << "{ ";
        for (size_t i = 0; i < vec.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << vec[i];
        }
        oss << " }";
        return oss.str();
    }
};

} // namespace ROL

#endif // ROL_XMLREADER_HPP
