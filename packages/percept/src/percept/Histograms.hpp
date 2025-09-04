// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Histograms_h
#define percept_Histograms_h

#include <percept/Percept.hpp>
#include <percept/HistogramV.hpp>
#include <sstream>
#include <fstream>
#include <map>
#include <set>

#if HAVE_YAML
#include <percept/YamlUtils.hpp>
#endif

  namespace percept {

    template<typename T>
    struct Histogram
    {
      Histogram()
        : m_data(), m_hist(m_data), m_bar_symbol('='), m_max_column_width(40)
      {
      }
      Histogram(const Histogram& other)
        : m_data(other.m_data), m_hist(m_data), m_bar_symbol(other.m_bar_symbol), m_max_column_width(other.m_max_column_width)
      {
        m_titles[0] = other.m_titles[0];
        m_titles[1] = other.m_titles[1];
        m_titles[2] = other.m_titles[2];
        m_titles[3] = other.m_titles[3];
      }

      Histogram& operator=(const Histogram& other)
      {
        m_data = other.m_data;
        m_bar_symbol = other.m_bar_symbol;
        m_max_column_width = other.m_max_column_width;
      }

      std::vector<T> m_data;
      HistogramV<T> m_hist;
      char m_bar_symbol;
      unsigned m_max_column_width;
      std::string m_titles[4];

      void set_titles(std::string title="Histogram",
                      std::string title1="Ranges / Mid-Point", std::string title2="Counts", std::string title3="Percentage")
      {
        m_titles[0]=title;
        m_titles[1]=title1;
        m_titles[2]=title2;
        m_titles[3]=title3;
      }

      void print_simple_table(std::ostream& stream)
      {
        HistogramV<T>::print_simple_table(m_hist, m_titles, stream);
      }
      void print_table(std::ostream& os, bool use_percentage = false)
      {
        HistogramV<T>::print_table(m_hist, m_titles, os, m_max_column_width, m_bar_symbol, use_percentage);
      }
      void print(std::ostream& stream)
      {
        HistogramV<T>::print_histogram(m_hist, stream);
      }
    };

    // a simple database of histograms
    template<typename T>
    class Histograms : public std::map<std::string, Histogram<T> >
    {
    public:
      std::string m_file_root;
      double m_database_time;
      int m_database_step;
      typedef std::map<std::string, Histogram<T> > HistogramMap;

      Histograms(std::string file_root="cout") : m_file_root(file_root), m_database_time(-1.0), m_database_step(-1) {}
      ~Histograms()
      {
      }

      void compute_uniform_bins(unsigned num_bins, bool use_log_scale=false)
      {
        for (typename HistogramMap::iterator iter=this->begin(); iter != this->end(); ++iter)
          {
            iter->second.m_hist.compute_uniform_bins(num_bins, use_log_scale);
          }
      }

      void print(bool print_bar_chart=false)
      {
        for (typename HistogramMap::iterator iter=this->begin(); iter != this->end(); ++iter)
          {
            if (iter->second.m_data.size())
              {
                std::string file_name = m_file_root+"."+iter->first+".hist";
                if (m_file_root == "cout")
                  {
                    iter->second.set_titles(iter->first);
                    iter->second.print_simple_table(std::cout);
                    if (print_bar_chart)
                      iter->second.print_table(std::cout, true);
                  }
                else
                  {
                    std::ofstream fout(file_name.c_str());
                    iter->second.set_titles(iter->first);
                    iter->second.print_simple_table(fout);
                    if (print_bar_chart)
                      {
                        file_name = m_file_root+"."+iter->first+".bar.hist";
                        std::ofstream fout1(file_name.c_str());
                        iter->second.print_table(fout1, true);
                      }
                  }
              }
          }
      }

    };

    template<typename T>
    class HistogramsParser
    {
      std::string m_root_string;
    public:


      HistogramsParser(std::string root_string) : m_root_string(root_string) {
        Util::replace(m_root_string, ":", ": ");
      }

      void create([[maybe_unused]] Histograms<T>& histograms)
      { 
#if HAVE_YAML
        std::stringstream ss(m_root_string);

        //YAML::Parser parser(ss);
        YAML::Node node, node1;

        try {

          if (1) {
            //while(parser.GetNextDocument(node)) {
            node = YAML::Load(ss);
            //std::cout << "\n read doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (node.Type() == YAML::NodeType::Scalar)
              {
                std::string file_name;
                file_name = node.as<std::string>();
                std::cout << "HistogramsParser:: redirecting input from "<< file_name << std::endl;
                std::ifstream file(file_name.c_str());
                //YAML::Parser parser1(file);
                //while(parser1.GetNextDocument(node1)) {
                node1 = YAML::Load(file);
                parse(node1, histograms);
              }
            else
              {
                parse(node, histograms);
              }
          }
        }
        catch(YAML::ParserException& e) {
          std::cout << e.what() << " input= " << m_root_string << "\n";
        }
#endif
      }

    private:

#if HAVE_YAML
      void parse(const YAML::Node& node, Histograms<T>& histograms)
      {
        set_if_present(node, "time", histograms.m_database_time, double(-1.0));

        set_if_present(node, "step", histograms.m_database_step, int(-1));

        const YAML::Node y_element_fields = node["fields"];
        if (y_element_fields)
          {
            for (size_t iter = 0; iter < y_element_fields.size(); iter++)
              {
                std::string element_field_name;
                element_field_name = y_element_fields[iter].as<std::string>();
                std::string title="Field "+element_field_name;
                histograms["field."+element_field_name].set_titles(title);
              }
          }
        std::string file_root;
        set_if_present_no_default(node, "file_root", file_root);
        if (file_root.length()) histograms.m_file_root = file_root;

        {
          std::string mesh_fields[4] =  {"edge_length", "quality_edge", "quality_vol_edge_ratio", "volume"};
          std::set<std::string> valid_values(mesh_fields, mesh_fields+4);

          const YAML::Node y_mesh_fields = node["mesh"];
          if (y_mesh_fields)
            {
              for (size_t iter = 0; iter < y_mesh_fields.size(); iter++)
                {
                  std::string mesh_field_name;
                  mesh_field_name = y_mesh_fields[iter].as<std::string>();
                  if (valid_values.find(mesh_field_name) == valid_values.end())
                    throw std::runtime_error("HistogramsParser:: unrecognized option: " + mesh_field_name);
                  std::string title="Mesh Field "+mesh_field_name;
                  //std::cout << "HistogramsParser::parse: adding " << mesh_field_name << std::endl;
                  histograms["mesh."+mesh_field_name].set_titles(title);
                }
            }
        }
      }
#endif

    };

  }

#endif
