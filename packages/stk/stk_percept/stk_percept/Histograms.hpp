#ifndef stk_percept_Histograms_h
#define stk_percept_Histograms_h

#include <stk_percept/Percept.hpp>
#include <stk_percept/Histogram.hpp>
#include <sstream>
#include <fstream>
#include <map>
#include <set>

#if STK_ADAPT_HAVE_YAML_CPP
#include <yaml-cpp/yaml.h>

#ifndef YAML_CHECK
#define YAML_CHECK(emitter) do { if (1 && !emitter.good()) { std::cout << "Emitter error: " << __FILE__ << ":" << __LINE__ << " emitter.good()= " \
                                                                       << emitter.good() << " Error Message: " << emitter.GetLastError() << std::endl; return;} } while(0)
#endif

#ifndef YAML_ERRCHECK
#define YAML_ERRCHECK YAML_CHECK(emitter)
#endif

#endif

namespace stk {
  namespace percept {

    // a simple database of histograms
    template<typename T>
    class Histograms : public std::map<std::string, Histogram<T> >
    {
    public:
      std::string m_file_root;
      typedef std::map<std::string, Histogram<T> > HistogramMap;

      Histograms(std::string file_root="cout") : m_file_root(file_root) {}
      ~Histograms()
      {
      }

      void compute_uniform_bins(unsigned num_bins, bool use_log_scale=false)
      {
        for (typename HistogramMap::iterator iter=this->begin(); iter != this->end(); ++iter)
          {
            if (iter->second.size())
              iter->second.compute_uniform_bins(num_bins, use_log_scale);
          }
      }

      void print(bool print_bar_chart=false)
      {
        for (typename HistogramMap::iterator iter=this->begin(); iter != this->end(); ++iter)
          {
            if (iter->second.size())
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
                      iter->second.print_table(fout, true);
                  }
              }
          }
      }

    };

#if STK_ADAPT_HAVE_YAML_CPP

    /// Set @param result if the @param key is present in the @param node, else set it to the given default value
    template<typename T>
    void set_if_present(const YAML::Node & node, const std::string& key, T& result, const T& default_if_not_present = T())
    {
      const YAML::Node *value = node.FindValue(key);
      if (value)
        *value >> result;
      else
        result = default_if_not_present;
    }

    /// this version doesn't change @param result unless the @param key is present in the @param node
    template<typename T>
    void set_if_present_no_default(const YAML::Node & node, const std::string& key, T& result)
    {
      const YAML::Node *value = node.FindValue(key);
      if (value)
        *value >> result;
    }

    template<typename T>
    class HistogramsParser
    {
      std::string m_root_string;
    public:
      HistogramsParser(std::string root_string) : m_root_string(root_string) {}

      void create(Histograms<T>& histograms)
      {
        std::stringstream ss(m_root_string);

        YAML::Parser parser(ss);
        YAML::Node node, node1;

        try {
          while(parser.GetNextDocument(node)) {
            //std::cout << "\n read doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (node.Type() == YAML::NodeType::Scalar)
              {
                std::string file_name;
                node >> file_name;
                std::cout << "HistogramsParser:: redirecting input from "<< file_name << std::endl;
                std::ifstream file(file_name.c_str());
                YAML::Parser parser1(file);
                while(parser1.GetNextDocument(node1)) {
                  parse(node1, histograms);
                }
              }
            else
              {
                parse(node, histograms);
              }
          }
        }
        catch(YAML::ParserException& e) {
          std::cout << e.what() << "\n";
        }
      }

    private:

      void parse(const YAML::Node& node, Histograms<T>& histograms)
      {
        const YAML::Node *y_element_fields = node.FindValue("fields");
        if (y_element_fields)
          {
            for (size_t iter = 0; iter < y_element_fields->size(); iter++)
              {
                std::string element_field_name;
                (*y_element_fields)[iter] >> element_field_name;
                std::string title="Field "+element_field_name;
                histograms["field."+element_field_name].set_titles(title);
              }
          }
        std::string file_root;
        set_if_present(node, "file_root", file_root, std::string("cout"));
        histograms.m_file_root = file_root;

        {
          std::string mesh_fields[4] =  {"edge_length", "quality_edge", "quality_vol_edge_ratio", "volume"};
          std::set<std::string> valid_values(mesh_fields, mesh_fields+4);

          const YAML::Node *y_mesh_fields = node.FindValue("mesh");
          if (y_mesh_fields)
            {
              for (size_t iter = 0; iter < y_mesh_fields->size(); iter++)
                {
                  std::string mesh_field_name;
                  (*y_mesh_fields)[iter] >> mesh_field_name;
                  if (valid_values.find(mesh_field_name) == valid_values.end())
                    throw std::runtime_error("HistogramsParser:: unrecognized option: " + mesh_field_name);
                  std::string title="Mesh Field "+mesh_field_name;
                  histograms["mesh."+mesh_field_name].set_titles(title);
                }
            }
        }
      }


    };
#endif


  }
}

#endif
