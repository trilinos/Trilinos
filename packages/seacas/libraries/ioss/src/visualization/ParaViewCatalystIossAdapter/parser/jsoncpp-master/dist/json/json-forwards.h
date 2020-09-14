/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
//////////////////////////////////////////////////////////////////////
// End of content of file: include/json/assertions.h
// //////////////////////////////////////////////////////////////////////

#endif // ifndef JSON_AMALGATED_H_INCLUDED
cause an exception.
      - This is a security issue (seg-faults caused by deeply nested JSON),
        so the default is low.
    - `"failIfExtra": false or true`
      - If true, `parse()` returns false when extra non-whitespace trails
        the JSON value in the input string.
    - `"rejectDupKeys": false or true`
      - If true, `parse()` returns false when a key is duplicated within an object.

    You can examine 'settings_` yourself
    to see the defaults. You can also write and read them just like any
    JSON Value.
    \sa setDefaults()
    */
  Json::Value settings_;

CharReaderBuilder();
virtual ~CharReaderBuilder();

virtual CharReader *newCharReader() const;

/** \return true if 'settings' are legal and consistent;
 *   otherwise, indicate bad settings via 'invalid'.
 */
bool validate(Json::Value *invalid) const;

/** A simple way to update a specific setting.
 */
Value &operator[](std::string key);

/** Called by ctor, but you can use this to reset settings_.
 * \pre 'settings' != NULL (but Json::null is fine)
 * \remark Defaults:
 * \snippet src/lib_json/json_reader.cpp CharReaderBuilderStrictMode
 */
static void setDefaults(Json::Value *settings);
/** Same as old Features::strictMode().
 * \pre 'settings' != NULL (but Json::null is fine)
 * \remark Defaults:
 * \snippet src/lib_json/json_reader.cpp CharReaderBuilderDefaults
 */
static void strictMode(Json::Value *settings);
}
;

/** Consume entire stream and use its begin/end.
 * Someday we might have a real StreamReader, but for now this
 * is convenient.
 */
bool JSON_API parseFromStream(CharReader::Factory const &, std::istream &, Value *root,
                              std::string *errs);

/** \brief Read from 'sin' into 'root'.

 Always keep comments from the input JSON.

 This can be used to read a file into a particular sub-object.
 For example:
 \code
 Json::Value root;
 cin >> root["dir"]["file"];
 cout << root;
 \endcode
 Result:
 \verbatim
 {
 "dir": {
     "file": {
     // The input stream JSON would be nested here.
     }
 }
 }
 \endverbatim
 \throw std::exception on parse error.
 \see Json::operator<<()
*/
JSON_API std::istream &operator>>(std::istream &, Value &);

} // namespace Json

#if defined(JSONCPP_DISABLE_DLL_INTERFACE_WARNING)
#pragma warning(pop)
#endif // if defined(JSONCPP_DISABLE_DLL_INTERFACE_WARNING)

#endif // CPPTL_JSON_READER_H_INCLUDED

// //////////////////////////////////////////////////////////////////////
// End of content of file: include/json/reader.h
// //////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////
// Beginning of content of file: include/json/writer.h
// //////////////////////////////////////////////////////////////////////

// Copyright 2007-2010 Baptiste Lepilleur
// Distributed under MIT license, or public domain if desired and
// recognized in your jurisdiction.
// See file LICENSE for detail or copy at http://jsoncpp.sourceforge.net/LICENSE

#ifndef JSON_WRITER_H_INCLUDED
#define JSON_WRITER_H_INCLUDED

#if !defined(JSON_IS_AMALGAMATION)
#include "value.h"
#endif // if !defined(JSON_IS_AMALGAMATION)
#include <ostream>
#include <string>
#include <vector>

// Disable warning C4251: <data member>: <type> needs to have dll-interface to
// be used by...
#if defined(JSONCPP_DISABLE_DLL_INTERFACE_WARNING)
#pragma warning(push)
#pragma warning(disable : 4251)
#endif // if defined(JSONCPP_DISABLE_DLL_INTERFACE_WARNING)

namespace Json {

  class Value;

  /**

  Usage:
  \code
    using namespace Json;
    void writeToStdout(StreamWriter::Factory const& factory, Value const& value) {
      std::unique_ptr<StreamWriter> const writer(
        factory.newStreamWriter());
      writer->write(value, &std::cout);
      std::cout << std::endl;  // add lf and flush
    }
  \endcode
  */
  class JSON_API StreamWriter
  {
  protected:
    std::ostream *sout_; // not owned; will not delete
  public:
    StreamWriter();
    virtual ~StreamWriter();
    /** Write Value into document as configured in sub-class.
        Do not take ownership of sout, but maintain a reference during function.
        \pre sout != NULL
        \return zero on success (For now, we always return zero, so check the
       stream instead.)
        \throw std::exception possibly, depending on configuration
     */
    virtual int write(Value const &root, std::ostream *sout) = 0;

    /** \brief A simple abstract factory.
     */
    class JSON_API Factory
    {
    public:
      virtual ~Factory();
      /** \brief Allocate a CharReader via operator new().
       * \throw std::exception if something goes wrong (e.g. invalid settings)
       */
      virtual StreamWriter *newStreamWriter() const = 0;
    }; // Factory
  };   // StreamWriter

  /** \brief Write into stringstream, then return string, for convenience.
   * A StreamWriter will be created from the factory, used, and then deleted.
   */
  std::string JSON_API writeString(StreamWriter::Factory const &factory, Value const &root);

  /** \brief Build a StreamWriter implementation.

  Usage:
  \code
    using namespace Json;
    Value value = ...;
    StreamWriterBuilder builder;
    builder["commentStyle"] = "None";
    builder["indentation"] = "   ";  // or whatever you like
    std::unique_ptr<Json::StreamWriter> writer(
        builder.newStreamWriter());
    writer->write(value, &std::cout);
    std::cout << std::endl;  // add lf and flush
  \endcode
  */
  class JSON_API StreamWriterBuilder : public StreamWriter::Factory
  {
  public:
    // Note: We use a Json::Value so that we can add data-members to this class
    // without a major version bump.
    /** Configuration of this builder.
      Available settings (case-sensitive):
      - "commentStyle": "None" or "All"
      - "indentation":  "<anything>"
      - "enableYAMLCompatibility": false or true
        - slightly change the whitespace around colons
      - "dropNullPlaceholders": false or true
        - Drop the "null" string from the writer's output for nullValues.
          Strictly speaking, this is not valid JSON. But when the output is being
          fed to a browser's Javascript, it makes for smaller output and the
          browser can handle the output just fine.

      You can examine 'settings_` yourself
      to see the defaults. You can also write and read them just like any
      JSON Value.
      \sa setDefaults()
      */
    Json::Value settings_;

    StreamWriterBuilder();
    virtual ~StreamWriterBuilder();

    /**
     * \throw std::exception if something goes wrong (e.g. invalid settings)
     */
    virtual StreamWriter *newStreamWriter() const;

    /** \return true if 'settings' are legal and consistent;
     *   otherwise, indicate bad settings via 'invalid'.
     */
    bool validate(Json::Value *invalid) const;
    /** A simple way to update a specific setting.
     */
    Value &operator[](std::string key);

    /** Called by ctor, but you can use this to reset settings_.
     * \pre 'settings' != NULL (but Json::null is fine)
     * \remark Defaults:
     * \snippet src/lib_json/json_writer.cpp StreamWriterBuilderDefaults
     */
    static void setDefaults(Json::Value *settings);
  };

  /** \brief Abstract class for writers.
   * \deprecated Use StreamWriter. (And really, this is an implementation detail.)
   */
  class JSON_API Writer
  {
  public:
    virtual ~Writer();

    virtual std::string write(const Value &root) = 0;
  };

  /** \brief Outputs a Value in <a HREF="http://www.json.org">JSON</a> format
   *without formatting (not human friendly).
   *
   * The JSON document is written in a single line. It is not intended for 'human'
   *consumption,
   * but may be useful to support feature such as RPC where bandwidth is limited.
   * \sa Reader, Value
   * \deprecated Use StreamWriterBuilder.
   */
  class JSON_API FastWriter : public Writer
  {

  public:
    FastWriter();
    virtual ~FastWriter() {}

    void enableYAMLCompatibility();

    /** \brief Drop the "null" string from the writer's output for nullValues.
     * Strictly speaking, this is not valid JSON. But when the output is being
     * fed to a browser's Javascript, it makes for smaller output and the
     * browser can handle the output just fine.
     */
    void dropNullPlaceholders();

    void omitEndingLineFeed();

  public: // overridden from Writer
    virtual std::string write(const Value &root);

  private:
    void writeValue(const Value &value);

    std::string document_;
    bool        yamlCompatiblityEnabled_;
    bool        dropNullPlaceholders_;
    bool        omitEndingLineFeed_;
  };

  /** \brief Writes a Value in <a HREF="http://www.json.org">JSON</a> format in a
 *human friendly way.
 *
 * The rules for line break and indent are as follow:
 * - Object value:
 *     - if empty then print {} without indent and line break
 *     - if not empty the print '{', line break & indent, print one value per

