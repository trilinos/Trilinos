// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "catch.hpp"
#include "ImageSetValidTests.h"

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCamera", "[ImageSetValid]")
{
  Json::Value dj = getDefaultImageSetJSON();

  Json::Value camera;
  camera["camera type"]            = "camera";
  dj["camera blocks"]["fooCamera"] = camera;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["camera"]    = "fooCamera";

  checkPhactoriStringValidParse(ImagesetCamera, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentation", "[ImageSetValid]")
{
  Json::Value dj                                   = getDefaultImageSetWithCameraJSON();
  dj["representation blocks"]["fooRepresentation"] = Json::objectValue;

  Json::Value &imageset      = dj["imageset blocks"]["fooImageset"];
  imageset["representation"] = "fooRepresentation";

  checkPhactoriStringValidParse(ImagesetRepresentation, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperation", "[ImageSetCameraValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value clip;
  clip["type"]                           = "clip";
  dj["operation blocks"]["fooOperation"] = clip;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["operation"] = "fooOperation";

  checkPhactoriStringValidParse(ImagesetOperation, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetBasename", "[ImageSetValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset      = dj["imageset blocks"]["fooImageset"];
  imageset["image basename"] = "foo";

  checkPhactoriStringValidParse(ImagesetBasename, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetBasedirectory", "[ImageSetValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset           = dj["imageset blocks"]["fooImageset"];
  imageset["image basedirectory"] = "foo";

  checkPhactoriStringValidParse(ImagesetBasedirectory, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetFormat", "[ImageSetValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset    = dj["imageset blocks"]["fooImageset"];
  imageset["image format"] = "jpg";

  checkPhactoriStringValidParse(ImagesetFormat, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetColorSettings", "[ImageSetValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset     = dj["imageset blocks"]["fooImageset"];
  imageset["axes color"][0] = 0.0;
  imageset["axes color"][1] = 1.0;
  imageset["axes color"][2] = 1.0;

  imageset["edge color"][0] = 1.0;
  imageset["edge color"][1] = 1.0;
  imageset["edge color"][2] = 0.0;

  imageset["text color"][0] = 1.0;
  imageset["text color"][1] = 0.0;
  imageset["text color"][2] = 1.0;

  imageset["background color"][0] = 0.0;
  imageset["background color"][1] = 0.0;
  imageset["background color"][2] = 0.0;

  imageset["show edges"]           = true;
  imageset["show axes"]            = true;
  imageset["show time annotation"] = true;
  imageset["color by scalar"]      = "VON_MISES";

  checkPhactoriStringValidParse(ImagesetColorSettings, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetDigitCount", "[ImageSetValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset         = dj["imageset blocks"]["fooImageset"];
  imageset["image digit count"] = 6;

  checkPhactoriStringValidParse(ImagesetDigitCount, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImageFilenameDatetimeSimtime", "[notWorking]")
{
  Json::Value dj = getDefaultPhactoriJSON();

  Json::Value camera;
  camera["camera type"]                   = "camera";
  camera["look direction"][0]             = -1.0;
  camera["look direction"][1]             = -1.0;
  camera["look direction"][2]             = -1.0;
  dj["camera blocks"]["singletestcamera"] = camera;

  Json::Value fooImageset1;
  fooImageset1["camera"]        = "singletestcamera";
  fooImageset1["image size"][0] = 800;
  fooImageset1["image size"][1] = 450;

  Json::Value fooImageset2(fooImageset1);
  Json::Value fooImageset3(fooImageset1);
  Json::Value fooImageset4(fooImageset1);
  Json::Value fooImageset5(fooImageset1);
  Json::Value fooImageset6(fooImageset1);

  fooImageset1["image basename"]        = "is_datetime.";
  fooImageset1["image digit count"]     = 6;
  fooImageset1["filename use datetime"] = true;

  fooImageset2["image basename"]                       = "is_iso_datetime.";
  fooImageset2["filename use datetime"]                = true;
  fooImageset2["filename datetime convert to integer"] = false;

  fooImageset3["image basename"]               = "is_simulation_time.";
  fooImageset3["image digit count"]            = 6;
  fooImageset3["filename use simulation time"] = true;

  fooImageset4["image basename"]               = "is_datetime_and_simulation_time.";
  fooImageset4["image digit count"]            = 6;
  fooImageset4["filename use simulation time"] = true;
  fooImageset4["filename use datetime"]        = true;

  fooImageset5["image basename"]               = "is_datetime_and_simulation_time_no_call_count.";
  fooImageset5["filename use simulation time"] = true;
  fooImageset5["filename use datetime"]        = true;
  fooImageset5["filename use call count"]      = false;

  fooImageset6["filename use simulation time"] = true;
  fooImageset6["filename use datetime"]        = true;
  fooImageset6["filename use call count"]      = true;

  dj["imageset blocks"]["fooImageset1"] = fooImageset1;
  dj["imageset blocks"]["fooImageset2"] = fooImageset2;
  dj["imageset blocks"]["fooImageset3"] = fooImageset3;
  dj["imageset blocks"]["fooImageset4"] = fooImageset4;
  dj["imageset blocks"]["fooImageset5"] = fooImageset5;
  dj["imageset blocks"]["fooImageset6"] = fooImageset6;

  checkPhactoriStringValidParse(ImageFilenameDatetimeSimtime, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetDefault", "[ImageSetValid]")
{
  Json::Value dj = getDefaultImageSetJSON();

  Json::Value camera;
  camera["camera type"]                   = "camera";
  camera["look direction"][0]             = -1.0;
  camera["look direction"][1]             = -1.0;
  camera["look direction"][2]             = -1.0;
  dj["camera blocks"]["singletestcamera"] = camera;

  Json::Value imageset;
  imageset["camera"]                   = "singletestcamera";
  dj["imageset blocks"]["fooImageset"] = imageset;

  checkPhactoriStringValidParse(ImagesetDefault, dj);
}
