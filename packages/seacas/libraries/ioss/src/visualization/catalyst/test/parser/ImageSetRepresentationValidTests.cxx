// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ImageSetRepresentationValidTests.h"
#include "CatalystTestFixture.h"
#include "catch.hpp"

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationSurfaces",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset     = dj["imageset blocks"]["fooImageset"];
  imageset["show surfaces"] = false;
  imageset["show edges"]    = true;

  checkPhactoriStringValidParse(ImagesetRepresentationSurfaces, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationEdges",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset  = dj["imageset blocks"]["fooImageset"];
  imageset["show edges"] = true;

  checkPhactoriStringValidParse(ImagesetRepresentationEdges, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationBoundingBox",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset         = dj["imageset blocks"]["fooImageset"];
  imageset["show surfaces"]     = false;
  imageset["show bounding box"] = true;

  checkPhactoriStringValidParse(ImagesetRepresentationBoundingBox, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorByScalar",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset       = dj["imageset blocks"]["fooImageset"];
  imageset["color by scalar"] = "VON_MISES";

  checkPhactoriStringValidParse(ImagesetRepresentationColorByScalar, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorByVectorMagnitude",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset                 = dj["imageset blocks"]["fooImageset"];
  imageset["color by vector magnitude"] = "displ";

  checkPhactoriStringValidParse(ImagesetRepresentationColorByVectorMagnitude, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorByVectorComponent",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset                 = dj["imageset blocks"]["fooImageset"];
  imageset["color by vector component"] = "displ_x";

  checkPhactoriStringValidParse(ImagesetRepresentationColorByVectorComponent, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorByTensorComponent",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset                 = dj["imageset blocks"]["fooImageset"];
  imageset["color by tensor component"] = "stress_xy";

  checkPhactoriStringValidParse(ImagesetRepresentationColorByTensorComponent, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorBySolidColor",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  color;
  color[0]                         = 1.0;
  color[1]                         = 0.5;
  color[2]                         = 0.0;
  imageset["color by solid color"] = color;

  checkPhactoriStringValidParse(ImagesetRepresentationColorBySolidColor, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorByBlockId",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset        = dj["imageset blocks"]["fooImageset"];
  imageset["color by blockid"] = true;

  checkPhactoriStringValidParse(ImagesetRepresentationColorByBlockId, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationTimeAnnotation",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset            = dj["imageset blocks"]["fooImageset"];
  imageset["show time annotation"] = true;

  checkPhactoriStringValidParse(ImagesetRepresentationTimeAnnotation, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorLegend",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset         = dj["imageset blocks"]["fooImageset"];
  imageset["show color legend"] = true;
  imageset["color by scalar"]   = "VON_MISES";

  checkPhactoriStringValidParse(ImagesetRepresentationColorLegend, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorLegendRange",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  range;
  range[0]                       = 22.0;
  range[1]                       = 32.0;
  imageset["color legend range"] = range;
  imageset["color by scalar"]    = "VON_MISES";

  checkPhactoriStringValidParse(ImagesetRepresentationColorLegendRange, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorLegendUseCurrentDataRange",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset                           = dj["imageset blocks"]["fooImageset"];
  imageset["color legend use current data range"] = true;
  imageset["color by tensor component"]           = "stress_yz";

  checkPhactoriStringValidParse(ImagesetRepresentationColorLegendUseCurrentDataRange, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorLegendUseCumulativeDataRange",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset                              = dj["imageset blocks"]["fooImageset"];
  imageset["color legend use cumulative data range"] = true;
  imageset["color by tensor component"]              = "stress_yz";

  checkPhactoriStringValidParse(ImagesetRepresentationColorLegendUseCumulativeDataRange, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorLegendMinimumRange",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  range;
  range[0]                               = -5.0;
  range[1]                               = 15.0;
  imageset["color legend minimum range"] = range;
  imageset["color by scalar"]            = "VON_MISES";

  checkPhactoriStringValidParse(ImagesetRepresentationColorLegendMinimumRange, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorLegendMaximumRange",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  range;
  range[0]                               = -5.0;
  range[1]                               = 15.0;
  imageset["color legend maximum range"] = range;
  imageset["color by scalar"]            = "VON_MISES";

  checkPhactoriStringValidParse(ImagesetRepresentationColorLegendMaximumRange, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationColorLegendPosition",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  position;
  position[0]                       = "right";
  position[1]                       = 0.5;
  imageset["color legend position"] = position;
  imageset["color by scalar"]       = "VON_MISES";

  checkPhactoriStringValidParse(ImagesetRepresentationColorLegendPosition, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationPresetColorScale",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();
  Json::Value rep;
  rep["color by tensor component"]                           = "stress_xy";
  rep["preset color scale"]                                  = "Blue_to_Red_Rainbow";
  rep["invert color scale"]                                  = true;
  dj["representation blocks"]["BlueToRedRainbowInvertedRep"] = rep;

  Json::Value &imageset      = dj["imageset blocks"]["fooImageset"];
  imageset["representation"] = "BlueToRedRainbowInvertedRep";
  imageset["image basename"] = "BlueToRedRainbow.";
  std::swap(dj["imageset blocks"]["BlueToRedRainbowInverted"], imageset);
  dj["imageset blocks"].removeMember("fooImageset");

  checkPhactoriStringValidParse(ImagesetRepresentationPresetColorScale, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationTimeAnnotationPosition",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  position;
  position[0]                          = "top left";
  position[1]                          = 0.5;
  imageset["time annotation position"] = position;
  imageset["show time annotation"]     = true;

  checkPhactoriStringValidParse(ImagesetRepresentationTimeAnnotationPosition, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetMultiColorBy", "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultPhactoriJSON();

  Json::Value camera;
  camera["camera type"]                   = "camera";
  camera["look direction"][0]             = -1.0;
  camera["look direction"][1]             = -1.0;
  camera["look direction"][2]             = -1.0;
  dj["camera blocks"]["singletestcamera"] = camera;

  Json::Value is1;
  is1["camera"]        = "singletestcamera";
  is1["image size"][0] = 800;
  is1["image size"][1] = 450;

  Json::Value is2(is1);
  Json::Value is3(is1);
  Json::Value is4(is1);
  Json::Value is5(is1);

  is1["image basename"]  = "is1.";
  is1["color by scalar"] = "VON_MISES";

  is2["image basename"]            = "is2.";
  is2["color by tensor component"] = "stress_xy";

  is3["image basename"]            = "is3.";
  is3["color by vector magnitude"] = "displ";

  is4["image basename"] = "is4.";
  Json::Value color;
  color[0]                    = 0.1;
  color[1]                    = 0.9;
  color[2]                    = 0.95;
  is4["color by solid color"] = color;

  is5["image basename"] = "is5.";

  dj["imageset blocks"]["is1"] = is1;
  dj["imageset blocks"]["is2"] = is2;
  dj["imageset blocks"]["is3"] = is3;
  dj["imageset blocks"]["is4"] = is4;
  dj["imageset blocks"]["is5"] = is5;

  checkPhactoriStringValidParse(ImagesetMultiColorBy, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetAxesLegendOnOffMultiImageset",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultPhactoriJSON();

  Json::Value camera;
  camera["camera type"]                   = "camera";
  camera["look direction"][0]             = -1.0;
  camera["look direction"][1]             = -1.0;
  camera["look direction"][2]             = -1.0;
  dj["camera blocks"]["singletestcamera"] = camera;

  Json::Value is1;
  is1["camera"]        = "singletestcamera";
  is1["image size"][0] = 800;
  is1["image size"][1] = 450;

  Json::Value is2(is1);
  Json::Value is3(is1);
  Json::Value is4(is1);
  Json::Value is5(is1);
  Json::Value is6(is1);
  Json::Value is7(is1);

  is1["image basename"]        = "is1.";
  is1["color by scalar"]       = "VON_MISES";
  is1["show axes"]             = true;
  is1["show orientation axes"] = true;
  is1["show time annotation"]  = true;
  is1["show edges"]            = true;

  is2["image basename"]            = "is2.";
  is2["color by tensor component"] = "stress_xy";
  is2["show color legend"]         = false;
  is2["show orientation axes"]     = false;

  is3["image basename"]    = "is3.";
  is3["color by scalar"]   = "VON_MISES";
  is3["show color legend"] = true;
  is3["show axes"]         = true;

  is4["image basename"] = "is4.";
  is4["show axes"]      = true;
  is4["show edges"]     = true;

  is5["image basename"] = "is5.";
  Json::Value color;
  color[0]                    = 1.0;
  color[1]                    = 0.5;
  color[2]                    = 0.0;
  is5["color by solid color"] = color;

  is6["image basename"] = "is6.";
  is6["show axes"]      = true;
  is6["show edges"]     = true;
  is6["show surfaces"]  = false;

  is7["image basename"] = "is7.";
  Json::Value color2;
  color2[0]                   = 0.0;
  color2[1]                   = 1.0;
  color2[2]                   = 1.0;
  is7["color by solid color"] = color2;
  is7["show axes"]            = true;

  dj["imageset blocks"]["is1"] = is1;
  dj["imageset blocks"]["is2"] = is2;
  dj["imageset blocks"]["is3"] = is3;
  dj["imageset blocks"]["is4"] = is4;
  dj["imageset blocks"]["is5"] = is5;
  dj["imageset blocks"]["is6"] = is6;
  dj["imageset blocks"]["is7"] = is7;

  checkPhactoriStringValidParse(ImagesetAxesLegendOnOffMultiImageset, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationAxes", "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["show axes"] = true;

  checkPhactoriStringValidParse(ImagesetRepresentationAxes, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationAxisLabelName",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset         = dj["imageset blocks"]["fooImageset"];
  imageset["show axes"]         = true;
  imageset["x axis label name"] = "xAxis";
  imageset["y axis label name"] = "yAxis";
  imageset["z axis label name"] = "zAxis";

  checkPhactoriStringValidParse(ImagesetRepresentationAxisLabelName, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationAxisLabel",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset         = dj["imageset blocks"]["fooImageset"];
  imageset["show axes"]         = true;
  imageset["show x axis label"] = false;
  imageset["show y axis label"] = false;
  imageset["show z axis label"] = false;

  checkPhactoriStringValidParse(ImagesetRepresentationAxisLabel, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationAxisTicMarks",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset             = dj["imageset blocks"]["fooImageset"];
  imageset["show axes"]             = true;
  imageset["show x axis tic marks"] = false;
  imageset["show y axis tic marks"] = false;
  imageset["show z axis tic marks"] = false;

  checkPhactoriStringValidParse(ImagesetRepresentationAxisTicMarks, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationAxisMinorTicMarks",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset                   = dj["imageset blocks"]["fooImageset"];
  imageset["show axes"]                   = true;
  imageset["show x axis minor tic marks"] = false;
  imageset["show y axis minor tic marks"] = false;
  imageset["show z axis minor tic marks"] = false;

  checkPhactoriStringValidParse(ImagesetRepresentationAxisMinorTicMarks, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationAxisTicMarkLabels",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset                   = dj["imageset blocks"]["fooImageset"];
  imageset["show axes"]                   = true;
  imageset["show x axis tic mark labels"] = false;
  imageset["show y axis tic mark labels"] = false;
  imageset["show z axis tic mark labels"] = false;

  checkPhactoriStringValidParse(ImagesetRepresentationAxisTicMarkLabels, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetRepresentationOrientationAxes",
                 "[ImageSetRepresentationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset             = dj["imageset blocks"]["fooImageset"];
  imageset["show orientation axes"] = true;

  Json::Value fooImageset2(imageset);
  fooImageset2["show orientation axes"] = false;
  dj["imageset blocks"]["fooImageset2"] = fooImageset2;

  checkPhactoriStringValidParse(ImagesetRepresentationOrientationAxes, dj);
}
