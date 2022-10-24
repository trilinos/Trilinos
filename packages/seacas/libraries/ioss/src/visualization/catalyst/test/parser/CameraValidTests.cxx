// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "catch.hpp"
#include "CameraValidTests.h"

TEST_CASE_METHOD(CatalystTestFixture, "Camera", "[cameraValid]")
{
  Json::Value dj = getDefaultCameraJSON();
  checkPhactoriStringValidParse(Camera, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraLookAtAbsolutePoint", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["look at absolute point"][0] = 1.1;
  camera["look at absolute point"][1] = 2.0;
  camera["look at absolute point"][2] = 3e-8;

  checkPhactoriStringValidParse(CameraLookAtAbsolutePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraLookDirection", "[cameraValid]")
{
  Json::Value  dj             = getDefaultCameraJSON();
  Json::Value &camera         = dj["camera blocks"]["fooCamera"];
  camera["look direction"][0] = 1.0;
  camera["look direction"][1] = 2.0;
  camera["look direction"][2] = 3.0;

  checkPhactoriStringValidParse(CameraLookDirection, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraLookAtRelativeDistance", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["look at relative distance"] = 2.0;
  checkPhactoriStringValidParse(CameraLookAtRelativeDistance, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraLookAtAbsoluteDistance", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["look at absolute distance"] = 15.0;
  checkPhactoriStringValidParse(CameraLookAtAbsoluteDistance, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraAtAbsolutePoint", "[cameraValid]")
{
  Json::Value  dj                       = getDefaultCameraJSON();
  Json::Value &camera                   = dj["camera blocks"]["fooCamera"];
  camera["camera at absolute point"][0] = -2.0;
  camera["camera at absolute point"][1] = 3.0;
  camera["camera at absolute point"][2] = 30.0;
  checkPhactoriStringValidParse(CameraAtAbsolutePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraAtRelativePoint", "[cameraValid]")
{
  Json::Value  dj                       = getDefaultCameraJSON();
  Json::Value &camera                   = dj["camera blocks"]["fooCamera"];
  camera["camera at relative point"][0] = -0.5;
  camera["camera at relative point"][1] = 1.5;
  camera["camera at relative point"][2] = 20.0;
  camera["look direction"][0]           = 0.1;
  camera["look direction"][1]           = -0.1;
  camera["look direction"][2]           = -1.0;
  checkPhactoriStringValidParse(CameraAtRelativePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraAtElement", "[cameraValid]")
{
  Json::Value  dj             = getDefaultCameraJSON();
  Json::Value &camera         = dj["camera blocks"]["fooCamera"];
  camera["camera at element"] = 1;
  checkPhactoriStringValidParse(CameraAtElement, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraAtNode", "[cameraValid]")
{
  Json::Value  dj          = getDefaultCameraJSON();
  Json::Value &camera      = dj["camera blocks"]["fooCamera"];
  camera["camera at node"] = 1;
  checkPhactoriStringValidParse(CameraAtNode, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraAtElementDisplaced", "[cameraValid]")
{
  Json::Value  dj                          = getDefaultCameraJSON();
  Json::Value &camera                      = dj["camera blocks"]["fooCamera"];
  camera["camera at element displaced"][0] = 1;
  camera["camera at element displaced"][1] = -3.0;
  camera["camera at element displaced"][2] = 10.0;
  camera["camera at element displaced"][3] = 20.0;
  camera["look at element"]                = 1;
  checkPhactoriStringValidParse(CameraAtElementDisplaced, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraAtNodeDisplaced", "[cameraValid]")
{
  Json::Value  dj                       = getDefaultCameraJSON();
  Json::Value &camera                   = dj["camera blocks"]["fooCamera"];
  camera["camera at node displaced"][0] = 1;
  camera["camera at node displaced"][1] = -3.0;
  camera["camera at node displaced"][2] = 10.0;
  camera["camera at node displaced"][3] = 20.0;
  camera["look at node"]                = 1;
  checkPhactoriStringValidParse(CameraAtNodeDisplaced, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraUpVector", "[cameraValid]")
{
  Json::Value  dj        = getDefaultCameraJSON();
  Json::Value &camera    = dj["camera blocks"]["fooCamera"];
  camera["up vector"][0] = 0.0;
  camera["up vector"][1] = 1.0;
  camera["up vector"][2] = 2.0;
  checkPhactoriStringValidParse(CameraUpVector, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraFOV", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["camera fov"]                = 45.0;
  camera["look at relative point"][0] = 0.0;
  camera["look at relative point"][1] = 0.0;
  camera["look at relative point"][2] = 0.0;
  camera["look at absolute distance"] = 10.0;
  checkPhactoriStringValidParse(CameraFOV, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraImageNameAddon", "[cameraValid]")
{
  Json::Value  dj            = getDefaultCameraJSON();
  Json::Value &camera        = dj["camera blocks"]["fooCamera"];
  camera["image name addon"] = "_foo_";
  checkPhactoriStringValidParse(CameraImageNameAddon, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraParallelProjection1", "[notWorking]")
{
  Json::Value dj = getDefaultCameraParallelProjectionJSON();
  checkPhactoriStringValidParse(CameraParallelProjection1, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8", "[cameraValid]")
{
  Json::Value  dj       = getDefaultCameraJSON();
  Json::Value &camera   = dj["camera blocks"]["fooCamera"];
  camera["camera type"] = "multicamera8";
  checkPhactoriStringValidParse(MultiCamera8, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8LookAtAbsolutePoint", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["camera type"]               = "multicamera8";
  camera["look at absolute point"][0] = 1.1;
  camera["look at absolute point"][1] = 2.0;
  camera["look at absolute point"][2] = 3e-8;
  checkPhactoriStringValidParse(MultiCamera8LookAtAbsolutePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8LookAtRelativePoint", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["camera type"]               = "multicamera8";
  camera["look at relative point"][0] = 0.5;
  camera["look at relative point"][1] = -0.5;
  camera["look at relative point"][2] = 0.5;
  checkPhactoriStringValidParse(MultiCamera8LookAtRelativePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8LookAtElement", "[cameraValid]")
{
  Json::Value  dj           = getDefaultCameraJSON();
  Json::Value &camera       = dj["camera blocks"]["fooCamera"];
  camera["camera type"]     = "multicamera8";
  camera["look at element"] = 17;
  checkPhactoriStringValidParse(MultiCamera8LookAtElement, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8LookAtNode", "[cameraValid]")
{
  Json::Value  dj        = getDefaultCameraJSON();
  Json::Value &camera    = dj["camera blocks"]["fooCamera"];
  camera["camera type"]  = "multicamera8";
  camera["look at node"] = 20;
  checkPhactoriStringValidParse(MultiCamera8LookAtNode, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8LookAtRelativeDistance", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["camera type"]               = "multicamera8";
  camera["look at relative distance"] = 2.0;
  checkPhactoriStringValidParse(MultiCamera8LookAtRelativeDistance, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8LookAtAbsoluteDistance", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["camera type"]               = "multicamera8";
  camera["look at absolute distance"] = 15.0;
  checkPhactoriStringValidParse(MultiCamera8LookAtAbsoluteDistance, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8UpVector", "[cameraValid]")
{
  Json::Value  dj        = getDefaultCameraJSON();
  Json::Value &camera    = dj["camera blocks"]["fooCamera"];
  camera["camera type"]  = "multicamera8";
  camera["up vector"][0] = 0.0;
  camera["up vector"][1] = 1.0;
  camera["up vector"][2] = 2.0;
  checkPhactoriStringValidParse(MultiCamera8UpVector, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8FOV", "[cameraValid]")
{
  Json::Value  dj                     = getDefaultCameraJSON();
  Json::Value &camera                 = dj["camera blocks"]["fooCamera"];
  camera["camera type"]               = "multicamera8";
  camera["camera fov"]                = 45.0;
  camera["look at relative point"][0] = 0.0;
  camera["look at relative point"][1] = 0.0;
  camera["look at relative point"][2] = 0.0;
  camera["look at absolute distance"] = 10.0;
  checkPhactoriStringValidParse(MultiCamera8FOV, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8ImageNameAddon", "[cameraValid]")
{
  Json::Value  dj            = getDefaultCameraJSON();
  Json::Value &camera        = dj["camera blocks"]["fooCamera"];
  camera["camera type"]      = "multicamera8";
  camera["image name addon"] = "_foo_";
  checkPhactoriStringValidParse(MultiCamera8ImageNameAddon, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8ParallelProjection1", "[notWorking]")
{
  Json::Value  dj            = getDefaultCameraParallelProjectionJSON();
  Json::Value &camParallel   = dj["camera blocks"]["parallel_projection_cam1"];
  camParallel["camera type"] = "multicamera8";
  camParallel.removeMember("look direction");

  Json::Value &camPerspective   = dj["camera blocks"]["perspective_projection_cam1"];
  camPerspective["camera type"] = "multicamera8";
  camPerspective.removeMember("look direction");
  checkPhactoriStringValidParse(MultiCamera8ParallelProjection1, dj);
}
