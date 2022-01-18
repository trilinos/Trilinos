// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "catch.hpp"

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidMultipleLookAt", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         look at absolute point = 0 0 0
         look at relative point = 0 0 0.5
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidLookAtInvalidElement", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         look at element = -1
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidLookAtInvalidNode", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         look at node = -1
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidLookDirectionZeros", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         look direction = 0 0 0
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidMultipleLookDistance", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         look at relative distance = 1
         look at absolute distance = 1
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidLookAtInvalidRelativeDistance",
                 "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         look at relative distance = -1
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidLookAtInvalidAbsoluteDistance",
                 "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         look at absolute distance = 0
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidMultipleCameraPosition", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         camera at absolute point = 0 0 0
         camera at element = 1
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidAtInvalidNode", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         camera at node = -1
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidUpVectorZeros", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         up vector = 0 0 0
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidUpVectorColinear", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         up vector = 1 0 0
         camera at = 2 0 0
         look at = 1 0 0
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidInvalidFOV", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         camera fov = 0
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidInvalidFOV2", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         camera fov = 180.1
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "CameraInvalidInvalidImageName", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin camera MyCamera
         image name addon = "$/%foo"
       end camera MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidMultipleLookAt", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         look at absolute point = 0 0 0
         look at relative point = 0 0 0.5
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidLookAtInvalidElement", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         look at element = -1
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidLookAtInvalidNode", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         look at node = -1
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidMultipleLookDistance", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         look at relative distance = 1
         look at absolute distance = 1
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidLookAtInvalidRelativeDistance",
                 "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         look at relative distance = -1
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidLookAtInvalidAbsoluteDistance",
                 "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         look at absolute distance = 0
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidUpVectorZeros", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         up vector = 0 0 0
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidInvalidFOV", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         camera fov = 0
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidInvalidFOV2", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         camera fov = 180.1
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "MultiCamera8InvalidInvalidImageName", "[cameraInvalid]")
{

  std::string ps = R"(
     begin catalyst
       begin multicamera8 MyCamera
         image name addon = "$/%foo"
       end multicamera8 MyCamera
     end
     )";

  checkPhactoriStringInvalidParse(ps);
}
