// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "catch.hpp"

TEST_CASE_METHOD(CatalystTestFixture,
    "Camera", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraLookAtAbsolutePoint", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         look at absolute point = 1.1 2 3e-8
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["look at absolute point"][0] = 1.1;
    camera["look at absolute point"][1] = 2.0;
    camera["look at absolute point"][2] = 3e-8;

    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraLookDirection", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         look direction = 1 2 3
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["look direction"][0] = 1.0;
    camera["look direction"][1] = 2.0;
    camera["look direction"][2] = 3.0;

    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraLookAtRelativeDistance", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         look at relative distance = 2.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["look at relative distance"] = 2.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraLookAtAbsoluteDistance", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         look at absolute distance = 15.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["look at absolute distance"] = 15.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraAtAbsolutePoint", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         camera at absolute point = -2.0 3.0 30.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera at absolute point"][0] = -2.0;
    camera["camera at absolute point"][1] = 3.0;
    camera["camera at absolute point"][2] = 30.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraAtRelativePoint", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         camera at relative point = -0.5 1.5 20.0
         look direction = 0.1 -0.1 -1.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera at relative point"][0] = -0.5;
    camera["camera at relative point"][1] = 1.5;
    camera["camera at relative point"][2] = 20.0;
    camera["look direction"][0] = 0.1;
    camera["look direction"][1] = -0.1;
    camera["look direction"][2] = -1.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraAtElement", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         camera at element = 1
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera at element"] = 1;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraAtNode", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         camera at node = 1
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera at node"] = 1;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraAtElementDisplaced", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         camera at element displaced = 1 -3.0 10.0 20.0
         look at element = 1
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera at element displaced"][0] = 1;
    camera["camera at element displaced"][1] = -3.0;
    camera["camera at element displaced"][2] = 10.0;
    camera["camera at element displaced"][3] = 20.0;
    camera["look at element"] = 1;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraAtNodeDisplaced", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         camera at node displaced = 1 -3.0 10.0 20.0
         look at node = 1
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera at node displaced"][0] = 1;
    camera["camera at node displaced"][1] = -3.0;
    camera["camera at node displaced"][2] = 10.0;
    camera["camera at node displaced"][3] = 20.0;
    camera["look at node"] = 1;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraUpVector", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         up vector = 0 1 2
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["up vector"][0] = 0.0;
    camera["up vector"][1] = 1.0;
    camera["up vector"][2] = 2.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraFOV", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         camera fov = 45
         look at relative point = 0.0 0.0 0.0
         look at absolute distance = 10.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera fov"] = 45.0;
    camera["look at relative point"][0] = 0.0;
    camera["look at relative point"][1] = 0.0;
    camera["look at relative point"][2] = 0.0;
    camera["look at absolute distance"] = 10.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraImageNameAddon", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin camera fooCamera
         image name addon = "_foo_"
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["image name addon"] = "_foo_";
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "CameraParallelProjection1", "[notWorking]") {

     std::string ps = R"(
     begin catalyst
       begin camera parallel_projection_cam1
         projection type = parallel
         look direction = -5 -1 -1
       end camera
       begin camera perspective_projection_cam1
         projection type = perspective
         look direction = -5 -1 -1
       end camera
       begin imageset parallel_projection_is1
         camera = parallel_projection_cam1
         image size = 800 450
         image basename = parallel_projection_is1.
       end imageset
       begin imageset perspective_projection_is1
         camera = perspective_projection_cam1
         image size = 800 450
         image basename = perspective_projection_is1.
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraParallelProjectionJSON();
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8LookAtAbsolutePoint", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at absolute point = 1.1 2 3e-8
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["look at absolute point"][0] = 1.1;
    camera["look at absolute point"][1] = 2.0;
    camera["look at absolute point"][2] = 3e-8;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8LookAtRelativePoint", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at relative point = 0.5 -0.5 0.5
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["look at relative point"][0] = 0.5;
    camera["look at relative point"][1] = -0.5;
    camera["look at relative point"][2] = 0.5;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8LookAtElement", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at element = 17
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["look at element"] = 17;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8LookAtNode", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at node = 20
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["look at node"] = 20;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8LookAtRelativeDistance", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at relative distance = 2.0
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["look at relative distance"] = 2.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8LookAtAbsoluteDistance", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at absolute distance = 15.0
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["look at absolute distance"] = 15.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8UpVector", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         up vector = 0 1 2
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["up vector"][0] = 0.0;
    camera["up vector"][1] = 1.0;
    camera["up vector"][2] = 2.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8FOV", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         camera fov = 45
         look at relative point = 0.0 0.0 0.0
         look at absolute distance = 10.0
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["camera fov"] = 45.0;
    camera["look at relative point"][0] = 0.0;
    camera["look at relative point"][1] = 0.0;
    camera["look at relative point"][2] = 0.0;
    camera["look at absolute distance"] = 10.0;
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8ImageNameAddon", "[cameraValid]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 fooCamera
         image name addon = "_foo_"
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

    Json::Value dj = getDefaultCameraJSON();
    Json::Value& camera = dj["camera blocks"]["fooCamera"];
    camera["camera type"] = "multicamera8";
    camera["image name addon"] = "_foo_";
    checkPhactoriStringValidParse(ps, dj);
}

TEST_CASE_METHOD(CatalystTestFixture,
    "MultiCamera8ParallelProjection1", "[notWorking]") {

     std::string ps = R"(
     begin catalyst
       begin multicamera8 parallel_projection_cam1
         projection type = parallel
       end multicamera8
       begin multicamera8 perspective_projection_cam1
         projection type = perspective
       end multicamera8
       begin imageset parallel_projection_is1
         camera = parallel_projection_cam1
         image size = 800 450
         image basename = parallel_projection_is1.
       end imageset
       begin imageset perspective_projection_is1
         camera = perspective_projection_cam1
         image size = 800 450
         image basename = perspective_projection_is1.
       end imageset
     end
     )";


    Json::Value dj = getDefaultCameraParallelProjectionJSON();
    Json::Value& camParallel = dj["camera blocks"]["parallel_projection_cam1"];
    camParallel["camera type"] = "multicamera8";
    camParallel.removeMember("look direction");

    Json::Value& camPerspective =
        dj["camera blocks"]["perspective_projection_cam1"];
    camPerspective["camera type"] = "multicamera8";
    camPerspective.removeMember("look direction");
    checkPhactoriStringValidParse(ps, dj);
}

