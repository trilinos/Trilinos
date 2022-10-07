// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __OPERATIONS_VALID_TESTS_H
#define __OPERATIONS_VALID_TESTS_H

#include <string>

std::string OperationBaseline = R"(
     begin catalyst
      begin clip clipNone
        relative point on plane = -0.5 0 0
        plane normal = 1 0 0
        side to keep = positive
      end clip
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = clipNone
        image size = 800 450
      end imageset
    end
     )";

std::string OperationClip = R"(
     begin catalyst
      begin clip clipNone
        relative point on plane = -0.5 0 0
        plane normal = 1 0 0
        side to keep = positive
      end clip
      begin clip fooOperation
        input = clipNone
      end clip
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = fooOperation
        image size = 800 450
      end imageset
    end
     )";

std::string OperationClipAbsolutePoint = R"(
     begin catalyst
      begin clip clipNone
        relative point on plane = -0.5 0 0
        plane normal = 1 0 0
        side to keep = positive
      end clip
      begin clip fooOperation
        input = clipNone
        absolute point on plane = 3.0 0.1 0.2
        plane normal = 1 1 0
      end clip
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = fooOperation
        image size = 800 450
      end imageset
    end
     )";

std::string OperationClipRelativePoint = R"(
    begin catalyst
      begin clip clipNone
        relative point on plane = -0.5 0 0
        plane normal = 1 0 0
        side to keep = positive
      end clip
      begin clip fooOperation
        input = clipNone
        #relative point on plane = 0.25 0.1 0.05
        relative point on plane = 0.25 0.3 0.5
        plane normal = 1 1 0
      end clip
            begin camera singletestcamera
        look direction = -1 -1 -1
        look at absolute distance = 10.0
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = fooOperation
        image size = 800 450
      end imageset
    end
     )";

std::string OperationClipNode = R"(
    begin catalyst
      begin clip clipNone
        relative point on plane = -0.5 0 0
        plane normal = 1 0 0
        side to keep = positive
      end clip
      begin clip fooOperation
        input = clipNone
        node on plane = 20
        plane normal = 1.0 1.0 0.0
      end clip
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = fooOperation
        image size = 800 450
      end imageset
    end
     )";

std::string OperationClipElement = R"(
    begin catalyst
      begin clip clipNone
        relative point on plane = -0.5 0 0
        plane normal = 1 0 0
        side to keep = positive
      end clip
      begin clip fooOperation
        input = clipNone
        element on plane = 17
        plane normal = 1.0 0.25 0.0
      end clip
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = fooOperation
        image size = 800 450
      end imageset
    end
     )";

std::string OperationClipPlaneNormal = R"(
    begin catalyst
      begin clip clipNone
        relative point on plane = -0.5 0 0
        plane normal = 1 0 0
        side to keep = positive
      end clip
      begin clip fooOperation
        input = clipNone
        plane normal = 1 2 3
        absolute point on plane = 3.0 0.1 0.2
      end clip
      begin camera singletestcamera
        look direction = -1 -1 -1
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = fooOperation
        image size = 800 450
      end imageset
    end
     )";

std::string OperationClipPlaneThreePoints1 = R"(
    begin catalyst
      begin extractblock eb1
        include blocks = block_1
      end
      begin clip fooOperation
        #input = clipNone
        input = eb1
        plane specification = three points
        data point on plane A = max scalar VON_MISES eb1
        data point on plane B = min scalar VON_MISES eb1
        relative point on plane C = 0.0 10.0 0.0
        side to keep = negative
      end clip
      begin clip clipop2
        plane specification = three points
        data point on plane A = max scalar VON_MISES eb1
        relative point on plane B = 0.0 0.0 10.0
        relative point on plane C = 0.0 10.0 0.0
      end
      begin camera singletestcamera
        #look direction = -1 -1 -1
        camera at absolute point = 15 10 10
        look at absolute point = 5 0 0
      end
      begin imageset fooImageset
        camera = singletestcamera
        operation = fooOperation
        image size = 800 450
        show edges = true
        color by scalar = VON_MISES
      end imageset
      begin imageset is2
        camera = singletestcamera
        operation = clipop2
        image size = 800 450
        show edges = true
        color by scalar = VON_MISES
      end imageset
      begin imageset is3
        camera = singletestcamera
        image size = 800 450
        show edges = true
        color by scalar = VON_MISES
      end imageset
    end
     )";

#endif /* __OPERATIONS_VALID_TESTS_H */
