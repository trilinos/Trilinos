/*
 * MultiImageReader.h
 *
 *  Created on: Sep 2, 2015
 *      Author: jonchu
 */

#ifndef _MULTIIMAGEREADER_HPP_
#define _MULTIIMAGEREADER_HPP_

#include <GameOfLife/PNGProcessor.hpp>
#include <MeshBuilder/MeshBuilder.hpp>

typedef std::pair<unsigned, unsigned> CoordinatePair;
typedef std::vector<CoordinatePair> CoordinateLayer;
typedef std::map<unsigned, CoordinateLayer> CoordinateMap;
typedef std::pair<unsigned,CoordinateLayer> DimensionLayer;

struct MultiImageReader
{
public:
    MultiImageReader(std::string baseName, unsigned numFiles)
    {
        std::vector<ColoredPNGProcessor> images;
        construct_images(baseName, images, numFiles);
        get_active_coordinates(images, numFiles);
    }

    ~MultiImageReader(){}

    void create_randomly_decomposed_mesh(HexMeshBuilder& mesh) const
    //requires mod_begin(obv)
    {
        unsigned procCounter = 0;
        unsigned numProcs = mesh.num_procs();

        for (const DimensionLayer& dimLayer : zToxy)
            for (const CoordinatePair& xy : dimLayer.second)
                mesh.create_element(xy.first, xy.second, dimLayer.first, procCounter++%numProcs);
    }
    void create_x_layered_decomposed_mesh(HexMeshBuilder& mesh) const
    {
        unsigned procCounter= 0;
        unsigned numProcs = mesh.num_procs();

        for (const DimensionLayer& dimLayer : xToyz)
        {
            for (const CoordinatePair& yz : dimLayer.second)
                mesh.create_element(dimLayer.first, yz.first, yz.second, procCounter%numProcs);
            procCounter++;
        }
    }
    void create_y_layered_decomposed_mesh(HexMeshBuilder& mesh) const
    {
        unsigned procCounter= 0;
        unsigned numProcs = mesh.num_procs();

        for (const DimensionLayer& dimLayer : yToxz)
        {
            for (const CoordinatePair& xz : dimLayer.second)
                mesh.create_element(xz.first, dimLayer.first, xz.second, procCounter%numProcs);
            procCounter++;
        }
    }
    void create_z_layered_decomposed_mesh(HexMeshBuilder& mesh) const
    {
        unsigned procCounter= 0;
        unsigned numProcs = mesh.num_procs();

        for (const DimensionLayer& dimLayer : zToxy)
        {
            for (const CoordinatePair& xy : dimLayer.second)
                mesh.create_element(xy.first, xy.second, dimLayer.first, procCounter%numProcs);
            procCounter++;
        }
    }
    void create_x_blocked_decomposed_mesh(HexMeshBuilder& mesh) const
    {
        unsigned numLayers = xToyz.size();
        unsigned numProcs = mesh.num_procs();
        unsigned basePerProc = numLayers/numProcs;
        unsigned leftOver = numLayers%numProcs;

        std::vector<unsigned> layersPerProc(numProcs);
        for (unsigned index = 0; index < numProcs; index++)
        {
            layersPerProc[index] = basePerProc;
            if (leftOver > index)
                layersPerProc[index]++;
        }
        CoordinateMap::const_iterator mapIter = xToyz.cbegin();

        for (unsigned procIndex = 0; procIndex < numProcs; procIndex++)
            for (unsigned layerIndex = 0; layerIndex < layersPerProc[procIndex]; layerIndex++, mapIter++)
                for (const CoordinatePair& yz : (*mapIter).second)
                    mesh.create_element((*mapIter).first, yz.first, yz.second, procIndex);
    }
    void create_y_blocked_decomposed_mesh(HexMeshBuilder& mesh)
    {
        unsigned numLayers = yToxz.size();
        unsigned numProcs = mesh.num_procs();
        unsigned basePerProc = numLayers/numProcs;
        unsigned leftOver = numLayers%numProcs;

        std::vector<unsigned> layersPerProc(numProcs);
        for (unsigned index = 0; index < numProcs; index++)
        {
            layersPerProc[index] = basePerProc;
            if (leftOver > index)
                layersPerProc[index]++;
        }
        CoordinateMap::const_iterator mapIter = yToxz.cbegin();

        for (unsigned procIndex = 0; procIndex < numProcs; procIndex++)
            for (unsigned layerIndex = 0; layerIndex < layersPerProc[procIndex]; layerIndex++, mapIter++)
                for (const CoordinatePair& xz : (*mapIter).second)
                    mesh.create_element(xz.first, (*mapIter).first, xz.second, procIndex);
    }
    void create_z_blocked_decomposed_mesh(HexMeshBuilder& mesh)
    {
        unsigned numLayers = zToxy.size();
        unsigned numProcs = mesh.num_procs();
        unsigned basePerProc = numLayers/numProcs;
        unsigned leftOver = numLayers%numProcs;

        std::vector<unsigned> layersPerProc(numProcs);
        for (unsigned index = 0; index < numProcs; index++)
        {
            layersPerProc[index] = basePerProc;
            if (leftOver > index)
                layersPerProc[index]++;
        }
        CoordinateMap::const_iterator mapIter = zToxy.cbegin();

        for (unsigned procIndex = 0; procIndex < numProcs; procIndex++)
            for (unsigned layerIndex = 0; layerIndex < layersPerProc[procIndex]; layerIndex++, mapIter++)
                for (const CoordinatePair& xy : (*mapIter).second)
                    mesh.create_element(xy.first, xy.second, (*mapIter).first, procIndex);
    }

private:
    //stuff
    CoordinateMap xToyz;

    CoordinateMap yToxz;

    CoordinateMap zToxy;

    //constructor
    void construct_images(std::string baseName, std::vector<ColoredPNGProcessor>& images,
                          unsigned numFiles)
    {

        for (unsigned index = 0; index < numFiles; index++)
        {
            std::string fileName = baseName + std::to_string(index) + ".png";
            images.emplace_back(ColoredPNGProcessor(fileName));
            images[index].commit_image_vector_to_pixel_vector_with_greyscale();
        }
    }
    void get_active_coordinates(std::vector<ColoredPNGProcessor>& images, unsigned numFiles)
    {
        for (unsigned index = 0; index < numFiles; index++)
        {
            CoordinateLayer layer;
            images[index].get_coordinates_of_inactive_pixels(layer);
            for (CoordinatePair xyPair : layer)
                add_triplet_to_maps(xyPair.first, xyPair.second, index+1);
        }
    }
    void add_triplet_to_maps(unsigned x, unsigned y, unsigned z)
    {
        xToyz[x].emplace_back(CoordinatePair(y,z));
        yToxz[y].emplace_back(CoordinatePair(x,z));
        zToxy[z].emplace_back(CoordinatePair(x,y));
    }
};

#endif /* _MULTIIMAGEREADER_HPP_ */

