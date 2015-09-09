/*
 * MultiImageReader.h
 *
 *  Created on: Sep 2, 2015
 *      Author: jonchu
 */

#ifndef _MULTIIMAGEREADER_HPP_
#define _MULTIIMAGEREADER_HPP_

#include <GameOfLife/PNGProcessor.hpp>

typedef std::pair<unsigned, unsigned> CoordinatePair;
typedef std::vector<std::pair<unsigned, unsigned>> CoordinateLayer;
typedef std::map<unsigned, CoordinateLayer> CoordinateMap;

struct MultiImageReader
{
public:
    MultiImageReader(std::string baseName, unsigned numFiles);

    ~MultiImageReader(){}

    void clease();

    //stuff
    std::vector<CoordinateLayer> coordinates;

    CoordinateMap xToyz;

    CoordinateMap yToxz;

    CoordinateMap zToxy;

private:

    //constructor
    void construct_images(std::string baseName, std::vector<ColoredPNGProcessor>& images,
                          unsigned numFiles);

    void get_active_coordinates(std::vector<ColoredPNGProcessor>& images, unsigned numFiles);

        void add_triplet_to_maps(unsigned x, unsigned y, unsigned z);

        void get_coordinate_depths();

};

//pulic
MultiImageReader::MultiImageReader(std::string baseName, unsigned numFiles)
{
    std::vector<ColoredPNGProcessor> images;
    construct_images(baseName, images, numFiles);
    get_active_coordinates(images, numFiles);
}
void MultiImageReader::clease()
{
}

//private
void MultiImageReader::construct_images(std::string baseName, std::vector<ColoredPNGProcessor>& images,
                                        unsigned numFiles)
{
    for (unsigned index = 0; index < numFiles; index++)
    {
        std::string fileName = baseName + std::to_string(index) + ".png";
        images.emplace_back(ColoredPNGProcessor(fileName));
        images[index].commit_image_vector_to_pixel_vector_with_greyscale();
    }
}
void MultiImageReader::get_active_coordinates(std::vector<ColoredPNGProcessor>& images, unsigned numFiles)
{
    coordinates.resize(numFiles);
    for (unsigned index = 0; index < numFiles; index++)
    {
        images[index].get_coordinates_of_active_pixels(coordinates[index]);
        for (CoordinatePair xyPair : coordinates[index])
            add_triplet_to_maps(xyPair.first, xyPair.second, index+1);
    }
}
void MultiImageReader::add_triplet_to_maps(unsigned x, unsigned y, unsigned z)
{
    xToyz[x].emplace_back(CoordinatePair(y,z));
    yToxz[y].emplace_back(CoordinatePair(x,z));
    zToxy[z].emplace_back(CoordinatePair(x,y));
}

#endif /* _MULTIIMAGEREADER_HPP_ */

