// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 // 
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 // 
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef _PNGPROCESSOR_HPP_
#define _PNGPROCESSOR_HPP_

#include <stk_mesh/base/Types.hpp>      // for EntityIdVector
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector


/*
 * How to use:
 * Pass the file name of the image (make sure it's in the same directory, else segfault),
 * and call the commit_image_vector function. There's another option if it's colored. Don't
 * use the bordered one unless the image pulled is specifically from ConwayLife, because that
 * is the only image type it really supports, and the regular one only takes in black and white,
 * however the colored one can take in literally anything. The padding is pretty self-explanatory,
 * and you can get the ids if you're using a quad mesh by using the fill_element_id_vector thing.
 */

class PNGProcessor
{
public:
    PNGProcessor(std::string fileName);

    virtual ~PNGProcessor() {}

    virtual void commit_image_vector_to_pixel_vector();

    void add_this_much_pixel_padding_to_right(unsigned amount);

    void add_this_much_pixel_padding_to_left(unsigned amount);

    void add_this_much_pixel_padding_to_top(unsigned amount);

    void add_this_much_pixel_padding_to_bottom(unsigned amount);

    //accessors and other stuff
    inline unsigned get_image_width() const;

    inline unsigned get_image_height() const;

    virtual void fill_id_vector_with_active_pixels(stk::mesh::EntityIdVector& elemIds) const;

protected:
    unsigned m_imageWidth;

    unsigned m_imageHeight;

    std::vector<unsigned char> m_byteVector;

    std::vector<std::vector<unsigned>> m_imageVector;

    std::vector<std::vector<bool>> m_pixelVector;

private:
    //constructor
    void process_bytes_into_image_vector();
    void process_image_vector_row(unsigned row, unsigned& byteVectorIndex);
    void smear_character_bytes_into_unsigned_int(unsigned row, unsigned col,
                                                 unsigned& byteVectorIndex);

    //convert image vector to pixel vector
    void commit_image_row_to_pixel_row(unsigned row);
};
inline unsigned PNGProcessor::get_image_width() const
{
   return m_imageWidth;
}
inline unsigned PNGProcessor::get_image_height() const
{
    return m_imageHeight;
}
class BorderedPNGProcessor : public PNGProcessor
{
public:
    BorderedPNGProcessor(std::string fileName);

    virtual ~BorderedPNGProcessor() {}

    void shrink_image();

private:
    //random variables
    unsigned m_numVerticalLines;
    unsigned m_numHorizontalLines;
    unsigned m_squareHeight;
    unsigned m_squareWidth;

    //compress_bordered_image
    void count_vertical_and_horizontal_lines();

    void declare_new_height_and_width();

    void shrink_image_vector();
    void shrink_image_row(unsigned row);
    void shrink_image_byte(unsigned row, unsigned col, unsigned verticalOffset);
};

struct Pixel
{
    unsigned x;
    unsigned y;
};

class SimpleColoredPng : public PNGProcessor
{
public:
    SimpleColoredPng(std::string fileName);
    virtual ~SimpleColoredPng() {}

    virtual void fill_id_vector_with_active_pixels(stk::mesh::EntityIdVector& elemIds) const;

    size_t get_number_other_pixels() { return mOtherPixelCount; }

    std::vector<Pixel> get_red_color_coords() { return mRedPixels; }
    std::vector<Pixel> get_green_color_coords() { return mGreenPixels; }
    std::vector<Pixel> get_blue_color_coords() { return mBluePixels; }

    stk::mesh::EntityIdVector get_elemIds_for_colored_pixels(const std::vector<Pixel> & coloredPixels);
private:
    void store_special_colors_with_coordinates(unsigned row, unsigned col);
    void update_image_value_ignoring_white(unsigned row, unsigned col);
private:
    size_t mOtherPixelCount;
    stk::mesh::EntityIdVector mRedElementIds;

    std::vector<Pixel> mRedPixels;
    std::vector<Pixel> mGreenPixels;
    std::vector<Pixel> mBluePixels;
};


#endif /*_PNGPROCESSOR_HPP_ */
