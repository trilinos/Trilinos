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

    unsigned get_ideal_rows_per_processor(unsigned numProcs) const;

    virtual void fill_id_vector_with_active_pixels(stk::mesh::EntityIdVector& elemIds) const;

    void get_coordinates_of_active_pixels(std::vector<std::pair<unsigned, unsigned>>& coordinates) const;

    void get_coordinates_of_inactive_pixels(std::vector<std::pair<unsigned, unsigned>>& coordinates) const;

    //data dump
    void print_image();

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

class ColoredPNGProcessor : public PNGProcessor
{
public:
    ColoredPNGProcessor(std::string fileName);

    virtual ~ColoredPNGProcessor() {}

    virtual void commit_image_vector_to_pixel_vector();

    void commit_image_vector_to_pixel_vector_with_exclusion();

    void commit_image_vector_to_pixel_vector_with_greyscale();

    //io
    void print_grey_bits();

private:
    // members
    unsigned char m_medianValue;
    unsigned char m_lowerBound;
    unsigned char m_upperBound;
    std::vector<unsigned char> m_greyBits;
    std::vector<Pixel> mRedPixels;
    std::vector<Pixel> mGreenPixels;
    std::vector<Pixel> mBluePixels;

    // constructor
    void convert_to_grey_bits();
    void process_unsigned_int_to_grey_bit(unsigned row, unsigned col);

    void find_median_grey_bit_value();

    // commit image
    void update_image_value_according_to_relation_with_median_value(unsigned row, unsigned col);
    void update_image_value_according_to_proximity_with_median_value(unsigned row, unsigned col);
    void update_image_value_according_to_grayscale(unsigned row, unsigned col);

    void find_upper_and_lower_bounds();
};

class SimpleColoredPng : public PNGProcessor
{
public:
    SimpleColoredPng(std::string fileName);
    virtual ~SimpleColoredPng() {}

    virtual void fill_id_vector_with_active_pixels(stk::mesh::EntityIdVector& elemIds) const;

    std::vector<Pixel> get_red_color_coords() { return mRedPixels; }
    std::vector<Pixel> get_green_color_coords() { return mGreenPixels; }
    std::vector<Pixel> get_blue_color_coords() { return mBluePixels; }

    stk::mesh::EntityIdVector get_elemIds_for_colored_pixels(const std::vector<Pixel> & coloredPixels);
private:
    void store_special_colors_with_coordinates(unsigned row, unsigned col);
    void update_image_value_ignoring_white(unsigned row, unsigned col);
private:
    stk::mesh::EntityIdVector mRedElementIds;

    std::vector<Pixel> mRedPixels;
    std::vector<Pixel> mGreenPixels;
    std::vector<Pixel> mBluePixels;
};


#endif /*_PNGPROCESSOR_HPP_ */
