#include <vector>
#include <string>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include "LodePNG.hpp"

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

    void fill_id_vector_with_active_pixels(stk::mesh::EntityIdVector& elemIds) const;

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
class ColoredPNGProcessor : public PNGProcessor
{
public:
    ColoredPNGProcessor(std::string fileName);

    virtual ~ColoredPNGProcessor() {}

    virtual void commit_image_vector_to_pixel_vector();

    void commit_image_vector_to_pixel_vector_with_exclusion();

private:
    // members
    unsigned char m_medianValue;
    unsigned char m_lowerBound;
    unsigned char m_upperBound;
    std::vector<unsigned char> m_greyBits;

    // constructor
    void convert_to_grey_bits();
    void process_unsigned_int_to_grey_bit(unsigned row, unsigned col);

    void find_medium_grey_bit_value();

    // commit image
    void update_image_value_according_to_relation_with_median_value(unsigned row, unsigned col);

    void update_image_value_according_to_proximity_with_median_value(unsigned row, unsigned col);

    void find_upper_and_lower_bounds();
};
