#include "PNGProcessor.hpp"
#include <stdio.h>                      // for printf
#include <algorithm>                    // for sort
#include "GameOfLife/LodePNG.hpp"       // for decode
#include "stk_mesh/base/Types.hpp"      // for EntityIdVector

unsigned char greyDivisor = 0x2f;

PNGProcessor::PNGProcessor(std::string fileName)
{
    lodepng::decode(m_byteVector, m_imageWidth, m_imageHeight, fileName);
    process_bytes_into_image_vector();
}
void PNGProcessor::commit_image_vector_to_pixel_vector()
{
    m_pixelVector.resize(m_imageHeight);
    for (unsigned rowIndex = 0; rowIndex < m_imageHeight; rowIndex++)
        commit_image_row_to_pixel_row(rowIndex);
    m_imageVector.clear(); // rest in pepperonis
}
void PNGProcessor::add_this_much_pixel_padding_to_right(unsigned amount)
{
    m_imageWidth += amount;
    for (unsigned rowIndex = 0; rowIndex < m_imageHeight; rowIndex++)
        m_pixelVector[rowIndex].insert(m_pixelVector[rowIndex].end(), amount, false);
}

void PNGProcessor::add_this_much_pixel_padding_to_left(unsigned amount)
{
    m_imageWidth += amount;
    for (unsigned rowIndex = 0; rowIndex < m_imageHeight; rowIndex++)
        m_pixelVector[rowIndex].insert(m_pixelVector[rowIndex].begin(), amount, false);
}

void PNGProcessor::add_this_much_pixel_padding_to_top(unsigned amount)
{
    m_imageHeight += amount;
    std::vector<bool> example(m_imageWidth, false);
    m_pixelVector.insert(m_pixelVector.begin(), amount, example);
}

void PNGProcessor::add_this_much_pixel_padding_to_bottom(unsigned amount)
{
    m_imageHeight += amount;
    std::vector<bool> example(m_imageWidth, false);
    m_pixelVector.insert(m_pixelVector.end(), amount, example);
}

//IO and Checking
//unsigned PNGProcessor::get_ideal_rows_per_processor(unsigned numProcs) const
//{
//    return m_imageHeight/numProcs + !!(m_imageHeight%numProcs);
//}

void PNGProcessor::fill_id_vector_with_active_pixels(stk::mesh::EntityIdVector& elemIds) const
{
    unsigned id = 1;
    for (int rowIndex = m_imageHeight-1; rowIndex >= 0; rowIndex--)
        for (unsigned colIndex = 0; colIndex < m_imageWidth; colIndex++, id++)
            if (m_pixelVector[rowIndex][colIndex])
                elemIds.push_back(id);
}

//void PNGProcessor::get_coordinates_of_active_pixels(std::vector<std::pair<unsigned, unsigned>>&
//                                                    coordinates) const
//{
//    coordinates.clear();
//    for (unsigned rowIndex = 0; rowIndex < m_imageHeight; rowIndex++)
//        for (unsigned colIndex = 0; colIndex < m_imageWidth; colIndex++)
//            if (m_pixelVector[rowIndex][colIndex])
//                coordinates.emplace_back(colIndex+1, m_imageHeight-rowIndex);
//}

void PNGProcessor::get_coordinates_of_inactive_pixels(std::vector<std::pair<unsigned,
                                                      unsigned>>& coordinates) const
{
    coordinates.clear();
    for (unsigned rowIndex = 0; rowIndex < m_imageHeight; rowIndex++)
        for (unsigned colIndex = 0; colIndex < m_imageWidth; colIndex++)
            if (!m_pixelVector[rowIndex][colIndex])
                coordinates.emplace_back(colIndex+1, m_imageHeight-rowIndex);
}

//void PNGProcessor::print_image()
//{
//    for (unsigned row = 0; row < m_imageHeight; row++)
//    {
//        for (unsigned col = 0; col < m_imageWidth;
//                col++)
//            printf("%8x ", m_imageVector[row][col]);
//        printf("\n");
//    }
//}

//private
void PNGProcessor::process_bytes_into_image_vector()
{
    m_imageVector.resize(m_imageHeight);
    unsigned byteVectorIndex = 0;
    for (unsigned row = 0; row < m_imageHeight; row++)
        process_image_vector_row(row, byteVectorIndex);
}
void PNGProcessor::process_image_vector_row(unsigned row, unsigned& byteVectorIndex)
{
    m_imageVector[row].resize(m_imageWidth);
    for (unsigned col = 0; col < m_imageWidth; col++)
        smear_character_bytes_into_unsigned_int(row, col, byteVectorIndex);
}
void PNGProcessor::smear_character_bytes_into_unsigned_int(unsigned row, unsigned col,
                                                           unsigned& byteVectorIndex)
{
    unsigned newValue = 0x00;
    for (int counter = 3; counter >= 0; counter --, byteVectorIndex++)
        newValue |= (m_byteVector[byteVectorIndex] << (counter*8));
    m_imageVector[row][col] = newValue;
}
void PNGProcessor::commit_image_row_to_pixel_row(unsigned row)
{
    m_pixelVector[row].resize(m_imageWidth);
    for (unsigned index = 0; index < m_imageWidth; index++)
        m_pixelVector[row][index] = (0xff == m_imageVector[row][index]);
}
/*
 * Bordered PNG Processor
 */
BorderedPNGProcessor::BorderedPNGProcessor(std::string fileName)
:PNGProcessor(fileName), m_numVerticalLines(0), m_numHorizontalLines(0)
{
    shrink_image();
}

void BorderedPNGProcessor::shrink_image()
{
    count_vertical_and_horizontal_lines();
    declare_new_height_and_width();
    shrink_image_vector();
}

void BorderedPNGProcessor::count_vertical_and_horizontal_lines()
{
    for (unsigned index = 0; index < m_imageWidth; index++)
        if (0xc6c6c6ff == m_imageVector[1][index])
            m_numVerticalLines++;
    for (unsigned index = 0; index < m_imageHeight; index++)
        if (0xc6c6c6ff == m_imageVector[index][1])
            m_numHorizontalLines++;
}
void BorderedPNGProcessor::declare_new_height_and_width()
{
    m_squareHeight = (m_imageHeight-m_numHorizontalLines)/(m_numHorizontalLines-1);
    m_squareWidth = (m_imageWidth-m_numVerticalLines)/(m_numVerticalLines-1);
    m_imageHeight = m_numHorizontalLines-1;
    m_imageWidth = m_numVerticalLines-1;
}

void BorderedPNGProcessor::shrink_image_vector()
{
    for (unsigned row = 0; row < m_imageHeight; row++)
        shrink_image_row(row);
    m_imageVector.resize(m_imageHeight);
}
void BorderedPNGProcessor::shrink_image_row(unsigned row)
{
    unsigned verticalOffset = m_squareHeight*row + row + 1;
    m_imageVector[row].resize(m_imageWidth);
    for (unsigned col = 0; col < m_imageWidth; col++)
        shrink_image_byte(row, col, verticalOffset);
}
void BorderedPNGProcessor::shrink_image_byte(unsigned row, unsigned col, unsigned verticalOffset)
{
    unsigned horizontalOffset = m_squareWidth*col + col + 1;
    m_imageVector[row][col] =
            m_imageVector[verticalOffset][horizontalOffset];
}
/*
 * Colored PNG processor
 */
ColoredPNGProcessor::ColoredPNGProcessor(std::string fileName)
:PNGProcessor(fileName)
{
    convert_to_grey_bits();
}
ColoredPNGProcessor::~ColoredPNGProcessor() { }

void ColoredPNGProcessor::commit_image_vector_to_pixel_vector()
{
    find_median_grey_bit_value();
    for (unsigned row = 0; row < m_imageHeight; row++)
        for (unsigned col = 0; col < m_imageWidth; col++)
            update_image_value_according_to_relation_with_median_value(row, col);
    PNGProcessor::commit_image_vector_to_pixel_vector();
}
void ColoredPNGProcessor::commit_image_vector_to_pixel_vector_with_exclusion()
{
    find_median_grey_bit_value();
    find_upper_and_lower_bounds();
    for (unsigned row = 0; row < m_imageHeight; row++)
        for (unsigned col = 0; col < m_imageWidth; col++)
            update_image_value_according_to_proximity_with_median_value(row, col);
    PNGProcessor::commit_image_vector_to_pixel_vector();
}

void ColoredPNGProcessor::commit_image_vector_to_pixel_vector_with_greyscale()
{
    for (unsigned row = 0; row < m_imageHeight; row++)
        for (unsigned col = 0; col < m_imageWidth; col++)
            update_image_value_according_to_grayscale(row, col);
    PNGProcessor::commit_image_vector_to_pixel_vector();
}

//io and checking
void ColoredPNGProcessor::print_grey_bits()
{
    for (unsigned row = 0; row < m_imageHeight; row++)
    {
        for (unsigned col = 0; col < m_imageWidth;
                col++)
            printf("%4x ", m_greyBits[row*m_imageWidth + col]);
        printf("\n");
    }
}
//private
void ColoredPNGProcessor::convert_to_grey_bits()
{
    for (unsigned row = 0; row < m_imageHeight; row++)
        for (unsigned col = 0; col < m_imageWidth; col++)
            process_unsigned_int_to_grey_bit(row, col);
}



void ColoredPNGProcessor::process_unsigned_int_to_grey_bit(unsigned row, unsigned col)
{
    unsigned mask = 0xff;
    unsigned imageVal = m_imageVector[row][col];
    unsigned char greyBit = (((imageVal>>24)&mask)*54 + ((imageVal>>16)&mask)*183 +((imageVal>>8)&mask)*19)/256;
    m_imageVector[row][col] = greyBit;
    m_greyBits.push_back(greyBit);
}

void ColoredPNGProcessor::find_median_grey_bit_value()
{
    std::sort(m_greyBits.begin(), m_greyBits.end());
    unsigned middleIndex = m_greyBits.size()/2;
    m_medianValue = m_greyBits[middleIndex];
}
void ColoredPNGProcessor::update_image_value_according_to_relation_with_median_value(unsigned row,
                                                                                     unsigned col)
{
    if (m_imageVector[row][col] <= m_medianValue)
        m_imageVector[row][col] = 0xff;
    else
        m_imageVector[row][col] = 0xffffffff;
}
void ColoredPNGProcessor::update_image_value_according_to_proximity_with_median_value(unsigned row,
                                                                                      unsigned col)
{
    if (m_imageVector[row][col] >= m_lowerBound && m_imageVector[row][col] <= m_upperBound)
        m_imageVector[row][col] = 0xff;
    else
        m_imageVector[row][col] = 0xffffffff;
}
void ColoredPNGProcessor::update_image_value_according_to_grayscale(unsigned row, unsigned col)
{
    if (m_greyBits[row*m_imageWidth + col] < greyDivisor)
        m_imageVector[row][col] = 0xff;
    else
        m_imageVector[row][col] = 0xffffffff;
}
void ColoredPNGProcessor::find_upper_and_lower_bounds()
{
    m_lowerBound = m_medianValue-4;
    m_upperBound = m_medianValue+4;
    if (m_lowerBound > m_medianValue)
        m_lowerBound=0x00;
    if (m_upperBound < m_medianValue)
        m_lowerBound = 0xff;
}

SimpleColoredPng::SimpleColoredPng(std::string fileName)
: PNGProcessor(fileName), mOtherPixelCount(0)
{
    for (unsigned row = 0; row < m_imageHeight; row++)
        for (unsigned col = 0; col < m_imageWidth; col++)
        {
            store_special_colors_with_coordinates(row, col);
            update_image_value_ignoring_white(row, col);
        }

    PNGProcessor::commit_image_vector_to_pixel_vector();
}

const unsigned RED   = 0xFF0000;
const unsigned GREEN = 0x00FF00;
const unsigned BLUE  = 0x0000FF;

void SimpleColoredPng::store_special_colors_with_coordinates(unsigned row, unsigned col)
{
    unsigned colorNoAlpha = (m_imageVector[row][col]>>8);
    switch(colorNoAlpha)
    {
        case RED:
            mRedPixels.push_back({row, col});
            break;
        case GREEN:
            mGreenPixels.push_back({row, col});
            break;
        case BLUE:
            mBluePixels.push_back({row, col});
            break;
        default:
            break;
    }
}

void SimpleColoredPng::update_image_value_ignoring_white(unsigned row, unsigned col)
{
    bool isWhite = m_imageVector[row][col] == 0xffffffff;
    bool isTransparent = (m_imageVector[row][col] & 0xff) == 0;
    if(!isWhite && !isTransparent)
    {
        m_imageVector[row][col] = 0xff;
        mOtherPixelCount++;
    }
}

void SimpleColoredPng::fill_id_vector_with_active_pixels(stk::mesh::EntityIdVector& elemIds) const
{
    unsigned id = 1;
    for (int rowIndex = m_imageHeight-1; rowIndex >= 0; rowIndex--)
        for (unsigned colIndex = 0; colIndex < m_imageWidth; colIndex++, id++)
            if (m_pixelVector[rowIndex][colIndex])
                elemIds.push_back(id);
}

stk::mesh::EntityIdVector SimpleColoredPng::get_elemIds_for_colored_pixels(const std::vector<Pixel> & coloredPixels)
{
    stk::mesh::EntityIdVector elemIds;
    for(Pixel pixel : coloredPixels)
    {
        stk::mesh::EntityId id = (m_imageHeight - pixel.x) * m_imageWidth - (m_imageWidth - (pixel.y+1));
        elemIds.push_back(id);
    }
    return elemIds;
}
