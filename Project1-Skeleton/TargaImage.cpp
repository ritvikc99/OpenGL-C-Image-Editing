///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <iterator>
#include <map>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel should be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    unsigned char* rgb = To_RGB();

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        unsigned char grayscale = (0.299 * rgb[j] + 0.587 * rgb[j + 1] + 0.114 * rgb[j + 2]);
        data[i] = data[i + 1] = data[i + 2] = grayscale;
    }

    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    for (int i = 0; i < width * height * 4; i += 4)
    { 
        data[i] = 32 * int(data[i] / 32.0);
        data[i + 1] = 32 * int(data[i + 1] / 32.0);
        data[i + 2] = 64 * int(data[i + 2] / 64.0);
    }

    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    map<vector<int>, int> histogram;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        int r, g, b;
        data[i] = r = 8 * int(data[i] / 8.0);
        data[i + 1] = g = 8 * int(data[i + 1] / 8.0);
        data[i + 2] = b = 8 * int(data[i + 2] / 8.0);

        vector<int> rgb{ r, g, b };
        histogram[rgb]++;
    }

    /*
    vector<pair<int, vector<int>>> sorted(histogram.size());

    transform(
        histogram.begin(), histogram.end(),
        sorted.begin(),
        [](pair<vector<int>, int> const& original) {
            return pair<int, vector<int>>(original.second, original.first);
        });
    sort(sorted.end(), sorted.begin(), greater<int>());

    for (int i = 0; i < width * height * 4; i += 4)
    {
        auto iter = sorted.begin();
        int j = 0, lowest = 0;
        vector<int> color;
        bool first = true;
        while (j < 255)
        {
            int r = (iter->second)[0], g = (iter->second)[1], b = (iter->second)[2];
            int current = sqrt((data[i] - r)^2 + (data[i + 1] - g)^2 + (data[i + 2] - b)^2);
            if (first)
            {
                lowest = current;
                color[0] = r, color[1] = g, color[2] = b;
            }
            else if (current < lowest)
            {
                lowest = current;
                color[0] = r, color[1] = g, color[2] = b;
            }
            ++iter;
            ++j;
        }
        data[i] = color[0];
        data[i + 1] = color[1];
        data[i + 2] = color[2];
    }
    */

    ClearToBlack();
    return false;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    To_Grayscale();
    float intensity;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        intensity = (float(data[i]) / 255.0);
        
        if (intensity < 0.5)
            data[i] = data[i + 1] = data[i + 2] = 0;
        else
            data[i] = data[i + 1] = data[i + 2] = 255;
    }

    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    srand(static_cast <unsigned> (time(0)));
    float intensity;
    for (int i = 0; i < width * height * 4; i += 4)
    {
        float random = (-0.2) + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 0.4));
      
        intensity = (float(data[i]) / 255.0);
        intensity += random;

        if (intensity < 0.5)
            data[i] = data[i + 1] = data[i + 2] = 0;
        else
            data[i] = data[i + 1] = data[i + 2] = 255;
    }

    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    To_Grayscale();

    float* intensities = new float[width * height];
    int current, index;
    float intensity, error, dataVal;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        intensities[i / 4] = float(data[i]) / 255.0;
    }

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            if (y % 2 == 0) // If this is an even row
                current = (y * width + x) * 4;
            else
                current = ((y * width + (width - x - 1)) * 4);

            if (intensities[current / 4] < 0.5)
            {
                error = float(intensities[current / 4] / 255.0);
                data[current] = data[current + 1] = data[current + 2] = 0;
            }
            else
            {
                error = float((intensities[current / 4] - 1.0) / 255.0);
                data[current] = data[current + 1] = data[current + 2] = 255;
            }

            // Very next pixel (x + 1, y) for even rows, (x - 1, y) for odd rows
            if (y % 2 == 0)
                index = (y * width + (x + 1));
            else
                index = (y * width + (width - x - 1));
            if (index < width * height) // If neighboring pixel is still within image borders
                intensities[index] += 255 * (float(7.0 / 16.0) * error);

            // Bottom next pixel (x + 1, y + 1) for even rows, (x - 1, y + 1) for odd rows
            if (y % 2 == 0)
                index = ((y + 1) * width + (x + 1));
            else
                index = ((y + 1) * width + (width - x - 1));
            if (index < width * height)
                intensities[index] += 255 * (float(1.0 / 16.0) * error);

            // Pixel directly underneath (x, y + 1)
            if (y % 2 == 0)
                index = ((y + 1) * width + x);
            else
                index = ((y + 1) * width + (width - x));
            if (index < width * height)
                intensities[index] += 255 * (float(5.0 / 16.0) * error);

            // Bottom previous pixel (x - 1, y + 1) for even rows, (x + 1, y + 1) for odd rows
            if (y % 2 == 0)
                index = ((y + 1) * width + (x - 1));
            else
                index = ((y + 1) * width + (width - x + 1));
            if (index < width * height)
                intensities[index] += 255 * (float(3.0 / 16.0) * error);
        }
    }

    delete[] intensities;

    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    To_Grayscale();

    unsigned int avgintensity = 0;
    unsigned char intensity;
    long int count = 0;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        avgintensity += (unsigned char)floor((data[i] / (float)255) * 255);
        count++;
    }

    avgintensity /= count;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        intensity = (unsigned char)floor((data[i] / (float)255) * 255);
        if (intensity < (255 - avgintensity))
            data[i] = data[i + 1] = data[i + 2] = 0;
        else
            data[i] = data[i + 1] = data[i + 2] = 255;
    }

    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    To_Grayscale();

    float matrix[4][4] = {
        {0.7500, 0.3750, 0.6250, 0.2500},
        {0.0625, 1.0000, 0.8750, 0.4375},
        {0.5000, 0.8125, 0.9375, 0.1250},
        {0.1875, 0.5625, 0.3125, 0.6875}
    };

    float intensity;
    int current = 0;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            current = (y * width + x) * 4;
            intensity = ((float)data[current] / 255.0);
            if (intensity < matrix[y % 4][x % 4])
                data[current] = data[current + 1] = data[current + 2] = 0;
            else
                data[current] = data[current + 1] = data[current + 2] = 255;
        }
    }

    return true;
}// Dither_Cluster

///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    ClearToBlack();
    return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    float fc, gc, af, ag;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        af = float(data[i + 3]);
        ag = float(pImage->data[i + 3]);

        for (int j = 0; j < 4; j++)
        {
            fc = float(data[i + j]) * af;
            gc = float(pImage->data[i + j]) * ag;

            data[i + j] = (unsigned char) fc + (1.0 - af) * gc;
        }
    }

    return true;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    float fc, gc, af, ag;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        af = float(pImage->data[i + 3]);
        ag = float(data[i + 3]);

        for (int j = 0; j < 4; j++)
        {
            fc = float(pImage->data[i + j]) * af;
            gc = float(data[i + j]) * ag;

            data[i + j] = (unsigned char) (ag * fc) + (0 * gc);
        }
    }

    return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }
    
    float fc, gc, af, ag;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        af = float(pImage->data[i + 3]);
        ag = float(data[i + 3]);

        for (int j = 0; j < 4; j++)
        {
            fc = float(pImage->data[i + j]) * af;
            gc = float(data[i + j]) * ag;

            data[i + j] = (unsigned char)(1 - ag) * fc + (0 * gc);
        }
    }

    return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this image and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference

///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    float filter[5][5] = {
        {0.04, 0.04, 0.04, 0.04, 0.04},
        {0.04, 0.04, 0.04, 0.04, 0.04},
        {0.04, 0.04, 0.04, 0.04, 0.04},
        {0.04, 0.04, 0.04, 0.04, 0.04},
        {0.04, 0.04, 0.04, 0.04, 0.04}
    };

    int yOffset, xOffset, dataIndex;
    float r, g, b;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            r = 0.0, g = 0.0, b = 0.0;
            for (int i = -2; i < 3; i++)
            {
                for (int j = -2; j < 3; j++)
                {
                    if (y < 2)
                        yOffset = 2;
                    else if (y > (height - 2))
                        yOffset = height - 2;
                    else
                        yOffset = y + i;

                    if (x < 2)
                        xOffset = 2;
                    else if (x > (width - 2))
                        xOffset = width - 2;
                    else
                        xOffset = x + j;

                    dataIndex = (yOffset * width + xOffset) * 4;

                    r += data[dataIndex] * filter[i + 2][j + 2];
                    g += data[dataIndex + 1] * filter[i + 2][j + 2];
                    b += data[dataIndex + 2] * filter[i + 2][j + 2];
                }
            }
            int current = (y * width + x) * 4;
            data[current] = (unsigned char) r;
            data[current + 1] = (unsigned char) g;
            data[current + 2] = (unsigned char) b;
        }
    }

    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    double filter[5][5] = {
        {0.012345679, 0.024691358, 0.037037037, 0.024691358, 0.012345679},
        {0.024691358, 0.049382716, 0.074074074, 0.049382716, 0.024691358},
        {0.037037037, 0.074074074, 0.111111111, 0.074074074, 0.037037037},
        {0.024691358, 0.049382716, 0.074074074, 0.049382716, 0.024691358},
        {0.012345679, 0.024691358, 0.037037037, 0.024691358, 0.012345679}
    };

    int yOffset, xOffset, dataIndex;
    double r, g, b;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            r = 0.0, g = 0.0, b = 0.0;
            for (int i = -2; i < 3; i++)
            {
                for (int j = -2; j < 3; j++)
                {
                    if (y < 2)
                        yOffset = 2;
                    else if (y > (height - 2))
                        yOffset = height - 2;
                    else
                        yOffset = y + i;

                    if (x < 2)
                        xOffset = 2;
                    else if (x > (width - 2))
                        xOffset = width - 2;
                    else
                        xOffset = x + j;

                    dataIndex = (yOffset * width + xOffset) * 4;

                    r += data[dataIndex] * filter[i + 2][j + 2];
                    g += data[dataIndex + 1] * filter[i + 2][j + 2];
                    b += data[dataIndex + 2] * filter[i + 2][j + 2];
                }
            }
            int current = (y * width + x) * 4;
            data[current] = (unsigned char) r;
            data[current + 1] = (unsigned char) g;
            data[current + 2] = (unsigned char) b;
        }
    }

    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    double filter[5][5] = {
        {0.00390625, 0.015625, 0.0234375, 0.015625, 0.00390625},
        {0.015625, 0.0625, 0.09375, 0.0625, 0.015625},
        {0.0234375, 0.09375, 0.140625, 0.09375, 0.0234375},
        {0.015625, 0.0625, 0.09375, 0.0625, 0.015625},
        {0.00390625, 0.015625, 0.0234375, 0.015625, 0.00390625}
    };

    int yOffset, xOffset, dataIndex;
    double r, g, b;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            r = 0.0, g = 0.0, b = 0.0;
            for (int i = -2; i < 3; i++)
            {
                for (int j = -2; j < 3; j++)
                {
                    if (y < 2)
                        yOffset = 2;
                    else if (y > (height - 2))
                        yOffset = height - 2;
                    else
                        yOffset = y + i;

                    if (x < 2)
                        xOffset = 2;
                    else if (x > (width - 2))
                        xOffset = width - 2;
                    else
                        xOffset = x + j;

                    dataIndex = (yOffset * width + xOffset) * 4;

                    r += data[dataIndex] * filter[i + 2][j + 2];
                    g += data[dataIndex + 1] * filter[i + 2][j + 2];
                    b += data[dataIndex + 2] * filter[i + 2][j + 2];
                }
            }
            int current = (y * width + x) * 4;
            data[current] = (unsigned char) r;
            data[current + 1] = (unsigned char) g;
            data[current + 2] = (unsigned char) b;
        }
    }

    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    double* oneD = new double[N];      // 1-D Gaussian row
    double** filter = new double* [N]; // 2-D Gaussian Filter

    for (int i = 0; i < N; i++)
    {
        oneD[i] = Binomial(N, i);
        filter[i] = new double[N]; // Allocating each inner array as I calc that row's Pascal triangle value xD
    }

    int sum = 0, value;

    for (int y = 0; y < N; y++)
    {
        for (int x = 0; x < N; x++)
        {
            value = oneD[y] * oneD[x];
            filter[y][x] = value;
            sum += value;
        }
    }

    for (int y = 0; y < N; y++)
    {
        for (int x = 0; x < N; x++)
        {
            filter[x][y] /= sum;
        }
    }

    int yOffset, xOffset, dataIndex;
    int lower = floor(N / 2), upper = floor(N / 2) + 1;
    double r, g, b;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            r = 0.0, g = 0.0, b = 0.0;
            for (int i = -1 * lower; i < upper; i++)
            {
                for (int j = -1 * lower; j < upper; j++)
                {
                    if (y < lower)
                        yOffset = lower;
                    else if (y > (height - lower))
                        yOffset = height - lower;
                    else
                        yOffset = y + i;

                    if (x < lower)
                        xOffset = lower;
                    else if (x > (width - lower))
                        xOffset = width - lower;
                    else
                        xOffset = x + j;

                    dataIndex = (yOffset * width + xOffset) * 4;

                    r += data[dataIndex] * filter[i + lower][j + lower];
                    g += data[dataIndex + 1] * filter[i + lower][j + lower];
                    b += data[dataIndex + 2] * filter[i + lower][j + lower];
                }
            }
            int current = (y * width + x) * 4;
            data[current] = (unsigned char)r;
            data[current + 1] = (unsigned char)g;
            data[current + 2] = (unsigned char)b;
        }
    }
    
    for (int i = 0; i < N; i++)
    {
        delete[] filter[i];
    }

    delete[] filter;
    delete[] oneD;

    return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    unsigned char* filtered = new unsigned char[width * height * 4];

    for (int i = 0; i < width * height * 4; i++)
        filtered[i] = data[i];

    TargaImage* pImage = new TargaImage(width, height, filtered);

    double filter[5][5] = {
        {0.012345679, 0.024691358, 0.037037037, 0.024691358, 0.012345679},
        {0.024691358, 0.049382716, 0.074074074, 0.049382716, 0.024691358},
        {0.037037037, 0.074074074, 0.111111111, 0.074074074, 0.037037037},
        {0.024691358, 0.049382716, 0.074074074, 0.049382716, 0.024691358},
        {0.012345679, 0.024691358, 0.037037037, 0.024691358, 0.012345679}
    };

    int yOffset, xOffset, dataIndex;
    double r, g, b;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            r = 0.0, g = 0.0, b = 0.0;
            for (int i = -2; i < 3; i++)
            {
                for (int j = -2; j < 3; j++)
                {
                    if (y < 2)
                        yOffset = 2;
                    else if (y > (height - 2))
                        yOffset = height - 2;
                    else
                        yOffset = y + i;

                    if (x < 2)
                        xOffset = 2;
                    else if (x > (width - 2))
                        xOffset = width - 2;
                    else
                        xOffset = x + j;

                    dataIndex = (yOffset * width + xOffset) * 4;

                    r += data[dataIndex] * filter[i + 2][j + 2];
                    g += data[dataIndex + 1] * filter[i + 2][j + 2];
                    b += data[dataIndex + 2] * filter[i + 2][j + 2];
                }
            }
            int current = (y * width + x) * 4;
            pImage->data[current] = (unsigned char)r;
            pImage->data[current + 1] = (unsigned char)g;
            pImage->data[current + 2] = (unsigned char)b;
        }
    }

    Difference(pImage);

    delete[] filtered;

    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    unsigned char* filtered = new unsigned char[width * height * 4];
    unsigned char* copy = new unsigned char[width * height * 4];

    for (int i = 0; i < width * height * 4; i++)
    {
        filtered[i] = data[i];
        copy[i] = data[i];
    }

    TargaImage* pImage = new TargaImage(width, height, filtered);
    TargaImage* qImage = new TargaImage(width, height, copy);

    double filter[5][5] = {
        {0.012345679, 0.024691358, 0.037037037, 0.024691358, 0.012345679},
        {0.024691358, 0.049382716, 0.074074074, 0.049382716, 0.024691358},
        {0.037037037, 0.074074074, 0.111111111, 0.074074074, 0.037037037},
        {0.024691358, 0.049382716, 0.074074074, 0.049382716, 0.024691358},
        {0.012345679, 0.024691358, 0.037037037, 0.024691358, 0.012345679}
    };

    int yOffset, xOffset, dataIndex;
    double r, g, b;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            r = 0.0, g = 0.0, b = 0.0;
            for (int i = -2; i < 3; i++)
            {
                for (int j = -2; j < 3; j++)
                {
                    if (y < 2)
                        yOffset = 2;
                    else if (y > (height - 2))
                        yOffset = height - 2;
                    else
                        yOffset = y + i;

                    if (x < 2)
                        xOffset = 2;
                    else if (x > (width - 2))
                        xOffset = width - 2;
                    else
                        xOffset = x + j;

                    dataIndex = (yOffset * width + xOffset) * 4;

                    r += data[dataIndex] * filter[i + 2][j + 2];
                    g += data[dataIndex + 1] * filter[i + 2][j + 2];
                    b += data[dataIndex + 2] * filter[i + 2][j + 2];
                }
            }
            int current = (y * width + x) * 4;
            pImage->data[current] = (unsigned char)r;
            pImage->data[current + 1] = (unsigned char)g;
            pImage->data[current + 2] = (unsigned char)b;
        }
    }
    
    qImage -> Difference(pImage);

    for (int i = 0; i < width * height * 4; i += 4)
    {
        data[i] += qImage->data[i];
        data[i + 1] += qImage->data[i + 1];
        data[i + 2] += qImage->data[i + 2];
    }

    delete[] filtered;
    delete[] copy;

    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}


///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

