#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "image.cpp"

using namespace std;

struct Point {
  int x,y;
};
void Correlation2D(ImageType &input, ImageType &output, Point mask_size,
                    float **mask, Point pad, int boundaries);
void Correlate(ImageType &input, ImageType &output, Point size, float **kernel,
                    Point pad, int boundaries);
void ImageToMask(ImageType image, float **mask);
void Normalize255(ImageType& image);
void Magnitude(ImageType& out, ImageType& part_x, ImageType& part_y);

int main(int argc, char *argv[]) {
  int rows, cols, bytes;
  bool type;
  readImageHeader(argv[2], rows, cols, bytes, type);
  ImageType image(rows, cols, bytes);
  readImage(argv[2], image);
  ImageType out_y(rows, cols, bytes);
  ImageType out_x(rows, cols, bytes);
  ImageType gradient_magnitude(rows, cols, bytes);
  ImageType out_image(rows, cols, bytes);
  readImageHeader(argv[1], rows, cols, bytes, type);
  ImageType mask_image(rows, cols, bytes);
  readImage(argv[1], mask_image);
/*  Point size;
  size.x = 3;
  size.y = 3;
  int partial_YY[3][3] =
  {
    {-1,-1,-1},
    {0,0,0},
    {1,1,1}
  };
  int partial_XX[3][3] =
  {
    {-1,0,1},
    {-1,0,1},
    {-1,0,1}
  };
  float ** partial_y = new float*[size.x];
  for (int i =0; i< size.x; ++i) {
      partial_y[i] = new float[size.y];
  }
  for (int i = 0; i<size.x; ++i) {
      for (int j =0; j < size.y; ++j) {
        partial_y[i][j] = partial_YY[i][j];
    }
  }
  float ** partial_x = new float*[size.x];
  for (int i =0; i< size.x; ++i) {
      partial_x[i] = new float[size.y];
  }
  for (int i = 0; i<size.x; ++i) {
      for (int j =0; j < size.y; ++j) {
        partial_x[i][j] = partial_XX[i][j];
    }
  }
*/
  Point size = {rows, cols};
  float **mask;
  mask = new float*[rows];
  for (int i=0; i<rows; ++i) {
    mask[i] = new float[cols];
  }

  ImageToMask(mask_image, mask);

  Point pad = {size.x/2, size.y/2};
  Correlation2D(image, out_image, size, mask, pad, 0);
  Normalize255(out_image);
  writeImage(argv[3], out_image);
  /*Correlation2D(image, out_y, size, partial_y, pad, 0);
  Correlation2D(image, out_x, size, partial_x, pad, 0);
  Magnitude(gradient_magnitude, out_x, out_y);
  Normalize255(out_y);
  Normalize255(out_x);
  Normalize255(gradient_magnitude);
  writeImage(argv[3], out_y);
  writeImage(argv[4], out_x);
  writeImage(argv[5], gradient_magnitude);*/
  return 0;
}


void Correlation2D(ImageType &input, ImageType &output, Point mask_size,
                   float **mask, Point pad, int boundaries)
{
  Correlate(input, output, mask_size, mask, pad, boundaries);
}



void Correlate(ImageType &input, ImageType &output, Point size,
            float **kernel, Point pad, int boundaries)
{
  int rows, cols, levels, dummy;
  input.getImageInfo(rows, cols, levels);
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      float sum = 0;
      Point sample;
      for (int k = 0; k < size.x; ++k)
      {
        for (int l = 0; l < size.y; ++l)
        {
          sample.x = abs(i + k - pad.x);
          sample.y = abs(j + l - pad.y);
          if (sample.x >= rows)
            sample.x = rows - 1 - sample.x % rows;
          if (sample.y >= cols)
            sample.y = cols - 1 - sample.y % cols;
          sum += static_cast<float>(input.getPixelVal(sample.x, sample.y, dummy)) * kernel[k][l];
        }
      }
      output.setPixelVal(i, j, sum);
    }
  }
}

void Normalize255(ImageType& image) {
  int rows, cols, levels, val, newval;
  int max = 0;
  int min = 0;
  image.getImageInfo(rows, cols, levels);
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      if (image.getPixelVal(i,j,val) > max) {
        max = val;
      }
      else if (image.getPixelVal(i,j,val) < min) {
        min = val;
      }
    }
  }
  int range = max - min;
  int newmin = 0;
  int newmax = 255;
  int newrange = newmax - newmin;
  cout << range;
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      newval = newmin + (((image.getPixelVal(i,j,val) - min) * newrange) / range);
      image.setPixelVal(i,j,newval);
    }
  }
}

void Magnitude(ImageType& out, ImageType& part_x, ImageType& part_y) {
  int rows, cols, levels, val, newval, part_X, part_Y;
  int max = 0;
  int min = 0;
  out.getImageInfo(rows, cols, levels);
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      part_x.getPixelVal(i,j,part_X);
      part_y.getPixelVal(i,j,part_Y);
      newval = floor(sqrt((part_X*part_X) + (part_Y*part_Y)));
      out.setPixelVal(i,j,newval);
    }
  }
}

void ImageToMask(ImageType image, float **mask)
{
  int rows, cols, levels, dummy;

  image.getImageInfo(rows, cols, levels);
  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      mask[i][j] = 1.0f / ((static_cast<float>(image.getPixelVal(i, j, dummy)) * static_cast<float>(rows * cols))) * 1000;
      if (image.getPixelVal(i, j, dummy) == 0) {
        mask[i][j] = 0;
      }
    }
  }
}
