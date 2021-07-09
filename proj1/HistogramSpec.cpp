#include "image.h"
#include <cstring>
#include <iostream>
#include <memory>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

void PrintOriginalImageValues(int N, int M, ImageType &image){
  int i,j,val;
  ofstream myfile;
  myfile.open("OriginalImage.csv");

  for (i=0; i<N;i++)
  {
    for (j=0; j<M; j++)
    {
      image.getPixelVal(i,j,val);
      myfile << val << ",";
    }
    myfile << endl;
  }
  myfile.close();
}

void PrintHistogramValues(int N, int M, ImageType &image){
  int i,j,val;
  ofstream myfile1;
  myfile1.open("hist.csv");

  for (i=0; i<N;i++)
  {
    for (j=0; j<M; j++)
    {
      image.getPixelVal(i,j,val);
      myfile1 << val << ",";
    }
    myfile1 << endl;
  }
  myfile1.close();
}

float * histogram(ImageType image) {
  int rows, cols, bytes;
  image.getImageInfo(rows, cols, bytes);
  int total_pixels = rows * cols;
  float * hist = new float[256];
  memset(hist, 0, sizeof(float)*256);
  int dummy = 0;
  for (int i=0; i<rows; ++i) {
    for (int j=0; j<cols; ++j) {
      if(image.getPixelVal(i,j, dummy) < 256 && image.getPixelVal(i,j, dummy) >= 0)
        hist[image.getPixelVal(i,j, dummy)] += 1.0f/(float)total_pixels;
    }
  }
  return hist;
}

float * cumulativeHist(const float * hist) {
  float * cumulativehist;
  float sum = 0;
  cumulativehist = new float[256];
  memset(cumulativehist, 0, sizeof(float)*256);
  for (int i=0; i<256; ++i) {
    sum += hist[i];
    cumulativehist[i] = sum;
  }
  return cumulativehist;
}

int * inverseMap(const float * in_hist, const float * out_hist) {
  int * map = new int[256];
  float var = 0, prev_var = 256;

  for (int i=0; i<256; ++i) {
    for(int j=0; j<256; ++j) {
      var = std::abs(out_hist[j] - in_hist[i]);
      if (var > prev_var) {
        map[i] = j-1;
        break;
      } else if (j==255) {
        map[i] = 255;
      }
      prev_var = var;
    }
    prev_var = 256;
  }
  return map;
}

ImageType * mapImage(const int * map, ImageType image) {
  int rows, cols, bytes;
  image.getImageInfo(rows, cols, bytes);
  ImageType * out_image = new ImageType(rows, cols, bytes);
  int dummy = 0;
  for (int i=0; i<rows; ++i) {
    for (int j=0; j<cols; ++j) {
      out_image->setPixelVal(i, j, map[image.getPixelVal(i,j, dummy)]);
      std::cout << out_image->getPixelVal(i,j, dummy) << std::endl;
    }
  }

  return out_image;
}


ImageType* histogramSpec(ImageType& im, ImageType& imSpec) {
  float * in_hist = histogram(im);
  float * out_hist = histogram(imSpec);
  float * in_cumulative_hist = cumulativeHist(in_hist);
  float * out_cumulative_hist = cumulativeHist(out_hist);
  int * map = inverseMap(in_cumulative_hist, out_cumulative_hist);
  ImageType * mapped_image = mapImage(map, im);
  //writeImage("test.pgm", *mapped_image);
  return mapped_image;
}
