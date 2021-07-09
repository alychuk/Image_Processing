// Adam Lychuk
// Oct 5
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "image.cpp"

using namespace std;

int main(int argc, char *argv[]) {
  // Read Images
  int i, j;
  int M, N, Q;
  bool type;
  vector<unique_ptr<ImageType>> images;
  vector<unique_ptr<ImageType>> images_out;
  // read images
  for (int i=1; i < argc; ++i) {
    cout << argv[i] << endl;
    readImageHeader(argv[i], N, M, Q, type);
    unique_ptr<ImageType> ptr(new ImageType(N,M,Q));
    readImage(argv[i], *ptr);
    images.push_back(move(ptr));
  }

  string fname = "out";
  ImageType * out_image = quantize((*(images[1])), 128);
  writeImage( const_cast<char*>((fname + to_string(i+1) + ".pgm").c_str()) , *out_image);

  // histogramspec
  /*string fname2 = "spec_out.pgm";
  ImageType *out;
  histogramSpec(*images[0], *images[1], out);
  writeImage(const_cast<char*>(fname2.c_str()), *out);*/
  //ImageType* changed_image = histogramSpec(*images[0], *images[1]);
  //Q3PrintOriginalImageValues(N,M, *images[0]);
  //Q3PrintEqualizedHistogramValues(N,M, *changed_image);
/*
  string fname = "out";
  // write images out
  for (int i=0; i < argc-1; ++i) {
    // cast to a c string with c_str and get rid of const with const_cast
    writeImage( const_cast<char*>((fname + to_string(i+1) + ".pgm").c_str()) , *(images[i]));
  }
*/
  return 0;
}
