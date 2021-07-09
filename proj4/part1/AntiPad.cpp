// Steven Hernandez
// CS 474
// 10/28/2020
//This removes the pad from the image
//execute like this: ./executable Paddedimage.pgm P/2xQ/2Image.pgm output.pgm
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#include "image.h"

int readImageHeader(char[], int&, int&, int&, bool&);
int readImage(char[], ImageType&);
int writeImage(char[], ImageType&);

int main(int argc, char *argv[])
{
 int M, N, Q;
 int O, P;
 bool type;
 int val;
 int sample;
 int resolution;

 // read image header
 readImageHeader(argv[1], N, M, Q, type);

 // allocate memory for the image array

 ImageType image(N, M, Q);

 // read image
 readImage(argv[1], image);

 O = N / 2;
 P = M / 2;

 ImageType image2(O, P, Q);

 // read image
 readImage(argv[2], image2);

 int newImage[O][P];

 for(int i = 0; i < O; i++)
 {
   for(int j = 0; j < P; j++)
   {
     image.getPixelVal(i,j,val);
     newImage[i][j] = val;
   }
 }

 for(int i = 0; i < O; i++)
 {
   for(int j = 0; j < P; j++)
   {
    image2.getPixelVal(i,j,val);
    image2.setPixelVal(i,j,newImage[i][j]);
   }
 }

 // write image
 writeImage(argv[3], image2);


 return (1);
}
