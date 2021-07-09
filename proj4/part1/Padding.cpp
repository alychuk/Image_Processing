// Steven Hernandez
// CS 474
// 10/28/2020
//This add the pad in the image PxQ
//When executing: ./exe paddedimage.pgm ReceivingImage.pgm
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
 bool type;
 int val;
 int sample;
 int resolution;

 // read image header
 readImageHeader(argv[1], N, M, Q, type);

 // allocate memory for the image array

 ImageType image(2*N, 2*M, Q);

 // read image
 readImage(argv[1], image);

 int newImage[2*N][2*M];

 for(int i = 0; i < 2*M; i++)
 {
   for(int j = 0; j < 2*N; j++)
   {
     image.getPixelVal(i,j,val);
     if(val < 1)
     {
       newImage[i][j] = 0;
     }
     else
     {
       newImage[i][j] = val;
     }
   }
 }


 for(int i = 0; i < 2*M; i++)
 {
   for(int j = 0; j < 2*N; j++)
   {
    image.getPixelVal(i,j,val);
    image.setPixelVal(i,j,newImage[i][j]);
   }
 }

 // write image
 writeImage(argv[2], image);


 return (1);
}
