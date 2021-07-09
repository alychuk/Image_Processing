// Steven Hernandez
// CS 474
// 10/28/2020
//This adds the pad P Q to the image
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
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

 ImageType image(N, M, Q);

 // read image
 readImage(argv[1], image);

 //Frequency butterworth equation
 //H(u,v) = 1 / (1+ (((DW)/ ((D^2)-(Do^2)))^2n))
 //W = 1, Do = 80, D= sqrt(u^2 + v^2), n = 1
 //change up W, Do and n as see fit
 float W = 15.0;
 float n = 1.0;
 float D;
 float Dzero = 71.0;

 for(int i = 0; i < N; i ++)
 {
   for(int j = 0; j < M; j++)
   {
     //There is no exact "middle" element in this 2D array so dividing into
     //4 quadrants to calculate distance
     //Upper left quadrant
     if((i <= 511) && (j <= 511))
     {
       float u = 511 - i;
       float v = 511 - j;
       //Get the distance from the current element to the middle of the image
       D = sqrt(pow(u,2) + pow(v,2));
       //DW / (D^2 - Do^2)
       float part1 = 0.0;
       part1 = (D*W) / (pow(D,2) - pow(Dzero, 2));
       //  (part1) ^ 2n
       float part2 = pow(part1,2*n);
       // 1 / (1 + part2)
       float part3 = 1 / (1 + part2);
       //part3 = 1 - part3;                  //Comment out
       part3 *= 255;
       if(val > 255)
       {
         val = 255;
       }else if(val < 0)
       {
         val = 0;
       }
       image.setPixelVal(i,j,part3);
     }


     //Upper right quadrant
     if((i >= 512) && (j <= 511))
     {
       float u = i - 512;
       float v = 511 - j;
       D = sqrt(pow(u,2) + pow(v,2));
       float part1 = 0.0;
       part1 = (D*W) / (pow(D,2) - pow(Dzero, 2));
       float part2 = pow(part1,2*n);
       float part3 = 1 / (1 + part2);
       //part3 = 1 - part3;                  //Comment out
       part3 *= 255;
       if(val > 255)
       {
         val = 255;
       }else if(val < 0)
       {
         val = 0;
       }
       image.setPixelVal(i,j,part3);
     }


     //Lower right quadrant
     if((i <= 511) && (j >= 512))
     {
       float u = 511 - i;
       float v = j - 512;
       D = sqrt(pow(u,2) + pow(v,2));
       float part1 = 0.0;
       part1 = (D*W) / (pow(D,2) - pow(Dzero, 2));
       float part2 = pow(part1,2*n);
       float part3 = 1 / (1 + part2);
       //part3 = 1 - part3;                  //Comment out
       part3 *= 255;
       if(val > 255)
       {
         val = 255;
       }else if(val < 0)
       {
         val = 0;
       }
       image.setPixelVal(i,j,part3);
     }

     //Lower left quadrant
     if((i >= 512) && (j >= 512))
     {
       float u = i - 512;
       float v = j - 512;
       D = sqrt(pow(u,2) + pow(v,2));
       float part1 = 0.0;
       part1 = (D*W) / (pow(D,2) - pow(Dzero, 2));
       float part2 = pow(part1,2*n);
       float part3 = 1 / (1 + part2);
      // part3 = 1 - part3;                  //Comment out
       part3 *= 255;
       if(val > 255)
       {
         val = 255;
       }else if(val < 0)
       {
         val = 0;
       }
       image.setPixelVal(i,j,part3);
     }
   }
 }

 // write image
 writeImage(argv[2], image);


 return (1);
}
