//Steven Hernandez and Adam Lychuk
//CS 474
//Project 4 experiment 3
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <complex>
using namespace std;

#include "image.h"

int readImageHeader(char[], int &, int &, int &, bool &);
int readImage(char[], ImageType &);
int writeImage(char[], ImageType &);


void fft2D(int N, int M, float **real_Fuv, float **imag_Fuv, int isign);
//This from the function the prof gave us
#define SWAP(a, b) \
  tempr = (a);     \
  (a) = (b);       \
  (b) = tempr
/*

   The real part is stored in the odd locations of the array
   (data[1], data[3], data[5]. etc) and the imaginary part
   in the even locations (data[2], data[4], data[6], etc.

   The elements in the array data must be stored from data[1]
   to data[2*nn] -  data[0] is not used!

   nn is the length of the input which should be power of 2.
   Warning: the program does not check this condition.

   isign: -1 Forward FFT, isign: 1  Inverse FFT (based on our definition)

   Warning: the FFT routine provided does not multiply by the normalization
   factor 1/N that appears in the forward DFT equation; you should do this
   yourself (see page 506 from the fft handout).

*/
void fft(float data[], unsigned long nn, int isign)
{
  unsigned long n, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta;
  float tempr, tempi;

  n = nn << 1;
  j = 1;
  for (i = 1; i < n; i += 2)
  {
    if (j > i)
    {
      SWAP(data[j], data[i]);
      SWAP(data[j + 1], data[i + 1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m)
    {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax)
  {
    istep = mmax << 1;
    theta = isign * (6.28318530717959 / mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m += 2)
    {
      for (i = m; i <= n; i += istep)
      {
        j = i + mmax;
        tempr = wr * data[j] - wi * data[j + 1];
        tempi = wr * data[j + 1] + wi * data[j];
        data[j] = data[i] - tempr;
        data[j + 1] = data[i + 1] - tempi;
        data[i] += tempr;
        data[i + 1] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}
#undef SWAP

void MapValues(float **image, int w, int h)
 {
   float min = 255, max = 0;
   for (int i = 0; i < h; ++i)
   {
     for (int j = 0; j < w; ++j)
     {
       float pixelval = image[i][j];
       if (min > pixelval)
       {
         min = pixelval;
       }
       if (max < pixelval)
       {
         max = pixelval;
       };
     }
   }
   float scale = 255.0f / (max - min);
   float shift = -scale * min;
   for (int i = 0; i < h; ++i)
   {
     for (int j = 0; j < w; ++j)
     {
       image[i][j] = scale * image[i][j] + shift;
     }
   }
 }

float** homomorphic_filter(float yh, float yl, float c, float d0, int rows,  int cols)
{
  float **filter = new float *[rows + 1];
  for (int i = 0; i < rows; ++i)
  {
    //add a new col for each row
    filter[i + 1] = new float[cols + 1];

    for (int j = 0; j < cols; ++j)
    {
      //insert value into real num array and have 0 on imaginary
      //take the natural log for homomorphic filtering
      float part0, part1, part2, part3, part4;
      part0 = float(pow(i-rows/2,2))+float(pow(j-cols/2,2));
      part1 = -c*(part0 / pow(d0,2));
      cout << part1 << endl;
      part2 = exp(part1);
      cout << part2 << endl;
      part3 = 1-part2;
      part3 *= (yh-yl);
      //part3 *= pow(-1, ((i+1)+(j+1)));
      filter[i + 1][j + 1] = part3 + yl;
      //filter[i + 1][j + 1] = (yh-yl)*(1-exp(-c*( (float(pow(i,2))+float(pow(j,2))) / pow(d0,2) ) ));
      //cout << filter[i+1][j+1] << endl;
      //filter[i+1][j+1] *= pow(-1, ((i+1)+(j+1)));
    }
  }
  return filter;
}

int main(int argc, char *argv[])
{
  int rows, cols, bytes;
  bool type;
  int val;
  readImageHeader(argv[1], rows, cols, bytes, type);
  ImageType image(rows, cols, bytes);
  readImage(argv[1], image);
  //2 double arrays using pointers to store real and imaginary numbers
  float **image_real = new float *[rows + 1];
  float **image_imag = new float *[rows + 1];
  float **image_imag_filter = new float *[rows + 1];
  //initialize the first element of the arrays
  image_real[0] = new float[cols + 1];
  image_imag[0] = new float[cols + 1];
  image_imag_filter[0] = new float[cols + 1];

  for (int i = 0; i < rows; ++i)
  {
    //add a new col for each row
    image_real[i + 1] = new float[cols + 1];
    image_imag[i + 1] = new float[cols + 1];
    image_imag_filter[i+1] = new float[cols+1];

    for (int j = 0; j < cols; ++j)
    {
      //gets the value of the pixel from image
      image.getPixelVal(i, j, val);

      //insert value into real num array and have 0 on imaginary
      //take the natural log for homomorphic filtering
      image_real[i + 1][j + 1] = log(val);
      image_imag[i + 1][j + 1] = 0;
      image_imag_filter[i + 1][j + 1] = 0;
      image_real[i + 1][j + 1] *= pow(-1, ((i+1)+(j+1)));
      // if multiplied by -1, it is centered
    } // if multiplied by 1, then it isn't centered
  }

  int N = rows;
  int M = cols;
  //call 2D fft and passes the size of row,col, the 2d arrays for
  //real and imaginary numbers along with the isign, -1 = forward,
  //1 = Inverse
  fft2D(M, N, image_real, image_imag, -1);

  float yh = 1.5;
  float yl = .5;
  float c = 1;
  float d0 = 150;
  float** filter = homomorphic_filter(yh, yl, c, d0, rows, cols);
  fft2D(M,N, filter, image_imag_filter, 1);
  for (int i = 1; i < N + 1; ++i)
  {
    for (int j = 1; j < N + 1; ++j)
    {
      filter[i + 1][j + 1] *= pow(-1, ((i+1)+(j+1)));
    }
  }
  fft2D(M,N, filter, image_imag_filter, -1);


  std::complex<float> complex_num;
  std::complex<float> complex_filter;
  for (int i = 1; i < N + 1; ++i)
  {
    for (int j = 1; j < N + 1; ++j)
    {
      complex_num = std::complex<float>(image_real[i][j], image_imag[i][j]);
      complex_filter = std::complex<float>(filter[i][j], image_imag_filter[i][j]);
      complex_num /= N * N;
      complex_filter /= N*N;

      complex_num = complex_num * complex_filter;
      image_real[i][j] = complex_num.real();
      image_imag[i][j] = complex_num.imag();
      //image_real[i][j] = complex_num.real() * filter[i][j];
      //image_imag[i][j] = complex_num.imag() * filter[i][j];
    }
  }

  fft2D(M,N, image_real, image_imag, 1);

  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      complex_num = std::complex<float>(image_real[i][j], image_imag[i][j]);
      complex_num = std::exp(complex_num);
      image_real[i][j] = complex_num.real();
      image_imag[i][j] = complex_num.imag();
    }
  }

  float **image_out = new float *[rows];
  for (int i = 0; i < rows; ++i)
  {
    image_out[i] = new float[cols];
    for (int j = 0; j < cols; ++j)
    {
      image_out[i][j] = log(1. + sqrt(pow(image_real[i + 1][j + 1], 2.) + pow(image_imag[i + 1][j + 1], 2.)));
      //image_out[i][j] *= pow(-1, ((i+1)+(j+1)));
    }
  }
  //MapValues(image_out,N,M);

  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      image.setPixelVal(i, j, image_out[i][j]);
    }
  }
  writeImage(argv[2], image);
  return 0;
}

void fft2D(int N, int M, float **real_fuv, float **imag_fuv, int isign)
{
  float *buffer = new float[2 * N + 1];
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      buffer[2 * j + 1] = real_fuv[i + 1][j + 1];
      buffer[2 * j + 2] = imag_fuv[i + 1][j + 1];
    }
    fft(buffer, N, isign);
    for (int j = 0; j < N; ++j)
    {
      real_fuv[i + 1][j + 1] = buffer[2 * j + 1];
      imag_fuv[i + 1][j + 1] = buffer[2 * j + 2];
    }
  }
  for (int j = 0; j < N; ++j)
  {
    for (int i = 0; i < N; ++i)
    {
      buffer[2 * i + 1] = real_fuv[i + 1][j + 1];
      buffer[2 * i + 2] = imag_fuv[i + 1][j + 1];
    }
    fft(buffer, N, isign);
    for (int i = 0; i < N; ++i)
    {
      real_fuv[i + 1][j + 1] = buffer[2 * i + 1];
      imag_fuv[i + 1][j + 1] = buffer[2 * i + 2];
    }
  }
}
