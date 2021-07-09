//Steven Hernandez and Adam Lychuk
//CS 474
//Project 3 part 1a,1b and 1c
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <complex>
using namespace std;

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

int main(int argc, char *argv[])
{
  /////////////////////////////   part a   /////////////////////////////////////

  //signal of 2,3,4,4
  float data[] = {0, 2, 0, 3, 0, 4, 0, 4, 0};

  //multiply real numbers by -1 to center for frequency domain
  for (int i = 0; i < 4; ++i)
  {
    if ((i + 1) % 2)
    {
      data[2 * i + 1] *= -1;
    }
  }

  //calls fft function, passes array, power 4 and signal for forward FT
  fft(data, 4, -1);

  //Complex float variable
  std::complex<float> complex_num;

  //i increases by 2 to only go over real numbers in array
  for (int i = 1; i < 9; i += 2)
  {
    //stores the absolute value of the complex value a+bi
    complex_num = std::complex<float>(data[i], data[i + 1]);
    //Normalize
    complex_num /= 4;
    //stores real value in odd elements of array
    data[i] = complex_num.real();
    //stores imaginary value in even elements of array
    data[i + 1] = complex_num.imag();
  }

  std::ofstream fout;
  fout.open("function_f.dat");
  for (int i = 0; i < 4; ++i)
  {
    //Magnitude from slides is sqrt((R^2(u))+(I^2(u)))
    //store real, imaginary and magnitude
    fout << data[2 * i + 1] << "," << data[2 * i + 2] << "," << sqrt(pow(data[2 * i + 1], 2) + pow(data[2 * i + 2], 2)) << "\n";
  }
  fout.close();
  fout.clear();

  //Compute Inverse FT with the current data
  fft(data, 4, 1);
  for (int i = 0; i < 4; ++i)
  {
    if ((i + 1) % 2)
    {
      data[2 * i + 1] *= -1;
    }
  }

  /////////////////////////////////// end part a /////////////////////////////////

  ////////////////////////////////   part b   //////////////////////////////////
  std::cout << " \n PART B: COSINE FUNCTION" << std::endl;
  //float array of size 257 since the first element is not used
  float f_cos[257];

  //Store the function result in odd elements of array and 0 in even(imaginary)
  for (int i = 0; i < 128; ++i)
  {
    f_cos[2 * i + 1] = cos(2.0 * 3.14159265359 * 8.0 * (float)(i + 1) / 128.0);
    f_cos[2 * i + 2] = 0;
  }
  fout.open("function_cos_128.dat");

  //store values
  for (int i = 0; i < 128; ++i)
  {
    fout << f_cos[2 * i + 1];
  }
  fout.close();
  fout.clear();
  std::cout << std::endl;

  //multiply real numbers by -1 to center frequency domain
  for (int i = 0; i < 128; ++i)
  {
    if ((i + 1) % 2)
    {
      f_cos[2 * i + 1] *= -1;
    }
  }

  //Forward FT of the cosine results
  fft(f_cos, 7, -1); //had 128 in the second parameter

  //Goes through the array for odd elements
  for (int i = 1; i < 257; i += 2)
  {
    //stores the absolute value of the complex value a+bi
    complex_num = std::complex<float>(f_cos[i], f_cos[i + 1]);
    //normalize by 128
    complex_num /= 128.0;
    //store real num
    f_cos[i] = complex_num.real();
    //store imag num
    f_cos[i + 1] = complex_num.imag();
  }

  //Experiment does not require  inverse FT
  fout.open("function_cos_128_fft.dat");

  //store values
  for (int i = 0; i < 128; ++i)
  {
    fout << f_cos[2 * i + 1] << endl;
  }
  fout << endl;
  for (int i = 0; i < 128; ++i)
  {
    fout << f_cos[2 * i + 2] << endl;
  }
  fout << endl;
  for (int i = 0; i < 128; ++i)
  {
    fout << sqrt(pow(f_cos[2 * i + 1], 2) + pow(f_cos[2 * i + 2], 2)) << endl;
  }
  fout << endl;
  for (int i = 0; i < 128; ++i)
  {
    fout << atan2(f_cos[2 * i + 2], f_cos[2 * i + 1]) << endl;
  }
  fout.close();
  fout.clear();
  std::cout << std::endl;

  ///////////////////////////////// end part b /////////////////////////////////

  ////////////////////////////////   part c    /////////////////////////////////
  std::ifstream fin;
  float rect[257];
  //Uses data from the given file
  fin.open("Rect_128.dat");
  //Stores data in the odd elements of the array (real numbers)
  if (fin.is_open())
  {
    for (int i = 0; i < 128; ++i)
    {
      fin >> rect[2 * i + 1];
      rect[2 * i + 2] = 0;
    }
    fin.close();
  }
  else
  {
    std::cout << "error" << std::endl;
  }

  //Centers values by using -1
  for (int i = 0; i < 128; ++i)
  {
    if ((i + 1) % 2)
      rect[2 * i + 1] *= -1;
  }

  //Applies Forward FT
  fft(rect, 7, 1);

  //Goes over odd elements of array
  for (int i = 1; i < 257; i += 2)
  {
    //stores the absolute value of the complex value a+bi
    complex_num = std::complex<float>(rect[i], rect[i + 1]);
    //normalize by 128 (num of values in file),don't normalize if doing inverse
    complex_num /= 128.0;
    //store real num
    rect[i] = complex_num.real();
    //store imag num
    rect[i + 1] = complex_num.imag();
  }

  std::cout << std::endl;

  fout.open("Rect_128_results.dat");
  //store values
  for (int i = 0; i < 128; ++i)
  {
    fout << rect[2 * i + 1] << endl;
  }
  fout << endl;
  for (int i = 0; i < 128; ++i)
  {
    fout << rect[2 * i + 2] << endl;
  }
  fout << endl;
  for (int i = 0; i < 128; ++i)
  {
    fout << sqrt(pow(rect[2 * i + 1], 2) + pow(rect[2 * i + 2], 2)) << endl;
  }
  fout << endl;
  for (int i = 0; i < 128; ++i)
  {
    fout << atan2(rect[2 * i + 2], rect[2 * i + 1]) << endl;
  }
  fout.close();
  fout.clear();

  /////////////////////////////// end part c /////////////////////////////////////
  return 0;
}
