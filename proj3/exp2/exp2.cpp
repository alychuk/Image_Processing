// modified exp2 for exp3 this is the code for exp2 copied over in case needed later

for (int i = 0; i < rows; ++i)
{
  //add a new col for each row
  image_real[i + 1] = new float[cols + 1];
  image_imag[i + 1] = new float[cols + 1];

  for (int j = 0; j < cols; ++j)
  {
    //gets the value of the pixel from image
    image.getPixelVal(i, j, val);

    //insert value into real num array and have 0 on imaginary
    image_real[i + 1][j + 1] = val;
    image_imag[i + 1][j + 1] = 0;

    // image_real[i+1][j+1] *= pow(-1, ((i+1)+(j+1)));
    // if multiplied by -1, it is centered
  } // if multiplied by 1, then it isn't centered
}

int N = rows;
int M = cols;
//call 2D fft and passes the size of row,col, the 2d arrays for
//real and imaginary numbers along with the isign, -1 = forward,
//1 = Inverse
fft2D(M, N, image_real, image_imag, -1);
std::complex<float> complex_num;
for (int i = 1; i < N + 1; ++i)
{
  for (int j = 1; j < N + 1; ++j)
  {
    complex_num = std::complex<float>(image_real[i][j], image_imag[i][j]);
    complex_num /= N * N;

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
  }
}

for (int i = 0; i < rows; ++i)
{
  for (int j = 0; j < cols; ++j)
  {
    image.setPixelVal(i, j, image_out[i][j]);
  }
}
writeImage(argv[2], image);
