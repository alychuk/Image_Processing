#include "image.h"
//quantize function
// assumes out image is the same size and type as in image
ImageType* quantize(ImageType& image, int type) {
  int rows, cols, levels, val;
  image.getImageInfo(rows, cols, levels);
  ImageType * out_image = new ImageType(rows, cols, levels);
  for (int i=0; i<rows; ++i) {
    for (int j=0; j<cols; ++j) {
      image.getPixelVal(i,j,val);
      out_image->setPixelVal(i, j, (256/type) * (val / (256/type)));
    }
  }
  return out_image;
}
