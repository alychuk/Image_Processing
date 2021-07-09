ImageType* padZeros(ImageType& im, int rpad, int cpad) {
  int new_rows, new_cols, rows, cols, bytes, val;
  im.getImageInfo(rows, cols, bytes);
  new_rows = (2 * rpad) + rows;
  new_cols = (2 * cpad) + cols;
  ImageType* zero_padded = new ImageType(rows, cols, bytes);
  for (int i=0; i < new_rows; ++i) {
    // if row is all zeros
    if ( (i < rpad) || (i > rpad + rows) ) {
      for (int j=0; j < new_cols; ++j) {
        zero_padded->setPixelVal(i,j,0);
      }
    // otherwise only pad for row start and end
    } else {
      for (int j=0; j < new_cols; ++j) {
        if ( (j < cpad) || (j > cpad + cols) ) {
          zero_padded->setPixelVal(i,j,0);
        } else {
          // subtract pad to re-align to original image
          im.getPixelVal(i-pad,j-pad,val);
          zero_padded->setPixelVal(i,j,val);
        }
      }
    }
  }
}

ImageType* correlate(ImageType& mask, ImageType& im) {
  int rows, cols, mrows, mcols, bytes, val, mval, sum;
  mask.getImageInfo(mrows, mcols, bytes);
  im->getImageInfo(rows, cols, bytes);
  ImageType* zero_padded = padZeros(im, floor(rows / 2), floor(cols / 2));
  zero_padded->getImageInfo(rows, cols, bytes);
  ImageType* out_image = new ImageType(rows, cols, bytes);
  // iterate through image
  for (int i=floor(mrows/2); i < (rows-floor(mrows/2)); ++i) {
    for (int j=floor(mcols/2); j < (rows-floor(mcols/2)); ++j) {
      sum = 0;
      //iterate through mask
      int m,n;
      for (int k=0; k < mrows; ++k) {
        m = -floor(mrows / 2);
        for (int l=0; l < mcols; ++l) {
          n = -floor(mcols / 2);
          zero_padded->getPixelVal(i+m,j+n,val);
          mask.getPixelVal(k,l,mval);
          sum += val * mval;
          n++;
        }
      }
      out_image->setPixelVal(i,j,sum);
    }
  }
}
