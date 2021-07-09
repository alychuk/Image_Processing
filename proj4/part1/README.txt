There are multiple programs:

//Padding.cpp: takes the original image and outputs the image with padding PxQ
//To execute: ./exe boy_noisy.pgm paddedImage.pgm

//FButterWFilter.cpp: Reads in an image that is PxQ and then computes the equation to make the
//filter on the image
//To execute: ./exe PXQimage.pgm filter.pgm

//Mix.cpp: Takes in the padded image and the mask, then outputs the padded image with the mask
//uncomment line 190 to extract the frequency.
//To execute: ./exe paddedImage.pgm filter.pgm outputImage.pgm

//AntiPad.cpp: takes in the padded image and an image that is NxM, the filtered image
// is outputted to the NxM image
//To execute: ./exe outputImage.pgm NxMImage.pgm FinalImage.pgm

//You will need the boy_noisy.pgm, an image that is NxM for the final filtered image
// and a PxQ image to make the filter

//To compile:
//g++ -Wall -std=c++11 -o nameofexecutable oneOfTheAboveFiles.cpp image.cpp ReadImage.cpp ReadImageHeader.cpp WriteImage.cpp
//To execute: read the above comments



//The exp2.cpp was from the previous project but it is used to help show the frequncy image
//of the padded image. That way I can visually see if the frequency was affected by the
//chosen filter 
