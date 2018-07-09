#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


struct Pixel {
	unsigned char R, G, B;  // Blue, Green, Red
};

class ColorImage {
	Pixel *pPixel;
	int xRes, yRes;
public:
	ColorImage();
	~ColorImage();
	void init(int xSize, int ySize);
	void clear(Pixel background);
	Pixel readPixel(int x, int y);
	void writePixel(int x, int y, Pixel p);
	void outputPPM(char *filename);
};


/*

ColorImage::ColorImage();

ColorImage::~ColorImage();

void ColorImage::init(int xSize, int ySize);

void ColorImage::clear(Pixel background);

Pixel ColorImage::readPixel(int x, int y);

void ColorImage::writePixel(int x, int y, Pixel p);

void ColorImage::outputPPM(char *filename);
*/
