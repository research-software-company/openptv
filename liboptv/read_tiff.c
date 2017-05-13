//
//  read_tiff.c
//  OpenPTV
//
//  Created by Alex Liberzon on 13/05/2017.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "tiffio.h"


//using namespace std;


void printArray(uint16 * array, uint16 width);

// int main()
// {
//     
//     
//     TIFF* tif = TIFFOpen("cam1.tif", "r");
//     if (tif) {
//         uint32 imagelength,height;
//         tdata_t buf;
//         uint32 row;
//         uint32 config;
//         
//         TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
//         TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
//         TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
//         buf = _TIFFmalloc(TIFFScanlineSize(tif));
//         
//         
//         uint16 s, nsamples;
//         uint16* data;
//         TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
//         for (s = 0; s < nsamples; s++)
//         {
//             for (row = 0; row < imagelength; row++)
//             {
//                 TIFFReadScanline(tif, buf, row, s);
//                 data=(uint16*)buf;
//                 printArray(data,imagelength);
//             }
//             // printArray(data,imagelength,height);
//         }
//         
//         
//         _TIFFfree(buf);
//         TIFFClose(tif);
//     }
//     exit(0);
// }

// int main()
// {
//     TIFF* tif = TIFFOpen("cam1.tif", "r");
//     if (tif) {
// 	uint32 imageWidth, imageLength;
// 	uint32 tileWidth, tileLength;
// 	uint32 x, y, z;
// 	tdata_t buf;
// 
// 	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
// 	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
// 	TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
// 	TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength);
// 	buf = _TIFFmalloc(TIFFTileSize(tif));
// 	for (y = 0; y < imageLength; y += tileLength)
// 	    for (x = 0; x < imageWidth; x += tileWidth)
// 		TIFFReadTile(tif, buf, x, y, 0, 0);
// 	_TIFFfree(buf);
// 	TIFFClose(tif);
//     }
// }


void printArray(uint16 * array, uint16 width)
{
    uint32 i;
    for (i=0;i<width;i++)
    {
        printf("%u ", array[i]);
    }
    printf("\n");
    
    
}

void printArray8(unsigned char * array, uint32 width)
{
    uint32 i;
    for (i=0;i<width;i++)
    {
        printf("%u ", array[i]);
    }
    printf("\n");
    
    
}

int main(void)
{
    TIFF* tif = TIFFOpen("32pix.tif", "r");
    if (tif) {
        uint32 imagelength;
        tsize_t scanline;
        tdata_t buf;
        uint32 row;
        uint32 col;

        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
        scanline = TIFFScanlineSize(tif);
        buf = _TIFFmalloc(scanline);
        for (row = 0; row < imagelength; row++)
        {
            TIFFReadScanline(tif, buf, row, 1);
            for (col = 0; col < scanline; col++)
                // printf("%u", buf[col]);
                printArray8(buf,imagelength);

            printf("\n");
        }
        _TIFFfree(buf);
        TIFFClose(tif);
    }
}

