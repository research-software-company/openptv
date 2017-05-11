//
//  openptv.c
//  OpenPTV
//
//  Created by Alex Liberzon on 11/05/2017.
//
//

#include "optv/tracking_frame_buf.h"
#include <stdio.h>

/* 
    Compile it using: 
       make 
    Run it using: 
       ./openptv test_cavity 10000 10004
*/
int main(int argc, char *argv[])
{
    
    // 1. process inputs: directory, first frame, last frame
    
      int count;

  printf ("This program was called with \"%s\".\n",argv[0]);
  
  if (argc != 4) 
    {
        printf("Wrong number of inputs, expecting: ./openptv test_cavity 10000 10004 \n");
        return 0;
    }

  if (argc > 1)
    {
      for (count = 1; count < argc; count++)
	    {
	        printf("argv[%d] = %s\n", count, argv[count]);
	    }
    }
    
    // 2. init_proc
    // 3. start_proc
    // 4. read parameters
    // 5. sequence (init, set images, loop)
    // 6. tracking (init, loop, finish)
        
    return 0;
}
