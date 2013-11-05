#include "trafo.h"

/*  transformation detection pixel coordinates -> geometric coordinates */
void pixel_to_metric (double * x_metric
		       , double * y_metric
		       , double x_pixel
		       , double y_pixel
		       , int im_size_x
		       , int im_size_y
		       , double pix_size_x
		       , double pix_size_y
		       , int y_remap_mode){
  
  switch (y_remap_mode){
  case NO_REMAP:
    break;

  case DOUBLED_PLUS_ONE:  
    y_pixel = 2. * y_pixel + 1.;  
    break;
    
  case DOUBLED:  
    y_pixel *= 2.;  
    break;  
  }
  
  *x_metric = (x_pixel - ((double)im_size_x) / 2.) * pix_size_x;
  *y_metric = ( ((double) im_size_y)/2. - y_pixel) * pix_size_y;

}

/*  wraps previous one, parameters are read directly from control_par* structure */
void pixel_to_metric_control_par(double * x_metric
				 , double * y_metric
				 , double x_pixel
				 , double y_pixel
				 , control_par* parameters				 
				 ){
  pixel_to_metric(x_metric
		  , y_metric
		  , x_pixel
		  , y_pixel
		  , parameters->imx
		  , parameters->imy
		  , parameters->pix_x
		  , parameters->pix_y
		  , parameters->chfield);


}





/*  transformation detection geometric coordinates -> pixel coordinates */
void metric_to_pixel (double * x_pixel
		      , double * y_pixel
		      , double x_metric
		      , double y_metric
		      , int im_size_x
		      , int im_size_y
		      , double pix_size_x
		      , double pix_size_y
		      , int y_remap_mode){
  
  
  *x_pixel = ( x_metric / pix_size_x ) + ( (double) im_size_x)/2.;
  *y_pixel = ((double)im_size_y) /2. - (y_metric / pix_size_y);

  switch (y_remap_mode){
  case NO_REMAP:
    break;

  case DOUBLED_PLUS_ONE:  
    *y_pixel = (*y_pixel - 1.)/2.;  
    break;
    
  case DOUBLED:  
    *y_pixel /= 2.;
    break;  
  }
}

