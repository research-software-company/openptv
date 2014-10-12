
#include "diagnostic.h"


void DIAG_PRINT_F(Exterior) (Exterior * Ex){
  //void diag_print_Exterior (Exterior * Ex){

  printf("Exterior:\n\tx0 = %f y0 = %f z0 = %f\n",Ex->x0,Ex->y0,Ex->z0);
  printf("\tomega = %f phi = %f kappa = %f\n",Ex->omega,Ex->phi,Ex->kappa);
  printf("\tdm = \n");
  int temp_idx_i,temp_idx_j;
  for (temp_idx_i = 0 ; temp_idx_i < 3 ;temp_idx_i++){
    printf("\t\t");
    for(temp_idx_j = 0 ; temp_idx_j < 3 ; temp_idx_j ++){
      printf("%f ", Ex->dm[temp_idx_i][temp_idx_j]);
    }
    printf("\n");
  }

}


void DIAG_PRINT_F(Interior) (Interior * I){

  printf("Interior:\n\txh = %f yh = %f cc = %f\n",I->xh,I->yh,I->cc);

}

void DIAG_PRINT_F(Glass) (Glass * G){

  printf("Glass:\n\tvec_x = %f vec_y = %f vec_z = %f\n",G->vec_x,G->vec_y,G->vec_z);

}

void DIAG_PRINT_F(mm_np) (mm_np * mm){
  
  printf("mm_np:\n\tnlay = %d n1 = %f\n",mm->nlay,mm->n1);
  printf("\tn2 = %f %f %f\n",mm->n2[0],mm->n2[1],mm->n2[2]);
  printf("\td = %f %f %f\n",mm->d[0],mm->d[1],mm->d[2]);
  printf("\tn3 = %f lut = %d\n",mm->n3,mm->lut);


}

