#include "full_s_2.h"




int main(void){
  
  printf("\n.... Programa full V 2.0 ....\n");
  
  presion=malloc(sizeof(double)*TAM_C*TAM_F);
  saturacion=malloc(sizeof(double)*TAM_C*TAM_F);
  Q_W=malloc(sizeof(double)*TAM_C*TAM_F);
  Q_N=malloc(sizeof(double)*TAM_C*TAM_F);
  K=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_n=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_w=malloc(sizeof(double)*TAM_C*TAM_F);
  krn=malloc(sizeof(double)*TAM_C*TAM_F);
  krw=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_t=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_d=malloc(sizeof(double)*TAM_C*TAM_F);
  
  lambda_d_1=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_d_2=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_d_3=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_d_4=malloc(sizeof(double)*TAM_C*TAM_F);
  
  inicializar_presion();
  inicializar_saturacion();
  inicializar_Q_N();
  inicializar_Q_W();
  inicializar_K();
  calcular_lambda_w();
  calcular_lambda_n();
  calcular_lambda_t();
  calcular_lambda_d();
  calcular_lambda_d_1_2_3_4();
  
  guardar(lambda_d);
  
  calculo_presion();
  
  
  
  printf("\n \n \n.... FIN ....\n \n");
  
  return 0;
}
