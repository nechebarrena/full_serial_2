#include "full_s_2.h"




int main(void){
  
  printf("\n.... Programa full V 2.0 ....\n");
  int t=0;
  int t_aux=0;
  double tiempo_total=0.0;
  int iteraciones=2;
  
  presion=malloc(sizeof(double)*TAM_C*TAM_F);
  presion_N=malloc(sizeof(double)*TAM_C*TAM_F);
  presion_capilar=malloc(sizeof(double)*TAM_C*TAM_F);
  
  saturacion=malloc(sizeof(double)*TAM_C*TAM_F);
  saturacion_N=malloc(sizeof(double)*TAM_C*TAM_F);
  
  Q_W=malloc(sizeof(double)*TAM_C*TAM_F);
  Q_N=malloc(sizeof(double)*TAM_C*TAM_F);
  
  K=malloc(sizeof(double)*TAM_C*TAM_F);
  
  lambda_n=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_w=malloc(sizeof(double)*TAM_C*TAM_F);
  
  krn=malloc(sizeof(double)*TAM_C*TAM_F);
  krw=malloc(sizeof(double)*TAM_C*TAM_F);
  
  lambda_t=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_d=malloc(sizeof(double)*TAM_C*TAM_F);
  
  termino_presion_capilar=malloc(sizeof(double)*TAM_C*TAM_F);
  
  lambda_d_1=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_d_2=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_d_3=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_d_4=malloc(sizeof(double)*TAM_C*TAM_F);

  lambda_t_1=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_t_2=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_t_3=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_t_4=malloc(sizeof(double)*TAM_C*TAM_F);

  presion_n=malloc(sizeof(double)*TAM_C*TAM_F);
  presion_w=malloc(sizeof(double)*TAM_C*TAM_F);
  
  
  lambda_n_1=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_n_2=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_n_3=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_n_4=malloc(sizeof(double)*TAM_C*TAM_F);

  lambda_w_1=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_w_2=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_w_3=malloc(sizeof(double)*TAM_C*TAM_F);
  lambda_w_4=malloc(sizeof(double)*TAM_C*TAM_F);
  
  U_w_1=malloc(sizeof(double)*TAM_C*TAM_F);
  U_w_2=malloc(sizeof(double)*TAM_C*TAM_F);
  U_w_3=malloc(sizeof(double)*TAM_C*TAM_F);
  U_w_4=malloc(sizeof(double)*TAM_C*TAM_F);
  
  U_n_1=malloc(sizeof(double)*TAM_C*TAM_F);
  U_n_2=malloc(sizeof(double)*TAM_C*TAM_F);
  U_n_3=malloc(sizeof(double)*TAM_C*TAM_F);
  U_n_4=malloc(sizeof(double)*TAM_C*TAM_F);
  
  U_t_1=malloc(sizeof(double)*TAM_C*TAM_F);
  U_t_2=malloc(sizeof(double)*TAM_C*TAM_F);
  U_t_3=malloc(sizeof(double)*TAM_C*TAM_F);
  U_t_4=malloc(sizeof(double)*TAM_C*TAM_F);
  
  f_w=malloc(sizeof(double)*TAM_C*TAM_F);
  
  velocidad=malloc(sizeof(double)*TAM_C*TAM_F);
  
  integral=malloc(sizeof(double));
  integral_total=malloc(sizeof(double)*iteraciones);
  
  inicializar_presion();
  inicializar_saturacion();
  inicializar_Q_N();
  inicializar_Q_W();
  inicializar_K();
  delta_tiempo();
  tiempo_total=delta_t*iteraciones;
  printf("\n \n \n Delta_T = %.15f    Tiempo Total=%.15f\n \n",delta_t,tiempo_total);
  
  for(t=1 ; t<= iteraciones ; t++){
  //printf("\n Iteracion = %i \n",t);
  calcular_lambda_w();
  calcular_lambda_n();
  calcular_lambda_t();
  calcular_lambda_d();
  
  calcular_lambda_n_1_2_3_4();
  calcular_lambda_w_1_2_3_4();
  
  calcular_lambda_d_1_2_3_4();
  calcular_lambda_t_1_2_3_4_A(); // Alternativa

  calcular_presion_capilar();
  calcular_termino_presion_capilar();
  t_aux=0;
  diferencia=0.001;
  if(t==1){
    
   while(diferencia>error && t_aux<40000){ 
    //calculo_presion();
    calculo_laplace();
   t_aux = t_aux +1;
   
   }
   
  }  
  else {  
  /*
  while(diferencia>error && t_aux<5000 ){
  //for(t=0 ; t<1 ; t++){
  //printf("\n Iteracion =%i    diferencia =%.16f \n",t,diferencia);  
  calculo_presion();
  t_aux=t_aux + 1;
   if(t_aux>60000 && t>1){
    printf("\n Iteracion =%i    diferencia =%.16f  tiempo=%f/%f  [%f %%] \n",t_aux,diferencia,t*delta_t,tiempo_total,(t*delta_t*100)/(tiempo_total));  
   }
  }
  */
  }
  
  /*
  int k=0;
  for(k=1 ; k<1000 ; k++){
  //printf("\n Iteracion =%i   Iteracion error=%i  diferencia =%.16f \n",t,k,diferencia);  
  calculo_presion();
  }
  */
  
  calcular_presion_w();
  calcular_presion_n();
  
  
  
  
  calcular_lambda_w();
  calcular_lambda_n();
  calcular_lambda_t();
  calcular_lambda_d();
  
  calcular_lambda_n_1_2_3_4();
  calcular_lambda_w_1_2_3_4();
  
  calcular_lambda_d_1_2_3_4();
  calcular_lambda_t_1_2_3_4_A(); // Alternativa
  
  
  
  
  

  calcular_U_w_1_2_3_4();
  calcular_U_n_1_2_3_4();
  calcular_U_t_1_2_3_4();
  calcular_f_w();
  
  calculo_saturacion();
  calculo_velocidad();
  
  guardar_1(lambda_t_3);
  guardar_2(lambda_t_4);
  calculo_integral();
  integral_total[t]=*integral;
  guardar_v(integral_total,iteraciones);
  printf("\n Tiempo= %f/%f [ %.3f %%]  iteraciones= %i \n",t*delta_t,tiempo_total,(t*delta_t*100)/(tiempo_total),t_aux);
  }
  
  //guardar(saturacion);
  
  
  
  
  printf("\n \n \n Delta_T = %.15f    Tiempo Total=%.15f\n \n",delta_t,tiempo_total);
  printf("\n \n \n.... FIN ....\n \n");
  
  return 0;
}
