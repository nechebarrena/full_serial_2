#include "full_s_2.h"




int main(void){
  
  printf("\n.... Programa full V 2.0 ....\n");
  int t=0;
  int t_aux=0;
  double tiempo_total=0.0;
  

  alocar_variables();
  inicializar_variables();
  delta_tiempo();
  tiempo_total=delta_t*iteraciones;
    
  t=1;
  //while(saturacion[TAM_C*55 + 55]<0.05){
  for(t=1 ; t<= iteraciones ; t++){ // Inicio del FOR
  //t=t+1;
  calcular_lambdas();
  t_aux=0;
  diferencia=0.001;
  
  if(t==1){  // Inicio del IF ...... Solo para la primera vez que calculo las presiones
   while(diferencia>error && t_aux<400000){ 
    calculo_presion();
    t_aux = t_aux +1;
   }
  }// Fin del IF
  
  
  if(t>1){  // Inicio del IF ....... Para todas las veces restantes
   while(diferencia>error_1 && t_aux<1000000){ 
    calculo_presion();
    t_aux = t_aux +1;
   }
  }  // Fin del IF
  
  calcular_presion_w();
  calcular_presion_n();
  calcular_lambdas();
  calcular_velocidades();
  calculo_velocidad();
  calcular_f_w();
  calculo_saturacion();
  
  printf("\n Tiempo= %f/%f [ %.3f %%]  iteraciones= %i \n",t*delta_t,tiempo_total,(t*delta_t*100)/(tiempo_total),t_aux);
  }  // Fin del FOR  
  
  
  
  printf("\n \n \n Delta_T = %.15f    Tiempo Total=%.15f\n \n",delta_t,tiempo_total);
  printf("\n \n \n.... FIN ....\n \n");
  
  guardar_1(saturacion);
  guardar_2(presion);
  
  return 0;
}
