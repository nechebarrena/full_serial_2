#include "full_s_2.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculo_saturacion(){
  
  recinto_interno_s();
  borde_izquierdo_s();
  borde_derecho_s();
  borde_superior_s();
  borde_inferior_s();
  esquina_superior_izquierda_s();
  esquina_superior_derecha_s();
  esquina_inferior_izquierda_s();
  esquina_inferior_derecha_s();
  copiar_s();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void recinto_interno_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
     
      if(U_t_1[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_1=f_w[TAM_C*f + (c-1)] ;
      }
      else if(U_t_1[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_1=f_w[TAM_C*f + (c)] ;
      }
      
      
      
      if(U_t_3[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_3=f_w[TAM_C*(f-1) + (c)] ;
      }
      else if(U_t_3[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_3=f_w[TAM_C*f + (c)] ;
      }

      
      
      if(U_t_2[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_2=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_2[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_2=f_w[TAM_C*f + (c+1)] ;
      }
      
      
      if(U_t_4[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_4=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_4[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_4=f_w[TAM_C*(f+1) + (c)] ;
      }
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + Q_W[TAM_C*f + (c)]);       
      
    }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void borde_izquierdo_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double f_1=0.0;
  double f_2=0.0;
  double f_3=0.0;
  double f_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    c=0;
     
     f_w_1=0.0;
      
      
      
      if(U_t_3[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_3=f_w[TAM_C*(f-1) + (c)] ;
      }
      else if(U_t_3[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_3=f_w[TAM_C*f + (c)] ;
      }

      
      
      if(U_t_2[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_2=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_2[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_2=f_w[TAM_C*f + (c+1)] ;
      }
      
      
      if(U_t_4[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_4=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_4[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_4=f_w[TAM_C*(f+1) + (c)] ;
      }
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + Q_W[TAM_C*f + (c)]);       
      //printf("\n  Saturacion = %f FILA=%i     COLUMNA=%i\n",saturacion_N[TAM_C*f + (c)],f,c);
      
      
      
      
      if(c==0 && f==1){
      f_1= f_w_1*U_t_1[TAM_C*f + (c)] + f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_2= f_w_2*U_t_2[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_3= f_w_3*U_t_3[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_4= f_w_4*U_t_4[TAM_C*f + (c)] + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      //printf("\n FILA=%i   COLUMNA=%i  flujo 1 =%f   flujo 2 =%f  flujo 3=%f  flujo 4=%f  Q_W =%f \n",f,c,f_1,f_2,f_3,f_4,Q_W[TAM_C*f + (c)]);
      }
      
      
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void borde_derecho_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C-1;
     
      if(U_t_1[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_1=f_w[TAM_C*f + (c-1)] ;
      }
      else if(U_t_1[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_1=f_w[TAM_C*f + (c)] ;
      }
      
      
      
      if(U_t_3[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_3=f_w[TAM_C*(f-1) + (c)] ;
      }
      else if(U_t_3[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_3=f_w[TAM_C*f + (c)] ;
      }

      
      
     f_w_2=0.0;
      
      
      if(U_t_4[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_4=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_4[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_4=f_w[TAM_C*(f+1) + (c)] ;
      }
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + Q_W[TAM_C*f + (c)]);       
      
    
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void borde_superior_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double f_1=0.0;
  double f_2=0.0;
  double f_3=0.0;
  double f_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
  
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=0;
      
      if(U_t_1[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_1=f_w[TAM_C*f + (c-1)] ;
      }
      else if(U_t_1[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_1=f_w[TAM_C*f + (c)] ;
      }
      
      
      
      f_w_3=0.0;

      
      
      if(U_t_2[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_2=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_2[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_2=f_w[TAM_C*f + (c+1)] ;
      }
      
      
      if(U_t_4[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_4=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_4[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_4=f_w[TAM_C*(f+1) + (c)] ;
      }
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + Q_W[TAM_C*f + (c)]);       
      
      if(c==1 && f==0){
      f_1= f_w_1*U_t_1[TAM_C*f + (c)] + f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_2= f_w_2*U_t_2[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_3= f_w_3*U_t_3[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_4= f_w_4*U_t_4[TAM_C*f + (c)] + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      //printf("\n FILA=%i   COLUMNA=%i  flujo 1 =%f   flujo 2 =%f  flujo 3=%f  flujo 4=%f  Q_W =%f \n",f,c,f_1,f_2,f_3,f_4,Q_W[TAM_C*f + (c)]);
      }
    
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void borde_inferior_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
  
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=TAM_F-1;
     
      if(U_t_1[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_1=f_w[TAM_C*f + (c-1)] ;
      }
      else if(U_t_1[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_1=f_w[TAM_C*f + (c)] ;
      }
      
      
      
      if(U_t_3[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_3=f_w[TAM_C*(f-1) + (c)] ;
      }
      else if(U_t_3[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_3=f_w[TAM_C*f + (c)] ;
      }

      
      
      if(U_t_2[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_2=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_2[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_2=f_w[TAM_C*f + (c+1)] ;
      }
      
      
     f_w_4=0.0;
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + Q_W[TAM_C*f + (c)]);       
      
    
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void esquina_superior_izquierda_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double f_1=0.0;
  double f_2=0.0;
  double f_3=0.0;
  double f_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
      c=0;
      f=0;
     
     f_w_1=0.0;
     
     f_w_3=0.0;

      
      
      if(U_t_2[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_2=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_2[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_2=f_w[TAM_C*f + (c+1)] ;
      }
      
      
      if(U_t_4[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_4=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_4[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_4=f_w[TAM_C*(f+1) + (c)] ;
      }
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + 4.0*Q_W[TAM_C*f + (c)]);  
      
      f_1= f_w_1*U_t_1[TAM_C*f + (c)] + f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_2= f_w_2*U_t_2[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_3= f_w_3*U_t_3[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      f_4= f_w_4*U_t_4[TAM_C*f + (c)] + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;
      
      /*
      printf("\n FILA=%i   COLUMNA=%i  flujo 1 =%.10f   flujo 2 =%.10f  flujo 3=%.10f  flujo 4=%.10f  Q_W =%.10f \n",f,c,f_1,f_2,f_3,f_4,4.0*Q_W[TAM_C*f + (c)]);
      printf("\n FILA=%i   COLUMNA=%i  flujo total = %.10f  epsilon=%.10f   gama=%.10f \n",f,c,1.0*courant*( -1.0*epsilon -1.0*gama + 4.0*Q_W[TAM_C*f + (c)]), -1.0*epsilon ,-1.0*gama);
      printf("\n U_t_2=%.10f \n",U_t_2[TAM_C*f + (c)]);
      */ 
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void esquina_superior_derecha_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
      c=TAM_C -1;
      f=0;
     
      if(U_t_1[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_1=f_w[TAM_C*f + (c-1)] ;
      }
      else if(U_t_1[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_1=f_w[TAM_C*f + (c)] ;
      }
      
      
      
      f_w_3=0.0;

      
      
      f_w_2=0.0;
      
      
      if(U_t_4[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_4=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_4[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_4=f_w[TAM_C*(f+1) + (c)] ;
      }
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + Q_W[TAM_C*f + (c)]);       
      

  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void esquina_inferior_izquierda_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
      c=0;
      f=TAM_F-1;
  
      f_w_1=0.0;
      
      
      
      if(U_t_3[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_3=f_w[TAM_C*(f-1) + (c)] ;
      }
      else if(U_t_3[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_3=f_w[TAM_C*f + (c)] ;
      }

      
      
      if(U_t_2[TAM_C*f + (c)] > 0.0  ){ //si es positiva va para "adelante"
        f_w_2=f_w[TAM_C*f + (c)] ;
      }
      else if(U_t_2[TAM_C*f + (c)] < 0.0){ //si es negativa va para "atras"
	f_w_2=f_w[TAM_C*f + (c+1)] ;
      }
      
      
     f_w_4=0.0;
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + Q_W[TAM_C*f + (c)]);       
      

  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void esquina_inferior_derecha_s(){

  int c=0;
  int f=0;
  
  double f_w_1=0.0;
  double f_w_2=0.0;
  double f_w_3=0.0;
  double f_w_4=0.0;
  
  double gama=0.0;
  double epsilon=0.0;
  
      c=TAM_C -1;
      f=TAM_F -1;
     
      if(U_t_1[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_1=f_w[TAM_C*f + (c-1)] ;
      }
      else if(U_t_1[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_1=f_w[TAM_C*f + (c)] ;
      }
      
      
      
      if(U_t_3[TAM_C*f + c] < 0.0){ //si es negativa va para "adelante"
	f_w_3=f_w[TAM_C*(f-1) + (c)] ;
      }
      else if(U_t_3[TAM_C*f + c] > 0.0){ // si es positiva va para "atras"
        f_w_3=f_w[TAM_C*f + (c)] ;
      }

      
      
     f_w_2=0.0;
     
     f_w_4=0.0;
      
      
      
      
      epsilon=  f_w_1*U_t_1[TAM_C*f + (c)] + f_w_2*U_t_2[TAM_C*f + (c)] + f_w_3*U_t_3[TAM_C*f + (c)] + f_w_4*U_t_4[TAM_C*f + (c)] ;
      gama= f_w_1*lambda_n_1[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_2*lambda_n_2[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] + f_w_3*lambda_n_3[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)]
      + f_w_4*lambda_n_4[TAM_C*f + (c)]*presion_capilar[TAM_C*f + (c)] ;

      saturacion_N[TAM_C*f + (c)] = saturacion[TAM_C*f + (c)] + 1.0*courant*( -1.0*epsilon -1.0*gama + Q_W[TAM_C*f + (c)]);       
      

  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void copiar_s(){

  int c=0;
  int f=0;
  

  
  for(f=0 ; f<(TAM_F) ; f++){
    for(c=0 ; c<(TAM_C) ; c++){
       saturacion[TAM_C*f + c]= saturacion_N[TAM_C*f + c];
    }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////