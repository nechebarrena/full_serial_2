#include "full_s_2.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void guardar_1(double *matriz){    // la funcion guardar graba los datos de matriz en un archivo llamado salia-1.dat
  
  FILE *archivo;   // creo el archivo
  archivo= fopen("salida-1.dat","w+"); //abro el archivo y le asigno el nombre
  int f=0; 
  int c=0;
  
  for(f=0 ; f<TAM_F ; f++){ //barro todas las filas de la matriz
    for(c=0 ; c<TAM_C ; c++){ //barro todas las columnas de la matriz
      fprintf(archivo,"%.15f\n",matriz[TAM_C*f + c]); //guardo un dato y le doy enter
    }
      fprintf(archivo,"\n"); //cuando termino una fila le doy un enter
  }
  
  
  fclose(archivo); //cierro el archivo
}
// El formato con el que guardo los datos es [ fila_i -renglon vacio- fila_i+1 -renglon vacio- fila_i+2 ..... ]   
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void guardar_2(double *matriz){    // la funcion guardar graba los datos de matriz en un archivo llamado salia-1.dat
  
  FILE *archivo;   // creo el archivo
  archivo= fopen("salida-2.dat","w+"); //abro el archivo y le asigno el nombre
  int f=0; 
  int c=0;
  
  for(f=0 ; f<TAM_F ; f++){ //barro todas las filas de la matriz
    for(c=0 ; c<TAM_C ; c++){ //barro todas las columnas de la matriz
      fprintf(archivo,"%.15f\n",matriz[TAM_C*f + c]); //guardo un dato y le doy enter
    }
      fprintf(archivo,"\n"); //cuando termino una fila le doy un enter
  }
  
  
  fclose(archivo); //cierro el archivo
}
// El formato con el que guardo los datos es [ fila_i -renglon vacio- fila_i+1 -renglon vacio- fila_i+2 ..... ]   
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void inicializar_presion(){
  
  int c=0;
  int f=0;
  
   
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      presion[TAM_C*f + c]= 0.0 ; 
    }
  }
  
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void inicializar_K(){
  
  int c=0;
  int f=0;
  
   
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      K[TAM_C*f + c]= 1.0 ; 
    }
  }
  
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void inicializar_Q_W(){
  
  int c=0;
  int f=0;
  
   
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      Q_W[TAM_C*f + c]= 0.0 ; 
    }
  }
  
  
  
   if(CB==0){
      for(f=0 ; f<TAM_F ; f++){
	c=0;
	  Q_W[TAM_C*f + c]= q_w; 
	
      }

  }
  
  
  
  if(CB==1){
      f=0 ;
	c=0 ;
	  Q_W[TAM_C*f + c]= q_w ; 
	
      

  }
  
  
  
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void inicializar_Q_N(){
  
  int c=0;
  int f=0;
  
   
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      Q_N[TAM_C*f + c]= 0.0 ; 
    }
  }
  
  
   if(CB==0){
      for(f=0 ; f<TAM_F ; f++){
	c=TAM_C -1;
	  Q_N[TAM_C*f + c]= -q_n; 
	
      }

  }
  
  
  
  if(CB==1){
      f=TAM_F -1 ;
	c=TAM_C-1 ;
	  Q_N[TAM_C*f + c]= -q_n ; 
	
      

  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void inicializar_saturacion(){
  
  int c=0;
  int f=0;
  double alfa=40;
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      saturacion[TAM_C*f + c]= 0.0 ; 
    }
  }
  
  
  if(CB==0){
      for(f=0 ; f<TAM_F ; f++){
	for(c=0 ; c<TAM_C ; c++){
	  saturacion[TAM_C*f + c]= 1.0*exp( -(c*c)/(alfa) ) ; 
	}
      }

  }
  
  
  
  if(CB==1){
      for(f=0 ; f<TAM_F ; f++){
	for(c=0 ; c<TAM_C ; c++){
	  saturacion[TAM_C*f + c]= 1.0*exp( -(c*c + f*f)/(alfa) ) ; 
	}
      }

  }
  
  
  
  
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_lambda_n(){
 
  int c=0;
  int f=0;
  
  calcular_krn();
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      lambda_n[TAM_C*f + c]= (K[TAM_C*f + c]*krn[TAM_C*f + c])/(mu_n) ; 
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_lambda_w(){
  
  int c=0;
  int f=0;
  
  calcular_krw();
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      lambda_w[TAM_C*f + c]= (K[TAM_C*f + c]*krw[TAM_C*f + c])/(mu_w) ; 
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_krn(){
  
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      krn[TAM_C*f + c]=  (1 - saturacion[TAM_C*f + c])*(1 - saturacion[TAM_C*f + c])  ; 
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_krw(){
  
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      krw[TAM_C*f + c]=  ( saturacion[TAM_C*f + c])*( saturacion[TAM_C*f + c])  ; 
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_lambda_t(){
  
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      lambda_t[TAM_C*f + c]= lambda_n[TAM_C*f + c] + lambda_w[TAM_C*f + c]   ; 
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_lambda_d(){
  
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      lambda_d[TAM_C*f + c]= lambda_n[TAM_C*f + c] - lambda_w[TAM_C*f + c]   ; 
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_lambda_d_1_2_3_4(){
  
  int c=0;
  int f=0;
  
/*  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
    }
  }
*/ 

  
  
  
   for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
 
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
    }
  }
 
  
  
  
  
  
  
 
/* Borde izquierdo   c=0 */ 
 for(f=1 ; f<(TAM_F-1) ; f++){
    c=0;
      lambda_d_1[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0);
    
  }
  
  
/* Borde derecho c=TAM_C -1*/  

for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C -1 ;
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
      lambda_d_2[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0);
    
  }
 
/* Borde superior f=0*/ 
 
 
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=0;
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_3[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0);
    }
  
/* Borde inferior f=TAM_F -1*/
 

    for(c=1 ; c<(TAM_C-1) ; c++){
      f=TAM_F -1 ;
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
      lambda_d_4[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0);
    }
  

/* esquina superior izquierda c=0 f=0*/  
 
 c=0;
 f=0;
 lambda_d_1[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
 lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
 lambda_d_3[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
 lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0);

/* esquina superior derecha c= TAM_C -1 f=0*/ 
 
 c=TAM_C -1;
 f=0;
 lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
 lambda_d_2[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
 lambda_d_3[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
 lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0);

 
/* esquina inferior izquierda  c=0  f=TAM_F -1*/ 
 
 c=0;
 f=TAM_F -1;
 lambda_d_1[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
 lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
 lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
 lambda_d_4[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0);

 
/* esquina inferior derecha c=TAM_C -1  f= TAM_F -1*/ 
 c=TAM_C -1;
 f=TAM_F -1;
 lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c])/(2.0);
 lambda_d_2[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c])/(2.0); 
 lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c])/(2.0); 
 lambda_d_4[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c])/(2.0);


}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_lambda_t_1_2_3_4(){
  
  int c=0;
  int f=0;
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
      lambda_t_1[TAM_C*f + c]=(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_2[TAM_C*f + c]=(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_3[TAM_C*f + c]=(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_4[TAM_C*f + c]=(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
    }
  }
 
 
/* Borde izquierdo   c=0 */ 
 for(f=1 ; f<(TAM_F-1) ; f++){
    c=0;
      lambda_t_1[TAM_C*f + c]= 0.0 ; //(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_2[TAM_C*f + c]=(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_3[TAM_C*f + c]=(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_4[TAM_C*f + c]=(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
    
  }
  
  
/* Borde derecho c=TAM_C -1*/  

for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C -1 ;
      lambda_t_1[TAM_C*f + c]=(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_2[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_3[TAM_C*f + c]=(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_4[TAM_C*f + c]=(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
    
  }
 
/* Borde superior f=0*/ 
 
 
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=0;
      lambda_t_1[TAM_C*f + c]=(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_2[TAM_C*f + c]=(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_3[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_4[TAM_C*f + c]=(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
    }
  
/* Borde inferior f=TAM_F -1*/
 

    for(c=1 ; c<(TAM_C-1) ; c++){
      f=TAM_F -1 ;
      lambda_t_1[TAM_C*f + c]=(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_2[TAM_C*f + c]=(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_3[TAM_C*f + c]=(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
      lambda_t_4[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
    }
  

/* esquina superior izquierda c=0 f=0*/  
 
 c=0;
 f=0;
 lambda_t_1[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_2[TAM_C*f + c]=(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_3[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_4[TAM_C*f + c]=(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 

/* esquina superior derecha c= TAM_C -1 f=0*/ 
 
 c=TAM_C -1;
 f=0;
 lambda_t_1[TAM_C*f + c]=(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_2[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_3[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_4[TAM_C*f + c]=(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 

 
/* esquina inferior izquierda  c=0  f=TAM_F -1*/ 
 
 c=0;
 f=TAM_F -1;
 lambda_t_1[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_2[TAM_C*f + c]=(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_3[TAM_C*f + c]=(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_4[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 

 
/* esquina inferior derecha c=TAM_C -1  f= TAM_F -1*/ 
 c=TAM_C -1;
 f=TAM_F -1;
 lambda_t_1[TAM_C*f + c]=(lambda_t[TAM_C*f + (c-1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c-1)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_2[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*f + (c+1)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*f + (c+1)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_3[TAM_C*f + c]=(lambda_t[TAM_C*(f-1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f-1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 
 lambda_t_4[TAM_C*f + c]=0.0;//(lambda_t[TAM_C*(f+1) + (c)]*lambda_t[TAM_C*f + c])/(lambda_t[TAM_C*(f+1) + (c)] + lambda_t[TAM_C*f + c] + cachito); 


}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_presion_capilar(){
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      presion_capilar[TAM_C*f + c]= -po*log(saturacion[TAM_C*f + c])   ; 
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_termino_presion_capilar(){
  int c=0;
  int f=0;
  
  double flujo_1=0;
  double flujo_2=0;
  double flujo_3=0;
  double flujo_4=0;
  
/* Recinto interior */  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
      flujo_1 = lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ; 
    }
  }



/* borde izquierdo c=0 */  
  for(f=1 ; f<(TAM_F-1) ; f++){
    c=0;
      flujo_1 = 0.0;//lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ; 
    }
  



/* borde derecho c=TAM_C -1 */  
  for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C -1;
      flujo_1 = lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = 0.0;//lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ; 
    }
  



/* borde superior f=0 */  
  
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=0;
      flujo_1 = lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = 0.0;//lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ; 
    }
    
    
/* borde inferior f=TAM_F -1 */  
  
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=TAM_F -1;
      flujo_1 = lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = 0.0;//lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ; 
    }
  
  
/* Esquina superior izquierda c=0 f=0 */

      c=0;
      f=0;
      flujo_1 = 0.0;//lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = 0.0;//lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ; 

      
/* Esquina superior derecha c=TAM_C-1 f=0 */

      c=TAM_C -1;
      f=0;
      flujo_1 = lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = 0.0;//lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = 0.0;//lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ;       
   
      
/* Esquina inferior izquierda c=0 f=TAM_F-1 */

      c=0;
      f=TAM_F -1;
      flujo_1 = 0.0;//lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = 0.0;//lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ;       
      
      
/* Esquina inferior derecha c=TAM_C-1 f=TAM_F-1 */

      c=TAM_C-1;
      f=TAM_F-1;
      flujo_1 = lambda_d_1[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c-1)] - presion_capilar[TAM_C*f + c]);
      flujo_2 = 0.0;//lambda_d_2[TAM_C*f + c]*(presion_capilar[TAM_C*f + (c+1)] - presion_capilar[TAM_C*f + c]);
      flujo_3 = lambda_d_3[TAM_C*f + c]*(presion_capilar[TAM_C*(f-1) + (c)] - presion_capilar[TAM_C*f + c]);
      flujo_4 = 0.0;//lambda_d_4[TAM_C*f + c]*(presion_capilar[TAM_C*(f+1) + (c)] - presion_capilar[TAM_C*f + c]);
      
      termino_presion_capilar[TAM_C*f + c]=  0.5*( flujo_1 + flujo_2 + flujo_3 +  flujo_4)  ;       
    
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculo_diferencia(){

  int c=0;
  int f=0;
  double aux=100.0;
  diferencia=0.0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      aux= fabs( presion[TAM_C*f + c] - presion_N[TAM_C*f + c] );
      
      if(aux>diferencia){
	diferencia=aux;
	//printf("\n .....diferencia = %f \n",*diferencia);
      }
      
    }
  }

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
void calcular_presion_n(){
 
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      presion_n[TAM_C*f + c]= presion[TAM_C*f + c] + (presion_capilar[TAM_C*f + c])/(2.0)  ; 
    }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
void calcular_presion_w(){
 
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      presion_w[TAM_C*f + c]= presion[TAM_C*f + c] - (presion_capilar[TAM_C*f + c])/(2.0)  ; 
    }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
void calcular_lambda_n_1_2_3_4(){
  
  
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      lambda_n_1[TAM_C*f + c]= (  lambda_t_1[TAM_C*f + c] + lambda_d_1[TAM_C*f + c] )/(2.0) ;
      lambda_n_2[TAM_C*f + c]= (  lambda_t_2[TAM_C*f + c] + lambda_d_2[TAM_C*f + c] )/(2.0) ;
      lambda_n_3[TAM_C*f + c]= (  lambda_t_3[TAM_C*f + c] + lambda_d_3[TAM_C*f + c] )/(2.0) ;
      lambda_n_4[TAM_C*f + c]= (  lambda_t_4[TAM_C*f + c] + lambda_d_4[TAM_C*f + c] )/(2.0) ;
    }
  }
  
  
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
void calcular_lambda_w_1_2_3_4(){
 
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      lambda_w_1[TAM_C*f + c]= (  lambda_t_1[TAM_C*f + c] - lambda_d_1[TAM_C*f + c] )/(2.0) ;
      lambda_w_2[TAM_C*f + c]= (  lambda_t_2[TAM_C*f + c] - lambda_d_2[TAM_C*f + c] )/(2.0) ;
      lambda_w_3[TAM_C*f + c]= (  lambda_t_3[TAM_C*f + c] - lambda_d_3[TAM_C*f + c] )/(2.0) ;
      lambda_w_4[TAM_C*f + c]= (  lambda_t_4[TAM_C*f + c] - lambda_d_4[TAM_C*f + c] )/(2.0) ;
    }
  }
  
  
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_U_w_1_2_3_4(){
  
  
  int c=0;
  int f=0;
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
      U_w_1[TAM_C*f + c]= -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
      U_w_2[TAM_C*f + c]= -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
      U_w_3[TAM_C*f + c]= -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
      U_w_4[TAM_C*f + c]= -1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );
    }
  }
 
 
/* Borde izquierdo   c=0 */ 
 for(f=1 ; f<(TAM_F-1) ; f++){
    c=0;
      U_w_1[TAM_C*f + c]=0.0;// -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
      U_w_2[TAM_C*f + c]= -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
      U_w_3[TAM_C*f + c]= -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
      U_w_4[TAM_C*f + c]= -1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );
 
  }
  
  
/* Borde derecho c=TAM_C -1*/  

for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C -1 ;
      U_w_1[TAM_C*f + c]= -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
      U_w_2[TAM_C*f + c]=0.0;// -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
      U_w_3[TAM_C*f + c]= -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
      U_w_4[TAM_C*f + c]= -1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );

  }
 
/* Borde superior f=0*/ 
 
 
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=0;
      U_w_1[TAM_C*f + c]= -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
      U_w_2[TAM_C*f + c]= -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
      U_w_3[TAM_C*f + c]=0.0;// -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
      U_w_4[TAM_C*f + c]= -1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );
    }
  
/* Borde inferior f=TAM_F -1*/
 

    for(c=1 ; c<(TAM_C-1) ; c++){
      f=TAM_F -1 ;
      U_w_1[TAM_C*f + c]= -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
      U_w_2[TAM_C*f + c]= -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
      U_w_3[TAM_C*f + c]= -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
      U_w_4[TAM_C*f + c]= 0.0;//-1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );
    }
  

/* esquina superior izquierda c=0 f=0*/  
 
 c=0;
 f=0;
 U_w_1[TAM_C*f + c]=0.0;// -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
 U_w_2[TAM_C*f + c]= -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
 U_w_3[TAM_C*f + c]=0.0;// -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
 U_w_4[TAM_C*f + c]= -1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );

/* esquina superior derecha c= TAM_C -1 f=0*/ 
 
 c=TAM_C -1;
 f=0;
U_w_1[TAM_C*f + c]= -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
U_w_2[TAM_C*f + c]=0.0;// -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
U_w_3[TAM_C*f + c]=0.0;// -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
U_w_4[TAM_C*f + c]= -1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );

 
/* esquina inferior izquierda  c=0  f=TAM_F -1*/ 
 
 c=0;
 f=TAM_F -1;
 U_w_1[TAM_C*f + c]=0.0;// -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
 U_w_2[TAM_C*f + c]= -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
 U_w_3[TAM_C*f + c]= -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
 U_w_4[TAM_C*f + c]=0.0;// -1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );

 
/* esquina inferior derecha c=TAM_C -1  f= TAM_F -1*/ 
 c=TAM_C -1;
 f=TAM_F -1;
 U_w_1[TAM_C*f + c]= -1.0*(lambda_w_1[TAM_C*f + c])*( presion_w[TAM_C*f + (c-1)] - presion_w[TAM_C*f + c]  );
 U_w_2[TAM_C*f + c]=0.0;// -1.0*(lambda_w_2[TAM_C*f + c])*( presion_w[TAM_C*f + (c+1)] - presion_w[TAM_C*f + c]  );
 U_w_3[TAM_C*f + c]= -1.0*(lambda_w_3[TAM_C*f + c])*( presion_w[TAM_C*(f-1) + c] - presion_w[TAM_C*f + c]  );
 U_w_4[TAM_C*f + c]=0.0;// -1.0*(lambda_w_4[TAM_C*f + c])*( presion_w[TAM_C*(f+1) + c] - presion_w[TAM_C*f + c]  );

  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_U_n_1_2_3_4(){
  
  
  
  int c=0;
  int f=0;
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
      U_n_1[TAM_C*f + c]= -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
      U_n_2[TAM_C*f + c]= -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
      U_n_3[TAM_C*f + c]= -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
      U_n_4[TAM_C*f + c]= -1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );
    }
  }
 
 
/* Borde izquierdo   c=0 */ 
 for(f=1 ; f<(TAM_F-1) ; f++){
    c=0;
      U_n_1[TAM_C*f + c]=0.0;// -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
      U_n_2[TAM_C*f + c]= -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
      U_n_3[TAM_C*f + c]= -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
      U_n_4[TAM_C*f + c]= -1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );
 
  }
  
  
/* Borde derecho c=TAM_C -1*/  

for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C -1 ;
      U_n_1[TAM_C*f + c]= -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
      U_n_2[TAM_C*f + c]=0.0;// -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
      U_n_3[TAM_C*f + c]= -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
      U_n_4[TAM_C*f + c]= -1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );

  }
 
/* Borde superior f=0*/ 
 
 
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=0;
      U_n_1[TAM_C*f + c]= -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
      U_n_2[TAM_C*f + c]= -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
      U_n_3[TAM_C*f + c]=0.0;// -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
      U_n_4[TAM_C*f + c]= -1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );
    }
  
/* Borde inferior f=TAM_F -1*/
 

    for(c=1 ; c<(TAM_C-1) ; c++){
      f=TAM_F -1 ;
      U_n_1[TAM_C*f + c]= -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
      U_n_2[TAM_C*f + c]= -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
      U_n_3[TAM_C*f + c]= -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
      U_n_4[TAM_C*f + c]= 0.0;//-1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );
    }
  

/* esquina superior izquierda c=0 f=0*/  
 
 c=0;
 f=0;
 U_n_1[TAM_C*f + c]=0.0;// -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
 U_n_2[TAM_C*f + c]= -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
 U_n_3[TAM_C*f + c]=0.0;// -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
 U_n_4[TAM_C*f + c]= -1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );

/* esquina superior derecha c= TAM_C -1 f=0*/ 
 
 c=TAM_C -1;
 f=0;
U_n_1[TAM_C*f + c]= -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
U_n_2[TAM_C*f + c]=0.0;// -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
U_n_3[TAM_C*f + c]=0.0;// -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
U_n_4[TAM_C*f + c]= -1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );

 
/* esquina inferior izquierda  c=0  f=TAM_F -1*/ 
 
 c=0;
 f=TAM_F -1;
 U_n_1[TAM_C*f + c]=0.0;// -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
 U_n_2[TAM_C*f + c]= -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
 U_n_3[TAM_C*f + c]= -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
 U_n_4[TAM_C*f + c]=0.0;// -1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );

 
/* esquina inferior derecha c=TAM_C -1  f= TAM_F -1*/ 
 c=TAM_C -1;
 f=TAM_F -1;
 U_n_1[TAM_C*f + c]= -1.0*(lambda_n_1[TAM_C*f + c])*( presion_n[TAM_C*f + (c-1)] - presion_n[TAM_C*f + c]  );
 U_n_2[TAM_C*f + c]=0.0;// -1.0*(lambda_n_2[TAM_C*f + c])*( presion_n[TAM_C*f + (c+1)] - presion_n[TAM_C*f + c]  );
 U_n_3[TAM_C*f + c]= -1.0*(lambda_n_3[TAM_C*f + c])*( presion_n[TAM_C*(f-1) + c] - presion_n[TAM_C*f + c]  );
 U_n_4[TAM_C*f + c]=0.0;// -1.0*(lambda_n_4[TAM_C*f + c])*( presion_n[TAM_C*(f+1) + c] - presion_n[TAM_C*f + c]  );

  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_U_t_1_2_3_4(){
 
  int c=0;
  int f=0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      U_t_1[TAM_C*f + c]= 1.0*(U_n_1[TAM_C*f + c] + U_w_1[TAM_C*f + c]) ;
      U_t_2[TAM_C*f + c]= 1.0*(U_n_2[TAM_C*f + c] + U_w_2[TAM_C*f + c]) ;
      U_t_3[TAM_C*f + c]= 1.0*(U_n_3[TAM_C*f + c] + U_w_3[TAM_C*f + c]) ;
      U_t_4[TAM_C*f + c]= 1.0*(U_n_4[TAM_C*f + c] + U_w_4[TAM_C*f + c]) ;
    }
  }
  
  //printf(" \n U_t_1 = %f    U_t_2 = %f    U_t_3 = %f    U_t_4 = %f",U_t_1[TAM_C*0 + 0],U_t_2[TAM_C*0 + 0],U_t_3[TAM_C*0 + 0],U_t_4[TAM_C*0 + 0]);
  
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_f_w(){
  
  int c=0;
  int f=0;
  double arriba=0.0;
  double abajo=0.0;
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      arriba= saturacion[TAM_C*f + c]*saturacion[TAM_C*f + c] ;
      abajo= arriba + (1.0*mu_w/1.0*mu_n)*(1-saturacion[TAM_C*f + c])*(1-saturacion[TAM_C*f + c]) ;
      f_w[TAM_C*f + c]= (1.0*arriba)/(1.0*abajo);
    }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void delta_tiempo(){
  
  double delta_x=1.0/TAM_C;
  delta_t=courant*delta_x*delta_x ;
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






















