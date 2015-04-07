#include "full_s_2.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void guardar(double *matriz){    // la funcion guardar graba los datos de matriz en un archivo llamado salia-1.dat
  
  FILE *archivo;   // creo el archivo
  archivo= fopen("salida-1.dat","w+"); //abro el archivo y le asigno el nombre
  int f=0; 
  int c=0;
  
  for(f=0 ; f<TAM_F ; f++){ //barro todas las filas de la matriz
    for(c=0 ; c<TAM_C ; c++){ //barro todas las columnas de la matriz
      fprintf(archivo,"%f\n",matriz[TAM_C*f + c]); //guardo un dato y le doy enter
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
  
  for(f=0 ; f<TAM_F ; f++){
    for(c=0 ; c<TAM_C ; c++){
      saturacion[TAM_C*f + c]= 0.0 ; 
    }
  }
  
  
  if(CB==0){
      for(f=0 ; f<TAM_F ; f++){
	for(c=0 ; c<TAM_C ; c++){
	  saturacion[TAM_C*f + c]= 1.0*exp( -(c*c)/(40.0) ) ; 
	}
      }

  }
  
  
  
  if(CB==1){
      for(f=0 ; f<TAM_F ; f++){
	for(c=0 ; c<TAM_C ; c++){
	  saturacion[TAM_C*f + c]= 1.0*exp( -(c*c + f*f)/(40.0) ) ; 
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
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
    }
  }
 
 
/* Borde izquierdo   c=0 */ 
 for(f=1 ; f<(TAM_F-1) ; f++){
    c=0;
      lambda_d_1[TAM_C*f + c]= 0.0 ; //(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
    
  }
  
  
/* Borde derecho c=TAM_C -1*/  

for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C -1 ;
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_2[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
    
  }
 
/* Borde superior f=0*/ 
 
 
    for(c=1 ; c<(TAM_C-1) ; c++){
      f=0;
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_3[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
    }
  
/* Borde inferior f=TAM_F -1*/
 

    for(c=1 ; c<(TAM_C-1) ; c++){
      f=TAM_F -1 ;
      lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
      lambda_d_4[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
    }
  

/* esquina superior izquierda c=0 f=0*/  
 
 c=0;
 f=0;
 lambda_d_1[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_3[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 

/* esquina superior derecha c= TAM_C -1 f=0*/ 
 
 c=TAM_C -1;
 f=0;
 lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_2[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_3[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_4[TAM_C*f + c]=(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 

 
/* esquina inferior izquierda  c=0  f=TAM_F -1*/ 
 
 c=0;
 f=TAM_F -1;
 lambda_d_1[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_2[TAM_C*f + c]=(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_4[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 

 
/* esquina inferior derecha c=TAM_C -1  f= TAM_F -1*/ 
 c=TAM_C -1;
 f=TAM_F -1;
 lambda_d_1[TAM_C*f + c]=(lambda_d[TAM_C*f + (c-1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c-1)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_2[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*f + (c+1)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*f + (c+1)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_3[TAM_C*f + c]=(lambda_d[TAM_C*(f-1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f-1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 
 lambda_d_4[TAM_C*f + c]=0.0;//(lambda_d[TAM_C*(f+1) + (c)]*lambda_d[TAM_C*f + c])/(lambda_d[TAM_C*(f+1) + (c)] + lambda_d[TAM_C*f + c] + cachito); 


}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



