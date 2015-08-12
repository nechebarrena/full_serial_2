#include "full_s_2.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculo_laplace(){
  
    recinto_interno_l();
    borde_izquierdo_l();
    borde_derecho_l();
    borde_superior_l();
    borde_inferior_l();
    esquina_superior_izquierda_l();
    esquina_superior_derecha_l();
    esquina_inferior_izquierda_l();
    esquina_inferior_derecha_l();
    //presion_N[TAM_C*50 + 50]=1.0;
    calculo_diferencia();
    copiar_l();
    //printf("\n HOLA \n");
  
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void recinto_interno_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
      flujo_1 =  1.0*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  1.0*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  1.0*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 4.0 ;
      presion_N[TAM_C*f + c]= (flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
    }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void borde_izquierdo_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
 
  for(f=1 ; f<(TAM_F-1) ; f++){
      c=0;
      flujo_1 =  0.0;//1.0*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  1.0*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  1.0*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 3.0 ;
      if(CB==1){
      presion_N[TAM_C*f + c]= (flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
      }
      else if(CB==0){
      presion_N[TAM_C*f + c]= PI;
      }
      
  }

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void borde_derecho_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
  for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C -1 ;
      flujo_1 =  1.0*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  0.0;//1.0*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  1.0*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 3.0 ;
      if(CB==1){
      presion_N[TAM_C*f + c]= (flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
      }
      else if(CB==0){
      presion_N[TAM_C*f + c]= -1.0*PI ;	
      }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void borde_superior_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
  f=0;
    for(c=1 ; c<(TAM_C-1) ; c++){
      flujo_1 =  1.0*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  1.0*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  0.0;//1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  1.0*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 3.0 ;
      presion_N[TAM_C*f + c]= (flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
    }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void borde_inferior_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
  f=TAM_F -1;
    for(c=1 ; c<(TAM_C-1) ; c++){
      flujo_1 =  1.0*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  1.0*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  0.0;//1.0*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 3.0 ;
      presion_N[TAM_C*f + c]= (flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
    }
  
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void esquina_superior_izquierda_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
      c=0;
      f=0;
      flujo_1 =  0.0;//1.0*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  1.0*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  0.0;//1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  1.0*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 2.0 ;
      
      presion_N[TAM_C*f + c]= PI ;//(flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
      
  
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void esquina_superior_derecha_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
      c=TAM_C -1;
      f=0;
      flujo_1 =  1.0*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  0.0;//1.0*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  0.0;//1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  1.0*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 2.0 ;
      if(CB==1){
      presion_N[TAM_C*f + c]= (flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
      }
      else if(CB==0){
      presion_N[TAM_C*f + c]= -1.0*PI ;
      }
  

  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void esquina_inferior_izquierda_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
      c=0;
      f=TAM_F-1;
      flujo_1 =  0.0;//*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  1.0*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  0.0;//*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 2.0 ;
      if(CB==1){
      presion_N[TAM_C*f + c]= (flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
      }
      else if(CB==0){
      presion_N[TAM_C*f + c]= PI;	
      }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void esquina_inferior_derecha_l(){
  
  int c=0;
  int f=0;
  
  double flujo_1=0.0;
  double flujo_2=0.0;
  double flujo_3=0.0;
  double flujo_4=0.0;
  double flujo_5=0.0;
  double suma_lambdas=0.0;
  
      c=TAM_C-1;
      f=TAM_F-1;
      flujo_1 =  1.0*presion[TAM_C*f + (c-1)] ;
      flujo_2 =  0.0;//*presion[TAM_C*f + (c+1)] ;
      flujo_3 =  1.0*presion[TAM_C*(f-1) + (c)] ;
      flujo_4 =  0.0;//*presion[TAM_C*(f+1) + (c)] ;
      flujo_5 =  0.0 ;
      suma_lambdas= 2.0 ;
      
      presion_N[TAM_C*f + c]=-1.0*PI;// (flujo_1 + flujo_2 + flujo_3 +  flujo_4 + flujo_5)/(suma_lambdas);
    
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void copiar_l(){

  int c=0;
  int f=0;
  

  
  for(f=0 ; f<(TAM_F) ; f++){
    for(c=0 ; c<(TAM_C) ; c++){
       presion[TAM_C*f + c]= presion_N[TAM_C*f + c];
    }
  }
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void calculo_velocidad_alternativa(){
  
  int c=0;
  int f=0;
  

  
  for(f=1 ; f<(TAM_F-1) ; f++){
    for(c=1 ; c<(TAM_C-1) ; c++){
    
      vel_1[TAM_C*f + c]=(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    }
  }
  
  
  
/* borde izquierdo*/

  for(f=1 ; f<(TAM_F-1) ; f++){
    c=0;
    
      vel_1[TAM_C*f + c]=0.0;//(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    
  }
  
/*borde derecho*/ 

  for(f=1 ; f<(TAM_F-1) ; f++){
    c=TAM_C -1;
    
      vel_1[TAM_C*f + c]=(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=0.0;//(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    
  }
 
/*borde superior*/

  f=0;
    for(c=1 ; c<(TAM_C-1) ; c++){
    
      vel_1[TAM_C*f + c]=(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=0.0;//(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    }
  
/*borde inferior*/

  f=TAM_C-1;
    for(c=1 ; c<(TAM_C-1) ; c++){
    
      vel_1[TAM_C*f + c]=(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=0.0;//(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    }
  

/*esquina superior izquierda*/  
  

c=0;
f=0;
      vel_1[TAM_C*f + c]=0.0;//(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=0.0;//(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    
  
/*esquina superior derecha*/  

c=TAM_C-1;
f=0;
      vel_1[TAM_C*f + c]=(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=0.0;//(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=0.0;//(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    
      
/*esquina inferior izquierda*/

c=0;
f=TAM_C-1;
      vel_1[TAM_C*f + c]=0.0;//(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=0.0;//(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    
/*esquina inferior derecha*/

c=TAM_C-1;
f=TAM_C-1;
      vel_1[TAM_C*f + c]=(presion[TAM_C*f + (c-1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_2[TAM_C*f + c]=0.0;//(presion[TAM_C*f + (c+1)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_3[TAM_C*f + c]=(presion[TAM_C*(f-1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
      vel_4[TAM_C*f + c]=0.0;//(presion[TAM_C*(f+1) + (c)]  - presion[TAM_C*f + (c)]  )*TAM_C;
    
  
  
  
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void copiar_velocidades(){
  
  int f=0;
  int c=0;
  
   for(f=0 ; f<(TAM_F) ; f++){
    for(c=0 ; c<(TAM_C) ; c++){
      
      U_t_1[TAM_C*f + c]=vel_1[TAM_C*f +c];
      U_t_2[TAM_C*f + c]=vel_2[TAM_C*f +c];
      U_t_3[TAM_C*f + c]=vel_3[TAM_C*f +c];
      U_t_4[TAM_C*f + c]=vel_4[TAM_C*f +c];
      
    }
   }
  
  
}


















