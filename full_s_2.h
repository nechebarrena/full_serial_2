#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>



///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/* Definicion de parametros globales                                             */
#define TAM_C 100
#define TAM_F 100
#define CB 1   // 0 para una condicion de borde tipo canal y 1 para una tipo FQS

#define mu_w 1
#define mu_n 2

#define q_n 0.0001
#define q_w 0.0001
#define po 0.1

#define cachito 0.0000000001

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/* Definicion de variables globales                                              */
double *presion;
double *presion_N;
double *presion_capilar;
double *saturacion;
double *K;
double *Q_W;
double *Q_N;
double *lambda_n;
double *lambda_w;
double *krn;
double *krw;
double *lambda_t;
double *lambda_d;
double *termino_presion_capilar;

double *lambda_t_1;
double *lambda_t_2;
double *lambda_t_3;
double *lambda_t_4;

double *lambda_d_1;
double *lambda_d_2;
double *lambda_d_3;
double *lambda_d_4;
double diferencia;
double error;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/* Funciones del archivo funciones_s_2.c                                         */
void guardar(double *matriz);
void inicializar_presion();
void inicializar_saturacion();
void inicializar_K();
void inicializar_Q_W();
void inicializar_Q_N();
void calcular_lambda_n();
void calcular_lambda_w();
void calcular_krw();
void calcular_krn();
void calcular_lambda_t();
void calcular_lambda_d();
void calcular_lambda_t_1_2_3_4();
void calcular_lambda_d_1_2_3_4();
void calcular_presion_capilar();
void calcular_termino_presion_capilar();
void calculo_diferencia();
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/*  Funciones del archivo calculo_presion.c                                      */
void calculo_presion();
void recinto_interno();
void borde_izquierdo();
void borde_derecho();
void borde_superior();
void borde_inferior();
void esquina_superior_izquierda();
void esquina_superior_derecha();
void esquina_inferior_izquierda();
void esquina_inferior_derecha();
void copiar();
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
