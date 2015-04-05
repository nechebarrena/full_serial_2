#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>



#define TAM_C 100
#define TAM_F 100

double *presion;
double *K;
double *Q_W;
double *Q_N;


void guardar(double *matriz);
void inicializar_presion();
void inicializar_K();
void inicializar_Q_W();
void inicializar_Q_N();