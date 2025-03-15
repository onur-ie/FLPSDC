#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <../../../../../encs/pkg/cplex-12.6.3/root/cplex/include/ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>


#include <sys/time.h>


int  CPXEnv, CounterRows = 0;
char FileName[200], InputData[100];

struct timeval start, stop;


float objval1;
double gap_p;
double slack[1];


FILE *Input, *Output;
CPXENVptr env = NULL;
CPXLPptr  lp = NULL;
CPXFILEptr logfile;



double OBJ1=0;
double WT1 = 0;
double TT1=0;
double StartTimeCheck, FinishTimeCheck;

int I, J, K;
double *t; //travel time
int *f; //fix cost of instalation
int *mu;
double *cv;
double *lambda;
int M = 999999;
double RMin = 1; 
double limitationTime = 7200;
double at, aw;
int B;
double epsilon = 0.000001; 
double *z, *w, *rho, *R;
double  *y, *x;
double StartTime, EndTime;
double Uhj[10000];
double Gap = 0;
int H; //for array UHj
const int Q = 200; //for itertion

double *LB, *UB;
double *TT, *WT;
int finishQ;


char OutName[100];
FILE *Out = NULL;

char		outfile[20];

time_t		ttt ;						//Variables for taking the time and knowing when the code was exectued
struct tm	*tm;	


//double second();
void Creat_Problem();
void ReadData();
void writeData ();
void ModelPopulate();


FILE   *Open_File(const char *, const char *);
double *create_double_vector(int);
int *create_int_vector(int);
void free_memory(void);
void i_vector(int **vector,int n,char *s);
void d_vector(double **vector,int n,char *s);
void c_vector(char **vector,int n,char *s);