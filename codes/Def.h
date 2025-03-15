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

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))



double objval1;
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
double timeDivider; 

int I, J, K;
double *t; 
int *f; 
int *mu;
double *cv;
double *lambda;
int M = 999999;
double RMin = 1; 
double limitationTime = 12*60*60;
double at, aw;
int B;
double epsilon = 0.0001; 
double *d,*lll;
double  *y, *x, *r, *xtemp;
double StartTime, EndTime;
double Uhj[100];
double Gap = 0;
int finishQ;


char OutName[100];
FILE *Out = NULL;

void Creat_Problem();
void ReadData (const char *name);
void ModelPopulate();


FILE   *Open_File(const char *, const char *);
double *create_double_vector(int);
int *create_int_vector(int);
void free_memory(void);
void i_vector(int **vector,int n,char *s);
void d_vector(double **vector,int n,char *s);
void c_vector(char **vector,int n,char *s);


char		outfile[20];
time_t		ttt ;					
struct tm	*tm;	

////coment out the rest for PO
//int CPXEnv1;
//CPXENVptr env1 = NULL;
//CPXLPptr  lp1 = NULL;
//struct timeval start, stop, start1, stop1;
//static int makeusercuts(CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);
//static int CPXPUBLIC user_cut_callback1(CPXCENVptr env, void* cbdata, int wherefrom, void* info, int* useraction_p);
//static int makelazyconstraint(CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);
//static int CPXPUBLIC lazy_feas_callback1(CPXCENVptr env, void* cbdata, int wherefrom, void* info, int* useraction_p);
//double MISOCP();
//int Comparevalue(const void* a, const void* b);
//int** create_int_matrix(int, int);
//double** create_double_matrix(int, int);
//void Initialize_memory();
//typedef struct AVAL
//{
//	int index;
//	double value;
//	int index1;
//	int index2;
//} AVAL;
//AVAL* Sorted;
//AVAL* Sorted2;
//struct cutinfo {
//	CPXLPptr lp;
//	int    numcols;
//	double bestlb;
//	long nodeid;
//	double nodeobjval;
//	int objsen;
//};
//typedef struct cutinfo CUTINFO, * CUTINFOptr;  
//int			vioptcuts;//