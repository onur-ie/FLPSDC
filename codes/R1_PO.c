#include "Def.h"

int        instancNumber;
clock_t	 StartT;
clock_t FinishT;
double	TotalT, checkTime;
double elapsed;
int Nnode = 0;
int q;
int totalQ;
int NoChangeBreake;



int	  			cur_numcols;
double			best_lbCall;
int				globaldepth = -1;
int				maxdepth = -1;
int				last_node;
int				usercutsnum;
int				rep_node;
int				flag_upd_av;
int				terminate;
int				size_av;	 
double			tol;
int				nodemult, depmult;
double* average;
int				flag_incumb = 0;
int				chgobj, rLimmtCBobj, rLimmtCB, tLimmtCBobj, tLimmtCB; 
int countbackLim;
int countbackObj;

int doLazy;
double	violation;
int numCutRelax, numCutBranchLazy, numCutBranchuser;

int checkSeq = 0;

int * pos_X, * pos_Y, * pos_r, * pos_d, *pos_d_r, * pos_tempX, * pos_Z, * pos_L;
double * pi, * sigma;



void ReadData(const char* name) {
	int i, j, k;

	Input = Open_File(name, "r");

	fscanf(Input, "%d", &I);


	fscanf(Input, "%d", &J);


	fscanf(Input, "%d", &K);

	Initialize_memory();

	for (i = 0; i < I; i++)
		fscanf(Input, "%lf", &lambda[i]);

	for (i = 0; i < I; i++)
		for (j = 0; j < J; j++)
			fscanf(Input, "%lf", &t[i * J + j]);

	for (j = 0; j < J; j++)
		for (k = 0; k < K; k++)
			fscanf(Input, "%d", &mu[j * K + k]);

	for (j = 0; j < J; j++)
		for (k = 0; k < K; k++)
			fscanf(Input, "%d", &f[j * K + k]);

	for (j = 0; j < J; j++)
		for (k = 0; k < K; k++)
			fscanf(Input, "%lf", &cv[j * K + k]);

	for (i = 0; i < 32; i++)
		fscanf(Input, "%lf", &Uhj[i]);
	fscanf(Input, "%lf", &aw);
	at = 1 - aw;
	fscanf(Input, "%d", &B);
	fclose(Input);



}


void Initialize_memory(void)
{

	lambda = create_double_vector(I);
	t = create_double_vector(I * J);
	mu = create_int_vector(J * K);
	f = create_int_vector(J * K);
	cv = create_double_vector(J * K);

	pi = create_double_vector(J * K * I);
	sigma = create_double_vector(J * K * I);

	x = create_double_vector(J * K * I);
	y = create_double_vector(J * K);

	Sorted = (AVAL*)calloc((J * K), sizeof(AVAL));
	Sorted2 = (AVAL*)calloc((I), sizeof(AVAL));
}



void free_memory(void)
{

	free(lambda);
	free(f);
	free(mu);
	free(cv);
	free(t);

	free(x);
	free(y);


	free(sigma);
	free(pi);

	free(Sorted);
	free(Sorted2);


}

double MISOCP(void)
{
	int	    i, j, k, l, status, status1, counter;
	int index, index1;  
	double best_upper_bound, best_lower_bound;
	int nodecount;
	
	CPXLPptr  lp;      
	CPXENVptr env;    
	int       numcols; 
	int       numrows; 
	int       numnz;   
	int       objsen;  
	double* obj;    
	double* rhs; 
	char* sense; 
	int* matbeg; 
	int* matind; 
	double* matval; 
	double* lbC;    
	double* ub;    
	char	probname[16]; 
	char* ctype;  
	double    value; 
	double	qrhs, * qmatval, * X_Cplex;
	int* qmatind, * qmatind1;

	int      num_x_var, num_Y_var, num_V_var, num_Z_var, num_d_var, num_TX_var, num_dv_var, num_s_var;
	static char const	qsense = 'L';


	CUTINFO lazyconinfo;
	CUTINFO usercutinfo;



	pos_X = create_int_vector(I * J * K);
	pos_Y = create_int_vector(J * K);
	pos_r = create_int_vector(J * K);
	pos_d = create_int_vector(J * K);
	pos_d_r = create_int_vector(J * K);
	pos_tempX = create_int_vector(J * K);
	pos_Z = create_int_vector(J * K);
	pos_L = create_int_vector(J * K);


	depmult = 1;
	nodemult = 100;
	size_av = 2;
	tol = 0.05;
	rLimmtCBobj = 7;
	rLimmtCB = 7;
	tLimmtCBobj = 1;
	tLimmtCB = 1;
	countbackLim = 0;
	countbackObj = 0;

	doLazy = 1;

	objsen = 1; 

	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "MISOCP");
	lp = CPXcreateprob(env, &status, probname);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	index1 = 0;  
	numcols = I * J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lbC, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++) {
				pos_X[(i * J + j) + (k * I * J)] = index1;		
				obj[index1] = (lambda[i] * t[i * J + j] * at) + (aw * (lambda[i] / mu[j * K + k]));
				ctype[index1] = 'C';
				lbC[index1] = 0;
				ub[index1] = 1;
				index1++;
			}
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lbC, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lbC);
	free(ub);
	free(ctype);
	num_x_var = index1;

	//Define x_ij variables
	index1 = 0;  // index of columns
	numcols = J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lbC, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			pos_Y[j * K + k] = num_x_var + index1;
			obj[index1] = 0;
			ctype[index1] = 'B';
			lbC[index1] = 0;
			ub[index1] = 1;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lbC, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lbC);
	free(ub);
	free(ctype);
	num_Y_var = index1;


	
	index1 = 0;  // index of columns
	numcols = J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lbC, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			pos_r[j * K + k] = num_Y_var + num_x_var + index1;
			obj[index1] = aw * (1 + (float)cv[j * K + k] * cv[j * K + k]) / (2 * (float)mu[j * K + k]);
			ctype[index1] = 'C';
			lbC[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lbC, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lbC);
	free(ub);
	free(ctype);
	num_s_var = index1;
	

	index1 = 0;  // index of columns
	numcols = J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lbC, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			pos_d[j * K + k] = num_Y_var + num_x_var + index1 + num_s_var;
			obj[index1] = 0;
			ctype[index1] = 'C';
			lbC[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lbC, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lbC);
	free(ub);
	free(ctype);
	num_d_var = index1;



	index1 = 0;  // index of columns
	numcols = J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lbC, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			pos_d_r[j * K + k] = num_Y_var + num_x_var + index1 + num_d_var + num_s_var;
			obj[index1] = 0;
			ctype[index1] = 'C';
			lbC[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lbC, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lbC);
	free(ub);
	free(ctype);
	num_dv_var = index1;


	index1 = 0;  // index of columns
	numcols = J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lbC, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			pos_tempX[j * K + k] = num_Y_var + num_x_var + index1 + num_d_var + num_s_var + num_dv_var;
			obj[index1] = 0;
			ctype[index1] = 'C';
			lbC[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lbC, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lbC);
	free(ub);
	free(ctype);
	num_TX_var = index1;


	index1 = 0;  // index of columns
	numcols = J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lbC, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			pos_Z[j * K + k] = num_Y_var + num_x_var + index1 + num_d_var + num_s_var + num_dv_var + num_TX_var;
			obj[index1] = 0;
			ctype[index1] = 'C';
			lbC[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lbC, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lbC);
	free(ub);
	free(ctype);
	num_Z_var = index1;



	index1 = 0;  // index of columns
	numcols = J * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lbC, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			pos_L[j * K + k] = num_Y_var + num_x_var + index1 + num_d_var + num_s_var + num_dv_var + num_TX_var, num_Z_var;
			obj[index1] = 0;
			ctype[index1] = 'C';
			lbC[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lbC, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lbC);
	free(ub);
	free(ctype);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	numrows = J * K;
	numnz = (J * K) * (I + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {

			sense[index1] = 'L';
			rhs[index1] = 0;
			matbeg[index1++] = index;

			for (i = 0; i < I; i++) {
				matind[index] = pos_X[(i * J + j) + (k * I * J)];
				matval[index++] = lambda[i];
			}

			matind[index] = pos_Y[j * K + k];
			matval[index++] = -1 * mu[j * K + k];
			

		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	////////////  Constraints set (2)	 //////////////	

	numrows = J;
	numnz = J * K;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		sense[index1] = 'L';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		for (k = 0; k < K; k++) {
			matind[index] = pos_Y[j * K + k];
			matval[index++] = 1;
		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	////////////  Constraints set (3)	 //////////////		
	numrows = I;
	numnz = I * (J * K);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;

	for (i = 0; i < I; i++) {
		sense[index1] = 'E';
		rhs[index1] = 1;
		matbeg[index1++] = index;

		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++) {
				matind[index] = pos_X[(i * J + j) + (k * I * J)];
				matval[index++] = 1;
			}
		}
	}

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	////////////  Constraints set (4)	 //////////////	

	numrows = J * K;
	numnz = (J * K) * (I + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {

			sense[index1] = 'G';
			rhs[index1] = 0;
			matbeg[index1++] = index;

			for (i = 0; i < I; i++) {
				matind[index] = pos_X[(i * J + j) + (k * I * J)];
				matval[index++] = lambda[i];
			}

			matind[index] = pos_Y[j * K + k];
			matval[index++] = -1 * RMin;


		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	////////////  Constraints set (5)	 //////////////	

	numrows = 1;
	numnz = J * K;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	

	sense[index1] = 'L';
	rhs[index1] = B;
	matbeg[index1++] = index;

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			matind[index] = pos_Y[j * K + k];
			matval[index++] = f[j * K + k];
		}
	}

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	////////////  Constraints set (6)	 //////////////	

	numrows = I * J * K;
	numnz = (I * J * K) * (2);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++) {

				sense[index1] = 'L';
				rhs[index1] = 0;
				matbeg[index1++] = index;

				matind[index] = pos_X[(i * J + j) + (k * I * J)];
				matval[index++] = 1;

				matind[index] = pos_Y[j * K + k];
				matval[index++] = -1;

			}
		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	////////////  Constraints set (d)	 //////////////	

	numrows = J * K;
	numnz = (J * K) * (I + 2);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {

			sense[index1] = 'L';
			rhs[index1] = 0;
			matbeg[index1++] = index;

			for (i = 0; i < I; i++) {
				matind[index] = pos_X[(i * J + j) + (k * I * J)];
				matval[index++] = lambda[i];
			}

			matind[index] = pos_Y[j * K + k];
			matval[index++] = -1 * mu[j * K + k];

			matind[index] = pos_d[j * K + k];
			matval[index++] = 1;

		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);



	////////////  Constraints set (dr)	 //////////////	

	numrows = J * K;
	numnz = (J * K) * (3);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {

			sense[index1] = 'G';
			rhs[index1] = 0;
			matbeg[index1++] = index;

			matind[index] = pos_d_r[j * K + k];
			matval[index++] = -1;

			matind[index] = pos_d[j * K + k];
			matval[index++] = 1;

			matind[index] = pos_r[j * K + k];
			matval[index++] = 1;

		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	////////////  Constraints set (tempX)	 //////////////	

	numrows = J * K;
	numnz = (J * K) * (I + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {

			sense[index1] = 'L';
			rhs[index1] = 0;
			matbeg[index1++] = index;

			for (i = 0; i < I; i++) {
				matind[index] = pos_X[(i * J + j) + (k * I * J)];
				matval[index++] = lambda[i];
			}

			matind[index] = pos_tempX[j * K + k];
			matval[index++] = -1;


		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);



	////////////  Constraints set (Redundant)	 //////////////	

	numrows = J * K;
	numnz = (J * K) * (I + 2);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {

			sense[index1] = 'G';
			rhs[index1] = 0;
			matbeg[index1++] = index;

			for (i = 0; i < I; i++) {
				matind[index] = pos_X[(i * J + j) + (k * I * J)];
				matval[index++] = -1 * lambda[i] * lambda[i]; // it can be: -1 * lambda[i]
			}

			matind[index] = pos_r[j * K + k];
			matval[index++] = mu[j * K + k];

			matind[index] = pos_Z[j * K + k];
			matval[index++] = -1;
		}
	}
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//////////////  Constraints set (Redundant_L)	 //////////////	

	numrows = J * K;
	numnz = (J * K) * (I + 2);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {

			sense[index1] = 'L';
			rhs[index1] = 0;
			matbeg[index1++] = index;

			for (i = 0; i < I; i++) {
				matind[index] = pos_X[(i * J + j) + (k * I * J)];
				matval[index++] = -1 * lambda[i] * lambda[i]; // it can be: -1 * lambda[i]
			}

			matind[index] = pos_L[j * K + k];
			matval[index++] = 1;

			matind[index] = pos_Z[j * K + k];
			matval[index++] = -1;
		}
	}
	status = CPXaddusercuts(env, lp, index1, index, rhs, sense, matbeg, matind, matval, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	////////////  Constraints set (SOC)	 //////////////		

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			if (i != j) {
				index = 0;
				index1 = 0;
				numrows = 1;
				numnz = 4;

				i_vector(&matbeg, numrows, "open_cplex:4");
				i_vector(&qmatind, numnz, "open_cplex:6");
				i_vector(&qmatind1, numnz, "open_cplex:6");
				d_vector(&qmatval, numnz, "open_cplex:7");

				qrhs = 0;
				matbeg[index1++] = index;

				qmatind[index] = pos_tempX[j * K + k];
				qmatind1[index] = pos_tempX[j * K + k];
				qmatval[index++] = 2;

				qmatind[index] = pos_d[j * K + k];
				qmatind1[index] = pos_d[j * K + k];
				qmatval[index++] = 1;

				qmatind[index] = pos_r[j * K + k];
				qmatind1[index] = pos_r[j * K + k];
				qmatval[index++] = 1;

				qmatind[index] = pos_d_r[j * K + k];
				qmatind1[index] = pos_d_r[j * K + k];
				qmatval[index++] = -1;

				status = CPXaddqconstr(env, lp, 0, index, qrhs, qsense, NULL, NULL, qmatind, qmatind1, qmatval, NULL);
				if (status)
					fprintf(stderr, "CPXaddrows failed.\n");
				free(matbeg);
				free(qmatind);
				free(qmatind1);
				free(qmatval);
			}
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	CPXsetdblparam(env, CPX_PARAM_TILIM, limitationTime); 
	CPXsetdblparam(env, CPX_PARAM_EPGAP, epsilon2); 
	CPXsetintparam(env, CPX_PARAM_THREADS, threads);
	status = CPXsetdblparam(env, CPX_PARAM_WORKMEM, 12288);	
	CPXsetintparam(env, CPX_PARAM_NODEFILEIND, 0);


	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, CPX_OFF);	
	status = CPXsetintparam(env,CPXPARAM_Preprocessing_Reduce,1);
	status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL); 	
	status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	status = CPXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_ON);


	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, CPX_OFF);
	status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);

	status = makelazyconstraint(env, lp, &lazyconinfo); 											
	status = CPXsetlazyconstraintcallbackfunc(env, lazy_feas_callback1, &lazyconinfo);





	StartTime = clock();
	gettimeofday(&start, NULL);
	status = CPXmipopt(env, lp);  //solve the integer program
	gettimeofday(&stop, NULL);
	EndTime = clock();
	

	i = CPXgetstat(env, lp);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 103)
		printf(" infeasible solution\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);

	
	CPXgetmipobjval(env, lp, &value);
	
	best_upper_bound = value;
	
	CPXgetbestobjval(env, lp, &value);  
	best_lower_bound = value;
	

	nodecount = CPXgetnodecnt(env, lp);
	
	Nnode = nodecount;

	numcols = CPXgetnumcols(env, lp);
	d_vector(&X_Cplex, numcols, "open_cplex:0");
	CPXgetmipx(env, lp, X_Cplex, 0, numcols - 1);  

	status = CPXgetmiprelgap(env, lp, &gap_p);
	

	numrows = CPXgetnumrows(env, lp);

	totalNumrows = numrows;
	totalNumcols = numcols;

	counter = 0;

	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++) {
				x[(i * J + j) + (k * I * J)] = X_Cplex[counter];
				if (x[(i * J + j) + (k * I * J)] < 0.5) {
					x[(i * J + j) + (k * I * J)] = 0;
				}
				else {
					x[(i * J + j) + (k * I * J)] = 1;
				}
				counter = counter + 1;
			}
		}
	}

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			y[j * K + k] = X_Cplex[counter];
			if (y[j * K + k] < 0.5) {
				y[j * K + k] = 0;
			}
			else y[j * K + k] = 1;
			counter = counter + 1;
		}
	}



	if (lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	free(pos_X);
	free(pos_Y);
	free(pos_Z);
	free(pos_d);
	free(pos_d_r);
	free(pos_r);
	free(pos_tempX);
	free(X_Cplex);

	return best_upper_bound;
}



void writeData() {
	int		i, j, k, l, m;
	double Objective = 0;
	Out = Open_File(outfile, "a+");


	TT = WT = wt1 = wt2 = 0;

	for (i = 0; i < I; i++)
	{
		for (j = 0; j < J; j++)
		{
			for (k = 0; k < K; k++) {
				Objective += at * lambda[i] * 10 * t[(i * J + j)] * x[(i * J + j) + (k * I * J)];
				TT += lambda[i] * 10 * t[i * J + j] * x[(i * J + j) + (k * I * J)];
			}
		}
	}

	for (j = 0; j < J; j++)
	{
		for (k = 0; k < K; k++)
		{
			temp = 0;
			for (i = 0; i < I; i++)
				temp += lambda[i] * x[(i * J + j) + (k * I * J)];
			Objective += aw * ((((1 + cv[j * K + k] * cv[j * K + k]) * temp * temp) / (2 * mu[j * K + k] * (mu[j * K + k] - temp))) + (temp / mu[j * K + k]));
			WT += ((((1 + cv[j * K + k] * cv[j * K + k]) * temp * temp) / (2 * mu[j * K + k] * (mu[j * K + k] - temp))) + (temp / mu[j * K + k]));
		}
	}

	fprintf(Out, "%d\t", I);
	fprintf(Out, "%d\t", J);
	fprintf(Out, "%.2f\t", aw);
	fprintf(Out, "%d\t", B);

	
	fprintf(Out, "%f\t", opt_value);
	fprintf(Out, "%f\t", gap_p);

	fprintf(Out, "%d\t", totalNumrows);
	fprintf(Out, "%d\t", totalNumcols);
	fprintf(Out, "%d\t", Nnode);
	fprintf(Out, "%d\t", numCutBranchLazy);

	fprintf(Out, "%f\t", TT);
	fprintf(Out, "%f\t", WT);
	

	fprintf(Out, "%f\t", (double)(EndTime - StartTime) / CLOCKS_PER_SEC);
	fprintf(Out, "%f\t", ((double) (stop.tv_sec - start.tv_sec) * 1000 + (double) (stop.tv_usec - start.tv_usec) / 1000)/1000);

	fprintf(Out, "|\t");
	for (j = 0; j < J; j++)
	{
		for (k = 0; k < K; k++)
		{
			if (y[j * K + k] > 0)
			{
				fprintf(Out, "%d\t%d\t|\t", j, k);
			}
		}
	}


	fclose(Out);
}







double* create_double_vector(int dim)
{
	double* ptr;

	if ((ptr = (double*)calloc(dim, sizeof(double))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	return ptr;
}

int* create_int_vector(int dim)
{
	int* ptr;

	if ((ptr = (int*)calloc(dim, sizeof(int))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	return ptr;
}


void i_vector(int** vector, int n, char* s)
{
	if ((*vector = (int*)calloc(n, sizeof(int))) == NULL)
	
		printf("Error: Insuficient memory \n");
	return;
}

void d_vector(double** vector, int n, char* s)
{
	if ((*vector = (double*)calloc(n, sizeof(double))) == NULL)
		// error(s);
		printf("Error: Insuficient memory \n");
	return;
}

void c_vector(char** vector, int n, char* s)
{
	if ((*vector = (char*)calloc(n, sizeof(char))) == NULL)
		//error(s);
		printf("Error: Insuficient memory \n");
	return;
}

FILE* Open_File(const char* name, const char* mode)  
{
	FILE* file;

	if ((file = fopen(name, mode)) == NULL) {
		printf("\nError: File cannot be opened \n");
		//OK=1;
	}
	return file;
}

int** create_int_matrix(int rows, int Columns)
{
	int i;
	int** ptr;

	if ((ptr = (int**)calloc(rows, sizeof(int*))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0;i < rows;i++)
		ptr[i] = create_int_vector(Columns);
	return ptr;
}

double** create_double_matrix(int rows, int Columns)
{
	int i;
	double** ptr;

	if ((ptr = (double**)calloc(rows, sizeof(double*))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0;i < rows;i++) {
		ptr[i] = create_double_vector(Columns);
	}
	return ptr;
}

int Comparevalue(const void* a, const void* b)
{
	if (((AVAL*)a)->value < ((AVAL*)b)->value)
		return 1;
	if (((AVAL*)a)->value > ((AVAL*)b)->value)
		return -1;
	return 0;
}





static int makelazyconstraint(CPXENVptr  env, CPXLPptr   lp, CUTINFOptr lazyconinfo)
{
	int status = 0;
	lazyconinfo->nodeid = -1;
	cur_numcols = CPXgetnumcols(env, lp);
	lazyconinfo->lp = lp;
	lazyconinfo->numcols = cur_numcols;
	lazyconinfo->nodeobjval = 0.0;
	status = CPXgetbestobjval(env, lp, &best_lbCall);
	return (status);
}

static int makeusercuts(CPXENVptr  env, CPXLPptr   lp, CUTINFOptr usercutinfo)
{
	int status = 0;
	usercutinfo->nodeid = -1;
	cur_numcols = CPXgetnumcols(env, lp);
	usercutinfo->lp = lp;
	usercutinfo->nodeobjval = 0.0;
	usercutinfo->numcols = cur_numcols;
	status = CPXgetbestobjval(env, lp, &best_lbCall);
	return (status);
}


static int CPXPUBLIC user_cut_callback1(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	int status = 0;
	int counter = 0;

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;
	int      numcols = cutinfo->numcols;

	int* cutind;
	double* cutval;
	int* cutind2;
	double* cutval2;

	double   cutvio;
	int      addcuts = 0;
	int		num_cuts = 0;
	int		cutnz;
	int		depth;
	int      optvio = 0;
	int      i, j, k, index = 0, SEQNUM;
	double* temp_x;
	double* objval;
	double   obval;
	double	change;
	int		Flag_feas = 1;
	int		Found_Vio = 0;
	int		com_check;
	double  value;
	
	double    rhs, rhs2;    
	char* sense;  
	int* matbeg, * matbeg2;
	int* matind, * matind2;
	double* matval, * matval2; 
	double		temp1;
	int			l, h;
	
	int       numrows; 
	int       numnz;  
	int		 index1, index2, index12; 
	int		flagVio = 0;
	int indextemp1, indextemp2;
	int oldnodeid;
	double oldnodeobjval;
	double objchg;
	checkSeq = 0;
	objval = &obval;
	
	d_vector(&temp_x, cur_numcols, "open_cplex:0");			
	objval = &obval;										
	numcols = cur_numcols;									
	*useraction_p = CPX_CALLBACK_DEFAULT;					

	
	status = CPXgetcallbacknodex(env, cbdata, wherefrom, temp_x, 0, numcols - 1);			
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, objval);
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &cutinfo->nodeobjval);		
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
	
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &cutinfo->nodeid);
	if (status) {
		fprintf(stderr, "Failed to get node id.\n");
		goto Terminate;
	}
	
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	if (status) {
		fprintf(stderr, "Failed to get depth.\n");
		goto Terminate;
	}



	globaldepth = depth;
	if (globaldepth > maxdepth) maxdepth = globaldepth;
	violation = 0.0001;
	rep_node++;

	flagVio = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			flagVio = 0;
			counter = 0;
			

			for (i = 0; i < I; i++) {
				Sorted2[counter].index = i;
				Sorted2[counter].value = temp_x[pos_X[(i * J + j) + (k * I * J)]];
				counter++;
			}

			qsort((AVAL*)Sorted2, counter, sizeof(Sorted[0]), Comparevalue);

			for (i = 0; i < I; i++) {
				if (flagVio > 10)
					break;
				if (Sorted2[i].value > 0) {
					temp = 0;
					temp1 = 0;
					for (l = 0; l <= i; l++) {
						temp += lambda[Sorted2[l].index];
						if (l != 0)
							temp1 += lambda[Sorted2[l - 1].index];
					}
					pi[Sorted2[i].index] = temp - temp1;
					if (pi[Sorted2[i].index] * Sorted2[i].value > temp_x[pos_Z[j * K + k]]) {
						numrows = 1;
						numnz = 2;
						i_vector(&matbeg, numrows, "open_cplex:4");
						i_vector(&matind, numnz, "open_cplex:6");
						d_vector(&matval, numnz, "open_cplex:7");
						index = 0;
						index1 = 0;


						rhs = 0;
						indextemp1 = index1;
						matbeg[index1++] = index;

				
						matind[index] = pos_X[(Sorted2[i].index1 * J + j) + (k * I * J)];
						matval[index++] = pi[Sorted2[i].index];

						matind[index] = pos_Z[j * K + k];
						matval[index++] = -1;


						status = CPXcutcallbackadd(env, cbdata, wherefrom, index, rhs, 'L', matind, matval, CPX_USECUT_PURGE);
						numCutBranchLazy++;

						if (status)
							fprintf(stderr, "CPXaddrows failed.\n");
						free(matbeg);
						free(matind);
						free(matval);
						flagVio++;
						if (flagVio > 1)
							break;
					}
				}
			}
		}
	}



	if (flagVio > 0) {
		*useraction_p = CPX_CALLBACK_SET;
		usercutsnum++;
		goto Terminate;
	}
	goto Terminate;

Terminate:
	if (*useraction_p != CPX_CALLBACK_SET) {  /*printf("Exited without any violation \n");*/ /*getchar();*/ terminate = 1; }

	free(temp_x);
	return (status);
}

static int CPXPUBLIC lazy_feas_callback1(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	int status = 0;
	int counter;

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;
	int      numcols = cutinfo->numcols;
	int cutnz;

	int* cutind;
	double* cutval;
	int* cutind2;
	double* cutval2;
	double   cutvio;
	int      addcuts = 0;
	int		num_cuts = 0;
	int		addcuts2 = 0;
	int      optvio = 0;
	int      i, j, k, index = 0, SEQNUM;
	double* temp_x;
	double* objval;
	double   obval;
	int			depth;
	int			com_check; 
	double    rhs, rhs2;    
	char* sense;  
	int* matbeg, * matbeg2;
	int* matind, * matind2;
	double* matval, * matval2; 
	double		temp1;
	int			l, h, jj, kk;
	
	int       numrows;
	int       numnz;  
	int		 index1, index2, index12; 
	int		flagVio = 0;
	int indextemp1, indextemp2;
	
	int    oldnodeid;
	double oldnodeobjval;
	int flag2=0;

	d_vector(&temp_x, cur_numcols, "open_cplex:0");

	objval = &obval;												
	numcols = cur_numcols;
	*useraction_p = CPX_CALLBACK_DEFAULT;
	
	oldnodeid = cutinfo->nodeid;
	oldnodeobjval = cutinfo->nodeobjval;
	

	
	status = CPXgetcallbacknodex(env, cbdata, wherefrom, temp_x, 0, numcols - 1);				
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, objval);
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &cutinfo->nodeobjval);		
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
	
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &cutinfo->nodeid);
	if (status) {
		fprintf(stderr, "Failed to get node id.\n");
		goto Terminate;
	}
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	if (status) {
		fprintf(stderr, "Failed to get depth.\n");
		goto Terminate;
	}

	counter = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			Sorted[counter].index = (j * K + k);
			Sorted[counter].value = temp_x[pos_Y[j * K + k]];
			counter++;
		}
	}
	qsort((AVAL*)Sorted, counter, sizeof(Sorted[0]), Comparevalue);
	for (i = 1; i < J * K; i++) {
		if (Sorted[i].value > 0) {
			temp = 0;
			temp1 = 0;
			for (j = 0; j <= i; j++) {
				temp += mu[Sorted[j].index];
				if (j != 0)
					temp1 += mu[Sorted[j-1].index];
			}
			pi[Sorted[i].index] = temp - temp1;
			if (pi[Sorted[i].index] * Sorted[i].value > temp_x[pos_L[Sorted[i].index]]) {
				numrows = 1;
				numnz = 2;
				i_vector(&matbeg, numrows, "open_cplex:4");
				i_vector(&matind, numnz, "open_cplex:6");
				d_vector(&matval, numnz, "open_cplex:7");
				index = 0;
				index1 = 0;
				rhs = 0;
				indextemp1 = index1;
				matbeg[index1++] = index;
				matind[index] = pos_L[Sorted[i].index];
				matval[index++] = -1;
				matind[index] = pos_Y[Sorted[i].index];
				matval[index++] = pi[Sorted[i].index];
				status = CPXcutcallbackadd(env, cbdata, wherefrom, index, rhs, 'L', matind, matval, CPX_USECUT_FORCE);
				numCutBranchLazy++;
				if (status)
					fprintf(stderr, "CPXaddrows failed.\n");
				free(matbeg);
				free(matind);
				free(matval);
				flagVio++;
				if (flagVio > 1)
					break;
			}
		}
	}


	flagVio = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			flagVio = 0;
			counter = 0;
		

			for (i = 0; i < I; i++) {
				Sorted2[counter].index = i;
				Sorted2[counter].value = temp_x[pos_X[(i * J + j) + (k * I * J)]];
				counter++;
			}

			qsort((AVAL*)Sorted2, counter, sizeof(Sorted[0]), Comparevalue);

			for (i = 0; i < I; i++) {
				if (flagVio > 10)
					break;
				if (Sorted2[i].value > 0) {
					temp = 0;
					temp1 = 0;
					for (l = 0; l <= i; l++) {
						temp += lambda[Sorted2[l].index];
						if (l != 0)
							temp1 += lambda[Sorted2[l - 1].index];
					}
					pi[Sorted2[i].index] = temp - temp1;
					if (pi[Sorted2[i].index] * Sorted2[i].value > temp_x[pos_Z[j * K + k]]) {
						numrows = 1;
						numnz = 2;
						i_vector(&matbeg, numrows, "open_cplex:4");
						i_vector(&matind, numnz, "open_cplex:6");
						d_vector(&matval, numnz, "open_cplex:7");
						index = 0;
						index1 = 0;


						rhs = 0;
						indextemp1 = index1;
						matbeg[index1++] = index;

						matind[index] = pos_X[(Sorted2[i].index1 * J + j) + (k * I * J)];
						matval[index++] = pi[Sorted2[i].index];

						matind[index] = pos_Z[j * K + k];
						matval[index++] = -1;


						status = CPXcutcallbackadd(env, cbdata, wherefrom, index, rhs, 'L', matind, matval, CPX_USECUT_PURGE);
						numCutBranchLazy++;

						if (status)
							fprintf(stderr, "CPXaddrows failed.\n");
						free(matbeg);
						free(matind);
						free(matval);
						flagVio++;
						flag2++;
						if (flagVio > 1)
							break;
					}
				}
			}
		}
	}

	
	if (flag2 > 0) {
		*useraction_p = CPX_CALLBACK_SET;
		
		goto Terminate;
	}
	goto Terminate;

Terminate:

	free(temp_x);

	return (status);
}



void main(int argc, char* argv[]) {
	int		i, j, k;
	int		numSample = 0;
	FILE* ini;
	int		MaxNumInst;
	char	instance[20];
	char    path[50];

	if (argc == 1) {
		printf("Error: Input file not specified \n");
		exit(8);
	}
	ini = Open_File(argv[1], "r");
	fscanf(ini, "%d", &MaxNumInst);
	fscanf(ini, "%s", &outfile);

	ttt = time(NULL);
	tm = localtime(&ttt);
	Out = Open_File(outfile, "a+");
	fprintf(Out, "\n %s\n", asctime(tm));
	fclose(Out);

	for (numSample = 1; numSample < MaxNumInst; numSample++) {


		fscanf(ini, "%s", &instance);
		sprintf(path, "./Data/");
		strcat(path, instance);


		ReadData(path);

		StartT = FinishT = 0;

		Out = Open_File(outfile, "a+");					
		fprintf(Out, "\n%s | %s ;\t", argv[1], instance);
		fclose(Out);



		for (i = 0; i < I; i++)
		{
			if (i == 0)
				sigma[i] = lambda[i];
			else
				sigma[i] = sigma[i - 1] + lambda[i];
		}



		opt_value = 0;
		opt_value = MISOCP();

		writeData();
		free_memory();
	}
	fclose(ini);

}