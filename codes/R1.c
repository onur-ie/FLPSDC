#include "Def.h"

clock_t	 StartT;
clock_t FinishT;
double	TotalT;
double elapsed;
int Nnode = 0;


void Creat_Problem() {
	int fstatus, status;
	if ( !CPXEnv )	
	{
		env = CPXopenCPLEX (&status);
		if ( env == NULL ) 
		{
			printf ("Could not open CPLEX environment.\n");
			exit(0);
		}
		CPXEnv = 1;
	}

	if (lp != NULL) 
	{
		fstatus = CPXfreeprob (env, &lp);
		if ( fstatus )
			printf("CPXfreeprob fail, error code %d \n", fstatus);
	}

	lp = CPXcreateprob (env, &status, "example");
	CPXchgobjsen (env, lp, CPX_MIN); 
}

void ReadData (const char *name){
	int i , j, k;

	
	Input = Open_File (name,"r");

	fscanf(Input, "%d", &I);


	fscanf(Input, "%d", &J);


	fscanf(Input, "%d", &K);
	
	lambda = create_double_vector(I);
	t = create_double_vector(I * J);
	f = create_int_vector(J * K);
	mu = create_int_vector(J * K);
	cv = create_double_vector(J * K);


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
	at = 1-aw;
	fscanf(Input, "%d", &B);
	fclose(Input);

}

void ModelPopulate (){
	int	    i, j, k, status, status1, counter;
	char    errmsg[1024];
	int		*rmatind, *rmatbeg;
	char	*sense, *coltype;
	long	constraints, nzcnt, numcols, numrows;
	double	*rhs, *rmatval, *obj, *lb, *ub, *X_Cplex, objval, nodecount;


	double	qrhs, *qmatval;
	int		*qmatind;
	
	static char const	qsense = 'L';
	static char const	qsense2 = 'G';


	
	double TT,WT,wt1,wt2;
	double Objective = 0;
	double partObj1 = 0;
	double temp = 0;

	d_vector(&rhs,2,"open_cplex:2");
	if (rhs == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	d_vector(&rmatval,100000,"open_cplex:7");
	if (rmatval == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	d_vector(&qmatval,100000,"open_cplex:7");
	if (qmatval == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	i_vector(&rmatind,100000,"open_cplex:6");
	if (rmatind == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	i_vector(&qmatind,100000,"open_cplex:6");
	if (qmatind == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	i_vector(&rmatbeg,2,"open_cplex:4");
	if (rmatbeg == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	c_vector(&sense,2,"open_cplex:3");
	if (sense == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	d_vector(&obj,100000,"open_cplex:1");
	if (obj == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	d_vector(&lb,100000,"open_cplex:8");
	if (lb == NULL) {
		printf("could not allocate memory for lb.\n");
		exit(-1);
	}

	d_vector(&ub,100000,"open_cplex:9");
	if (ub == NULL) {
		printf("could not allocate memory for ub.\n");
		exit(-1);
	}

	c_vector(&coltype,100000,"open_cplex:01");
	if (coltype == NULL) {
		printf("could not allocate memory for coltype.\n");
		exit(-1);
	}




	/////////////////Variables allocate: /////////////////
	x = create_double_vector(I * J * K);
	y = create_double_vector(J * K);
	r = create_double_vector(J * K);
	d = create_double_vector(J * K);
	lll = create_double_vector(J * K);
	xtemp= create_double_vector(J * K);

	numcols = 0;
	numrows = 0;
	
	///////////////// Initialize ptr-vectors /////////////////
	for (i = 0; i < I; i++){
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++){
				x[(i * J + j) + ( k * I * J )] = -1;
			}
		}
	}

	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++) {
			y[j * K + k] = -1;
		}
	}
	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++) {
			r[j * K + k] = -1;
		}
	}
	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++) {
			d[j * K + k] = -1;
		}
	}
	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++) {
			lll[j * K + k] = -1;
		}
	}
	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++) {
			xtemp[j * K + k] = -1;
		}
	}


	//////////////	 Define variables in Obj Function	 //////////////

	
	counter = 0;
	for (i = 0; i <I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++){
				lb[counter] = 0;
				ub[counter] = 1;
				obj[counter] = (lambda[i]  * t[i * J + j] * at) + (aw * (lambda[i]/ mu[j * K + k]));
				coltype[counter++] = 'B';
				x[(i * J + j) + ( k * I * J )] = numcols++;
			}
		}
	}

	for (j = 0; j <J ; j++) {
		for (k = 0; k < K; k++) {
			lb[counter] = 0;
			ub[counter] = 1;
			obj[counter] = 0;
			coltype[counter++] = 'B';
			y[j * K + k] = numcols++;
		}
	}

	for (j = 0; j <J; j++) {
		for (k = 0; k < K; k++) {
			lb[counter] = 0;
			ub[counter] = CPX_INFBOUND;
			obj[counter] = aw * (1 + (float) cv[j * K + k] * cv[j * K + k])/(2 * (float)mu[j * K + k]);
			coltype[counter++] = 'C';
			r[j * K + k] = numcols++;
		}
	}

	for (j = 0; j <J; j++) {
		for (k = 0; k < K; k++) {
			lb[counter] = 0;
			ub[counter] = CPX_INFBOUND;
			obj[counter] =0;
			coltype[counter++] = 'C';
			d[j * K + k] = numcols++;
		}
	}
	for (j = 0; j <J; j++) {
		for (k = 0; k < K; k++) {
			lb[counter] = 0;
			ub[counter] = CPX_INFBOUND;
			obj[counter] =0;
			coltype[counter++] = 'C';
			lll[j * K + k] = numcols++;
		}
	}
	for (j = 0; j <J; j++) {
		for (k = 0; k < K; k++) {
			lb[counter] = 0;
			ub[counter] = CPX_INFBOUND;
			obj[counter] =0;
			coltype[counter++] = 'C';
			xtemp[j * K + k] = numcols++;
		}
	}

	status = CPXnewcols(env, lp, counter, obj, lb, ub, coltype, NULL);
	if (status)	{
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
		exit(-1);
	}

	constraints = 0;		
	////////////////  Consraints
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = 0; 
			sense[0] = 'L';
			rmatbeg[0] = 0;
			for (i = 0; i < I; i++) {
				if (x[(i * J + j) + ( k * I * J )] >= 0){
					rmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
					rmatval[nzcnt++] = lambda[i];
				}
			}
			if (y[j * K + k] >= 0){
				rmatind[nzcnt] = y[j * K + k];
				rmatval[nzcnt++] = -1 * mu[j * K + k];
			}


			if (nzcnt >= 100000) {
				sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
				printf("%s", errmsg);
				exit(-1);
			}
			rmatbeg[1] = nzcnt;
			if (nzcnt > 0) {
				status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status)	{
					CPXgeterrorstring(env, status, errmsg);
					printf("%s", errmsg);
					exit(-1);
				}
			}
		}
	}


	for (j = 0; j < J; j++) {
		constraints = constraints + 1;
		nzcnt = 0;
		rhs[0] = 1;  
		sense[0] = 'L';
		rmatbeg[0] = 0;
		for (k = 0; k < K; k++) {
			if (y[j * K + k] >= 0){
				rmatind[nzcnt] = y[j * K + k];
				rmatval[nzcnt++] = 1;
			}
		}
		if (nzcnt >= 100000) {
			sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
			printf("%s", errmsg);
			exit(-1);
		}
		rmatbeg[1] = nzcnt;
		if (nzcnt > 0) {
			status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
			if (status)	{
				CPXgeterrorstring(env, status, errmsg);
				printf("%s", errmsg);
				exit(-1);
			}
		}
	}

	for(i = 0; i < I; i++) {
		constraints = constraints + 1;
		nzcnt = 0;
		rhs[0] = 1;  
		sense[0] = 'E';
		rmatbeg[0] = 0;
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++) {
				if (x[(i * J + j) + ( k * I * J )] >= 0){
					rmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
					rmatval[nzcnt++] = 1;
				}
			}
		}
		if (nzcnt >= 100000) {
			sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
			printf("%s", errmsg);
			exit(-1);
		}
		rmatbeg[1] = nzcnt;
		if (nzcnt > 0) {
			status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
			if (status)	{
				CPXgeterrorstring(env, status, errmsg);
				printf("%s", errmsg);
				exit(-1);
			}
		}
	}

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = 0;  
			sense[0] = 'G';
			rmatbeg[0] = 0;
			for (i = 0; i < I; i++) {
					if (x[(i * J + j) + ( k * I * J )] >= 0){
						rmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
						rmatval[nzcnt++] = lambda[i];
					}
			}
				if (y[j * K + k] >= 0){
					rmatind[nzcnt] = y[j * K + k];
					rmatval[nzcnt++] = -1 * RMin;
				}
			if (nzcnt >= 100000) {
				sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
				printf("%s", errmsg);
				exit(-1);
			}
			rmatbeg[1] = nzcnt;
			if (nzcnt > 0) {
				status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status)	{
					CPXgeterrorstring(env, status, errmsg);
					printf("%s", errmsg);
					exit(-1);
				}
			}
		}
	}

	constraints = constraints + 1;
	nzcnt = 0;
	rhs[0] = B;  
	sense[0] = 'L';
	rmatbeg[0] = 0;
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			if (y[j * K + k] >= 0){
				rmatind[nzcnt] = y[j * K + k];
				rmatval[nzcnt++] = f[j * K + k];
			}
		}
	}
	if (nzcnt >= 100000) {
		sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
		printf("%s", errmsg);
		exit(-1);
	}
	rmatbeg[1] = nzcnt;
	if (nzcnt > 0) {
		status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
		if (status)	{
			CPXgeterrorstring(env, status, errmsg);
			printf("%s", errmsg);
			exit(-1);
		}
	}

	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++) {
				constraints = constraints + 1;
				nzcnt = 0;
				rhs[0] = 0;  
				sense[0] = 'L';
				rmatbeg[0] = 0;
				if (x[(i * J + j) + ( k * I * J )] >= 0){
					rmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
					rmatval[nzcnt++] = 1;
				}
				if (y[j * K + k] >= 0){
					rmatind[nzcnt] = y[j * K + k];
					rmatval[nzcnt++] = -1 ;
				}
				if (nzcnt >= 100000) {
					sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
					printf("%s", errmsg);
					exit(-1);
				}
				rmatbeg[1] = nzcnt;
				if (nzcnt > 0) {
					status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status)	{
						CPXgeterrorstring(env, status, errmsg);
						printf("%s", errmsg);
						exit(-1);
					}
				}
			}
		}
	}
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] =  0; 
			sense[0] = 'L';
			rmatbeg[0] = 0;
			for (i = 0; i < I; i++) {
				if (x[(i * J + j) + ( k * I * J )] >= 0){
					rmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
					rmatval[nzcnt++] = lambda[i];
				}
			}
			if (d[j * K + k] >= 0){
				rmatind[nzcnt] = d[j * K + k];
				rmatval[nzcnt++] = 1;
			}

			if (y[j * K + k] >= 0){
				rmatind[nzcnt] = y[j * K + k];
				rmatval[nzcnt++] = -1 *  mu[j * K + k];
			}

			if (nzcnt >= 100000) {
				sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
				printf("%s", errmsg);
				exit(-1);
			}
			rmatbeg[1] = nzcnt;
			if (nzcnt > 0) {
				status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status)	{
					CPXgeterrorstring(env, status, errmsg);
					printf("%s", errmsg);
					exit(-1);
				}
			}
		}
	}

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = 0;  
			sense[0] = 'E';
			rmatbeg[0] = 0;
			if (d[j * K + k] >= 0){
				rmatind[nzcnt] = d[j * K + k];
				rmatval[nzcnt++] = 1;
			}
			if (r[j * K + k] >= 0){
				rmatind[nzcnt] = r[j * K + k];
				rmatval[nzcnt++] = 1;
			}
			if (lll[j * K + k] >= 0){
				rmatind[nzcnt] = lll[j * K + k];
				rmatval[nzcnt++] = -1;
			}
			if (nzcnt >= 100000) {
				sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
				printf("%s", errmsg);
				exit(-1);
			}
			rmatbeg[1] = nzcnt;
			if (nzcnt > 0) {
				status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status)	{
					CPXgeterrorstring(env, status, errmsg);
					printf("%s", errmsg);
					exit(-1);
				}
			}
		}
	}
	
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = 0;  
			sense[0] = 'E';
			rmatbeg[0] = 0;
			for (i = 0; i < I; i++) {
				if (x[(i * J + j) + ( k * I * J )] >= 0){
					rmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
					rmatval[nzcnt++] = lambda[i];
				}
			}
			if (xtemp[j * K + k] >= 0){
				rmatind[nzcnt] = xtemp[j * K + k];
				rmatval[nzcnt++] = -1;
			}
			if (nzcnt >= 100000) {
				sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
				printf("%s", errmsg);
				exit(-1);
			}
			rmatbeg[1] = nzcnt;
			if (nzcnt > 0) {
				status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status)	{
					CPXgeterrorstring(env, status, errmsg);
					printf("%s", errmsg);
					exit(-1);
				}
			}
		}
	}
		
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			qrhs = 0;  
			rmatbeg[0] = 0;
			if (xtemp[j * K + k] >= 0){
				qmatind[nzcnt] = xtemp[j * K + k];
				qmatval[nzcnt++] = 2;
			}
			if (d[j * K + k] >= 0){
				qmatind[nzcnt] = d[j * K + k];
				qmatval[nzcnt++] = 1;
			}
			if (r[j * K + k] >= 0){
				qmatind[nzcnt] = r[j * K + k];
				qmatval[nzcnt++] = 1;
			}
			if (lll[j * K + k] >= 0){
				qmatind[nzcnt] =lll[j * K + k];
				qmatval[nzcnt++] = -1;
			}
			if (nzcnt >= 100000) {
				sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
				printf("%s", errmsg);
				exit(-1);
			}
			rmatbeg[1] = nzcnt;
			if (nzcnt > 0) {
				status =  CPXaddqconstr(env,lp,0,nzcnt,qrhs,qsense,NULL,NULL,qmatind,qmatind,qmatval,NULL);
				if (status)	{
					CPXgeterrorstring(env, status, errmsg);
					printf("%s", errmsg);
					exit(-1);
				}
			}
		}
	}

	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			qrhs = 0;  
			rmatbeg[0] = 0;
			if (xtemp[j * K + k] >= 0){
				qmatind[nzcnt] = xtemp[j * K + k];
				qmatval[nzcnt++] = 2;
			}
			if (d[j * K + k] >= 0){
				qmatind[nzcnt] = d[j * K + k];
				qmatval[nzcnt++] = 1;
			}
			if (r[j * K + k] >= 0){
				qmatind[nzcnt] = r[j * K + k];
				qmatval[nzcnt++] = 1;
			}
			if (lll[j * K + k] >= 0){
				qmatind[nzcnt] =lll[j * K + k];
				qmatval[nzcnt++] = -1;
			}
			if (nzcnt >= 100000) {
				sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
				printf("%s", errmsg);
				exit(-1);
			}
			rmatbeg[1] = nzcnt;
			if (nzcnt > 0) {
				status =  CPXaddqconstr(env,lp,0,nzcnt,qrhs,qsense2,NULL,NULL,qmatind,qmatind,qmatval,NULL);
				if (status)	{
					CPXgeterrorstring(env, status, errmsg);
					printf("%s", errmsg);
					exit(-1);
				}
			}
		}
	}


	/////////////////////////////////////////////////////////////////////////

	d_vector(&X_Cplex,numcols,"open_cplex:1");
	numcols = CPXgetnumcols(env,lp);
	numrows = CPXgetnumrows(env,lp);
	CPXcopyctype (env, lp, coltype);
	CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
	CPXsetintparam(env,CPX_PARAM_THREADS, 4);
	status = CPXsetdblparam(env , CPX_PARAM_TILIM, 7200);		
	status = CPXsetintparam(env , CPX_PARAM_MEMORYEMPHASIS , 1);	
	status = CPXsetdblparam(env , CPX_PARAM_WORKMEM , 1);			
	CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0001); 

CPXsetintparam(env , CPX_PARAM_NODEFILEIND , 0);


	//cout<<"Solving . . .";
	StartTime =  clock();
gettimeofday(&start, NULL);
	status = CPXmipopt(env, lp);
	EndTime =  clock();
gettimeofday(&stop, NULL);

	if (status) {
		printf ("\nFailed to optimize LP.\n");
	}

	status = CPXgetmipx (env, lp, X_Cplex, 0, numcols - 1);
	if ( status ) {
		printf ("\nFailed to obtain solution.\n");
	}

	CPXgetslack (env, lp, slack, 0, 0);
	status = CPXgetmiprelgap(env, lp, &gap_p);
	printf("\t%f\t",gap_p);
	
	counter = 0;
	for (i = 0; i < I; i++) {
		for (j = 0; j <J; j++) {
			for (k=0; k < K ; k++){
				x[(i * J + j) + ( k * I * J )] = X_Cplex[counter];
				if (x[(i * J + j) + ( k * I * J )] <0.5 ){				
					x[(i * J + j) + ( k * I * J )] = 0;
				}
				counter = counter + 1;
			}
		}
	}
	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			y[j * K + k] = X_Cplex[counter];
			if (y[j * K + k] < 0.5 ){
				y[j * K + k] = 0;
			}
			counter = counter + 1;
		}
	}


	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			r[j * K + k] = X_Cplex[counter];
			counter = counter + 1;
		}
	}

	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			d[j * K + k] = X_Cplex[counter];
			counter = counter + 1;
		}
	}

	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			lll[j * K + k] = X_Cplex[counter];
			counter = counter + 1;
		}
	}
	
	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			xtemp[j * K + k] = X_Cplex[counter];
			counter = counter + 1;
		}
	}
	objval=0;
	status1 = CPXgetmipobjval (env, lp, &objval);
	if ( status1 ) {
		printf ("Failed to get the MIP objective value.\n");
	}
	objval1 = 0;
	status1 = CPXgetbestobjval (env, lp, &objval1);
	if ( status1 ) {
		printf ("Failed to get the MIP objective value.\n");
	}
	nodecount = 0;
	nodecount = CPXgetnodecnt (env, lp);
	numcols = CPXgetnumcols (env, lp);
	


	TT = WT = wt1 = wt2 = 0;
	
	for(i = 0; i < I; i++)
	{
		for(j = 0; j < J; j++)
		{
			for(k = 0; k < K; k++){
				Objective += at * lambda[i] * 10 * t[(i * J + j)] * x[(i * J + j) + ( k * I * J )];
				TT += lambda[i] * 10 * t[i * J + j] * x[(i * J + j) + ( k * I * J )];
			}
		}
	}
	partObj1 = Objective;
	for(j = 0; j < J; j++)
	{
		for(k = 0; k < K; k++)
		{
			temp = 0;
			for(i = 0; i < I; i++)
				temp += lambda[i] * x[(i * J + j) + ( k * I * J )];
			Objective += aw * ((((1+cv[j * K + k]*cv[j * K + k])*temp*temp) / (2*mu[j * K + k] * (mu[j * K + k] - temp))) + (temp / mu[j * K + k]));
			WT += ((((1+cv[j * K + k]*cv[j * K + k])*temp*temp) / (2*mu[j * K + k] * (mu[j * K + k] - temp))) + (temp / mu[j * K + k]));

		}
	}


	printf("\nObj: %f\n", objval); 
	printf("Obj_changed: %f\n", Objective); 


	Out = Open_File(outfile,"a+");

	fprintf(Out,"%d\t",I);
	fprintf(Out,"%d\t", J);	
	fprintf(Out,"%.2f\t", aw);
	fprintf(Out,"%d\t", B);
	fprintf(Out, "%f\t", objval);
	fprintf(Out, "%f\t", objval1);
	fprintf(Out, "%f\t", gap_p);
	fprintf(Out, "%f\t", nodecount );
	fprintf(Out, "%d\t", numcols);
	fprintf(Out, "%f\t", Objective);
	fprintf(Out, "%f\t", TT);
	fprintf(Out, "%f\t", WT);
	fprintf(Out, "%f\t", (double) (EndTime-StartTime)/ CLOCKS_PER_SEC);

fprintf(Out, "%f\t", ((double) (stop.tv_sec - start.tv_sec) * 1000 + (double) (stop.tv_usec - start.tv_usec) / 1000)/1000);
	fprintf(Out, "|\t");

	for(j = 0; j < J; j++)
	{
		for(k = 0; k < K; k++)
		{
			if (y[j * K + k] > 0) 
			{
				fprintf(Out, "%d\t%d\t|\t", j, k);
			}
		}
	}

	fprintf(Out,"X:\t|\t");
	for(i = 0; i < I; i++){
		for(j = 0; j < J; j++)
		{
			for(k = 0; k < K; k++)
			{
				if (x[(i * J + j) + ( k * I * J )] > 0) 
				{
					fprintf(Out, "%d\t%d\t%d\t|\t",i, j, k); 
				}
			}
		}
	}

	fprintf(Out,"r:\t|\t");
	for(j = 0; j < J; j++)
	{
		for(k = 0; k < K; k++)
		{
			if (r[j * K + k] > 0) 
			{
				fprintf(Out, "%d\t%d\t=\t%f|\t", j, k, r[j * K + k]);
			}
		}
	}
fprintf(Out,"d:\t|\t");
	for(j = 0; j < J; j++)
	{
		for(k = 0; k < K; k++)
		{
			if (d[j * K + k] > 0) 
			{
				fprintf(Out, "%d\t%d\t=\t%f|\t", j, k, d[j * K + k]);
			}
		}
	}
fprintf(Out,"lll:\t|\t");
	for(j = 0; j < J; j++)
	{
		for(k = 0; k < K; k++)
		{
			if (lll[j * K + k] > 0) 
			{
				fprintf(Out, "%d\t%d\t=\t%f|\t", j, k, lll[j * K + k]);
			}
		}
	}
fprintf(Out,"xtemp:\t|\t");
	for(j = 0; j < J; j++)
	{
		for(k = 0; k < K; k++)
		{
			if (xtemp[j * K + k] > 0) 
			{
				fprintf(Out, "%d\t%d\t=\t%f|\t", j, k, xtemp[j * K + k]);
			}
		}
	}

	fclose(Out);

	free(rmatind);
	free(qmatind);
	free(rmatbeg);
	free(sense);
	free(coltype);
	free(rhs);
	free(rmatval);
	free(qmatval);
	free(obj);
	free(lb);
	free(ub);
	free(X_Cplex);
	objval = 0;
	CPXfreeprob(env, &lp);
	
}



void main (int argc, char *argv[]){
	int		i , j , k;
	int		numSample = 0;
	FILE	*ini;
	int		MaxNumInst;
	char	instance[20];  
	char    path[50];

	if(argc == 1) {
		printf("Error: Input file not specified \n");
		exit(8);
	}
	ini = Open_File (argv[1],"r");
	fscanf(ini, "%d", &MaxNumInst);
	fscanf(ini, "%s", &outfile);

	ttt= time(NULL);
	tm = localtime(&ttt);
	Out = Open_File(outfile,"a+");	
	fprintf(Out, "\n %s\n", asctime(tm));
	fclose(Out);

	for (numSample=1; numSample<MaxNumInst; numSample++){
		

		fscanf(ini,"%s",&instance);
		sprintf(path,"./New Data/");
		strcat(path,instance);

		ReadData(path);
		StartT = FinishT = 0;

		Out = Open_File(outfile,"a+");				
		fprintf(Out,"\n%s | %s ;\t", argv[1], instance);
		fclose(Out);

		Creat_Problem();
		ModelPopulate();

		free_memory();
	}
	fclose(ini);
}





void free_memory(void)
{
  free(lambda);
  free(t);
  free(f);
  free(mu);
  free(cv);
}

double *create_double_vector (int dim)
{
 double *ptr;

 if((ptr=(double *) calloc (dim, sizeof(double)))==NULL) {
    printf("\nError: Insuficient memory \n");
    exit(8);
 }
 return ptr;
}

int *create_int_vector (int dim)
{
 int *ptr;

 if((ptr=(int *) calloc (dim, sizeof(int)))==NULL) {
    printf("\nError: Insuficient memory \n");
    exit(8);
 }
 return ptr;
}

void i_vector(int **vector,int n,char *s)
{
if((*vector=(int *)calloc(n,sizeof(int)))==NULL)
  printf("Error: Insuficient memory \n");
return;
}

void d_vector(double **vector,int n,char *s)
{
if((*vector=(double *)calloc(n,sizeof(double)))==NULL)
 printf("Error: Insuficient memory \n");
return;
}

void c_vector(char **vector,int n,char *s)
{
if((*vector=(char *)calloc(n,sizeof(char)))==NULL)
  printf("Error: Insuficient memory \n");
return;
}

FILE *Open_File (const char *name, const char *mode)  
{
 FILE *file;

 if((file=fopen(name,mode))==NULL) {
    printf("\nError: File cannot be opened \n");
 }
 return file;
}
