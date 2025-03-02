#include "Def.h"

clock_t	 StartT;
clock_t FinishT;
double	TotalT;
double elapsed;


void Creat_Problem() {
	int fstatus, status;
	//CPXENVptr env = NULL;
	//CPXLPptr  lp = NULL;
	if ( !CPXEnv )	// if CPXEnv == 0
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

	//FILE *in;
	Input = Open_File (name,"r");

	fscanf(Input, "%d", &I);
	//cout << I << "\t";

	fscanf(Input, "%d", &J);
	//cout << J << "\t";

	fscanf(Input, "%d", &K);
	//cout << K << "\t";
	
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
	int		*qmatind, *qmatind1;
	
	static char const	qsense = 'L';
	long nzcnt2; 

	
	double TT,WT,wt1,wt2;
	double Objective = 0;
	double partObj1 = 0;
	double temp = 0;

	// Allocate temporary CPLEX data structures
	d_vector(&rhs,2,"open_cplex:2");
	if (rhs == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	d_vector(&rmatval,1000000,"open_cplex:7");
	if (rmatval == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	d_vector(&qmatval,1000000,"open_cplex:7");
	if (qmatval == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	i_vector(&rmatind,1000000,"open_cplex:6");
	if (rmatind == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	i_vector(&qmatind,1000000,"open_cplex:6");
	if (qmatind == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	i_vector(&qmatind1,1000000,"open_cplex:6");
	if (qmatind1 == NULL) {
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

	d_vector(&obj,1000000,"open_cplex:1");
	if (obj == NULL) {
		printf("could not allocate memory for obj.\n");
		exit(-1);
	}

	d_vector(&lb,1000000,"open_cplex:8");
	if (lb == NULL) {
		printf("could not allocate memory for lb.\n");
		exit(-1);
	}

	d_vector(&ub,1000000,"open_cplex:9");
	if (ub == NULL) {
		printf("could not allocate memory for ub.\n");
		exit(-1);
	}

	c_vector(&coltype,1000000,"open_cplex:01");
	if (coltype == NULL) {
		printf("could not allocate memory for coltype.\n");
		exit(-1);
	}

	/*colname = new char  *[100000];
	for (i = 0; i< 100000; i++)
		colname[i] = new char[100];
	if (colname == NULL) {
		printf("could not allocate memory for colname.\n");
		exit(-1);
	}*/


	/////////////////Variables allocate: /////////////////
	x = create_double_vector(I * J * K);
	y = create_double_vector(J * K);
	s = create_double_vector(J * K);
	d = create_double_vector(J * K);
	w = create_double_vector(I * J * K);
	
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
			s[j * K + k] = -1;
		}
	}
	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++) {
			d[j * K + k] = -1;
		}
	}
	/*for (j = 0; j < J; j++){
		for (k = 0; k < K; k++) {
			w[j * K + k] = -1;
		}
	}*/
	for (i = 0; i < I; i++){
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++){
				w[(i * J + j) + ( k * I * J )] = -1;
			}
		}
	}

	////////////		Note
	//Variable Types In Cplex
	//'C' continuous variable j (Default)
	//'B' binary variable j 
	//'I' general integer variable j 
	//'S' semi-continuous variable j 
	//'N' semi-integer variable j 
	///////////

	//////////////	 Define variables in Obj Function	 //////////////

	
	counter = 0;
	for (i = 0; i <I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++){
				//sprintf(colname[counter], "x%3d_%3d_%3d", i, j, k);
				lb[counter] = 0;
				ub[counter] = 1;
				obj[counter] =  (lambda[i] * t[i * J + j] * at) + ( (aw * lambda[i] * (1 - cv[j * K + k] * cv[j * K + k])) / (2 * mu[j * K + k]) );
				coltype[counter++] = 'B';
				x[(i * J + j) + ( k * I * J )] = numcols++;
			}
		}
	}

	for (j = 0; j <J ; j++) {
		for (k = 0; k < K; k++) {
			//sprintf(colname[counter], "y%3d_%3d", j, k);
			lb[counter] = 0;
			ub[counter] = 1;
			obj[counter] = 0;
			coltype[counter++] = 'B';
			y[j * K + k] = numcols++;
		}
	}

	for (j = 0; j <J; j++) {
		for (k = 0; k < K; k++) {
			//sprintf(colname[counter], "r%3d_%3d", j, k);
			lb[counter] = 0;
			ub[counter] = M;
			obj[counter] = (aw * (1 + (float) cv[j * K + k] * cv[j * K + k])) / 2 ;
			coltype[counter++] = 'C';
			s[j * K + k] = numcols++;
		}
	}

	for (j = 0; j <J; j++) {
		for (k = 0; k < K; k++) {
			//sprintf(colname[counter], "d%3d_%3d", j, k);
			lb[counter] = 0;
			ub[counter] = M;
			obj[counter] =0;
			coltype[counter++] = 'C';
			d[j * K + k] = numcols++;
		}
	}
	//for (j = 0; j <J; j++) {
	//	for (k = 0; k < K; k++) {
	//		//sprintf(colname[counter], "lll%3d_%3d", j, k);
	//		lb[counter] = 0;
	//		ub[counter] = M;
	//		obj[counter] =0;
	//		coltype[counter++] = 'C';
	//		w[j * K + k] = numcols++;
	//	}
	//}
	for (i = 0; i <I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++){
				//sprintf(colname[counter], "x%3d_%3d_%3d", i, j, k);
				lb[counter] = 0;
				ub[counter] = M;
				obj[counter] =  0;
				coltype[counter++] = 'C';
				w[(i * J + j) + ( k * I * J )] = numcols++;
			}
		}
	}
	status = CPXnewcols(env, lp, counter, obj, lb, ub, coltype, NULL);
	if (status)	{
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
		exit(-1);
	}

	//status = CPXwriteprob (env, lp, "lpex1.lp", NULL);
	constraints = 0;		//The initial value for "constraints" must be 1	,So here It's OK (Take a look at some lines below => constraints = constraints + 1 )
	
//	Note:	 'G', 'L', 'E', and 'R', indicating greater-than, less-than, equality, and ranged constraints, respectively.

	////////////////  Constraints set (5)	 //////////////		∀j		Σ (lamdai * Xij) <= Σ (mujk * yjk)
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = mu[j * K + k];  // meghdar un vare moadele
			sense[0] = 'L';
			rmatbeg[0] = 0;
			for (i = 0; i < I; i++) {
				if (x[(i * J + j) + ( k * I * J )] >= 0){
					rmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
					rmatval[nzcnt++] = lambda[i];
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
	}

	////////////  Constraints set (6)	 //////////////		∀j		Σ yjk <= 1
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

	////////////  Constraints set (7)	 //////////////		∀i		 Σ Xijk = 1
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

	////////////  Constraints set (8)	 //////////////		∀j		 Σ (lamdai * Xij) >= Σ (RMin * yjk)
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = 0;  
			sense[0] = 'G';
			rmatbeg[0] = 0;
			for (i = 0; i < I; i++) {
				//for (k = 0; k < K; k++) {
					if (x[(i * J + j) + ( k * I * J )] >= 0){
						rmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
						rmatval[nzcnt++] = lambda[i];
					}
				//}
			}
		//	for (k = 0; k < K; k++) {
				if (y[j * K + k] >= 0){
					rmatind[nzcnt] = y[j * K + k];
					rmatval[nzcnt++] = -1 * RMin;
				}
		//	}
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

	////////////  Constraints set (9)	 //////////////		Σ Σ (fjk * yjk) <= B
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

	////////////  Constraints set (16)	 //////////////		∀j,k	 xijk - yjk <= 0
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

	////////////  Constraints set (d)	 //////////////		∀j,k	 lambda * xijk + djk <= mujk
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = 0; //mu[j * K + k];  
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


	////////////  Constraints set (d2)	 //////////////		∀j,k	 lambda * xijk + djk <= mujk
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = 0;  
			sense[0] = 'L';
			rmatbeg[0] = 0;
			for (i = 0; i < I; i++) {
				if (w[(i * J + j) + ( k * I * J )] >= 0){
					rmatind[nzcnt] = w[(i * J + j) + ( k * I * J )];
					rmatval[nzcnt++] = 1;
				}
			}
			if (s[j * K + k] >= 0){
				rmatind[nzcnt] = s[j * K + k];
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
	////////////  Constraints set (SOCP2)	 //////////////		
	for (i = 0; i < I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++) {
				constraints = constraints + 1;
				nzcnt = 0;
				nzcnt2 = 0;
				qrhs = 0;  
				rmatbeg[0] = 0;
			
				if (x[(i * J + j) + ( k * I * J )] >= 0){
					qmatind[nzcnt] = x[(i * J + j) + ( k * I * J )];
					qmatind1[nzcnt] = x[(i * J + j) + ( k * I * J )];
					qmatval[nzcnt++] = lambda[i];
				}

				qmatind[nzcnt] = w[(i * J + j) + ( k * I * J )];
				qmatind1[nzcnt] = d[j * K + k];
				qmatval[nzcnt++] = -1;
				/*if (sd[j * K + k] >= 0){
					rmatind[nzcnt2] = sd[j * K + k];
					rmatval[nzcnt2++] = -1;
				}*/
				if (nzcnt >= 100000) {
					sprintf(errmsg, "MAXNZ exceeded, constraint n. %d \n", constraints);
					printf("%s", errmsg);
					exit(-1);
				}
				rmatbeg[1] = nzcnt;
				if (nzcnt > 0) {
					//status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
					status =  CPXaddqconstr(env,lp,0,nzcnt,qrhs,qsense,NULL,NULL,qmatind,qmatind1,qmatval,NULL);
					//status =  CPXaddqconstr(env,lp,nzcnt2,nzcnt,qrhs,qsense,rmatind,rmatval,qmatind,qmatind,qmatval,NULL);
					if (status)	{
						CPXgeterrorstring(env, status, errmsg);
						printf("%s", errmsg);
						exit(-1);
					}
				}
			}
		}
	}
	


	/////////////////////////////////////////////////////////////////////////

	d_vector(&X_Cplex,numcols,"open_cplex:1");
	numcols = CPXgetnumcols(env,lp);
	numrows = CPXgetnumrows(env,lp);
	CPXcopyctype (env, lp, coltype);

	//sprintf(FileName, "%s", "logfileinfo.txt");
	//strcat (FileName, ".txt");
	//cplex.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
	//cplex.extract(model);
	//cplex.setParam(IloCplex::EpInt, 0);
	//logfile = CPXfopen (FileName, "w");
	//status = CPXsetlogfile (env, logfile);
	//status = CPXwriteprob (env, lp, "lpex1.lp", NULL);
	CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
	CPXsetintparam(env,CPX_PARAM_THREADS, 4);
	//elapsed = ((double)(StartTimeCheck) / CLOCKS_PER_SEC);
	status = CPXsetdblparam(env , CPX_PARAM_TILIM, 7200);		//Run Time = 3 Hours
	status = CPXsetintparam(env , CPX_PARAM_MEMORYEMPHASIS , 1);	//conserve memory where possible
	//status = CPXsetdblparam(env , CPX_PARAM_WORKMEM , 1);			//trigger tree size in MB to create node file (Default 128 MB)
	//status = CPXsetintparam(env , CPX_PARAM_NODEFILEIND , 2);		//store node file on memory or disk?
	CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0001); // e-optimal solution (%gap)
	CPXsetintparam(env , CPX_PARAM_NODEFILEIND , 0);
	//status = CPXsetintparam(env , CPX_PARAM_NODESEL, 0);			//Switch to depth-first only
	//status = CPXsetdblparam(env , CPX_PARAM_BTTOL , 1.0);			//Always choose Depth-first in backtracking
	//status = CPXsetintparam(env , CPX_PARAM_BNDSTRENIND , 1);		//Apply bound strengthening
	//status = CPXsetintparam(env , CPX_PARAM_MIPEMPHASIS	, 1);	//Emphasize optimality over feasibility
	//CPXsetintparam(env, CPX_PARAM_MIQCPSTRAT, 0);


	//cout<<"Solving . . .";
	StartTime =  clock();
gettimeofday(&start, NULL);
	status = CPXmipopt(env, lp);
gettimeofday(&stop, NULL);

	EndTime =  clock();

	if (status) {
		printf ("\nFailed to optimize LP.\n");
		//EndTime = second();
		//cout<<"\nTime:\t"<<EndTime-StartTime;
		//exit(-11);
	}

	status = CPXgetmipx (env, lp, X_Cplex, 0, numcols - 1);
	if ( status ) {
		printf ("\nFailed to obtain solution.\n");
		//EndTime = second();
		//cout<<"\n\nTime:\t"<<EndTime-StartTime;
		//exit(-12);
	}

	CPXgetslack (env, lp, slack, 0, 0);
	status = CPXgetmiprelgap(env, lp, &gap_p);
	printf("\t%f\t",gap_p);
	

	/*EndTime = second();
	float time;
	time = EndTime - StartTime;*/

	//cout<<endl<<"X:__________"<<endl;
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
	//cout<<endl<<"Y:__________"<<endl;
	//////////////	 Print Y	 //////////////
	//cout<<"Y:__________"<<endl;
	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			y[j * K + k] = X_Cplex[counter];
			if (y[j * K + k] < 0.5 ){
				y[j * K + k] = 0;
			}
			counter = counter + 1;
		}
	}

	//////////////	 Print rho	 //////////////
	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			s[j * K + k] = X_Cplex[counter];
			counter = counter + 1;
		}
	}
	//////////////	 Print R	 //////////////
	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			d[j * K + k] = X_Cplex[counter];
			counter = counter + 1;
		}
	}

	//////////////	 Print R	 //////////////
	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			w[j * K + k] = X_Cplex[counter];
			counter = counter + 1;
		}
	}

	objval=0;
	status1 = CPXgetmipobjval (env, lp, &objval);
	if ( status1 ) {
		printf ("Failed to get the MIP objective value.\n");
		//exit(-13);
	}
	objval1 = 0;
	status1 = CPXgetbestobjval (env, lp, &objval1);
	if ( status1 ) {
		printf ("Failed to get the MIP objective value.\n");
		//exit(-13);
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
				Objective += at * lambda[i] * t[(i * J + j)] * x[(i * J + j) + ( k * I * J )];
				TT += lambda[i] * t[i * J + j] * x[(i * J + j) + ( k * I * J )];
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

			//wt1 += ((1 + cv[j * K + k]*cv[j * K + k])*temp*temp)/(2*mu[j * K + k]*(mu[j * K + k]-temp)) + (temp / mu[j * K + k]);
			//wt2 +=  ((1 + cv[j * K + k]*cv[j * K + k]) * (temp / (mu[j * K + k] - temp))) + ((1 - cv[j * K + k]*cv[j * K + k]) * (temp / mu[j * K + k]));
		}
	}


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
				fprintf(Out, "%d\t%d\t|\t", j, k);;
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

fprintf(Out,"s:\t|\t");
	for(j = 0; j < J; j++)
	{
		for(k = 0; k < K; k++)
		{
			if (s[j * K + k] > 0) 
			{
				fprintf(Out, "%d\t%d\t=\t%f|\t", j, k, s[j * K + k]);
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

	fprintf(Out,"W:\t|\t");
	for(i = 0; i < I; i++){
		for(j = 0; j < J; j++)
		{
			for(k = 0; k < K; k++)
			{
				if (w[(i * J + j) + ( k * I * J )] > 0) 
				{
					fprintf(Out, "%d\t%d\t%d\t=\t%f|\t",i, j, k,w[(i * J + j) + ( k * I * J )]); 
				}
			}
		}
	}

	fclose(Out);

	free(rmatind);
	free(qmatind);
	free(rmatbeg);
	free(sense);
	free(coltype);
	//free(colname);
	free(rhs);
	free(rmatval);
	free(qmatval);
	free(obj);
	free(lb);
	free(ub);
	free(X_Cplex);
	objval = 0;
	CPXfreeprob(env, &lp);
	//cout<<"Objective Value: "<<objval<<"\t\t"<<"Time Taken: "<<EndTime-StartTime<<endl;
	
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

		//H = 32;
		ReadData(path);
		/*double travelTime, waitingTime, temp1, temp2, temp3, UhjNewTemp1, UhjNewTemp2, tempUB;
		tempUB = M;*/
		StartT = FinishT = 0;
		//StartTimeCheck = FinishTimeCheck = 0;

		////sprintf(OutName, "%s","Out");
		////sprintf(OutName,"OutFinal2_%.2f_%d_%.1f",aw, B, countCV);
		//sprintf(OutName,"Out_%d",numSample);
		//strcat(OutName, ".txt");
		//Out = fopen(OutName,"w");
		//StartT = second();
		//StartTimeCheck = second();
		//
		Out = Open_File(outfile,"a+");					//Writing what instance we are solving
		fprintf(Out,"\n%s | %s ;\t", argv[1], instance);
		fclose(Out);

		/*CPXENVptr env = NULL;
		CPXLPptr  lp = NULL;*/
		Creat_Problem();
		ModelPopulate();
		
		//writeData();
		free_memory();
	}
	fclose(ini);
//	}
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


// CPLEX functions to allocate memeory to arrays

void i_vector(int **vector,int n,char *s)
{
if((*vector=(int *)calloc(n,sizeof(int)))==NULL)
  //error(s);
  printf("Error: Insuficient memory \n");
return;
}

void d_vector(double **vector,int n,char *s)
{
if((*vector=(double *)calloc(n,sizeof(double)))==NULL)
 // error(s);
 printf("Error: Insuficient memory \n");
return;
}

void c_vector(char **vector,int n,char *s)
{
if((*vector=(char *)calloc(n,sizeof(char)))==NULL)
  //error(s);
  printf("Error: Insuficient memory \n");
return;
}

FILE *Open_File (const char *name, const char *mode)  // This function opens the file "name" under mode "mode" (reading, writing, etc)
{
 FILE *file;

 if((file=fopen(name,mode))==NULL) {
    printf("\nError: File cannot be opened \n");
    //OK=1;
 }
 return file;
}
