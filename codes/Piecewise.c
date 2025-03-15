#include "PI-Def.h"

clock_t	 StartT;
clock_t FinishT;
double	TotalT, checkTime;
double elapsed;
int q;
int totalQ;
int nodecount = 0;


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
	int i , j , k;

	FILE *in;
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

	for (i = 0; i < I; i++){
		for (j = 0; j < J; j++){
			fscanf(Input, "%lf", &t[i * J + j]);
		}
	}

	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++)
			fscanf(Input, "%d", &mu[j * K + k]);
	}

	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++)
			fscanf(Input, "%d", &f[j * K + k]);
	}

	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++)
			fscanf(Input, "%lf", &cv[j * K + k]);
	}

	for (i = 0; i < H; i++)
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
	char	*sense, *coltype, **colname;
	long	constraints, nzcnt, numcols, numrows;
	double	*rhs, *rmatval, *obj, *lb, *ub, *X_Cplex, objval;
	double temp1 =0;
	int h;
	double time;


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

	i_vector(&rmatind,100000,"open_cplex:6");
	
	if (rmatind == NULL) {
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


	x = create_double_vector(I * J * K);
	y = create_double_vector(J * K);
	rho = create_double_vector(J * K);
	R = create_double_vector(J * K);
	

	numcols = 0;
	numrows = 0;
	
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
			rho[j * K + k] = -1;
		}
	}
	for (j = 0; j < J; j++){
		for (k = 0; k < K; k++) {
			R[j * K + k] = -1;
		}
	}


	counter = 0;
	for (i = 0; i <I; i++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++){
				//sprintf(colname[counter], "x%3d_%3d_%3d", i, j, k);
				lb[counter] = 0;
				ub[counter] = 1;
				obj[counter] = lambda[i] * t[i * J + j] * at;
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
			//sprintf(colname[counter], "rho%3d_%3d", j, k);
			lb[counter] = 0;
			ub[counter] = 1;
			obj[counter] = aw * 0.5 * (1 - cv[j * K + k] * cv[j * K + k]);
			coltype[counter++] = 'C';
			rho[j * K + k] = numcols++;
		}
	}

	for (j = 0; j <J; j++) {
		for (k = 0; k < K; k++) {
			//sprintf(colname[counter], "R%3d_%3d", j, k);
			lb[counter] = 0;
			ub[counter] = M;
			obj[counter] = aw * 0.5 * (1 + cv[j * K + k] * cv[j * K + k]);
			coltype[counter++] = 'C';
			R[j * K + k] = numcols++;
		}
	}

	
	status = CPXnewcols(env, lp, counter, obj, lb, ub, coltype, NULL);
	if (status)	{
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
		exit(-1);
	}

	
	constraints = 0;	
	for (j = 0; j < J; j++) {
		for (k = 0; k < K; k++) {
			constraints = constraints + 1;
			nzcnt = 0;
			rhs[0] = mu[j * K + k];  
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
			if (rho[j * K + k] >= 0){
				rmatind[nzcnt] = rho[j * K + k];
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

	
	for (h = 0; h < H; h++) {
		for (j = 0; j < J; j++) {
			for (k = 0; k < K; k++) {
				constraints = constraints + 1;
				nzcnt = 0;
				rhs[0] =  (((Uhj[h] * Uhj[h])) / ((1 + Uhj[h]) * (1 + Uhj[h]))) ;
				sense[0] = 'L';
				rmatbeg[0] = 0;
				if (R[j * K + k] >= 0){
					rmatind[nzcnt] = R[j * K + k];
					rmatval[nzcnt++] = -1 /  ((1 + Uhj[h]) * (1 + Uhj[h]));
				}
				if (rho[j * K + k] >= 0){
					rmatind[nzcnt] = rho[j * K + k];
					rmatval[nzcnt++] = 1;
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


	d_vector(&X_Cplex,numcols,"open_cplex:1");
	numcols = CPXgetnumcols(env,lp);
	numrows = CPXgetnumrows(env,lp);
	CPXcopyctype (env, lp, coltype);

	CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
	CPXsetintparam(env,CPX_PARAM_THREADS, 4);

	gettimeofday(&stop, NULL);
	elapsed = ((double) (stop.tv_sec - start.tv_sec) * 1000 + (double) (stop.tv_usec - start.tv_usec) / 1000)/1000;

	status = CPXsetdblparam(env , CPX_PARAM_TILIM, limitationTime - elapsed);	
	status = CPXsetintparam(env , CPX_PARAM_MEMORYEMPHASIS , 1);
	CPXsetdblparam(env,CPX_PARAM_EPGAP,  0.000001); 
CPXsetintparam(env , CPX_PARAM_NODEFILEIND , 0);

	StartTime = clock();
	status = CPXmipopt(env, lp);
	EndTime = clock();

	if (status) {
		printf ("\nFailed to optimize LP.\n");
		EndTime = clock();
	}

	status = CPXgetmipx (env, lp, X_Cplex, 0, numcols - 1);
	if ( status ) {
		printf ("\nFailed to obtain solution.\n");
		EndTime = clock();
	}

	CPXgetslack (env, lp, slack, 0, 0);
	status = CPXgetmiprelgap(env, lp, &gap_p);
	printf("\t%f\t",gap_p);

	nodecount += CPXgetnodecnt (env, lp);


	EndTime = clock();
	
	time = EndTime - StartTime;


	counter = 0;
	for (i = 0; i < I; i++) {
		for (j = 0; j <J; j++) {
			for (k=0; k < K ; k++){
				x[(i * J + j) + ( k * I * J )] = X_Cplex[counter];
				if (x[(i * J + j) + ( k * I * J )] <0.5 )
					x[(i * J + j) + ( k * I * J )] = 0;
				counter = counter + 1;
			}
		}
	}

	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			y[j * K + k] = X_Cplex[counter];
			if (y[j * K + k] < 0.5 )
					y[j * K + k] = 0;
			counter = counter + 1;
		}
	}

	objval=0;
	status1 = CPXgetmipobjval (env, lp, &objval);
	if ( status1 ) {
		printf ("Failed to get the MIP objective value.\n");
		//exit(-13);
	}
	objval1 = objval;

	WT1 = 0;
	TT1=0;
	for (i = 0; i < I; i++) {
		for (j = 0; j <J; j++) {
			for (k = 0; k <K; k++) 
				if (lambda[i] * t[i * J + j] * x[(i * J + j) + ( k * I * J )] > 0)
					TT1 += lambda[i] * t[i * J + j] * x[(i * J + j) + ( k * I * J )];
		}
	}
	
	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			temp1 = 0;
			for (i = 0; i < I; i++) {
				temp1 += lambda[i] * x[(i * J + j) + ( k * I * J )];
			}
			if ( mu[j * K + k] - temp1 != 0 && mu[j * K + k] !=0){
				if (((1 + cv[j * K + k]* cv[j * K + k] ) * R[j * K + k]) +  ((1 - cv[j * K + k]* cv[j * K + k] ) * rho[j * K + k]) > 0)
					WT1 +=  ((1 + cv[j * K + k]* cv[j * K + k] ) * R[j * K + k]) +  ((1 - cv[j * K + k]* cv[j * K + k] ) * rho[j * K + k]);
			}
		}
	}
	OBJ1=0;
 	OBJ1 = at * TT1 + 0.5 * aw * WT1;
	
	OBJ1= objval;
	free(rmatind);
	free(rmatbeg);
	free(sense);
	free(coltype);
	free(rhs);
	free(rmatval);
	free(obj);
	free(lb);
	free(ub);
	free(X_Cplex);
	objval = 0;

}

void writeData (){
	int  i, j , k;

	printf("%s", "X:\n");
	for (i = 0; i < I; i++) {
		for (j = 0; j <J; j++) {\
			for (k = 0; k <K; k++) {
				if (x[(i * J + j) + ( k * I * J )]>0) {
					printf("X[ %d ][ %d ][ %d ]= %f\n", i, j,k , x[(i * J + j) + ( k * I * J )]);
				}
			}
		}
	}


	printf("%s", "\nY:\n");
	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			if (y[j * K + k]>0) {
				printf("Y[ %d ][ %d ]= %f\n", j, k, y[j * K + k]);
			}
		}
	}


	printf("%s", "\nrho:\n");
	for (j = 0; j <J; j++) {
		if (rho[j * K + k]>0) {
			printf("rho[ %d ]= %f\n", j, rho[j * K + k]);
		}
	}

	printf("%s", "\nR:\n");
	for (j = 0; j <J; j++) {
		if (R[j * K + k]>0) {
			printf("R[ %d ]= %f\n", j, (double) R[j * K + k]);
		}
	}

	printf("%s%f%s", "\nObjective Value: ", LB[finishQ], "\n");
	printf("%s%f%s", "\nTravel Time: ", TT[finishQ], "\n");
	printf("%s%f%s", "\nWating Time: ", WT[finishQ]/2, "\n");
	printf("%s%f%s", "\nGAP: ", Gap, "\n");
	printf("%s%f%s", "\nTime Taken: ", EndTime - StartTime, "\n");
	printf("\nTime= %f\n", TotalT);



	Out = Open_File(outfile,"a+");
	fprintf(Out,"%d\t",I);
	fprintf(Out,"%d\t", J);
	fprintf(Out,"%.2f\t",aw);
	fprintf(Out,"%d\t", B);
	fprintf(Out,"%.1f\t", cv[1]);
	fprintf(Out,"%d\t", H);
	totalQ = q+1; 
	fprintf(Out,"%d\t", totalQ);
	fprintf(Out,"%d\t", nodecount);
	fprintf(Out, "%f\t", LB[finishQ]);
	fprintf(Out, "%f\t", UB[finishQ]);
	fprintf(Out, "%f\t", TT[finishQ]);
	fprintf(Out, "%f\t", WT[finishQ]/2);

gettimeofday(&stop, NULL);


	fprintf(Out, "%f\t", ((double) (stop.tv_sec - start.tv_sec) * 1000 + (double) (stop.tv_usec - start.tv_usec) / 1000)/1000);
	fprintf(Out, "|\t");

	for (j = 0; j <J; j++) {
		for (k = 0; k <K; k++) {
			if (y[j * K + k]>0) {
				fprintf(Out, "%d\t%d\t|\t", j, k);
			}
		}
	}


	fprintf(Out,"\n");
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

	
	fclose(Out);
}



FILE *Open_File (const char *name, const char *mode) 
{
 FILE *file;

 if((file=fopen(name,mode))==NULL) {
    printf("\nError: File cannot be opened \n");
    //OK=1;
 }
 return file;
}



void main (int argc, char *argv[]){
	int		i , j , k;
	char	instance[20];  
	char    path[50];
	int		numSample;
	FILE	*ini;
	int		MaxNumInst;
	int hh;
	double tempmu =0;
	double travelTime, waitingTime, temp1, temp2, temp3, UhjNewTemp1, UhjNewTemp2, tempUB;
	int jj;
	double waitingTime2;

	if(argc == 1) {
		printf("Error: Input file not specified \n");
		exit(8);
	}
	ini=Open_File (argv[1],"r");

	fscanf(ini, "%d", &MaxNumInst);
	fscanf(ini, "%s", &outfile);


	LB = create_double_vector(Q);
	UB = create_double_vector(Q);
	TT = create_double_vector(Q);
	WT = create_double_vector(Q);

	ttt= time(NULL);
	tm = localtime(&ttt);
	Out = Open_File(outfile,"a+");	
	fprintf(Out, "\n %s\n", asctime(tm));
	fclose(Out);

	for (numSample=1 ; numSample<MaxNumInst ; numSample++){
		


		fscanf(ini,"%s",&instance);
		sprintf(path,"./New Data/");
		strcat(path,instance);

		H = 32;
		ReadData(path);
		
		tempUB = M;
		StartT = FinishT = 0;
		StartTimeCheck = FinishTimeCheck = 0;


		Out = Open_File(outfile,"a+");			
		fprintf(Out,"\n%s | %s ;\t", argv[1], instance);
		fclose(Out);



		StartT = clock();
		StartTimeCheck = clock();

		gettimeofday(&start, NULL);
		
		for (q=0 ; q<Q ; q++) {
			CPXENVptr env = NULL;
			CPXLPptr  lp = NULL;
			Creat_Problem();
			ModelPopulate();
			travelTime = 0;
			waitingTime = 0;
			waitingTime2 =0;
			UhjNewTemp1 = UhjNewTemp2 =  0;
			LB[q]= OBJ1;
			
			for (i = 0; i < I; i++) {
				for (jj = 0; jj <J; jj++) {
					for (k = 0; k <K; k++) {
						if (lambda[i] * t[i * J + jj] * x[(i * J + jj) + ( k * I * J )] > 0 )
							travelTime += lambda[i] * t[i * J + jj] * x[(i * J + jj) + ( k * I * J )];
					}
				}
			}
			for (j = 0; j <J; j++) {
				for (k = 0; k <K; k++) {
					temp1 = 0;
					for (i = 0; i < I; i++) {
						temp1 += lambda[i] * x[(i * J + j) + ( k * I * J )];
					}
					if ( mu[j * K + k] - temp1 != 0 && mu[j * K + k] !=0){
						waitingTime +=  ((1 + cv[j * K + k]*cv[j * K + k]) * (temp1 / (mu[j * K + k] - temp1))) + ((1 - cv[j * K + k]*cv[j * K + k]) * (temp1 / mu[j * K + k]));
						
					}
				}
			}
			tempUB = (at * travelTime) + (0.5 * aw * waitingTime);
			if (q>0) {
				if (tempUB > UB[q-1]){
					UB[q] = UB[q-1];
					TT[q] = TT[q-1];
					WT[q] = WT[q-1];
				} else{
					UB[q] = tempUB;
					TT[q] = travelTime;
					WT[q] = waitingTime;
				}
			} else {
				UB[q] = tempUB;
				TT[q] = travelTime;
				WT[q] = waitingTime;
			}
			Gap = (UB[q] - LB[q]) / LB[q];
			printf("\n|||||||||||q: %d | LB: %lf | Gap: %lf\n", q, LB[q], Gap);
					
			
			for (i = 0; i < I; i++) {
				for (j = 0; j <J; j++) {
					for (k = 0; k <K; k++) {
						if (x[(i * J + j) + ( k * I * J )] > 0) {
							
							printf("X %d [ %d ][ %d ] [ %d ] = %f\n", q, i, j, k, x[(i * J + j) + ( k * I * J )]);
						}
					}
				}
			}

			printf("\nY%d :\n", q); 
			for (j = 0; j <J; j++) {
				for (k = 0; k <K; k++) {
					if (y[j * K + k]>0) {
						
						printf("Y %d [ %d ][ %d ]= %f\n", q, j, k, y[j * K + k]);
					}
				}
			}
			printf("Lower Bound %d :  %f\n", q, LB[q]);
			printf("Upper Bound %d :  %f\n", q, UB[q]);
			printf("Gap %d :  %f\n", q, Gap);
			printf("epsilon :  %f\n", epsilon);
		
			printf("Number of H in iteration %d :  %d\n", q, H);


			if (Gap < epsilon){
				printf("______________vared shode ke dar Q:  %d \n", q);
				finishQ = q;
				break;
			}
			
			for (j = 0; j <J; j++) {
				UhjNewTemp1 = 0;
				tempmu = 0;
				for (k = 0; k <K; k++) {
					tempmu += mu[j * K + k] * y[j * K + k];
					for (i = 0; i < I; i++) {
						UhjNewTemp1 += lambda[i] * x[(i * J + j) + ( k * I * J )];
					}
				}
				if ((tempmu - UhjNewTemp1) > 0){
					Uhj[H] = UhjNewTemp1 / (tempmu - UhjNewTemp1);
					H++;
				}
			}
			finishQ = q;

			gettimeofday(&stop, NULL);

			checkTime = ((double) (stop.tv_sec - start.tv_sec) * 1000 + (double) (stop.tv_usec - start.tv_usec) / 1000)/1000;
			if (checkTime > limitationTime){
				printf("______________vared shode baraie zaman dar Q:  %d \n", q);
				finishQ = q;
				break;
			}
		}

		FinishT = clock();
		TotalT = (double)(FinishT - StartT)/ CLOCKS_PER_SEC;
		printf("\n\n\nQQQQQQQQQ:		%d\n\n\n", q);
		writeData();

		for ( q=0 ; q<Q ; q++) {
			UB[q] = 0;
			LB[q] = 0;
			TT[q] = WT[q] = 0;
		}
		for (hh = 32; hh <5000; hh++) {
			Uhj[hh] = 0 ;
		}
		H =32;

		free_memory();
	
				
	}

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
