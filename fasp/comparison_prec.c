/**
 *		Regression test for FASP: CSR
 *
 *------------------------------------------------------
 *
 *		Created  by Chensong Zhang on 03/20/2010
 *      Modified by Chensong Zhang on 09/02/2011
 *      Modified by Chenosng Zhang on 12/03/2011
 *      Modified by Chensong Zhang on 03/20/2012
 *
 *------------------------------------------------------
 *
 */

/*! \file regression.c
 *  \brief Regression testing for CSR iterative solvers
 */

#include <time.h>

#include "fasp.h"
#include "fasp_functs.h"

unsigned INT  ntest;    /**< number of tests all together */
unsigned INT  nfail;    /**< number of failed tests */

/**
 * \fn static void check_solu(dvector *x, dvector *sol, double tol)
 *
 * \brief This function compares x and sol to a given tolerance tol. 
 */
static void check_solu(dvector *x, dvector *sol, double tol)
{
    double diff_u = fasp_dvec_maxdiff(x, sol);
    ntest++;
    
    if ( diff_u < tol ) {
        printf("Max diff %.4e smaller than tolerance................. [PASS]\n", diff_u);
    }
    else {
        nfail++;
        printf("### WARNING: Max diff %.4e BIGGER than tolerance..... [ATTENTION!!!]\n", diff_u);
    }	
}

/**
 * \fn int main (int argc, const char * argv[])
 * 
 * \brief This is the main function for regression test. 
 *
 * \author Chensong Zhang
 * \date   03/20/2009
 * 
 * Modified by Chensong Zhang on 09/02/2011 
 * Modified by Chensong Zhang on 12/03/2011
 * Modified by Chensong Zhang on 03/20/2012
 */
int main (int argc, const char * argv[]) 
{

     /* Local Variables */
    dCSRmat        A;            // coefficient matrix
    dvector        b, x, sol;    // rhs, numerical sol, exact sol 
    INT            indp;         // index for test problems
    
	// Step 0. Set parameters
	input_param     inpar;  // parameters from input files
	itsolver_param  itpar;  // parameters for itsolver
	AMG_param       amgpar; // parameters for AMG
	ILU_param       ilupar; // parameters for ILU
   
    // Set solver parameters: use ./ini/bamg.dat
    fasp_param_set(argc, argv, &inpar);
    fasp_param_init(&inpar, &itpar, &amgpar, &ilupar, NULL);

    const int output_type   = inpar.output_type;   

    // Set output device
    if (output_type) {
		char *outputfile = "out/test.out";
		printf("Redirecting outputs to file: %s ...\n", outputfile);
		freopen(outputfile,"w",stdout); // open a file for stdout
	}


    
	// Step 1. Input stiffness matrix and right-hand side
	char filename1[512], *datafile1;
	char filename2[512], *datafile2;


	strcpy(inpar.workdir,"../../data/");
	
	strncpy(filename1,inpar.workdir,128);
	strncpy(filename2,inpar.workdir,128);
    
   
    int problem_num   = inpar.problem_num;    

    do {

    printf("Test Problem %d\n", problem_num);
    // Default test problem from black-oil benchmark: SPE01
	if (problem_num == 10) {				
        // Read the stiffness matrix from bsrmat_SPE01.dat
        datafile1="FDM/fdm_mat_csr_127X127X127.dat"; 
        //datafile1="FDM/fdm_mat_csr_1023X1023.dat"; 
        strcat(filename1,datafile1);
        
        // Read the RHS from rhs_SPE01.dat
        datafile2="FDM/fdm_rhs_127X127X127.dat"; 
        //datafile2="FDM/fdm_rhs_1023X1023.dat"; 
        strcat(filename2,datafile2);

        fasp_dcsrvec2_read(filename1, filename2, &A, &b);
    }
    else if (problem_num == 20){
		
		datafile1="spe10-uncong/SPE1020.amg";
	    strncpy(filename1,inpar.workdir,128);
		strcat(filename1,datafile1);
	
		datafile2="spe10-uncong/SPE1020.rhs";
	    strncpy(filename2,inpar.workdir,128);
		strcat(filename2,datafile2);
		
		fasp_matrix_read_bin(filename1, &A);
		
		FILE *fp = fopen(filename2,"r");
		INT  i, n;
		REAL value;
		fscanf(fp,"%d",&n);
		fasp_dvec_alloc(n,&b);
		for (i=0;i<n;++i) {
			fscanf(fp, "%le", &value);
			b.val[i]=value;
		}
		fclose(fp);

        //fasp_dcsrvec2_write ("../../data/SPE1020.amg.dat", "../../data/SPE1020.rhs.dat", &A, &b);                       
	}
	else if (problem_num == 40){
		
		datafile1="spe10-uncong/SPE1040.amg";
	    strncpy(filename1,inpar.workdir,128);
		strcat(filename1,datafile1);
	
		datafile2="spe10-uncong/SPE1040.rhs";
	    strncpy(filename2,inpar.workdir,128);
		strcat(filename2,datafile2);
		
		fasp_matrix_read_bin(filename1, &A);
		//fasp_matrix_write("../data/spe10-data/SPE1040.amg", &A, 001);
		FILE *fp = fopen(filename2,"r");
		INT  i, n;
		REAL value;
		fscanf(fp,"%d",&n);
		fasp_dvec_alloc(n,&b);
		for (i=0;i<n;++i) {
			fscanf(fp, "%le", &value);
			b.val[i]=value;
		}
		fclose(fp);
        fasp_dcsrvec2_write ("../../data/SPE1040.amg.dat", "../../data/SPE1040.rhs.dat", &A, &b);                       
	}
	else {
		printf("### ERROR: Unrecognized problem number %d\n", problem_num);
		return ERROR_INPUT_PAR;
	}
	
        
    /************************************/
    /* Step 2. Check matrix properties  */
    /************************************/
    fasp_check_symm(&A);     // check symmetry
    fasp_check_diagpos(&A);  // check sign of diagonal entries
    fasp_check_diagdom(&A);  // check diagonal dominance

    /*****************************/
    /* Step 3. Solve the system  */ 
    /*****************************/

    // allocate mem for numerical solution
    fasp_dvec_alloc(b.row, &x);  
#if 1
    /* AMG V-cycle (Direct interpolation) with GS smoother as a solver */			
    printf("------------------------------------------------------------------\n");
    printf("Classical(RS) AMG (direct interp) V-cycle as preconditioned iterative solver ...\n");	
            
    // reset initial guess
    fasp_dvec_set(b.row, &x, 0.0); 

    amgpar.interpolation_type = INTERP_DIR;
	amgpar.coarsening_type = COARSE_RS;

    //fasp_solver_amg(&A, &b, &x, &amgpar);
	
	itpar.itsolver_type = 4;
	if (problem_num == 10) itpar.itsolver_type = 1;
    fasp_solver_dcsr_krylov_amg (&A, &b, &x, &itpar, &amgpar);
            
    /* AMG V-cycle (Standard interpolation) with GS smoother as a solver */			
    printf("------------------------------------------------------------------\n");
    printf("Classical(RS) AMG (standard interp) V-cycle as preconditioned iterative solver ...\n");	
            
    fasp_dvec_set(b.row, &x, 0.0); // reset initial guess

    amgpar.interpolation_type = INTERP_STD;
	amgpar.coarsening_type = COARSE_RS;

    //fasp_solver_amg(&A, &b, &x, &amgpar);
	itpar.itsolver_type = 4;
	if (problem_num == 10) itpar.itsolver_type = 1;
    fasp_solver_dcsr_krylov_amg (&A, &b, &x, &itpar, &amgpar);

    /* AMG V-cycle (Direct interpolation) with GS smoother as a solver */			
    printf("------------------------------------------------------------------\n");
    printf("Classic(AC1) AMG (direct interp) V-cycle as preconditioned iterative solver ...\n");	
            
    // reset initial guess
    fasp_dvec_set(b.row, &x, 0.0); 

    amgpar.interpolation_type = INTERP_DIR;
	amgpar.coarsening_type = COARSE_AC;
    amgpar.aggressive_path = 1;

    //fasp_solver_amg(&A, &b, &x, &amgpar);
	itpar.itsolver_type = 4;
	if (problem_num == 10) itpar.itsolver_type = 1;
    fasp_solver_dcsr_krylov_amg (&A, &b, &x, &itpar, &amgpar);
            
    /* AMG V-cycle (Standard interpolation) with GS smoother as a solver */			
    printf("------------------------------------------------------------------\n");
    printf("Classic(AC1) AMG (standard interp) V-cycle as preconditioned iterative solver ...\n");	
            
    fasp_dvec_set(b.row, &x, 0.0); // reset initial guess

    amgpar.interpolation_type = INTERP_STD;
	amgpar.coarsening_type = COARSE_AC;
    amgpar.aggressive_path = 1;

    //fasp_solver_amg(&A, &b, &x, &amgpar);
	itpar.itsolver_type = 4;
	if (problem_num == 10) itpar.itsolver_type = 1;
    fasp_solver_dcsr_krylov_amg (&A, &b, &x, &itpar, &amgpar);

    /* AMG V-cycle (Direct interpolation) with GS smoother as a solver */			
    printf("------------------------------------------------------------------\n");
    printf("Classic(AC2) AMG (direct interp) V-cycle as preconditioned iterative solver ...\n");	
            
    // reset initial guess
    fasp_dvec_set(b.row, &x, 0.0); 

    amgpar.interpolation_type = INTERP_DIR;
	amgpar.coarsening_type = COARSE_AC;
    amgpar.aggressive_path = 2;

    //fasp_solver_amg(&A, &b, &x, &amgpar);
	itpar.itsolver_type = 4;
	if (problem_num == 10) itpar.itsolver_type = 1;
    fasp_solver_dcsr_krylov_amg (&A, &b, &x, &itpar, &amgpar);
            
    /* AMG V-cycle (Standard interpolation) with GS smoother as a solver */			
    printf("------------------------------------------------------------------\n");
    printf("Classic(AC2) AMG (standard interp) V-cycle as preconditioned iterative solver ...\n");	
            
    fasp_dvec_set(b.row, &x, 0.0); // reset initial guess

    amgpar.interpolation_type = INTERP_STD;
	amgpar.coarsening_type = COARSE_AC;
    amgpar.aggressive_path = 2;

    //fasp_solver_amg(&A, &b, &x, &amgpar);
	itpar.itsolver_type = 4;
	if (problem_num == 10) itpar.itsolver_type = 1;
    fasp_solver_dcsr_krylov_amg (&A, &b, &x, &itpar, &amgpar);

#else	
    /* AMG V-cycle (Standard interpolation) with GS smoother as a solver */			
    printf("------------------------------------------------------------------\n");
    printf("UA NA-cycle as preconditioned iterative solver ...\n");	
            
    fasp_dvec_set(b.row, &x, 0.0); // reset initial guess

	if(problem_num == 10) {

    amgpar.AMG_type    = UA_AMG;
	amgpar.cycle_type  = NL_AMLI_CYCLE;

    //fasp_solver_amg(&A, &b, &x, &amgpar);
	itpar.itsolver_type = 4;
	if (problem_num == 10) itpar.itsolver_type = 1;
    fasp_solver_dcsr_krylov_amg (&A, &b, &x, &itpar, &amgpar);

	}

    /* AMG V-cycle (Standard interpolation) with GS smoother as a solver */			
    printf("------------------------------------------------------------------\n");
    printf("SA AMG V-cycle as preconditioned iterative solver ...\n");	
            
    fasp_dvec_set(b.row, &x, 0.0); // reset initial guess
    
    amgpar.AMG_type    = SA_AMG;
	amgpar.cycle_type  = V_CYCLE;

    //fasp_solver_amg(&A, &b, &x, &amgpar);
	itpar.itsolver_type = 4;
	if (problem_num == 10) itpar.itsolver_type = 1;
    fasp_solver_dcsr_krylov_amg (&A, &b, &x, &itpar, &amgpar);
         
#endif
    /* clean up memory */
    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&x);

    problem_num = 2*problem_num;

	} while(problem_num < 80);
	
	return SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
