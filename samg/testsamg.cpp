/***********************************************************/
/* Demo C++ driver which reads demo.* files and calls SAMG */
/*                      Modified by ZhengLi on 2.27.2014   */
/***********************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>

#include "fasp.h"
#include "fasp_functs.h"


// the old C++ compiler on Alpha cannot handle using namespace std !
#ifndef SAMG_OLDALPHA
using namespace std;
#endif

#include "samg.h"
//  USE THE HEADER FILE SAMG.H TO DECLARE THE AMG SUBROUTINE SAMG.
//  ADDITIONALLY, SAMG.H DEMONSTRATES HOW TO DEFINE INTERFACES TO ALL
//  ACCESS ROUTINES FOR PARAMETERS WHICH CONTROL THE BEHAVIOUR OF SAMG. 
//  THIS WAY, THE CALLER FROM C(++) CAN CONTROL ALL SAMG FEATURES
//  DESCRIBED IN THE MANUAL!

int main()
{
    cout << " *** Demo C++ driver which reads demo.* files and calls SAMG ***" << endl;

    // Set primary parameters. Others can be set by access functions as shown below.

    int    ndiu      = 1;        // dimension of (dummy) vector iu
    int    ndip      = 1;        // dimension of (dummy) vector ip
    int    nsolve    = 10;       // results in scalar approach (current system is scalar)
	                             // Plain relaxation, variable©\wise
    int    ifirst    = 1;        // first approximation = zero
    double eps       = 1.0e-6;   // required (relative) residual reduction 
    //int    ncyc      = 11050;   // V-cycle as pre-conditioner for CG; at most 50 iterations
	int    ncyc      = 109200;   // V-cycle as pre-conditioner for CG; at most 200iterations

    double a_cmplx   = 2.2;      // estimated dimensioning define the dimensions used to allocate SAMG¡¯s initial memory
    double g_cmplx   = 1.7;      // estimated dimensioning
    double w_avrge   = 2.4;      // estimated dimensioning
    double p_cmplx   = 0.0;      // estimated dimensioning (irrelevant for scalar case)

    double chktol    = -1.0;     // input checking de-activated (we know it's ok!)
    int    idump     = 0;        // minimum output during setup
    int    iout      = 2;        // display residuals per iteration and work statistics

    int    n_default = 0;        // select default settings for secondary parameters
                                 // CURRENTLY AVAILABLE: 10-13, 15-18, 20-23, 25-28
                                 // NOTE: the higher the respective SECOND digit, the
                                 // more aggressive the coarsening (--> lower memory at
                                 // the expense of slower convergence)
    int iswtch = 5100+n_default; // complete SAMG run ....
                                 // ... memory de-allocation upon return ....
                                 // ... memory extension feature activated ....
                                 // ... residuals measured in the L2-norm .....
                                 // ... secondary parameter default setting # n_default

    // Secondary parameters which have to be set if n_default=0
    //(at the same time a demonstration of how to access secondary or hidden parameters)

    int intin;
    double dblin;
 
    // NCGRAD_DEFAULT=NCGRAD=0 AMG as iteration methods  1:CG 2,BICGstab 3,FGMRES 4 FCG
	intin=3;     SAMG_SET_NCGRAD_DEFAULT(&intin);


    if (n_default == 0) {
         intin=25;    SAMG_SET_LEVELX(&intin);

         intin=1;    SAMG_SET_MODE_MESS(&intin);

         intin=100;   SAMG_SET_NPTMN(&intin); // coarsest dof 

         intin=10;     SAMG_SET_NCG(&intin);  // Specifies the treatment of positive matrix 
		                                      // entries in the coarsening process.
		                                      // ncg[2] 0:standrad coarsening 1-4:aggresive coarsening

         intin=2;     SAMG_SET_NWT(&intin);   // types of interpolation:1 direc 2 standard 4 multi-pass
         intin=1;     SAMG_SET_NTR(&intin);   // Type of truncation 1: absolute value
         intin=131;   SAMG_SET_NRD(&intin);   // pre-smoothing steps
         intin=131;   SAMG_SET_NRU(&intin);   // post-smoothing  steps

         intin=6;     SAMG_SET_NRC(&intin);   // coarsest-level solver 6:FGE,7:SGE
 

		 intin=0;     SAMG_SET_NP_OPT(&intin); // selects the level of optimization

         dblin=21.25; SAMG_SET_ECG(&dblin);  // Two threshold for strong connectivity 
		                                     // and strong (local)diagonal dominance.

         dblin=0.20;  SAMG_SET_EWT(&dblin);  // Two threshold for strong violation
                                             // of diagonal dominance, 
		                                     // and large positive couplings which should 
		                                     //not be ignored in coarsening.

         dblin=12.20; SAMG_SET_ETR(&dblin);  // Two truncation parameters for the basis interpolation: default
    }

    // amg declarations 

    int npnt,nsys,matrix,nnu,nna;

    int * ia, * ja;
    int * iu     = new int[1];
    int * ip     = new int[1];
    int * iscale = new int[1];

	iscale[0] = 0;

    double * a, * u, * f;

    // output:
    int ierr,ncyc_done;
    double res_out,res_in;

	dCSRmat A;
	dvector b, sol;

	// Read A and b from two files in IJ format. 
    
	//fasp_dcsrvec2_read("../../data/fdm_mat_csr_127X127X127.dat", "../../data/fdm_rhs_127X127X127.dat", &A, &b);
    fasp_dcsrvec2_read("../../data/SPE1020.amg.dat", "../../data/SPE1020.rhs.dat", &A, &b);
	//fasp_dcsrvec2_read("../../data/SPE1040.amg.dat", "../../data/SPE1040.amg.dat", &A, &b);
	
	// Read ref. sol. from an indexed vec file.
    dvector2SAMGInput(&b,"demo.rhs");

    dCSRmat2SAMGInput(&A,"demo.frm", "demo.amg");

    //Read data from stdin

    ifstream frmfile("demo.frm", ios::in);    // open demo.frm

    int iversion;
    char ch;

    frmfile >> ch >> iversion;
    frmfile >> nna >> nnu >> matrix >> nsys >> npnt;

    frmfile.close();

    if (ch!='f' || iversion != 4) {
        cout << "invalid file format for this test driver " << endl;
        ierr=1; return ierr;
    }
      
    // matrix
    int i;

    ia = new int[nnu+1]; 
    ja = new int[nna];
    a  = new double[nna];

    if (!(ia && ja && a)) {
        cout << " allocation failed (ia,ja,a) " << endl;
        ierr=1; return ierr;
     }

     // open demo.amg
     ifstream amgfile("demo.amg", ios::in);    
     for (i=0;i<nnu+1;i++) amgfile >> ia[i];
     for (i=0;i<nna;i++)   amgfile >> ja[i]; 
     for (i=0;i<nna;i++)   amgfile >> a[i]; 

     amgfile.close();

     // right hand side
     f = new double[nnu]; 
	 u = new double[nnu]; 

     if (!(f && u)) {
         cout << " allocation failed (f,u) " << endl;
         ierr=1; return ierr;
     }
     for (i=0;i<nnu;i++) u[i]=0.0;

     // open demo.rhs
     ifstream rhsfile("demo.rhs", ios::in);    

     for (i=0;i<nnu;i++) rhsfile >> f[i];

     rhsfile.close();

     float told,tnew,tamg;
     SAMG_CTIME(&told);

     SAMG(&nnu,&nna,&nsys,
           &ia[0],&ja[0],&a[0],&f[0],&u[0],&iu[0],&ndiu,&ip[0],&ndip,&matrix,&iscale[0],
           &res_in,&res_out,&ncyc_done,&ierr,
           &nsolve,&ifirst,&eps,&ncyc,&iswtch,
           &a_cmplx,&g_cmplx,&p_cmplx,&w_avrge,
           &chktol,&idump,&iout);

     if (ierr > 0) {
          cout << endl << " SAMG terminated with error code " 
               << ierr << " **** " << endl;
     }
      else if (ierr < 0) {
          cout << endl << " SAMG terminated with warning code " 
               << ierr << " **** " << endl;
      }

      SAMG_CTIME(&tnew);
      tamg=tnew-told; 
      cout << endl << " ***** total run time: " << tamg << " ***** " << endl; 

      delete[] ia,ja,a,f,u,iu,ip,iscale;

      return ierr;
}
