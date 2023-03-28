// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// GalBar.cc                                                                   |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// defines Dehnen & Binney's (1998) Galaxy potential as NEMO potential         |
// with a bar                                                                  |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// The units differ from those used in class GalPot: we have                   |
//                                                                             |
//  quantity              | here          | GalPot        | unit               |
// -----------------------+---------------+---------------+--------------      |
// unit of length         | 1             | 1             | kpc                |
// unit of time           | 1 E+9         | 1 E+6         | yr                 |
// unit of mass           | 2.2228847 E+5 | 1             | solar masses       |
// unit of velocity       | 0.9777753     | 9.777753  E+2 | km/s               |
// unit of potential      | 0.95604457    | 9.5604457 E+5 | (km/s)^2           |
// unit of acceleration   | 0.9777753 E-9 | 0.977753  E-6 | km/s/yr            |
//                                                                             |
// The implied size of Newton's constant G are                                 |
// G = 1                 here                                                  |
// G = 4.49865897 E-12   in GalPot units                                       |
//                                                                             |
//-----------------------------------------------------------------------------+
//                                                                             |
// Versions                                                                    |
// 0.0   08-jul-2020    created                                            GT  |
//-----------------------------------------------------------------------------+

// Compilation
//g++ -o acc/GalBar.so src/public/acc/GalBar.cc -Iinc/ -Iinc/utils/ -I/Users/guillaume/nemo/inc -I/Users/guillaume/nemo/lib -D_FILE_OFFSET_BITS=64  -mfpmath=sse -mpreferred-stack-boundary=4 --param inline-unit-growth=50 -ggdb3 -Wall -Wextra -Winit-self  -Wshadow -Wno-format-security -O2 -fPIC -funroll-loops -fforce-addr   -Woverloaded-virtual   -Llib/ -lfalcON -L../utils/lib -lWDutils -L/Users/guillaume/nemo/lib -lnemo -ldl -shared -w; cp acc/GalBar.so /Users/guillaume/nemo/obj/acc/.
// mknemo gyrfalcON

//g++ -o acc/GalBar.so src/public/acc/GalBar.cc -Iinc/ -Iinc/utils/ -I/home/gthomas/nemo/inc -I/home/gthomas/nemo/lib -D_FILE_OFFSET_BITS=64 -std=c++03  -mfpmath=sse -mpreferred-stack-boundary=4 --param inline-unit-growth=50 -ggdb3 -Wall -Wextra -Winit-self  -Wshadow -Wno-format-security -Wno-misleading-indentation -O2 -fPIC -funroll-loops -fforce-addr   -rdynamic -march=native -Woverloaded-virtual   -Llib/ -lfalcON -L../utils/lib -lWDutils -L/home/gthomas/nemo/lib -lnemo -ldl -shared -w; cp acc/GalBar.so /home/gthomas/nemo/obj/acc/.

//-----------------------------------------------------------------------------+

#include <iostream>
#include <fstream>
#define POT_DEF
#include <defacc.h>
////////////////////////////////////////////////////////////////////////////////
using namespace std;
#include "utils/exception.h"   // so that tupel.h uses static assertion
#include "utils/spline.h"
#include "acc/GalPot.cc"
using namespace WDutils; // To use the spline

////////////////////////////////////////////////////////////////////////////////
namespace {

  using GalPot::GalaxyPotential;

  class GalaxyFile {
  protected:
    ifstream from;
    GalaxyFile(const char*file)
    {
      if(file==0 || file[0]==0) 
	::error("Need data file to initialize GalPot");
      from.open(file);
      if(!from.is_open())
	::error("GalPot: cannot open file \"%s\"",file);
      nemo_dprintf(4,"file \"%s\" opened\n",file);
    }
  };
  //////////////////////////////////////////////////////////////////////////////
class Bar{
	public:
	
	char fileAmpl[256];
	bool isBar=false;
	double omega_b; // Pattern speed
	double theta0; // Angle of origin of the bar compare to the X axis
	long nb_lines=0;
	double *rs;
	double *A20,*A40,*A60;
	double *zh2,*zh4,*zh6;
	double *phase2,*phase4,*phase6;
	double *derive_A20,*derive_A40,*derive_A60;
	double *derive_zh2,*derive_zh4,*derive_zh6;
	double *derive_phase2,*derive_phase4,*derive_phase6;
	cubic_splines spline_A20,spline_A40,spline_A60;
	cubic_splines spline_zh2,spline_zh4,spline_zh6;
	cubic_splines spline_phase2,spline_phase4,spline_phase6;
		
	//Method to get than initial angle of the bar and its pattern speed
	void loadBar(const char*namefile, int nb_disk, int nb_spher){
		long i,nline=0;
		char line[10000];

		// Get number line of the file containing the potential parameters
		FILE* file = NULL;
		file=fopen(namefile,"r");
		if(file!=NULL){
			while(fgets(line,10000,file)!=NULL){
				nline++;
			}
			fclose(file);
		}else{::error("Need data file to initialize GalBar"); return;}
		
		if(nline>nb_disk+nb_spher){ // if parameters of the bar
			file=fopen(namefile,"r");
			for(i=0;i<nb_disk+nb_spher+2;i++){fgets(line,10000,file);} // Skip line where are define the differents componants

			fscanf(file,"%lf %lf %s\n",&omega_b,&theta0,fileAmpl); // km/s/kpc, degree, file for spline
			fclose(file);
			omega_b/=0.9777753; //km/s/kpc-> gyr-1
			theta0=theta0/180.0*M_PI; // deg to rad
			isBar=true;
		}
	}
	
	// Method load the amplitude and phase file
	void loadFile(){  // Method/function defined inside the class
	
		// Compte the number of lines
		long i;
		char line[10000];
		FILE* file = NULL;
		file=fopen(fileAmpl,"r");
		nb_lines=0;
		if(file!=NULL){
			while(fgets(line,10000,file)!=NULL){
				nb_lines++;
			}
			fclose(file);
		}else{::error("No file containing the Amplitude and phase of the modes of the bar"); return;}
		nb_lines-=1; // remove the line of the header for the count
		
		// Load the file
		rs=(double*)malloc(nb_lines*sizeof(double));
		A20=(double*)malloc(nb_lines*sizeof(double));
		A40=(double*)malloc(nb_lines*sizeof(double));
		A60=(double*)malloc(nb_lines*sizeof(double));
		zh2=(double*)malloc(nb_lines*sizeof(double));
		zh4=(double*)malloc(nb_lines*sizeof(double));
		zh6=(double*)malloc(nb_lines*sizeof(double));
		phase2=(double*)malloc(nb_lines*sizeof(double));
		phase4=(double*)malloc(nb_lines*sizeof(double));
		phase6=(double*)malloc(nb_lines*sizeof(double));

		file=fopen(fileAmpl,"r");
		if(file!=NULL){
			fgets(line,10000,file); // Skip first line because is header
			for(i=0;i<(nb_lines);i++){
		fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",&rs[i],&A20[i],&zh2[i],&phase2[i],&A40[i],&zh4[i],&phase4[i],&A60[i],&zh6[i],&phase6[i]);
			}
		}
		fclose(file);
		return;
	}
	
	
	// Method to make the spline on the different variables
	void getSpline(void){
		// Amplitude
		derive_A20=(double*)malloc(nb_lines*sizeof(double));
		spline_A20.construct(nb_lines, rs, A20,derive_A20);
		derive_A40=(double*)malloc(nb_lines*sizeof(double));
		spline_A40.construct(nb_lines, rs, A40,derive_A40);
		derive_A60=(double*)malloc(nb_lines*sizeof(double));
		spline_A60.construct(nb_lines, rs, A60,derive_A60);

		// Scale height
		derive_zh2=(double*)malloc(nb_lines*sizeof(double));
		spline_zh2.construct(nb_lines, rs, zh2,derive_zh2);
		derive_zh4=(double*)malloc(nb_lines*sizeof(double));
		spline_zh4.construct(nb_lines, rs, zh4,derive_zh4);
		derive_zh6=(double*)malloc(nb_lines*sizeof(double));
		spline_zh6.construct(nb_lines, rs, zh6,derive_zh6);
		
		// Phase
		derive_phase2=(double*)malloc(nb_lines*sizeof(double));
		spline_phase2.construct(nb_lines, rs, phase2,derive_phase2);
		derive_phase4=(double*)malloc(nb_lines*sizeof(double));
		spline_phase4.construct(nb_lines, rs, phase4,derive_phase4);
		derive_phase6=(double*)malloc(nb_lines*sizeof(double));
		spline_phase6.construct(nb_lines, rs, phase6,derive_phase6);
	}
	
	
	// Method to get the amplitude of at the position
	void getAmp(double R,double *A20_r,double *A40_r,double *A60_r,double *zh2_r,double *zh4_r,double *zh6_r,double *phase2_r,double *phase4_r,double *phase6_r){
		long idxm,idxp;
		
		if(R<=rs[nb_lines-2]){
			// get the index of the radius boundaries around the particle; -2 to avoid boundary problems
			idxp=0;
			while(R>=rs[idxp]) idxp++;
			idxm=idxp-2;
			if(R==rs[idxm]){ // if at the position of one of the point of the core point of the spline
				*A20_r=A20[idxm];
				*A40_r=A40[idxm];
				*A60_r=A60[idxm];
				*zh2_r=zh2[idxm];
				*zh4_r=zh4[idxm];
				*zh6_r=zh6[idxm];
				*phase2_r=phase2[idxm];
				*phase4_r=phase4[idxm];
				*phase6_r=phase6[idxm];
				
			}else{
			spline_A20.evaluate(R,rs[idxm],rs[idxp],A20[idxm],A20[idxp],derive_A20[idxm],derive_A20[idxp],A20_r);
			spline_A40.evaluate(R,rs[idxm],rs[idxp],A40[idxm],A40[idxp],derive_A40[idxm],derive_A40[idxp],A40_r);
			spline_A60.evaluate(R,rs[idxm],rs[idxp],A60[idxm],A60[idxp],derive_A60[idxm],derive_A60[idxp],A60_r);
			spline_zh2.evaluate(R,rs[idxm],rs[idxp],zh2[idxm],zh2[idxp],derive_zh2[idxm],derive_zh2[idxp],zh2_r);
			spline_zh4.evaluate(R,rs[idxm],rs[idxp],zh4[idxm],zh4[idxp],derive_zh4[idxm],derive_zh4[idxp],zh4_r);
			spline_zh6.evaluate(R,rs[idxm],rs[idxp],zh6[idxm],zh6[idxp],derive_zh6[idxm],derive_zh6[idxp],zh6_r);			
			spline_phase2.evaluate(R,rs[idxm],rs[idxp],phase2[idxm],phase2[idxp],derive_phase2[idxm],derive_phase2[idxp],phase2_r);
			spline_phase4.evaluate(R,rs[idxm],rs[idxp],phase4[idxm],phase4[idxp],derive_phase4[idxm],derive_phase4[idxp],phase4_r);
			spline_phase6.evaluate(R,rs[idxm],rs[idxp],phase6[idxm],phase6[idxp],derive_phase6[idxm],derive_phase6[idxp],phase6_r);			
			}
		}else{
			*A20_r=0.0;
			*A40_r=0.0;
			*A60_r=0.0;
			*zh2_r=1.0;
			*zh4_r=1.0;
			*zh6_r=1.0;
			*phase2_r=0.0;
			*phase4_r=0.0;
			*phase6_r=0.0;
		}
	}
	
	// Method to check the spline interpolation
	void print(){
		double R,A20,A40,A60,zh2,zh4,zh6,phase2,phase4,phase6;
		long i;
		FILE* file = NULL;
		double inter=25.0/100000.0;

		file=fopen("spline.asc","w");
		for(i=0;i<100000.0;i++){
			R=i*inter;
			getAmp( R,&A20,&A40,&A60,&zh2,&zh4,&zh6,&phase2,&phase4,&phase6);
			fprintf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",R,A20,A40,A60,zh2,zh4,zh6,phase2,phase4,phase6);
		}
		fclose(file);
		
	}

	
};


double amp_bar(double Am0, double zhm, int m, double theta,double theta0,double t,double omega_b,double R,double Z){
	return Am0*cos(m*(theta+theta0+t*omega_b))/(1.0+pow(Z/zhm,2)); //Note here that the omega_b is positive because we are in simulation time where tcurrent = tend *(R*R/(R*R+Z*Z));
}

Bar bar;
double tbar;


  //////////////////////////////////////////////////////////////////////////////
  class GalaxyBar :
    private GalaxyFile,
    private GalaxyPotential
  {
  public:
    //--------------------------------------------------------------------------
    static const char* name() { return "GalBar"; }
    bool NeedMass() const { return false; }
    bool NeedVels() const { return false; }
    //--------------------------------------------------------------------------
    GalaxyBar(const double*, int, const char*file)
      : GalaxyFile     ( file ),
	GalaxyPotential( from )
    {
      if(nemo_debug(2) )
	std::cerr<<
	  " falcON debug info:\n"
	  " external potential \"GalPot\" requires data file in GalPot format.\n";
		// Load the Bar parameters
		bar.loadBar(file,GalaxyPotential::NumberofDisks(),GalaxyPotential::NumberofSpheroids());
		if(bar.isBar==true){
			// Load amplitude file and made the spline model
			bar.loadFile();
			bar.getSpline();
			bar.print(); // save the spline data in a file
		}
		from.close();
    }
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void set_time(double    s ,
		  int          ,
		  const scalar*,
		  const scalar*,
						const scalar*) const {tbar=s;}
    //--------------------------------------------------------------------------
    template<int NDIM, typename scalar>
    void acc(const scalar*,
	     const scalar*X,
	     const scalar*,
	     scalar      &P,
	     scalar      *A) const{
		 
			double delta=1.0e-7; // Step for compute acceleration (eqivalent to 20 UA)
			double P1,P2,P0;
			double fR1,fz1,R1;
			double fR2,fz2,R2;
			double theta,theta1,theta2;
		 //cout<<"Time "<<tbar<<"\n";
			// Compute axisymetric forces from GalPot
			double fR,fz,R=hypot(X[0],X[1]);
			theta=atan2(X[1],X[0]);
			double A20,A40,A60;
			double A20_1,A40_1,A60_1;
			double A20_2,A40_2,A60_2;
			double zh2,zh4,zh6;
			double zh2_1,zh4_1,zh6_1;
			double zh2_2,zh4_2,zh6_2;
			double phase2,phase4,phase6;
			double phase2_1,phase4_1,phase6_1;
			double phase2_2,phase4_2,phase6_2;

			bar.getAmp(R,&A20,&A40,&A60,&zh2,&zh4,&zh6,&phase2,&phase4,&phase6); // Get the amplitude of the Bar ¡
			//printf("Amp %lf %lf %lf %lf %lf \n",R,A20,A40,A60,theta/M_PI*180.0);
	
			if(NDIM > 2) {
				// Axisymetric GalPot potential at the position
				P    = 1.e6 * GalaxyPotential::operator()(R, X[2], fR, fz);
				if(bar.isBar==true){
//printf("A\n");
					P=P*(1.0+amp_bar(A20,zh2,2,theta,bar.theta0,tbar,bar.omega_b,R,X[2])+amp_bar(A40,zh4,4,theta,bar.theta0,tbar,bar.omega_b,R,X[2])+amp_bar(A60,zh6,6,theta,bar.theta0,tbar,bar.omega_b,R,X[2]));
				}
				if(R) {
					// Compute Ax
					fR1,fz1,R1=hypot(X[0]-delta,X[1]);
					fR2,fz2,R2=hypot(X[0]+delta,X[1]);
					theta1=atan2(X[1],X[0]-delta);
					theta2=atan2(X[1],X[0]+delta);
					bar.getAmp(sqrt(R1*R1),&A20_1,&A40_1,&A60_1,&zh2_1,&zh4_1,&zh6_1,&phase2_1,&phase4_1,&phase6_1); // Get the amplitude of
					bar.getAmp(sqrt(R2*R2),&A20_2,&A40_2,&A60_2,&zh2_2,&zh4_2,&zh6_2,&phase2_2,&phase4_2,&phase6_2); // Get the amplitude of
					P1    = 1.e6 * GalaxyPotential::operator()(R1, X[2], fR1, fz1);
					P2    = 1.e6 * GalaxyPotential::operator()(R2, X[2], fR2, fz2);
					if(bar.isBar==true){
					P1=P1*(1.0+amp_bar(A20_1,zh2_1,2,theta1,bar.theta0,tbar,bar.omega_b,R1,X[2])+amp_bar(A40_1,zh4_1,4,theta1,bar.theta0,tbar,bar.omega_b,R1,X[2])+amp_bar(A60_1,zh6_1,6,theta1,bar.theta0,tbar,bar.omega_b,R1,X[2]));
					P2=P2*(1.0+amp_bar(A20_2,zh2_2,2,theta2,bar.theta0,tbar,bar.omega_b,R2,X[2])+amp_bar(A40_2,zh4_2,4,theta2,bar.theta0,tbar,bar.omega_b,R2,X[2])+amp_bar(A60_2,zh6_2,6,theta2,bar.theta0,tbar,bar.omega_b,R2,X[2]));
					}
					A[0]=-(P2-P1)/(2.0*delta);
					
					// Compute Ay
					fR1,fz1,R1=hypot(X[0],X[1]-delta);
					fR2,fz2,R2=hypot(X[0],X[1]+delta);
					theta1=atan2(X[1]-delta,X[0]);
					theta2=atan2(X[1]+delta,X[0]);
					bar.getAmp(sqrt(R1*R1),&A20_1,&A40_1,&A60_1,&zh2_1,&zh4_1,&zh6_1,&phase2_1,&phase4_1,&phase6_1); // Get the amplitude of
					bar.getAmp(sqrt(R2*R2),&A20_2,&A40_2,&A60_2,&zh2_2,&zh4_2,&zh6_2,&phase2_2,&phase4_2,&phase6_2); // Get the amplitude of

					P1    = 1.e6 * GalaxyPotential::operator()(R1, X[2], fR1, fz1);
					P2    = 1.e6 * GalaxyPotential::operator()(R2, X[2], fR2, fz2);
					if(bar.isBar==true){
					P1=P1*(1.0+amp_bar(A20_1,zh2_1,2,theta1,bar.theta0,tbar,bar.omega_b,R1,X[2])+amp_bar(A40_1,zh4_1,4,theta1,bar.theta0,tbar,bar.omega_b,R1,X[2])+amp_bar(A60_1,zh6_1,6,theta1,bar.theta0,tbar,bar.omega_b,R1,X[2]));
					P2=P2*(1.0+amp_bar(A20_2,zh2_2,2,theta2,bar.theta0,tbar,bar.omega_b,R2,X[2])+amp_bar(A40_2,zh4_2,4,theta2,bar.theta0,tbar,bar.omega_b,R2,X[2])+amp_bar(A60_2,zh6_2,6,theta2,bar.theta0,tbar,bar.omega_b,R2,X[2]));
					}
					A[1]=-(P2-P1)/(2.0*delta);

				} else {
				  A[0] = 0.;
				  A[1] = 0.;
				}
				// Compute Az
				P1    = 1.e6 * GalaxyPotential::operator()(R, X[2]-delta, fR, fz);
				P2    = 1.e6 * GalaxyPotential::operator()(R, X[2]+delta, fR, fz);
				if(bar.isBar==true){
					P1=P1*(1.0+amp_bar(A20,zh2,2,theta,bar.theta0,tbar,bar.omega_b,R,X[2]-delta)+amp_bar(A40,zh4,4,theta,bar.theta0,tbar,bar.omega_b,R,X[2]-delta)+amp_bar(A60,zh6,6,theta,bar.theta0,tbar,bar.omega_b,R,X[2]-delta));
					P2=P2*(1.0+amp_bar(A20,zh2,2,theta,bar.theta0,tbar,bar.omega_b,R,X[2]+delta)+amp_bar(A40,zh4,4,theta,bar.theta0,tbar,bar.omega_b,R,X[2]+delta)+amp_bar(A60,zh6,6,theta,bar.theta0,tbar,bar.omega_b,R,X[2]+delta));
				}
				A[2]=-(P2-P1)/(2.0*delta);
				
				fR  *=-1.e6/R;
				
			}else{
				P    = 1.e6 * GalaxyPotential::operator()(R, 0.,   fR, fz);
				if(R) {
				  fR  *=-1.e6/R;
				  A[0] = fR * X[0];
				  A[1] = fR * X[1];
				} else {
				  A[0] = 0.;
				  A[1] = 0.;
				}
			}
    }
  };
} // namespace {
//------------------------------------------------------------------------------

__DEF__ACC(GalaxyBar)
__DEF__POT(GalaxyBar)

//------------------------------------------------------------------------------

