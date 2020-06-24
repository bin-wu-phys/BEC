/*
	In this version, the spacing in p grid is refined for i< nmin
	Modified on Jan 18, 2014.
*/
using namespace std;

#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>     /* atof */
#include "DoubleVector.h"
#include "DoubleMatrix.h"
#include "fVector.h"
#include "fMatrix.h"

#define Nav	1

double 	Qs = 0.5, Nc = 3.0, Nf = 3.0, fzero = 0.1;

double getn0(double f0d){
	return 0.03377372788077926*f0d*(-1.0 + Nc*Nc)*Qs*Qs*Qs;
}

double f0c(double nf){
	return 0.2315750139367751*pow(5.333333333333333 + 
   3.0*nf,4.0)/pow(10.666666666666666 + 7.0*nf,3.0);
}

double Teq(double fz, double nf){
	return 0.8005069837721376*pow(fz/(10.666666666666666 + 7.0*nf),0.25);
}

double ngeq(double temp){
	return 1.9487012517371693*pow(temp,3.0);
}

double nqeq(double temp, double nf){
	return 1.0961444541021579*nf*pow(temp,3.0);
}

double sgeq(double temp){
	return 7.018385351014531*pow(temp,3.0);
}

double sqeq(double temp, double nf){
	return 4.60581538723398*nf*pow(temp,3.0);
}

#define a 10.0

DoubleVector f0(DoubleVector p){	
	DoubleVector f0(p.size());
	for(int i =0; i<p.size();i++)
		if(p[i]<Qs)
			f0[i]=fzero;
		else
			f0[i] = 0;//fzero*exp(-a*pow((p[i]/Qs-1.0),2.0));
		
	return f0;
}

DoubleVector F0(DoubleVector p){
	return 0*p;//exp(-p*p);//1.0/(exp(p)-1.0);//
}

DoubleVector entropyG(DoubleVector f){
	return (1.0+f)*log(1.0+f) - f*log(f);
}

DoubleVector entropyQ(DoubleVector F){
	return (-1.0+F)*log(1.0-F) - F*log(F);
}

int main(int argc, char **argv){

	int 	n=800,
			nmin=20; 	//nmin of refined grids
	double t = 0, dt = 1.0e-6, tMax=50.0;	//setup time variables
	double tSave = 0.1, tCheck=0.001;	//tCheck is used to determine the time interval to check thermalization. 

	//Take the parameters
	switch (argc) {
		case 2:
			fzero = atof(argv[1]);
			break;
		case 3:
			fzero = atof(argv[1]);
			Nf = atof(argv[2]);
			break;
		case 4:
			fzero = atof(argv[1]);
			Nf = atof(argv[2]);
			n = atoi(argv[3]);
			break;
		case 5:
			fzero = atof(argv[1]);
			Nf = atof(argv[2]);
			n = atoi(argv[3]);
			dt = atof(argv[4]);
			break;
		case 6:
			fzero = atof(argv[1]);
			Nf = atof(argv[2]);
			n = atoi(argv[3]);
			dt = atof(argv[4]);
			tSave = atof(argv[5]);
			break;
		case 7:
			fzero = atof(argv[1]);
			Nf = atof(argv[2]);
			n = atoi(argv[3]);
			dt = atof(argv[4]);
			tSave = atof(argv[5]);
			tCheck = atof(argv[6]);
			break;
		default:
		cout << "Usage: ./condensate fzeoro (Nf) (n) (dt) (tSave) (tCheck)\n";exit(0);
	}

	double n0 = getn0(fzero);	

	//Global varia(bles
	double 	dp = 4.0/n,	//spacing in p grids
			dpr = 0.1*dp;	//refined dp
	DoubleVector 	p(n+1),		//p grids in p
			ps(n+1),			//p on each boundary
			ps2(n+1),			//p^2 on each boundary
	 		np(n+1),			//number density of gluons with momentum = p[i]
			jf(n+1),			//gluon current
	 		nFp(n+1),			//number density of quarks with momentum = p[i]
			jF(n+1),			//quark current
			vol(n+1),			//Keeps the volume of each momentum shell
			volInv(n+1),		//Keeps the inverse the volume of each momentum shell
			area(n+1),			//Keeps the area of each momentum shell
			len(n+1),			//Keeps the length of each momentum shell
			lend(n+1),			//Keeps the length of each momentum shell for differentiation
			ntoe(n+1),			//to get energy density from number density
			conversionTerm;		//The conversion term between quarks and gluons

	double			stot, sq, sg,		//entropy 
					etot, eg, eq,		//energy
					ntot, ng, nq;		//number density

	unsigned int nSave = (unsigned int)(tSave/dt), ns = nSave;//time to save the data
	unsigned int nCheck = (unsigned int)(tCheck/dt), nchk = 0;//nCheck;//time to save the data

 	//Temperary variables
	int i,j;
	double 	mD2, mD2G, mD2Q, qhat, qhatG, qhatQ,
		cm1=0,	//cm1 = c_{-1} the coefficient of f = c_{-1}/p-1/2+...
		Ic, 	//The coefficent of Ic
		pfi,Fi,			//f & F at p = i*dp
 		PIvol = 0.5/(PI*PI);	//volume element of phase space

	double 	CF = (Nc*Nc-1.0)/(2.0*Nc), CFONc = CF/Nc, NfCFONc = Nf*CFONc, dtQ=CFONc*dt;//TimeStep for quarks	
	double numg = 2.0*(Nc*Nc-1.0)*PIvol, numq = 4.0*Nf*Nc*PIvol;
	double	integf, integF;

	//Output
	ofstream out, outm, outt;
	ostringstream ostr;
	out.precision(15);
	outm.precision(15);

	for(i=0;i<=n;i++){
		if(i<nmin)
			ps[i]=(i+1)*dpr;
		else
			ps[i]=ps[i-1]+dp;
		if(i==0) p[0] = 0.5*ps[0];
		else p[i]=0.5*(ps[i]+ps[i-1]);
		
		ps2[i]=ps[i]*ps[i];
	}
		
	DoubleVector f = f0(p), pf = p*f, F=F0(p);
	
	for(i=0;i<=n;i++){
		if(i==0){
			vol[0]=pow(ps[0],3.0)/3.0;
			ntoe[0]=0.25*pow(ps[0],4.0)/vol[0];
			area[0]=2.0/pow(ps[0],2.0);
			len[0]=ps[0];
		}
		else{
			vol[i]=(pow(ps[i],3.0)-pow(ps[i-1],3.0))/3.0;
			ntoe[i]=0.25*( pow(ps[i],4.0)-pow(ps[i-1],4.0))/vol[i];
			area[i]=2.0/( pow(ps[i],2.0)-pow(ps[i-1],2.0) );
			len[i]=ps[i]-ps[i-1];
		}
		lend[i]=1.0/(p[i+1]-p[i]);
		np[i]=f[i]*vol[i];
		nFp[i]=F[i]*vol[i];
		volInv[i]=1.0/vol[i];	//Multiplication is cheaper than division.
	}

	ostr << "f0"  << setfill('0') << setw(6) << int(1000000.0*fzero);
	ostr << "Qs"  << setfill('0') << setw(2) << int(10.0*Qs);
	ostr << "Nc"  << setfill('0') << setw(1) << int(Nc);
	ostr << "Nf"  << setfill('0') << setw(1) << int(Nf);
	ostr << "N"  << setfill('0') << n;
	ostr << "dtm"  << setfill('0') << setw(1) << int(-log10(dt)) << "V3.2.1";

	out.open(string(ostr.str()+".dat").c_str(),ios::app);
	out << "# nz = " << n << ", dt = " << dt  << ", dp = " << dp << "\n"; 
	out << "# Nc = " << Nc  << ", Nf = " << Nf  << ", Qs = " << Qs   << ", f0 = " << fzero << "\n";	
	out << "# {p, p*f, F, j_f, j_F, ConversionTerm}\n"; 
	out.close();

	outm.open(string(ostr.str()+"mac.dat").c_str(),ios::app);
	outm << "# nz = " << n << ", dt = " << dt  << ", dp = " << dp << "\n"; 
	outm << "# Nc = " << Nc  << ", Nf = " << Nf  << ", Qs = " << Qs   << ", f0 = " << fzero << "\n";	
	outm << "#{t, qhat, m_d^2, n_g, n_q, n_tot, e_g, e_q, e_tot, s_g, s_q, s_tot, jf0, Ic, c_{-1}, c'_{-1}, c_{0}, c'_{0}, d_{0}, d'_{0}}\n";		
	outm.close();
	
	
	double running_time=clock()/CLOCKS_PER_SEC;	//Start the stopwatch!

	double twoNc = 2.0*Nc, twoNf = 2.0*Nf, mD2dp, integfnm, qhatGnm, jf0=0, pf0;	//Temperary variables to avoid repeated computation
	DoubleVector psp=(ps-p)*lend, plen=p*len, pslend=ps*lend;

	double cm1nm1 = 0,	//cm1 at t-dt to get the time derivative of c_{-1}
		c0nm1 = fzero,      //cm1 at t-dt to get the time derivative of c_{0}
		d0nm1 = 0;      //cm1 at t-dt to get the time derivative of d_{0}
	double tc = 0; 		// time for condensate to set in

	bool	condensate = 0,	// 1 after the condenstate to set in 
		isUnder = fzero<f0c(Nf),	//f0 is under-populated or over-populated.
		isThermal = 0;	//equilibrated?
	double Te, sge, sqe, nge,nqe;
	if(!isUnder){
		Te = Teq(fzero, Nf)*Qs; sge = sgeq(Te);nge = ngeq(Te);sqe = sqeq(Te,Nf);nqe = nqeq(Te,Nf);
	}
	double erreq = 0.05;

	//Time loops
	do{
		integfnm=0;qhatGnm=0;
		for(i=0;i<nmin;i++){
			integfnm+=pf[i];qhatGnm+=pf[i]*(p[i]+pf[i]);
		}

		integf=0;qhatG=0;
		for(i=nmin;i<=n;i++){
			integf+=pf[i];qhatG+=pf[i]*(p[i]+pf[i]);
		}

		integf=integf*dp+integfnm*dpr;integF=(plen*F).sum();
		qhatG = qhatG*dp+qhatGnm*dpr;qhatQ=(nFp*(1.0-F)).sum();

		mD2 = twoNc*integf+twoNf*integF;qhat = Nf*qhatQ + Nc*qhatG;Ic = (integf+integF);
		if(isnan(qhat)) break;

		mD2dp=mD2*dp;
		for(i=0;i<n;i++){
			pfi=pf[i]+ (pf[i+1]-pf[i])*psp[i];jf[i]=qhat*( pslend[i]*(pf[i+1]-pf[i]) - pfi ) + mD2*pfi*( ps[i] + pfi);
			Fi=F[i]+ (F[i+1]-F[i])*psp[i];jF[i]= ( qhat*(F[i+1]-F[i])*lend[i] + mD2*Fi*( 1.0 -Fi) )*ps2[i];
		}
		jf[n]=0;jF[n]=0;	//Boundary condtion jf & jF = 0 at p = pmax
		conversionTerm=Ic*(len*( F*(p+pf)-pf*(1.0-F) ));
		
		// Checking thermalization ...
		if((nchk==nCheck)&&(!isUnder)&&(!isThermal)){

			f=pf/p;ng=(numg*(np.sum())); nq=(numq*(nFp.sum()));
			sg=(numg*( (vol*entropyG(f)).sum())); sq=(numq*( (vol*entropyQ(F)).sum() ));
			if(Nf>0){
				isThermal = (abs(qhat/(mD2*Te)-1.0)<=erreq)&&(abs(ng/nge-1.0)<=erreq)&&(abs(nq/nqe-1.0)<=erreq)&&(abs(sg/sge-1.0)<=erreq)&&(abs(sq/sqe-1.0)<=erreq);
			}
			else
				isThermal = (abs(qhat/(mD2*Te)-1.0)<=erreq)&&(abs(ng/nge-1.0)<=erreq)&&(abs(sg/sge-1.0)<=erreq);
				
			if(isThermal){
				outt.open(string(ostr.str()+"time.dat").c_str(),ios::app);
				outt << "The equilibration time is " << t << endl;
				outt.close();
				if(ns!=0||ns!=nSave){
					ntot=ng+nq;stot=sg+sq;
					eg=(numg*( (ntoe*np).sum() )); eq=(numq*( (ntoe*nFp).sum())); etot=eg+eq;
			
					outm.open(string(ostr.str()+"mac.dat").c_str(),ios::app);
					outm << t << "   " << qhat << "   " << mD2 << "   " << ng << "   " << nq  << "   " << ntot << "   ";
					outm << eg << "   " << eq << "   " << etot << "   " << sg  << "   " << sq << "   " << stot  << "   ";
					outm  << jf0 << "   " << Ic  << "   " << cm1 << "   " << (cm1-cm1nm1)/dt;
					outm  << "   "  << f[0] << "   " << (f[0]-c0nm1)/dt << "   "  << F[0] << "   " << (F[0]-d0nm1)/dt;
					outm  << "   "  << -log(1.0+p[0]/pf[0]) << "   " << -log(1.0/F[0]-1.0) <<  endl;
					outm.close();
			
					out.open(string(ostr.str()+".dat").c_str(),ios::app);
					out << "#t = " << t << "\n";
					out << ((DoubleMatrix)p&pf&F&jf&jF&(conversionTerm/len)&ps) << endl;
					out.close();
				}

			}
			nchk=0;
		}
		//End of checking thermalization

		if(ns==nSave){
			f=pf/p;ng=(numg*(np.sum())); nq=(numq*(nFp.sum()));
			sg=(numg*( (vol*entropyG(f)).sum())); sq=(numq*( (vol*entropyQ(F)).sum() ));
			ntot=ng+nq;stot=sg+sq;
			eg=(numg*( (ntoe*np).sum() )); eq=(numq*( (ntoe*nFp).sum())); etot=eg+eq;
			
			outm.open(string(ostr.str()+"mac.dat").c_str(),ios::app);
			outm << t << "   " << qhat << "   " << mD2 << "   " << ng << "   " << nq  << "   " << ntot << "   ";
			outm << eg << "   " << eq << "   " << etot << "   " << sg  << "   " << sq << "   " << stot  << "   ";
			outm  << jf0 << "   " << Ic  << "   " << cm1 << "   " << (cm1-cm1nm1)/dt;
			outm  << "   "  << f[0] << "   " << (f[0]-c0nm1)/dt << "   "  << F[0] << "   " << (F[0]-d0nm1)/dt;
			outm  << "   "  << -log(1.0+p[0]/pf[0]) << "   " << -log(1.0/F[0]-1.0) <<  endl;
			outm.close();
			
			out.open(string(ostr.str()+".dat").c_str(),ios::app);
			out << "#t = " << t << "\n";
			out << ((DoubleMatrix)p&pf&F&jf&jF&(conversionTerm/len)&ps) << endl;
			out.close();
			ns=0;
		}

		cm1nm1=cm1; c0nm1=pf[0]/p[0]; d0nm1=F[0];
		for(i =0;i<=n;i++){
				if(i==0){
					if(!condensate){
						condensate = (abs(pf[0]*mD2/qhat-1.0)<1.0e-4);//&&(abs(pf[0]/pf[1]-1.0)<1.0e-3);//abs(log(1.0+p[0]/pf[0]))<1.0e-3)&&abs(pf[0]/pf[1]-1.0)<1.0e-3||jf0!=0;
					}else{
						ng=(numg*(np.sum())); nq=(numq*(nFp.sum()));ntot=ng+nq;
						pf0=0;for(j=0;j<Nav;j++) pf0+=pf[j];pf0/=Nav;
						jf0=-qhat*pf0+mD2*pf0*pf0;
						if(abs(ntot/n0-1.0) < 1.0e-8&&jf0<0&&(abs(tc/t-1.0)>1e-2)&&t>0){
							condensate = 0;tc=0;
							outt.open(string(ostr.str()+"time.dat").c_str(),ios::app);
							outt << "The condenstation stops at " << t << endl;
							outt.close();

						}
						
					}

					if(condensate||tc!=0){
						pf0=0;for(j=0;j<Nav;j++) pf0+=pf[j];pf0/=Nav;
						jf0=-qhat*pf0+mD2*pf0*pf0;cm1=pf0;						if(tc==0){ 
						tc = t;
						outt.open(string(ostr.str()+"time.dat").c_str(),ios::app);
						outt << "The condenstate time is " << tc << endl;
						outt.close();
							if(ns!=0||ns!=nSave){
								ntot=ng+nq;stot=sg+sq;
								eg=(numg*( (ntoe*np).sum() )); eq=(numq*( (ntoe*nFp).sum())); etot=eg+eq;
			
								outm.open(string(ostr.str()+"mac.dat").c_str(),ios::app);
								outm << t << "   " << qhat << "   " << mD2 << "   " << ng << "   " << nq  << "   " << ntot << "   ";
								outm << eg << "   " << eq << "   " << etot << "   " << sg  << "   " << sq << "   " << stot  << "   ";
								outm  << jf0 << "   " << Ic  << "   " << cm1 << "   " << (cm1-cm1nm1)/dt;
								outm  << "   "  << f[0] << "   " << (f[0]-c0nm1)/dt << "   "  << F[0] << "   " << (F[0]-d0nm1)/dt;
								outm  << "   "  << -log(1.0+p[0]/pf[0]) << "   " << -log(1.0/F[0]-1.0) <<  endl;
								outm.close();
			
								out.open(string(ostr.str()+".dat").c_str(),ios::app);
								out << "#t = " << t << "\n";
								out << ((DoubleMatrix)p&pf&F&jf&jF&(conversionTerm/len)&ps) << endl;
								out.close();

							}
						}

					}
					else{
						jf0=0;cm1=0;
					} 
					np[i]+= dt*( jf[i]-jf0 + NfCFONc*conversionTerm[i]);		//Boundary condition for jf =0 at p = 0
					nFp[i]+= dtQ*(jF[i] - CF*conversionTerm[i]);		//Boundary condition for jF =0 at p = 0
				}
				else{ 
					np[i]+=dt*(jf[i]-jf[i-1] + NfCFONc*conversionTerm[i]);
					nFp[i]+=dtQ*(jF[i]-jF[i-1] - CF*conversionTerm[i]);
				}

				pf[i] = np[i]*area[i];//f[i]=np[i]*volInv[i];
				F[i]=nFp[i]*volInv[i];
		}
		
		t+=dt;ns++;nchk++;
				
	}while(t<=tMax);
	
	//Compuation time
	out.open(string(ostr.str()+".dat").c_str(),ios::app);

	running_time=clock()/CLOCKS_PER_SEC - running_time;

	if(running_time>=60.0){
		if(running_time>=3600){
			out << "#Calculation is done and the running time is " << running_time/3600.0 << " hours.\n\n\n";
		}
		else {
			out << "#Calculation is done and the running time is " << running_time/60.0 << " minutes.\n\n\n";
		}
	}
	else {
		out << "#Calculation is done and the running time is " << running_time << " seconds.\n\n\n";
	}
	out.close();

}
