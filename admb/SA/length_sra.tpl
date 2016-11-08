//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Brett van Poorten
//Modified by: Roberto Licandeo and Catarina Wor
//Date:	June 21, 2013;  Update: March 2015
//Purpose: length-based SRA based on Carl's spreadsheet
//Notes: 	basic code structure taken from Rob Ahrens - thanks for that			 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
	
DATA_SECTION

	// model dimensions 	
	init_int syr;						// start year of survey data
	init_int eyr;						// end year of survey data
	init_int sage;						// first age
	init_int nage;						// last age
	init_int nlen;						// number of length-bins (assume start at 0)
	init_int lstp;						// width of length-bins
	init_int SR;						// stock-recruit relationship: 1==Beverton-Holt; 2==Ricker
	
	// model known parameters
	init_number m;						// natural mortality
	init_number alw;					// multiplier for length-weight relationship
	init_number blw;					// multiplier for length-weight relationship
	

	// population values
	init_vector vul(sage,nage);				// survey vulnerability
	init_vector fec(sage,nage); 			// iyr(syr,nyr)
	
	// fisheries data
	init_int nyt;						// number of survey observations
	init_vector iyr(1,nyt); 			// iyr(syr,nyr)
	init_vector survB(1,nyt);       	// yt(syr,nyr)
	init_matrix Clt(syr,eyr,1,nlen);	// catch at length and year

	//known growth and recruitment parameters	
	init_number ilinf; 					// Linf for VB growth curve
	init_number ik; 					// growth rate parameter fro VB growth curve
	init_number to; 					// time of length 0 for VB growth curve
	init_number icvl;					// coefficient of variantion for age at length curve
	init_number ireck;					// recruitment compensation ratio		
	init_number iRo; 					// Avearge unfished recruitment
 	init_vector iwt(syr,eyr-1);         // Recruitment deviations

	init_number cv_it;					// sd for survey

	init_number sigR;					// sigma for recruitment deviations, fixed? RL: yes, it is better fix this parameter than using the Error in Variable (EIV) approach. otherwise the likelihoods are going to be more complicated and messy 
	init_number sigVul;					// RL: it is the parameter control the variability for the vulnerability. low values means that the vul doesn't change much over time...I guess
	
	init_int phz_reck;					// phase for estimating reck
	init_int phz_growth;				// phase for growth parameters
	init_int use_prior;					// add priors to the likehoods ? (1 or 0)

	init_int dend;
	
	
	LOCAL_CALCS
		
		if( dend != 999 )
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			cout<< "dend is:"<<dend<<endl;
			ad_exit(1);
		}

	END_CALCS
	
	vector age(sage,nage);					// ages
	vector len(1,nlen);					// middle length of each length bin
	int Am1;							// maximum age minus 1
	number tiny;						// very small number to be used in the fpen function
	
	LOC_CALCS
		// FILL IN SEQUENCE VECTORS
		age.fill_seqadd( sage, 1 );
		len.fill_seqadd( lstp, lstp );
		Am1=nage-1;
		tiny=1.e-20;
	END_CALCS


PARAMETER_SECTION
	init_number log_Ro;			    	//Log of average unfished recruitment
	init_number log_Linf(phz_growth);	//log of l infinity
	init_number log_k(phz_growth);		//log of k from VB
	init_number log_cvl(phz_growth);	// log of coefficient of variantion for age at length curve
	init_number log_reck(phz_reck);		// log of recruitment compensation ratio

	// set growth parameters to true values
	!! log_Linf=log(ilinf);
	!! log_k=log(ik);
	!! log_cvl=log(icvl);
	!! log_reck=log(ireck);
	!! log_Ro=log(iRo);

	// log of recruitment deviation
	init_bounded_dev_vector log_wt(syr,eyr-1,-10.,10.,2);  
 	!! log_wt = log(iwt);

	objective_function_value nll;
	
	number fpen;						// penalty to be added to likelihood when posfun is used
	number Linf;						// von Bertalanffy asymptotic length (estimated - based on log_linf)
	number k;							// von Bertalanffy metabolic parameter (estimated - based on log_k)
	number cvl;							// coefficient of variation in length at age (estimated - based on log_cvl)
	number reck;						// Goodyear recruitment compensation parameter (estimated - based on log_reck)
	number Ro;							// unfished recruitment (estimated - based on log_Ro)
	
	
	number ssvul; 						// it is the sum sq devs for the length vul deviations (mean va(L) - va(L,t))^2
	
	number Eo;							// unfished egg deposition
	number reca;						// alpha of stock-recruit relationship
	number recb;						// beta of stock-recruit relationship
	number phie;
	number q;							// catchability coefficient (based on zstat)
	number Sa;							// survival-at-age (assume constant across ages)

	vector zstat(1,nyt);				// MLE of q
	vector wt(syr,eyr-1);				// recruitment anomalies
	
 	//vector vul(1,nage);					// age-specific vulnerabilities
	vector la(sage,nage);					// length-at-age
	vector std(sage,nage);
	vector wa(sage,nage);					// weight-at-age
	vector lxo(sage,nage);					// unfished survivorship at age
	
	
	vector maxUy(syr,eyr);				// maximum U over length classes for each year?
	vector muUl(1,nlen);				// RL; It is just the mean for the vul(L) (ie. integrated across t) to be used in the penalty for the vulnerability 
	
 	vector psurvB(syr,eyr);				// predicted survey biomass
	
	matrix Nlt(syr,eyr,1,nlen); 		// Matrix of numbers at length class
	matrix P_al(sage,nage,1,nlen);			// proportion of individual of each age at a given length class
	matrix P_la(1,nlen,sage,nage);			// transpose of above
	matrix Nat(syr,eyr,1,nage);			// Numbers of individuals at age
	matrix Ulength(syr,eyr,1,nlen); 	// U (explitation rate) for each length class
	matrix Uage(syr,eyr,1,nage);		// U (explitation rate) for each age

PRELIMINARY_CALCS_SECTION


PROCEDURE_SECTION
    
    trans_parms();
	incidence_functions();
	propAgeAtLengh();
	initialYear();   
	SRA();
	observation_model();
	objective_function();

FUNCTION trans_parms
	
	//Bring parameters from log to normal space
	Linf = exp( log_Linf );
	k = exp( log_k );
	cvl = exp( log_cvl );
	reck = exp( log_reck );
	Ro = exp( log_Ro );
	wt = exp( log_wt );


	
FUNCTION incidence_functions
	
	
	la.initialize();
 	std.initialize();
	double zn;
	

	la = Linf * ( 1. - exp( -k * ( age - to )));
	std = la * cvl;
	
	Sa = exp(-m);
	
	lxo( sage ) = 1.;
	for( int a = sage+1; a <= nage; a++ ) 
	{
		lxo( a ) = lxo(a-1)*Sa;
	}	
		//lxo( nage ) /= 1. - Sa;
	
	
	wa = alw * pow( la, blw );

	phie = lxo * fec;

	reca = reck/phie;
	recb = (reck - 1.)/(Ro*phie); 



 	
FUNCTION propAgeAtLengh
	// Calculate proportion of length at age class


 	dvector z1( 1, nlen );
	dvector z2( 1, nlen );

 	for( int a = 1; a <= nage; a++ )
	{
		
		// Calculate the integral for proportion age at each length
		z1 = (( len - lstp * 0.5 )-value( la( a )))/value( std( a ));
		z2 = (( len + lstp * 0.5 )-value( la( a )))/value( std( a ));
		
		for( int b=1; b<= nlen; b++ )
		{
			P_al( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b ));
		}
		
	}
	
	P_la = trans( P_al );	

	
	
FUNCTION initialYear

	
	// INITIAL YEAR (no fishing assumed)
	Nat( syr, 1 )= Ro;
	for( int a = 2; a <= nage; a++ )
	{
		Nat( syr, a ) = Nat( syr, a - 1 ) * Sa;	// initial age-structure
	}		
	//Nat( syr, nage ) /= 1. - Sa;
	
	// length-structure in year-1
	Nlt( syr) = Nat( syr ) * P_al;	
	
	// exploitation by length		
	Ulength( syr ) = elem_div(Clt( syr ), Nlt( syr ));

	// exploitation rate for fully recruited age class
	maxUy( syr ) = max( Ulength( syr ));
	
	// exploitation by age
	Uage( syr ) = Ulength( syr ) * P_la;
	


FUNCTION SRA

	// RL: these are the devs for the vul penalty  
	dvar_matrix Upen( syr, eyr, 1, nlen );
			

	// SUBSEQUENT YEARS
	for( int y = syr + 1; y <= eyr; y++ )
	{
	  
		dvariable sbt = fec * Nat(y - 1);				// eggs in year-y
	    	
		Nat( y, sage ) = reca * sbt / ( 1. + recb * sbt) * wt( y - 1 );	// B-H recruitment
		
		
		// age-distribution post-recruitment

		for( int a = sage +1; a <= nage; a++ ){

			Nat( y )( a ) =  Nat( y - 1 )( a - 1 )* Sa * (1. - Uage( y - 1 )( a -1 ));
			Nat( y )( a ) =  posfun( Nat( y )( a ), tiny, fpen);	
		}
		
		
		//Nat( y, nage ) /= 1. - Sa * ( 1. - Uage( y - 1, nage));


		//=====================================================================================
		//Numbers at lengh

		Nlt( y ) = Nat( y ) * P_al;	

		for( int b = 1; b <= nlen; b++ )
		{
			// length-distribution by year
			

			// exploitation by length										
			//Ulength( y, b ) = Clt( y, b ) / posfun( Nlt( y, b ), Clt( y, b ), fpen);	
			Ulength( y, b ) = Clt( y, b ) / Nlt( y, b );	
		}
			
		// exploitation by age
		Uage( y) = Ulength( y ) * P_la;			
		
		// max exploitation (fully selected) across lengths
		maxUy( y ) = max( Ulength( y ));

		//cout<<"Uage"<<endl<<Uage( syr)<<endl;
		//cout<<"Nat"<<endl<<Nat( syr)<<endl;
	

	}
	



		for( int b = 1; b <= nlen; b++ )
		{ 
			//  exploitation rate relative to fully recruited U(expected value?) at length over al years
			// RL: It is just the mean for the vul(L) (ie. integrated across t) to be used in the penalty for the vulnerability 
			muUl(b) = sum(elem_div( column(Ulength,b),maxUy ) ) / size_count(maxUy);		
		}
		
		for( int y = syr; y <= eyr; y++ )
		{
			// penalty against dramatic changes in vulnerability?? 
			Upen( y ) = pow( Ulength( y ) / posfun( maxUy( y ), tiny, fpen ) - muUl, 2. );	
		}
		ssvul = sum( Upen );				// vulnerability penalty


	 	psurvB = Nat * elem_prod(wa,vul);
	 	//cout<<"Nat "<<endl<<Nat <<endl;
	 	//cout<<"psurvB "<<endl<<psurvB <<endl;
	 	//exit(1);
	 	
	//psurvB = Nat * wa;


	
	

FUNCTION observation_model


	dvar_vector datry(sage,nage);
	for( int i = 1; i <= nyt ; i++ )
	{
		zstat(i)=log(survB(i))-log(psurvB(iyr(i)));
	}
			
	q=exp(mean(zstat));
	zstat -= mean(zstat);					// z-statistic used for calculating MLE of q
	
	//cout<<"survB "<<endl<<survB <<endl;
	//cout<<"psurvB "<<endl<<psurvB <<endl;
	//cout<<"zstat "<<endl<<zstat <<endl;
	// 	exit(1);
	
FUNCTION objective_function 

	dvar_vector lvec(1,2);
	lvec.initialize();

	lvec(1)=dnorm(zstat,cv_it);
	lvec(2)=dnorm(log_wt,sigR);

	dvar_vector pvec(1,2);
	pvec.initialize();
	
	//if(active(log_reck))
	//{  
 	//	dvariable h=reck/(4.+reck);	
 	//	// beta a=1 and b=1 --> flat
	//	pvec(1)=dbeta((h-0.2)/0.8,1.,1.);
	//	cout<<" (h-0.2)/0.8 is "<< (h-0.2)/0.8<< endl;
	//	cout<<" (h is "<< h << endl;
	//	
   	//}
	
	if(last_phase())
	{		
		pvec(1)=dnorm(log_wt,2.); 	// estimate recruitment deviations with dnorm function
	}
	else
	{
		pvec(1)=100.*norm2(log_wt); 	
	}
	

	//if(last_phase())
	//{
	//	cout<<"log_wt is: "<<log_wt<<endl;
	//	cout<<"lvec is: "<<lvec<<endl;
	//	exit(1);	
	//}

	//=====================================================================================
	//pergunta:
	// I am not sure about what kind of penalty this is
	pvec(2) = ssvul/sigVul;
	// RL: see commment above
	//=====================================================================================
	
	nll = sum(lvec) + sum(pvec)*use_prior;


REPORT_SECTION
	
	REPORT(Ro);
	REPORT(reck);
	REPORT(reca);
	REPORT(recb);
	report<<"bo\n"<<Eo<<endl;
	REPORT(Linf);
	REPORT(k);
	REPORT(cvl);
 	REPORT(wt);
	REPORT(psurvB);
	REPORT(survB);
	
	report<<"bt\n"<<Nat*wa<<endl;
        dvector yield=value(elem_prod(Nat,Uage)*wa);
 	REPORT(yield);
	REPORT(muUl);
	REPORT(maxUy);
 	REPORT(Nat);
	REPORT(Ulength);
	REPORT(Uage);
	REPORT(Clt);
 	dvector sbt=value(Nat.sub(syr,eyr)*fec);
 	REPORT(sbt);
	double depletion = sbt(eyr)/sbt(syr);
	REPORT(depletion);
//	REPORT(depletion);
 	ivector yr(syr,eyr);
	yr.fill_seqadd(syr,1); 
	report<<"yr\n"<<yr<<endl;
	report<<"iyr\n"<<iyr<<endl;
	REPORT(age);
	REPORT(len);
	REPORT(q);
	REPORT(Sa);
	REPORT(vul);
	REPORT(la);
	REPORT(vul);
	REPORT(Nlt)

TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	//#include<MyLikelihoods.cpp>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;

		cout<<"*******************************************"<<endl;
		cout<<"--Start time: "<<ctime(&start)<<endl;
		cout<<"--Finish time: "<<ctime(&finish)<<endl;
		cout<<"--Runtime: ";
		cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
		cout<<"*******************************************"<<endl;
