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
	init_int SR;						// stock-recruit relationship: 1==Beverton-Holt; 2==Ricker  BvP: not used
	
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
	//init_number ilinf; 					// Linf for VB growth curve
	//init_number ik; 					// growth rate parameter fro VB growth curve
	//init_number to; 					// time of length 0 for VB growth curve
	//init_number icvl;					// coefficient of variantion for age at length curve
	//init_number ireck;					// recruitment compensation ratio		
	//init_number iRo; 					// Avearge unfished recruitment

	init_number cv_it;					// sd for survey

	init_number sigR;					// sigma for recruitment deviations, fixed? RL: yes, it is better fix this parameter than using the Error in Variable (EIV) approach. otherwise the likelihoods are going to be more complicated and messy 
	init_number sigVul;					// RL: it is the parameter control the variability for the vulnerability. low values means that the vul doesn't change much over time...I guess
	
	//init_int phz_reck;					// phase for estimating reck
	//init_int phz_growth;				// phase for growth parameters
	//init_int use_prior;					// add priors to the likehoods ? (1 or 0)

	init_number u_init;					
	init_int dend;

	
	LOCAL_CALCS
		
		if( dend != 999 )
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			cout<< "dend is:"<<dend<<endl;
			cout<< "sigR is:"<<sigR<<endl;
			ad_exit(1);
		}

	END_CALCS
	
	int nag;
	vector age(sage,nage);					// ages
	vector len(1,nlen);					// middle length of each length bin
	int Am1;							// maximum age minus 1
	number tiny;						// very small number to be used in the fpen function
	
	LOC_CALCS
		// FILL IN SEQUENCE VECTORS
		age.fill_seqadd( sage, 1 );
		len.fill_seqadd( lstp, lstp );

		nag = nage-sage;
		Am1=nage-1;
		tiny=1.e-40;
	END_CALCS

		!! ad_comm::change_datafile_name("length_sra.ctl");
	
	// |---------------------------------------------------------------------------------|
	// | Control File - parameter intitial values and prior info
	// |---------------------------------------------------------------------------------|
	// | ilvec[1] 		-> fisheries catge data (1,fisharea)
	// | ilvec[2]       -> number of survey years       (1,nItNobs)
		

	init_int npar;
	init_matrix theta_control(1,npar,1,7);

	vector   theta_ival(1,npar);
	vector     theta_lb(1,npar);
	vector     theta_ub(1,npar);
	ivector   theta_phz(1,npar);
	ivector theta_prior(1,npar);

	init_vector iwt(syr+1,eyr);         // Recruitment deviations
	//init_vector iwt_init(sage+1,nage);         // Recruitment deviations in initial year
	
	LOC_CALCS
		

		theta_ival  = column(theta_control,1);
		theta_lb    = column(theta_control,2);
		theta_ub    = column(theta_control,3);
		theta_phz   = ivector(column(theta_control,4));
		theta_prior = ivector(column(theta_control,5));


	END_CALCS

INITIALIZATION_SECTION
 	 theta theta_ival;


PARAMETER_SECTION

	init_bounded_vector_vector theta(1,npar,1,1,theta_lb,theta_ub,theta_phz);

	//init_number log_Ro;			    	//Log of average unfished recruitment
	//init_number log_Linf(phz_growth);	//log of l infinity
	//init_number log_k(phz_growth);		//log of k from VB
	//init_number log_cvl(phz_growth);	// log of coefficient of variantion for age at length curve
	//init_number log_reck(phz_reck);		// log of recruitment compensation ratio

	// set growth parameters to true values
	//!! log_Linf=log(ilinf);
	//!! log_k=log(ik);
	//!! log_cvl=log(icvl);
	//!! log_reck=log(ireck);
	//!! log_Ro=log(iRo);

	// log of recruitment deviation
		
	init_bounded_dev_vector log_wt(syr+1,eyr,-5.,5.,2); 
	//!!cout<< "chegou aqui"<<endl;
	//init_bounded_dev_vector log_wt_init(sage+1,nage,-10.,10.,2); 



 	!! log_wt = log(iwt);


 	//!! log_wt_init = log(iwt_init);


	objective_function_value nll;
	
	number fpen;						// penalty to be added to likelihood when posfun is used
	number Linf;						// von Bertalanffy asymptotic length (estimated - based on log_linf)
	number k;							// von Bertalanffy metabolic parameter (estimated - based on log_k)
	number to;
	number cvl;							// coefficient of variation in length at age (estimated - based on log_cvl)
	number reck;						// Goodyear recruitment compensation parameter (estimated - based on log_reck)
	number Ro;							// unfished recruitment (estimated - based on log_Ro)
	//number Rinit;						// recruitment in the first year (estimated - based on log_Rinit)
	number sbo;
	


	number ssvul; 						// it is the sum sq devs for the length vul deviations (mean va(L) - va(L,t))^2
	
	//number Eo;							// unfished egg deposition
	number reca;						// alpha of stock-recruit relationship
	number recb;						// beta of stock-recruit relationship
	number phie;
	number q;							// catchability coefficient (based on zstat)
	number Sa;							// survival-at-age (assume constant across ages)

	vector zstat(1,nyt);				// MLE of q
	vector wt(syr+1,eyr);					// recruitment anomalies
	//vector wt_init(sage+1,nage);				// recruitment anomalies for initial year
	
 	//vector vul(1,nage);					// age-specific vulnerabilities
	vector la(sage,nage);					// length-at-age
	vector std(sage,nage);
	vector wa(sage,nage);					// weight-at-age
	vector lxo(sage,nage);					// unfished survivorship at age
	//vector lz(sage,nage);					// unfished survivorship at age

	//PB added after bias correction
	vector sbt(syr,eyr);					
	
	
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

	//output_runone();

	//cout<<"maxUy"<<endl<<maxUy<<endl;
	//exit(1);

FUNCTION trans_parms
	
	//Bring parameters from log to normal space
	Ro = exp( theta(1,1) );
	//Rinit = exp( theta(2,1) );
	//reck = exp( theta(3,1) );
	//Linf = exp( theta(4,1) );
	//k = exp( theta(5,1) );
	//to =  theta(6,1) ;
	//cvl = exp( theta(7,1) );

	reck = exp( theta(2,1) );
	Linf = exp( theta(3,1) );
	k = exp( theta(4,1) );
	to =  theta(5,1) ;
	cvl = exp( theta(6,1) );
	
	wt = exp( log_wt);
	//wt_init = exp( log_wt_init );

	//cout<<"ok after trans_parms"<<endl;


	
FUNCTION incidence_functions
	
	
	fpen=0.;

	la.initialize();
 	std.initialize();
 	Nat.initialize();
	
	

	la = Linf * ( 1. - exp( -k * ( age - to )));
	std = la * cvl;
	
	Sa = exp(-m);
	
	lxo( sage ) = 1.;
	
	for( int a = sage+1; a <= nage; a++ ) 
	{
		lxo( a ) = lxo(a-1)*Sa;
	}	
	lxo( nage ) /= 1. - Sa;

	
	
	wa = alw * pow( la, blw );

	phie = lxo * fec;

	reca = reck/phie;
	recb = (reck - 1.)/(Ro*phie); 
	sbo  = Ro*phie;


 	
FUNCTION propAgeAtLengh
	// Calculate proportion of length at age class


 	dvar_vector z1( 1, nlen );
	dvar_vector z2( 1, nlen );

 	for( int a = 1; a <= nage; a++ )
	{
		
		// Calculate the integral for proportion age at each length
		z1 = (( len - lstp * 0.5 )-( la( a )))/( std( a ));
		z2 = (( len + lstp * 0.5 )-( la( a )))/( std( a ));
		
		for( int b=1; b<= nlen; b++ )
		{
			P_al( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b ));
		}
		
	}
	
	P_la = trans( P_al );	

	
	
FUNCTION initialYear

	
	// INITIAL YEAR (no fishing assumed)


	Nat( syr, sage )= Ro;
	for( int a = 2; a <= nage; a++ )
	{
		Nat( syr, a ) = Nat( syr, a - 1 ) * Sa;	// initial age-structure
	}		
	Nat( syr, nage ) /= 1. - Sa;

	//Nat(syr,sage)= Rinit * (wt(syr));
	//for( int a = 2; a <= nage; a++ )
	//{
	//	Nat( syr, a ) = Nat( syr, a - 1 ) * Sa * (1. - u_init);	// initial age-structure
	//}		
	//Nat( syr, nage ) /= 1. - (Sa*(1.-u_init));

	//Nat(syr)(sage+1,nage) = Rinit* wt_init;
	//Nat(syr)(sage+1,nage) = elem_prod(Nat(syr)(sage+1,nage), lxo(sage+1,nage));

	
	// length-structure in year-1
	Nlt( syr) = Nat( syr ) * P_al;	
	
	// exploitation by length		
	Ulength( syr ) = elem_div(Clt( syr ), Nlt( syr ));

	// exploitation rate for fully recruited age class
	maxUy( syr ) = max( Ulength( syr ));
	
	// exploitation by age
	Uage( syr ) = Ulength( syr ) * P_la;

	sbt(syr) = fec * Nat(syr);
	
	//cout<<"wt_init"<<endl<<wt_init<<endl;
	//cout<<"Nat"<<endl<<Nat( syr)<<endl;
	//exit(1);

FUNCTION SRA

	// RL: these are the devs for the vul penalty  
	dvar_matrix Upen( syr, eyr, 1, nlen );
			

	// SUBSEQUENT YEARS
	for( int y = syr + 1; y <= eyr; y++ )
	{	
	    dvariable sbtm = fec * Nat(y - 1);	

		Nat( y, sage ) = (reca * sbtm / ( 1. + recb * sbtm)) * (wt( y ) );	///mfexp( sigR*sigR/2.) B-H recruitment
						
		
		// age-distribution post-recruitment
		//Nat(y)(sage+1,nage) = ++elem_prod(Nat(y-1)(sage,nage-1)*Sa,1.-Uage(y - 1)(sage,nage-1));
		
		for( int a = sage +1; a <= nage; a++ ){
			Nat( y )( a ) =  Nat( y - 1 )( a - 1 )* Sa * (1. - Uage( y - 1 )( a -1 ));
				
		}
			
		Nat( y, nage ) /= (1. - Sa * ( 1. - Uage( y - 1, nage)));
		
		for( int aa = sage; aa <= nage; aa++ ){
			Nat( y )( aa ) =  posfun( Nat( y )( aa ), tiny, fpen);
		}


		//=====================================================================================
		//Numbers at lengh

		Nlt( y ) = Nat( y ) * P_al;	
		for( int b = 1; b <= nlen; b++ )
		{
			// length-distribution by year
			
			// exploitation by length										
			Ulength( y, b ) = Clt( y, b ) / posfun( Nlt( y, b ), Clt( y, b ), fpen);	// BvP put this back in to keep values from going over one
			//Ulength( y, b ) = Clt( y, b ) / Nlt( y, b );	
		}


		//Upow( y ) = sum(pow(Ulength( y ),10));
			
		// exploitation by age
		Uage( y) = Ulength( y ) * P_la;			
		
		// max exploitation (fully selected) across lengths
		maxUy( y ) = max( Ulength( y ));
		

		//PB: added
		sbt(y) = fec * Nat(y);				// eggs in year-y

		//cout<<"Uage"<<endl<<Uage( syr)<<endl;
		//cout<<"Nat"<<endl<<Nat( y)<<endl;
	

	}
	
	//cout<<"maxUy"<<endl<<maxUy<<endl;
	//cout<<"Ro"<<endl<<Ro<<endl;
	//cout<<"reck"<<endl<<reck<<endl;
	//cout<<"wt"<<endl<<wt<<endl;
	//cout<<"Nat"<<endl<<Nat<<endl;
	//
	//exit(1);

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
			//Upen( y ) = pow( Ulength( y ) /  maxUy( y ) - muUl, 2. );	
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

	dvar_vector lvec(1,1);
	lvec.initialize();

	lvec(1)=dnorm(zstat,cv_it);
	//lvec(2)=dnorm(log_wt,sigR);
	//lvec(2)=0.;

	dvar_vector npvec(1,npar);
	npvec.initialize();
	
	dvar_vector pvec(1,npar);
	pvec.initialize();

	
	//Priors 
	//prior for h
	for(int i=1;i<=npar;i++)
	{
		switch(theta_prior(i))
		{
				case 1:		//normal
					npvec(i) = dnorm(theta(i,1),theta_control(i,6),theta_control(i,7));
					break;
					
				case 2:		//lognormal CHANGED RF found an error in dlnorm prior. rev 116
					npvec(i) = dlnorm(theta(i,1),theta_control(i,6),theta_control(i,7));
					break;
					
				case 3:		//beta distribution (0-1 scale)
					double lb,ub;
					lb=theta_lb(i);
					ub=theta_ub(i);
					npvec(i) = dbeta((theta(i,1)-lb)/(ub-lb),theta_control(i,6),theta_control(i,7));
					break;
					
				case 4:		//gamma distribution
					npvec(i) = dgamma(theta(i,1),theta_control(i,6),theta_control(i,7));
					break;
					
				default:	//uniform density
					dvariable dummy;
					dummy = posfun(theta(i,1),theta_control(i,2),fpen);   // BvP added to keep parameter within bounds
					dummy = posfun(theta_control(i,3)-theta(i,1),0.0,fpen);   // BvP added to keep parameter within bounds
					npvec(i) = log(1./(theta_control(i,3)-theta_control(i,2)));
					//npvec(i) = (1./(theta_control(i,3)-theta_control(i,2)));
					break;
			}
	}




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
		pvec(1)=dnorm(log_wt,sigR); 	// estimate recruitment deviations with dnorm function
	}
	else
	{
		pvec(1)=norm2(log_wt);//1000;///1000.0; 	
	}
	//pvec(1)=0;


	//if(last_phase())
	//{
	//	cout<<"log_wt is: "<<log_wt<<endl;
	//	cout<<"lvec is: "<<lvec<<endl;
	//	exit(1);	
	//}

	//=====================================================================================
	//pergunta:
	// 
	pvec(2) = ssvul/sigVul;


	//cout<<"sum(npvec) "<<endl<<sum(npvec)<<endl;

	// RL: see commment above
	//=====================================================================================
	
	//nll = sum(lvec) + sum(npvec)+ sum(pvec);
	//nll = sum(lvec) + sum(npvec);//+ sum(pvec);
	nll = sum(lvec) + sum(npvec)+ sum(pvec)+fpen;
	//nll = sum(lvec) + sum(npvec)+ sum(pvec)+fpen*1000;
	
	//nll = sum(lvec) +  sum(pvec);

	//cout<<"nll is "<< nll<<endl;

	//nll = sum(lvec) +  sum(pvec);

	

FUNCTION output_runone
	

	ofstream ofs("runone.rep");
	ofs<<"Nat" << endl << Nat <<endl;
	ofs<<"Nlt " << endl << Nlt <<endl;
	ofs<<"maxUy "<< endl << maxUy <<endl;
	ofs<<"Ulength "<< endl << Ulength <<endl;
	ofs<<"Ro "<< endl << Ro <<endl;
	ofs<<"reck "<< endl << reck <<endl;
	ofs<<"wt "<< endl << wt <<endl;
	ofs<<"reca "<< endl << reca <<endl;
	ofs<<"recb "<< endl << recb <<endl;
	ofs<<"phie "<< endl << phie <<endl;
	ofs<<"sbt "<< endl << sbt <<endl;
	ofs<<"pit "<< endl << psurvB <<endl;
		
		

	cout<<"OK after otput_runone"<<endl;
		

	
	

REPORT_SECTION
	
	
	REPORT(Ro);
	//REPORT(Rinit);
	REPORT(reck);
	REPORT(reca);
	REPORT(recb);
	REPORT(sbo);
	REPORT(Linf);
	REPORT(k);
	REPORT(to);
	REPORT(cvl);
 	REPORT(wt);
 	//REPORT(wt_init);
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
 	REPORT(sbt);
	double depletion = value(sbt(eyr)/sbt(syr));
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
	REPORT(la);
	REPORT(vul);
	REPORT(Nlt);
	REPORT(log_wt);
	

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
