//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Brett van Poorten
//Modified by: Roberto Licandeo and Catarina Wor
//Date:	June 21, 2013;  Update: March 2015
//Purpose: length-based SRA based on Carl's spreadsheet
//Notes: 	basic code structure taken from Rob Ahrens - thanks for that			 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
	
DATA_SECTION


	init_adstring DataFile;      ///< String for the input datafile name.
	/// | ControlFile.ctl        : controls for phases, selectivity options 
	init_adstring ControlFile;	 ///< String for the control file.
	
	!! ad_comm::change_datafile_name(DataFile);
	// model dimensions 	
	init_int syr;						// start year of survey data
	init_int eyr;						// end year of survey data
	init_int sage;						// first age
	init_int nage;						// last age
	init_int slen;						// start length bin
	init_int nlen;						// number of length-bins (assume start at 0)
	init_int lstp;						// width of length-bins
	init_int SR;						// stock-recruit relationship: 1==Beverton-Holt; 2==Ricker  BvP: not used
	
	// model known parameters
	init_number m;						// natural mortality
	init_number alw;					// multiplier for length-weight relationship
	init_number blw;					// multiplier for length-weight relationship
	init_vector waobs(sage,nage);
	init_int wasw;

	// population values
	init_vector vul(sage,nage);				// survey vulnerability
	init_vector fec(sage,nage); 			// iyr(syr,nyr)
	
	// fisheries data
	init_int nyt;						// number of survey observations
	init_vector iyr(1,nyt); 			// iyr(syr,nyr)
	init_vector survB(1,nyt);       	// yt(syr,nyr)
	init_matrix Clt(syr,eyr,1,nlen);	// catch at length and year

	//known growth and recruitment parameters	
	init_number Linf; 					// Linf for VB growth curve
	init_number k; 					// growth rate parameter fro VB growth curve
	init_number to; 					// time of length 0 for VB growth curve
	init_number cvl;					// coefficient of variantion for age at length curve
	//init_number ireck;					// recruitment compensation ratio		
	//init_number iRo; 					// Avearge unfished recruitment

	init_number cv_it;					// sd for survey
	//init_number sigR;					// sigma for recruitment deviations, fixed? RL: yes, it is better fix this parameter than using the Error in Variable (EIV) approach. otherwise the likelihoods are going to be more complicated and messy 
	
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
			cout<< "Linf is:"<<Linf<<endl;
			cout<< "to is:"<<to<<endl;
			ad_exit(1);
		}

	END_CALCS
	
	int nag;
	vector age(sage,nage);					// ages
	vector len(1,nlen);					// middle length of each length bin
	int Am1;							// maximum age minus 1
	number tiny;						// very small number to be used in the fpen function
	number tinyyy;


	LOC_CALCS
		// FILL IN SEQUENCE VECTORS
		age.fill_seqadd( sage, 1 );
		len.fill_seqadd( slen, lstp );

		nag = nage-sage;
		Am1=nage-1;
		tiny=1.e-7;
		tinyyy=1.e-4;
		
	END_CALCS

	!! ad_comm::change_datafile_name(ControlFile);
	
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

	//init_vector iwt(syr+1,eyr);         // Recruitment deviations
	//init_vector iwt_init(syr-nag,syr-1);         // Recruitment deviations in initial year
	
	init_vector iwt(syr-nag,eyr);

	init_vector plogq(1,3); 

	init_int dend2;  

	LOCAL_CALCS
		
		if( dend2 != 999 )
		{
			cout<<"Error reading control.\n Fix it."<<endl;
			cout<< "iwt is:"<<iwt<<endl;
			cout<< "dend is:"<<dend2<<endl;
			ad_exit(1);
		}

	END_CALCS
	

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
 	
		
	//init_bounded_dev_vector wt(syr+1,eyr,-5.,5.,3); 

	//init_vector wt(syr+1,eyr,2); 
	//init_bounded_vector wt(syr+1,eyr,-5.,5.,3); 
	//init_bounded_dev_vector wt(syr+1,eyr,-5.,5.,3); 
	
	//init_bounded_dev_vector log_wt_init(syr-nag,syr-1,-5.,5.,3); 
	
	//init_bounded_dev_vector wt(syr-nag,eyr,-5.,5.,3); 
	init_bounded_dev_vector wt(syr-nag,eyr,-5.,5.,2); 
	
	//!!cout<< "chegou aqui"<<endl;


	
 	!! wt = (iwt);
 	//!! log_wt_init = log(iwt_init);

 	vector lvec(1,2);
 	vector pvec(1,2);
	objective_function_value nll;
	
	number fpen;						// penalty to be added to likelihood when posfun is used
	number ffpen;	
	number fffpen;	
	//number Linf;						// von Bertalanffy asymptotic length (estimated - based on log_linf)
	//number k;							// von Bertalanffy metabolic parameter (estimated - based on log_k)
	//number to;
	//number cvl;							// coefficient of variation in length at age (estimated - based on log_cvl)
	number reck;						// Goodyear recruitment compensation parameter (estimated - based on log_reck)
	number Ro;							// unfished recruitment (estimated - based on log_Ro)
	//number Rbar;	
	number Rinit;						// recruitment in the first year (estimated - based on log_Rinit)
	number sbo;
	number sigR;
	//number sigVul;	
	//number cv_it;
	


	number ssvul; 						// it is the sum sq devs for the length vul deviations (mean va(L) - va(L,t))^2
	
	//number Eo;							// unfished egg deposition
	number reca;						// alpha of stock-recruit relationship
	number recb;						// beta of stock-recruit relationship
	number phie;
	number q;							// catchability coefficient (based on zstat)
	number Sa;							// survival-at-age (assume constant across ages)

	vector zstat(1,nyt);				// MLE of q

 	//vector vul(1,nage);					// age-specific vulnerabilities
	vector la(sage,nage);					// length-at-age
	vector std(sage,nage);
	vector wa(sage,nage);					// weight-at-age
	vector lxo(sage,nage);					// unfished survivorship at age
	//vector lz(sage,nage);					// unfished survivorship at age

	//PB added after bias correction
	vector sbt(syr,eyr);					
	
	//vector Ulpen(syr,eyr);
	vector maxUy(syr,eyr);				// maximum U over length classes for each year?
	vector penmaxUy(syr,eyr);
	//vector penUm(syr,eyr);

	vector avgUy(syr,eyr);
	vector avgUylast(syr,eyr);
	vector muUl(1,nlen);				// RL; It is just the mean for the vul(L) (ie. integrated across t) to be used in the penalty for the vulnerability 
	//vector Recs(syr+1,eyr);
	//vector Rdevz(syr+1,eyr);
 	//vector delta(syr,eyr);
 	vector psurvB(syr,eyr);				// predicted survey biomass
 	vector umsy(syr,eyr);
 	vector msy(syr,eyr);

	
	matrix Nlt(syr,eyr,1,nlen); 		// Matrix of numbers at length class
	matrix P_al(sage,nage,1,nlen);			// proportion of individual of each age at a given length class
	matrix P_la(1,nlen,sage,nage);			// transpose of above
	matrix Nat(syr,eyr,1,nage);			// Numbers of individuals at age
	matrix Ulength(syr,eyr,1,nlen); 	// U (explitation rate) for each length class
	matrix Uage(syr,eyr,1,nage);		// U (explitation rate) for each age
	matrix penUl(syr,eyr,1,nlen);

PRELIMINARY_CALCS_SECTION
	

PROCEDURE_SECTION
    
    trans_parms();

	incidence_functions();
	propAgeAtLengh();
	initialYear();   
	SRA();
	observation_model();
	objective_function();

	if(last_phase())
	{
		calc_msy();
		output_runone();
	}

	//cout<<"maxUy"<<endl<<maxUy<<endl;
	//exit(1);

FUNCTION trans_parms
	
	//Bring parameters from log to normal space
	Ro = mfexp( theta(1,1) );
	//reck = mfexp( theta(2,1) ); 
	//sigR = mfexp( theta(3,1) );
	
	Rinit = mfexp( theta(2,1) );
	reck = mfexp( theta(3,1) ); 
	sigR = mfexp( theta(4,1) );



	
	//sigVul = mfexp( theta(5,1) );

	//cout<<"Ro is "<< Ro<<endl;
	//cout<<"Rinit is "<< Rinit<<endl;
	//cout<<"sigR is "<< sigR<<endl;
	//cout<<"wt is "<< wt<<endl;
	
	//exit(1);
	//cout<<"ok after trans_parms"<<endl;


	
FUNCTION incidence_functions
	
	
	//fpen.initialize();
	ffpen.initialize();
	fffpen.initialize();
	penmaxUy.initialize();


	la.initialize();
 	std.initialize();
 	Nat.initialize();
 	//Ulpen.initialize();
	
	

	la = Linf * ( 1. - mfexp( -k * ( age - to )));
	std = la * cvl;
	
	Sa = exp(-m);
	
	lxo( sage ) = 1.;
	
	for( int a = sage+1; a <= nage; a++ ) 
	{
		lxo( a ) = lxo(a-1)*Sa;
	}	
	lxo( nage ) /= 1. - Sa;

	
	switch(wasw){
		case 1: 
			wa = waobs;
		break;

		default:
			wa = alw * pow( la, blw );
		break;

	}
	wa = alw * pow( la, blw );

	phie = lxo * fec;

	reca = reck/phie;
	recb = (reck - 1.)/(Ro*phie); 
	sbo  = Ro*phie;





 	
FUNCTION propAgeAtLengh
	// Calculate proportion of length at age class


 	dvar_vector z1( 1, nlen );
	dvar_vector z2( 1, nlen );

 	for( int a = sage; a <= nage; a++ )
	{
		
		// Calculate the integral for proportion age at each length
		z1 = (( len - lstp * 0.5 )- la( a ))/( std( a ));
		z2 = (( len + lstp * 0.5 )- la( a ))/( std( a ));
		
		for( int b=1; b< nlen; b++ )
		{
			P_al( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b ));
		}

		// plus length group
		P_al( a, nlen )= 1.-cumd_norm( z1( nlen ));

	}
	
	P_la = trans( P_al );	

	//cout<<"P_al"<<endl<<column(P_al,nlen )<<endl;
	//cout<<"P_al"<<endl<<P_al<<endl;
	//exit(1);
	
FUNCTION initialYear

	
	// INITIAL YEAR (no fishing assumed)

	

	//Nat( syr, sage )= Ro;
	//for( int a = 2; a <= nage; a++ )
	//{
	//	Nat( syr, a ) = Nat( syr, a - 1 ) * Sa;	// initial age-structure
	//}		
	//Nat( syr, nage ) /= 1. - Sa;




	Nat(syr,sage)= Rinit;  
	
	
	for( int a = sage+1; a <= nage; a++ )
	{
		Nat( syr, a ) = Nat( syr, a - 1 ) * Sa * (1.  - u_init) ;	//  initial age-structure
	}		
	//Nat( syr, nage ) /= 1. - (Sa*(1.-u_init));
	Nat(syr,nage)+=  Nat( syr, nage ) *Sa*  (1. - u_init);
	
	for( int a = sage; a <= nage; a++ ){
		Nat( syr, a ) *=  mfexp(wt(syr-a+1.) );//- sigR*sigR/2.
		Nat( syr )( a ) =  posfun( Nat( syr )( a ), tiny, fpen);
	}
	
	

	//Nat(syr)(sage+1,nage) = Rinit* wt_init;
	//Nat(syr)(sage+1,nage) = elem_prod(Nat(syr)(sage+1,nage), lxo(sage+1,nage));

	
	// length-structure in year-1
	Nlt( syr) = Nat( syr ) * P_al;	
	
	// exploitation by lengt
	for( int b = 1; b <= nlen; b++ )
	{	
		//Ulength( syr )(b) = 1./posfun2(  Nlt( syr, b )/Clt( syr, b ) ,1.0, fpen);	
		//Ulength( syr )(b) = 1.-posfun(1.-(Clt( syr, b ) /  Nlt( syr, b )),tinyyy,fffpen);
		Ulength( syr )(b) = 1.-posfun2(1.-(Clt( syr, b ) /  Nlt( syr, b )),tinyyy,fffpen);
		

		//Ulength( syr )(b) = Clt( syr, b ) /  Nlt( syr, b )t
		//penUl( syr )(b) = posfun(1.-Ulength( syr )(b),tinyyy,fffpen);
	}

	// exploitation rate for fully recruited age class
	maxUy( syr ) = max( Ulength( syr ));
	avgUy( syr ) = mean( Ulength( syr )(1,nlen));
	avgUylast( syr ) = mean( Ulength( syr )(nlen-5,nlen));
	
	penmaxUy(syr) = sum(pow( Ulength( syr )(1,nlen),10));
	

	// exploitation by age CW i think this is wrong
	//Uage( syr ) = Ulength( syr ) * P_la;


	for( int au = sage; au <= nage; au++ ){
		Uage( syr )(au) = Ulength( syr ) * P_al(au);
	}


	sbt(syr) = fec * Nat(syr);
	psurvB(syr) = Nat(syr) * elem_prod(wa,vul);
	
	
	//cout<<"Nat"<<endl<<Nat( syr)<<endl;
	//exit(1);

FUNCTION SRA

	// RL: these are the devs for the vul penalty  
	dvar_matrix Upen( syr, eyr, 1, nlen );
	Upen.initialize();
			

	// SUBSEQUENT YEARS
	for( int y = syr + 1; y <= eyr; y++ )
	{	
	    dvariable sbtm = fec * Nat(y - 1);	
		
		Nat( y, sage ) = (reca * sbtm / ( 1. + recb * sbtm)) * mfexp( wt(y)   );//- sigR*sigR/2.	///mfexp( sigR*sigR/2.) B-H recruitment
						
		
		// age-distribution post-recruitment
		//Nat(y)(sage+1,nage) = ++elem_prod(Nat(y-1)(sage,nage-1)*Sa,1.-Uage(y - 1)(sage,nage-1));
		
		for( int a = sage +1; a <= nage; a++ ){
			Nat( y )( a ) =  Nat( y - 1 )( a - 1 )* Sa * (1. - Uage( y - 1 )( a -1 ));
				
		}
		
		Nat( y, nage ) +=  Nat( y-1, nage ) *Sa*  (1. - Uage( y-1)( nage ));
		//Nat( y, nage ) /= (1. - Sa * ( 1. - Uage( y - 1, nage)));
		
		for( int aa = sage; aa <= nage; aa++ ){
			Nat( y )( aa ) =  posfun( Nat( y )( aa ), tiny, fpen);
		}

		//cout<<"Nat "<<y<<endl<<Nat( y)<<endl;

		//=====================================================================================
		//Numbers at lengh
		Nlt( y ) = Nat( y ) * P_al;	
		for( int b = 1; b <= nlen; b++ )
		{
			// length-distribution by year
			
			// exploitation by length										
			Ulength( y, b ) = 1.-posfun(1.-(Clt( y, b ) /  Nlt( y, b )),tinyyy,fffpen);	// BvP put this back in to keep values from going over one
			//Ulength( y, b ) = 1./posfun2( Nlt( y, b )/Clt( y, b ),1.0, fpen);	// CW inserted this back in to keep values from going over one
			//Ulength( y, b ) = Clt( y, b ) /  Nlt( y, b );	// BvP put this back in to keep values from going over one
			//penUl( y )(b) = posfun(1.-Ulength( y )(b),tinyyy,fffpen);
		}
		//Upow( y ) = sum(pow(Ulength( y ),10));
			
		// exploitation by age -- is this wrong??
		//Uage( y) = Ulength( y ) * P_la;	
		for( int au = sage; au <= nage; au++ ){
			Uage( y )(au) = Ulength( y ) * P_al(au);
		}
			
		
		// max exploitation (fully selected) across lengths
		maxUy( y ) = max( Ulength( y ));
		avgUy( y ) = mean( Ulength( y ) (1,nlen));
		avgUylast( y ) = mean( Ulength( y ) (nlen-5,nlen));
		penmaxUy(y) = sum(pow( Ulength( y )(1,nlen),10));
		//
		
		sbt(y) = fec * Nat(y);				// eggs in year-y
		//cout<<"Uage"<<endl<<Uage( syr)<<endl;
		//cout<<"Nat"<<endl<<Nat( y)<<endl;
	
		psurvB(y) =  Nat(y) * elem_prod(wa,vul);
	}
	
	//cout<<"maxUy"<<endl<<maxUy<<endl;
	//cout<<"Ro"<<endl<<Ro<<endl;
	//cout<<"reck"<<endl<<reck<<endl;
	//cout<<"wt"<<endl<<wt<<endl;
	//cout<<"Nat"<<endl<<Nat<<endl;
	//
	//exit(1);

		//tradition
		//for( int b = 1; b <= nlen; b++ )
		//{ 
		//	//  exploitation rate relative to fully recruited U(expected value?) at length over al years
		//	// RL: It is just the mean for the vul(L) (ie. integrated across t) to be used in the penalty for the vulnerability 
		//	muUl(b) = sum(elem_div( column(Ulength,b),avgUy ) ) / size_count(avgUy);	
		//	//muUl(b) = sum(elem_div( column(Ulength,b),avgUy ) ) / size_count(avgUy);		
		//}
		
		


		
		for( int y = syr+2; y <= eyr; y++ )
		{
			for( int b = 1; b <= nlen; b++ )
			{ 
				muUl(b) = sum(elem_div(  column(Ulength,b)(y-2,y),avgUy(y-2,y) ) ) / 4.;		
			}
			
			// penalty against dramatic changes in vulnerability?? 
			//Upen( y ) = pow( Ulength( y ) / posfun( avgUy( y ), tiny, fpen ) - muUl, 2. );	
			//Upen( y ) = pow( Ulength( y )(nlen-5,nlen) /  maxUy( y ) - avgUylast( y ), 2. );
			
			Upen( y ) = pow( Ulength( y )(1,nlen) /  avgUy( y ) - muUl, 2. );		
		}
		ssvul = sum( Upen );				// vulnerability penalty

		//for( int y = syr+3; y <= eyr; y++ )
		//{
		//	for( int b = 4; b <= nlen; b++ )
		//	{ 
		//		muUl(b) = Ulength(y)(b)-mean(Ulength(y)(b-3,b));
		//	}
		//	// penalty against dramatic changes in vulnerability?? 
		//	//Upen( y ) = pow( Ulength( y ) / posfun( avgUy( y ), tiny, fpen ) - muUl, 2. );	
		//	//Upen( y ) = pow( Ulength( y )(nlen-5,nlen) /  maxUy( y ) - avgUylast( y ), 2. );
		//	
		//	Upen( y ) = (pow( muUl, 2. ));		
		//}
		//ssvul = sum( Upen );				// vulnerability penalty


		//exit(1);


	//FUNCTION stock_recruit_residuals


	//delta = log(Recs)-log(Recz);//+0.5*sigR*sigR;
	//for( int y = syr + 1; y <= eyr; y++ )
	//{	
	 //   dvariable sbtm = fec * Nat(y - 1);
 	//Recs = (reca * sbtm / ( 1. + recb * sbtm)) * (wt( y ));
	
	
FUNCTION observation_model
	
	
	for( int i = 1; i <= nyt ; i++ )
	{
		zstat(i)=(log(survB(i))-log(psurvB(iyr(i))));//-cv_it*cv_it/2.
	}
	
	//cout<<"survB is "<<survB<<endl;
	//cout<<"psurvB is "<<endl<<psurvB<<endl;

	///exit(1);
			
	q=mfexp(mean(zstat));
	
	
	zstat -= mean(zstat);
	//cout<<"q is "<<q<<endl;					// z-statistic used for calculating MLE of q
	//cout<<"zstat is "<<zstat<<endl;
	//zstat -= cv_it*cv_it/2.;
	
	 	//exit(1);
	
FUNCTION objective_function 
	//dvar_vector lvec(1,1);
	lvec.initialize();
	//lvec(1)=dnorm(zstat,cv_it);///100;
	if(last_phase())
	{
		lvec(1)=dnorm(zstat,cv_it);
		dvariable op = 0.;
		op = mean(zstat);
		lvec(2)=1.e5 * op*op;
		
	}else{
		lvec(1)=norm2(zstat)/cv_it;//
		dvariable op = 0.;
		op = mean(zstat);
		lvec(2)=1.e5 * op*op;

		//lvec(1)=norm2(zstat);
	}
	//lvec(2)=dnorm(delta,sigR);
	//lvec(2)=dnorm(log_wt,sigR);
	//lvec(1)=norm2(zstat)/cv_it;
	//lvec(2)=0;
	dvar_vector npvec(1,npar);
	npvec.initialize();
	
	//dvar_vector pvec(1,1);
	pvec.initialize();

	dvar_vector qpvec(1,1);
	qpvec.initialize();
	
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
					
				case 5:		//no prior
					npvec(i) = 0.0;

				default:	//uniform density
					dvariable dummy;
					dummy = posfun(theta(i,1),theta_control(i,2),ffpen);   // BvP added to keep parameter within bounds
					dummy = posfun(theta_control(i,3)-theta(i,1),0.0,ffpen);   // BvP added to keep parameter within bounds
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
	
	//if(last_phase())
	//{		
	//pvec(1)=dnorm(wt,sigR); 	//  estimate recruitment deviations with dnorm function
	//pvec(2)=dnorm(log_wt_init,sigR); 	//  estimate recruitment deviations with dnorm function
	//}
	//else
	//{
	if(last_phase()){
		pvec(1)=dnorm(wt,2.);
		//pvec(1)=dnorm(wt,2.0);
		dvariable s = 0.;
		s = mean(wt);
		pvec(2)=1.e5 * s*s;
	}else{
		pvec(1)=100.*norm2(wt)/sigR;
		dvariable s = 0.;
		s = mean(wt);
		pvec(2)=1.e5 * s*s;
	}
	//sigR;///1000.0; 	
	//pvec(2)=(norm2(log_wt_init));///1000.0; 	
	//}
	//pvec(2)=0.;
	//if(last_phase())
	//{
	//	cout<<"log_wt is: "<<log_wt<<endl;
	//	cout<<"lvec is: "<<lvec<<endl;
	//	exit(1);	
	//}

	int qpr;
	 qpr=plogq(1);

	switch(qpr)
		{
				case 1:		//normal
					qpvec(1) = dnorm(log(q),plogq(2),plogq(3));
					break;
					
				case 2:		
					qpvec(1) = dlnorm(log(q),plogq(2),plogq(3));
					break;
					
				case 3:		//beta distribution (0-1 scale)
					double lb,ub;
					
					cout<<"beta distribution not implemented"<<endl;
					qpvec(1) =0.0;
					break;
					
				case 4:		//gamma distribution
					qpvec(1) = dgamma(log(q),plogq(2),plogq(3));
					break;
					
				case 5:		
					qpvec(1) = 0.0;
					break;

				default:	//uniform density
					dvariable dummy;
					dummy = posfun(log(q),plogq(2),ffpen);   // BvP added to keep parameter within bounds
					dummy = posfun(plogq(3)-log(q),0.0,ffpen);   // BvP added to keep parameter within bounds
					qpvec(1) = log(1./(plogq(3)-plogq(2)));
					//npvec(i) = (1./(the3a_control(i,3)-theta_control(i,2)));
					//cout<<"qpvec(1)"<<qpvec(1)<<endl;
					break;
			}



	//=====================================================================================
	//pergunta:
	//
	
	//pvec(2) = (ssvul/(sigVul));
	//pvec(2) = 0.0;
	//cout<<"sum(npvec) "<<endl<<sum(npvec)<<endl;
	// RL: see commment above
	//==================================================================================
	
	//scenario 1
	nll = sum(lvec) + sum(pvec)+sum(npvec)+ sum(qpvec)+ ssvul/(sigVul)+fffpen+fpen+ffpen +sum(penmaxUy);// +sum(penmaxUy)+ ssvul/(sigVul);// + ssvul/(sigVul)+ffpen; // sum(npvec)+ ssvul/(sigVul);// 
	
	//scenario 2
	//nll = sum(lvec) + sum(pvec)+sum(npvec)+ sum(qpvec)+ ssvul/(sigVul) +sum(penmaxUy);// +sum(penmaxUy)+ ssvul/(sigVul);// + ssvul/(sigVul)+ffpen; // sum(npvec)+ ssvul/(sigVul);// 
	
	//scenario 3 - much worse
	//nll = sum(lvec) + sum(pvec)+sum(npvec)+ sum(qpvec)+ ssvul/(sigVul)+fffpen+fpen+ffpen ;// +sum(penmaxUy) +sum(penmaxUy)+ ssvul/(sigVul);// + ssvul/(sigVul)+ffpen; // sum(npvec)+ ssvul/(sigVul);// 
	
	//scenario 4
	//nll = sum(lvec) + sum(pvec)+sum(npvec)+ sum(qpvec)+ fffpen+fpen+ffpen +sum(penmaxUy);// ssvul/(sigVul)
	
	//scenario 5 -worse of them all
	//nll = sum(lvec) + sum(pvec)+sum(npvec)+ sum(qpvec);// 
	

	//nll = sum(lvec) + sum(npvec);//+ sum(pvec);
	//nll = sum(lvec) + sum(npvec)+ sum(pvec)+fpen+sum(Ulpen);
	//nll = sum(lvec)+ sum(npvec)  + sum(pvec)  +  ffpen +   sum(qpvec)   ;//+ ssvul/(sigVul) + fpen +fffpen+ + sum(penmaxUy) //mean(penmaxUy)  + fffpen*100000 +
	
	//nll = sum(lvec) +  sum(pvec);
	//cout<<"mean(penmaxUy) is "<< mean(penmaxUy)<<endl;
	//cout<<"lvec is "<< sum(lvec)<<endl;
	//cout<<"pvec is "<< sum(pvec)<<endl;
	//cout<<"npvec is "<< sum(npvec)<<endl;
	//cout<<"ffpen is "<< ffpen<<endl;
	//cout<<"fffpen is "<< fffpen<<endl;
	
	//cout<<"q is "<< (q)<<endl;
	//cout<<"qpvec is "<< sum(qpvec)<<endl;

	//nll = sum(lvec) +  sum(pvec);


FUNCTION calc_msy
	dvector utest(1,1001);
	utest.fill_seqadd(0,0.001);
 //This function calculates MSY in the lazy and slow way. 
 	int k, kk ;
	int NF=size_count(utest);
	
	dmatrix selage(syr,eyr,sage,nage);
	
	
	
	for(int y=syr;y<=eyr;y++){
		selage(y)= value(Uage(y)/max(Uage(y)));
		dvector ye(1,NF);
		ye.initialize();
		
		for(k=1; k<=NF; k++)
		{
			dvector lz(sage,nage);
			lz.initialize();
			dvariable phieq;
			dvariable phiz;
			dvariable req;
			phieq.initialize();
			phiz.initialize();
			req.initialize();
			
			lz(sage) = 1.; //first age	
			for(int a = sage+1 ; a <= nage ; a++)
			{
				lz(a) = value(lz(a-1)*Sa*(1.-utest(k))); // proportion of individuals at age surviving M only
			}
			lz(nage) /= value(1.-Sa*(1.-utest(k))); // age plus group
			phiz= lz*fec;
			
			phieq = elem_prod(lz,selage(y))*wa;
			req = Ro*(reck-phie/phiz)/(reck-1.);
			
			ye(k)= value(utest(k)*req*phieq);
		}
		
		msy(y)= max(ye);
		double mtest;	
		for(kk=1; kk<=NF; kk++)
		{
			mtest=ye(kk);
				
			if(mtest==msy(y)){
				umsy(y)=utest(kk);
			} 
		}
	}
	
FUNCTION output_runone
	
	dvector predSurvB(1,nyt);
 	for( int i = 1; i <= nyt ; i++ )
	{
		predSurvB(i)=value(q*psurvB(iyr(i)));
	}
	ofstream ofs("runone.rep");
	ofs<<"Nat" << endl << Nat <<endl;
	ofs<<"Nlt " << endl << Nlt <<endl;
	ofs<<"maxUy "<< endl << maxUy <<endl;
	ofs<<"avgUy "<< endl << avgUy <<endl;
	ofs<<"Ulength "<< endl << Ulength <<endl;
	ofs<<"Ro "<< endl << Ro <<endl;
	ofs<<"Rinit "<< endl << Rinit <<endl;
	ofs<<"reck "<< endl << reck <<endl;
	//ofs<<"wt_init "<< endl << log_wt_init <<endl;
	ofs<<"wt "<< endl << wt <<endl;
	ofs<<"reca "<< endl << reca <<endl;
	ofs<<"recb "<< endl << recb <<endl;
	ofs<<"phie "<< endl << phie <<endl;
	ofs<<"sbt "<< endl << sbt <<endl;
	ofs<<"pit "<< endl << predSurvB <<endl;
	ofs<<"P_al "<< endl << P_al <<endl;
	ofs<<"zstat "<< endl << zstat <<endl;
	ofs<<"fpen "<< endl << fpen <<endl;ofs<<"fpen "<< endl << fpen <<endl;
	ofs<<"lvec "<< endl << sum(lvec) <<endl;
	ofs<<"pvec"<< endl << sum(pvec) <<endl;
	//ofs<<"Ulpen "<< endl << Ulpen <<endl;
		
		
	//cout<<"OK after otput_runone"<<endl;
		
	
	
REPORT_SECTION
	
	output_runone();
	
	REPORT(Ro);
	REPORT(Rinit);
	REPORT(reck);
	REPORT(cv_it);
	REPORT(sigR);
	REPORT(sigVul);
	REPORT(reca);
	REPORT(recb);
	REPORT(sbo);
	REPORT(Linf);
	REPORT(k);
	REPORT(to);
	REPORT(cvl);
 	REPORT(wt);
 	
 
 	dvector predSurvB(1,nyt);
 	for( int i = 1; i <= nyt ; i++ )
	{
		predSurvB(i)=value(psurvB(iyr(i)));
	}
 	//predSurvB=psurvB(iyr(1,nyt));
	
	REPORT(predSurvB);
	REPORT(survB);
	
	report<<"bt\n"<<Nat*wa<<endl;
        dvector yield=value(elem_prod(Nat,Uage)*wa);
 	REPORT(yield);
	REPORT(muUl);
	REPORT(maxUy);
	REPORT(avgUy);
 	REPORT(Nat);
	REPORT(Ulength);
	REPORT(Uage);
 	REPORT(sbt);
	double depletion = value(sbt(eyr)/sbo);
	REPORT(depletion);
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
	REPORT(Clt);
	REPORT(lvec);
	REPORT(pvec);
	REPORT(umsy);
	REPORT(msy);
	REPORT(phie);

RUNTIME_SECTION
convergence_criteria .00001, .0001, .00001
maximum_function_evaluations 100, 1000, 10000

	
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