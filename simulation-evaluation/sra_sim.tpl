//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Brett van Poorten
// Modified by: Roberto Licandeo
//				Catarina Wor
//Date:	June 21, 2013														 
//Purpose: length-based SRA / VPA based on Carl's spreadsheet
//Notes: 	basic code structure taken from Rob Ahrens - thanks for that			 
//							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	// given in directions.txt
	
	int phz_growth;									// phase for growth parameters
	int debug;										// if debug==1, only calculates the initial model run and spits out lots of numbers	
	int simeval;						  		    // if simeval== 1 print true and predicted recruitment deviations
	int varon;										// if varon==1, includes variation in simulation parameters, if not, no variation
	int autooff;									//if autooff quits after calculating objective function			
	
	int seed;										// seed to be used in the simulations

	LOCAL_CALCS
		ifstream ifsin("directions.txt");       //input stream :  read from this file
		ifsin>>phz_growth>>debug>>simeval>>varon>>autooff; // read in the following numbers
		//mod=1;          // set SRA mode on
		if( debug )	{ varon=0.;	autooff=1; 
			cout<<" phz_growth\t"<<phz_growth<<"\n debug\t"<<debug<<"\n simeval\t"<<simeval<<"\n varon\t"<<varon<<"\n autooff\t"<<autooff<<endl; // print all these variables to the screen
		}
		
		
		if( simeval ) 
		{	
			//this procedure changes the seed for every simulation run
			ifstream ifs( "seed.txt" ); // if this file is available
			
			//read in the seed
			ifs>>seed;

			// add 12 to the seed				
			seed += 12; 

			//put new seed in to seed.txt
			ofstream ofs( "seed.txt" ); 
			ofs<<seed<<endl; 

			if( debug )		cout<<"model initiated"<<endl;
		}
		
	END_CALCS
	

	//variables to be read in	
	init_int syr;								// start year of survey data
	init_int eyr;								// end year of survey data
	init_int nage;								// number of age-classes
	init_int nlen;								// number of length-bins (assume start at 0)
	init_int lstp;								// width of length-bins

	init_int SR;								// stock-recruit relationship: 1==Beverton-Holt; 2==Ricker
	init_number m;								// natural mortality
	init_number alw;							// multiplier for length-weight relationship
	init_number blw;							// exponent for length-weight relationship
	init_number wmat;							// initial weight at maturity
	
	init_number mat50;							// maturity parameter
	init_number matsd;							// maturity parameter
	init_number ahat;							// vulnerability parameter
	init_number ghat;							// vulnerability parameter
	init_vector va(1,nage);						// survey vulnerability
	init_int nyt;								// number of years for which survey is available
	init_ivector iyr(1,nyt); 					// survey counter
	init_vector survB(1,nyt); 					// survey biomass

	init_matrix Clt(syr,eyr,1,nlen);			// catch at length and year	
	
	init_number cv_it;			// CV for cpue
	init_number sigR;			// sigma R
	init_number sigVul;			// sd for vulnerability
	init_int 	phz_reck;		// phase for reck estimation
	init_int 	use_prior;		// add priors to the likehoods ? (1 or 0)


	init_int dend;							// eof
	
	LOCAL_CALCS	
		if( dend != 999 )
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}

	END_CALCS
	

	int iter;					//simulation counters
	int Am1;					// maximum age minus 1
	vector age(1,nage);			// ages
	vector len(1,nlen);		    // initial length of each length bin	
	number tiny;				//small number used for penalize likelihoods	
	
	LOC_CALCS
		age.fill_seqadd( 1, 1 ); 		 // fill in the age vector
		len.fill_seqadd( lstp, lstp ); 	 //fill in length bins in  this case by lstp unit intervals
		Am1=nage-1;        				 // set AM1 to nage-1
		tiny=1.e-20; 
	 	iter=0;                    // set tiny
	END_CALCS
	

	// import "true" parameter values in log space
	!! ad_comm::change_datafile_name( "Base_pars.dat" ); // change file name to "Base_pars.dat"
	init_number lLinf;              // log of Linf
	init_number lk;				    // log of vonb k
	init_number lto;				// vonb to
	init_number lcvl;				// log  for vonb curve
	init_number lreck;				// recruitment compensation ratio
	init_number lRo;				// average recruitment at virgin levels
	init_vector lwt(syr,eyr-1);		//recruitment deviations
	init_number dend2;

	LOCAL_CALCS
		if( dend2 != 999 )
		{
			cout<<"Error reading base parameters.\n Fix it."<<endl;
			cout<<dend2<<endl;
			ad_exit(1);
		}
		if( debug )		cout<<"data read"<<endl;

	END_CALCS


PARAMETER_SECTION
	
	init_number log_Linf(phz_growth); 		// log(linf) to be estimated in phase grow
	init_number log_k(phz_growth); 			// log(k), estimated in phase grow
	init_number to(phz_growth);				// to to be estimated in phase grow
	init_number log_cvl(phz_growth);  		// cv at length
	init_bounded_number log_reck(1,6,phz_reck); 		// Goodyear compensation ratio only estimated in SRA
	init_number log_Ro;  					// average unfished recruitment - only estimated in SRA
		
	init_bounded_dev_vector log_wt(syr,eyr-1,-10.,10.,2); // recruitment deviations. only estimated in SRA
	
	objective_function_value nll;
	
	number fpen;
	number Linf;								// von Bertalanffy asymptotic length
	number k;									// von Bertalanffy metabolic parameter
	number cvl;									// coefficient of variation in length at age
	number reck;								// Goodyear recruitment compensation parameter
	number Ro;									// unfished recruitment
	number ssvul								//stores penalty for dramatic changes in vulnerability
	vector zstat(1,nyt);						// posterior probability distribution of q
	vector wt(syr,eyr-1);						// recruitment anomalies
	
	vector Sa(1,nage);		// survival-at-age (assume constant across ages)
 	vector vul(1,nage);		// age-specific vulnerabilities
	
	number Eo;									// unfished egg deposition
	number reca;								// alpha of stock-recruit relationship
	number recb;								// beta of stock-recruit relationship
	vector la(1,nage);							// length-at-age
	vector wa(1,nage);							// weight-at-age
	vector lxo(1,nage);							// unfished survivorship at age
	vector fec(1,nage);							// age-specific fecundity - used for stock-recruit function
	
	matrix P_la(1,nage,1,nlen);					// probability of being in a length bin given that you are at a given age age
	matrix P_al(1,nlen,1,nage);					// transpose of above	
	
	matrix Nat(syr,eyr,1,nage);				//numbers at age
	matrix Ulength(syr,eyr,1,nlen);			//harvest rate by length classes 
	matrix Uage(syr,eyr,1,nage);			// harvest rate by age
	matrix Nlt(syr,eyr,1,nlen);				
	vector maxUy(syr,eyr);
	vector muUl(1,nlen);
 	vector psurvB(1,nyt);
	number q;		

	// performance measure storage objects
	vector T_parm(1,nage+eyr-syr);				// true initial age-structure and recruitment
	vector E_parm(1,nage+eyr-syr);				// estimated initial age-structure and recruitment
	number T_Ro;
	number E_Ro;
	number T_kappa;
	number E_kappa;
	vector T_Ulength_fin(1,nlen);
	vector E_Ulength_fin(1,nlen);

	vector bias(1,nage+eyr-syr);				// proportional bias in initial age-structure and recruitment estimates
	vector biasUlength_fin(1,nlen);
	number biasRo;
	number biaskappa;


PRELIMINARY_CALCS_SECTION
    
	// if you are simulating-estimating data
	if( simeval )
    {	
        gen_parms();			//generate parameters with or without variation	
        if( debug )				cout<<"parameters generated"<<endl;  
        trans_parms();			//transform parameters from logspace to normal space
        if( debug )				cout<<"parameters transformed"<<endl;       
        incidence_functions();	// biological processes that will not be directly estimated in the assessment
        if( debug )				cout<<"incidence functions generated"<<endl;
        sim_dynamics();			//simulate new data
        if( debug )				cout<<"dynamics simulated\n START ESTIMATING!!"<<endl;
    }

PROCEDURE_SECTION
	
	fpen = 0.;
	ssvul = 0.;
	trans_parms();
		
	incidence_functions();	
	if( debug ){cout<<"incidence functions calculated"<<endl;	}
		
	SRA();
	if( debug ){cout<<"SRA calculated"<<endl;	}
	
	observation_model();
	if( debug ){cout<<"observation model calculated"<<endl;	}


	objective_function();
	if( debug )			{cout<<"objective function calculated"<<endl;}		
	if( autooff )		exit(1); // stops all the procedure before the estimation starts
	

//_________________________________________
FUNCTION gen_parms

	// add variation to parameters
	random_number_generator rnd( seed );	//cout<<"rnd"<<endl;  create random number genetrators based on seed provided
	
	// variation in growth and recruitment parameters
	dvector r( 1, 4 ); 
	
	//fill the vector with normal random numbers  											
	r.fill_randn(rnd);					

	// vector to store recruitment deviations
	dvector rec_dev( syr, eyr-1 );	  							
	rec_dev.fill_randn(rnd);        			
	
	double cv1 = 0.1*varon,	cv2 = 0.2*varon, cv3 = 0.4*varon; 	// cv* variation switch, if 0 no variation
	log_Linf = lLinf + r( 3 ) * cv1;  					// variation in  linfcomes from r(3)
	log_k = lk + r( 4 ) * cv1 ;		   					// variation in k comes from r(4)
	log_cvl = lcvl;                            					// cv remains the same
	log_reck = lreck + r( 1 ) * cv1 ;  					// recruitment compensation 
	log_Ro = lRo + r( 2 ) * cv2 ;						// unexploited average recruitment
	log_wt += rec_dev * cv3;									// log_lwt= normally distributed error *cv
	
	ofstream ofs3( "sra_sim.pin" ); //print real numbers to a pin file, initial guesses
	ofs3<<"#log_Linf\n"<<log_Linf<<endl;
	ofs3<<"#log_k\n"<<log_k<<endl;
	ofs3<<"#to\n"<<to<<endl;
	ofs3<<"#log_cvl\n"<<log_cvl<<endl;
	ofs3<<"#log_reck\n"<<log_reck<<endl;
	ofs3<<"#log_Ro\n"<<log_Ro<<endl;
	ofs3<<"#log_wt\n"<<log_wt<<endl;
	if( debug )	cout<<"parameters output"<<endl;

//_________________________________________
FUNCTION trans_parms
	
	//bring the parameters back to normal space
	Linf = mfexp( log_Linf );
	k 	 = mfexp( log_k );
	cvl  = mfexp( log_cvl );
	reck = mfexp( log_reck );
	Ro 	 = mfexp( log_Ro );
	wt 	 = mfexp( log_wt );
	
	
	if( debug )	cout<<"reck "<<reck<<"\n Ro "<<Ro<<"\n wt "<<wt<<endl;


FUNCTION incidence_functions
	
	//  biological processes that will not be directly estimated in the assessment
	dvector z1( 1, nlen ); 
	dvector z2( 1, nlen );
	double zn;
	dvar_vector std( 1, nage );

	Sa = exp(-m);

	la.initialize();
	la = Linf * ( 1. - mfexp( -k * ( age - to ))); //average length at age
	
	std = la * cvl; // vector of proportional cds according to length at age
	
	lxo( 1 ) = 1.; //first age
	
	for( int a = 2; a <= nage; a++ )
	{
		lxo( a ) = lxo( 1 ) * pow( Sa( a - 1 ), age( a - 1 )); // proportion of individuals at age surviving M only
	}
	
	lxo( nage ) /= 1. - Sa( Am1 ); // age plus group
	
	wa = alw * pow( la, blw ); //weight at age
	fec = elem_prod(wa,plogis(age,mat50,matsd));; // wight at age-weight at 50% maturity
 	vul = plogis(age,ahat,ghat);

 	for( int a = 1; a <= nage; a++ )
	{
		
		z1 = (( len - lstp * 0.5 )-value( la( a )))/value( std( a ));//calculate for each length bin the probs tof being of specific ages

		z2 = (( len + lstp * 0.5 )-value( la( a )))/value( std( a ));

		for( int b=1; b<= nlen; b++ )
		{
			P_la( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b )); // calculates the proportion of being of a given age given your length
		}
	}
	
	P_al = trans( P_la ); //transpose matrix to length by age


	if( debug )	cout<<"la\t"<<la<<"\n wa\t"<<wa<<"\n lxo\t"<<lxo<<"\n fec\t"<<fec<<"\n len\t"<<len<<"\n P_la\n"<<P_la<<"\n P_al\n"<<P_al<<endl;

	Eo = Ro * sum( elem_prod( lxo, fec )); // calculate the number of eggs at a virgin level
	reca = reck * Ro / Eo; // parameter a of the BH recruitment function	
	
	if( SR )	recb = ( reck - 1. ) / Eo; // 
	else		recb = log( reck ) / Eo;
	
	if( debug )	cout<<"Eo\t"<<Eo<<"\n reca\t"<<reca<<"\n recb\t"<<recb<<"\n reck\t"<<reck<<endl;

//_________________________________________
FUNCTION sim_dynamics
	random_number_generator rnd3( seed );	//cout<<"rnd"<<endl; // overwrite rnd3 and 4 withdifferent seeds
	
	dvector r2( 1, 2 );		
	r2.fill_randn(rnd3);

	double r_fishery = mfexp( log( 0.8 ) * ( 1. + r2( 1 ) * 0.1 * varon ));				//  growth rate of exploitation of 0.8 with cv=0.1
	double k_fishery = -log(value( mean( Sa )));										//carrying capacity of exploitation: M=F
	
	dvector T_Uyear( syr, eyr ); 														// vector of exploitation rates
	dvector T_L50year( syr, eyr ); 														// length at 50% vulnerability
	dvector T_vuln( 1, nage ); 															// vulnerability at age vector
	dvector T_eggs( syr, eyr ); 														//no of eggs produced, derived from spawning stock biomass* fecundity
	dmatrix T_Ulength( syr, eyr, 1, nlen ); 											// matrix of exploitation rate at length
	dmatrix T_Uage( syr, eyr, 1, nage );  												// matrix of exploitation rate at age
	dmatrix T_Nat( syr, eyr, 1, nage ); 												// matrix to store numbers at age
	dmatrix T_Nlt( syr, eyr, 1, nlen ); 												//matrix of numbers at length classes
	dmatrix T_Clt(syr,eyr,1,nlen);														// matrix o store numbers at length
	T_Nat.initialize(); 																// make sure that Nat is empty
	
	// exploitation rate in the first year is set to very small number
	T_Uyear( syr ) = 0.001; 
	
	//logistic curve for development of fishery	3-step average being used as a smoothing function																		
	T_Uyear( syr + 1 ) = T_Uyear( syr ) + T_Uyear( syr ) * r_fishery * ( 1. - T_Uyear( syr ) / k_fishery ); 	
	T_Uyear( syr + 2 ) = T_Uyear( syr + 1 ) + T_Uyear( syr + 1) * r_fishery * ( 1. - mean( T_Uyear( syr, syr + 1)) / k_fishery );
	
	// L50  (selectivity) set to half of Linf
	T_L50year( syr ) = value( Linf ) * 0.5;

	// exploitation rate at length for the first year	u*vul ate length -- cumulative normal vulnerability with sd=1.5												
	T_Ulength( syr ) = T_Uyear( syr ) / ( 1. + mfexp( -1.7 * ( len - T_L50year( syr )) / 1.5 )); 
	
	
	//logistic curve for development of fishery : use running 3 yr average to smooth out trajectory 
	for( int y = syr + 3; y <= eyr; y++ )
	{
		T_Uyear( y ) = T_Uyear( y - 1 ) + T_Uyear( y - 1 ) * r_fishery * ( 1. - mean( T_Uyear( y - 3, y - 1 )) / k_fishery ); 
	}
	
	// length at 50% vulnerability declines until that last 10 years according to a logistic curve with up side on the right
	//exploitation at length is cumulative normal approximation with sd 1.5
	for( int y = syr + 1; y <= eyr - 10; y++ )
	{
		T_L50year( y ) = T_L50year( y - 1 ) - T_L50year( y - 1 ) * 0.1 * ( 1. - T_L50year( y - 1 ) / value( Linf * 0.6 ));	
		T_Ulength( y ) = T_Uyear( y ) / ( 1. +mfexp( -1.7 * ( len - T_L50year( y )) / 1.5 ));
	}
	
	// length at 50% vulnerability increases quickly in the last 10 years;
	//exploitation at length is cumulative normal approximation with sd 1.5
	for( int y = eyr - 9; y <= eyr; y++ )
	{
		T_L50year( y ) = T_L50year( y - 1 ) + T_L50year( y - 1 ) * 0.5 * ( 1. - T_L50year( y - 1 ) / value( Linf * 0.6 ));	
		T_Ulength( y ) = T_Uyear( y ) / ( 1. + mfexp( -1.7 * ( len - T_L50year( y )) / 1.5 ));
	}
	
	// exploitation at age= exploitation at length *probability of being of a given length given you are from a age class.
	for( int y = syr; y <= eyr; y++ )
	{
		for( int a = 1; a <= nage; a++ )
		{
			T_Uage( y, a ) = T_Ulength( y ) * value( P_la( a ));				
		}
	}

	//numbers at age
	for( int a = 1; a <= nage; a++ ) 
	{
		T_Nat( syr, a ) = value( Ro  * pow( Sa( a ), age( a ) - 1 ));// in the first year 
	}

	T_Nat( syr, nage ) /= value(1. - Sa( nage ));// plus group in the first year
	T_eggs( syr ) = value( fec * T_Nat( syr ));//eggs or ssb in the first year
	
	for( int y = syr + 1; y <= eyr; y++ )
	{
		//recruitment
		T_Nat( y, 1 ) = value( reca * T_eggs( y - 1 )/( 1. + value( recb ) * T_eggs( y - 1 )) *  wt( y - 1 ));
		//numbers at age
		T_Nat( y )( 2, nage ) =++ value(elem_prod( elem_prod( T_Nat( y - 1 )( 1, Am1 ), Sa( 1, Am1 )), 1. - T_Uage( y -  1 )( 1, Am1 )));
		// plus group
		T_Nat( y, nage ) += value(Nat( y - 1, nage ) * Sa( nage ) * ( 1. - T_Uage( y - 1, nage )));	
		T_eggs( y ) = value( fec * T_Nat( y ));
		
		for( int b = 1; b <= nlen; b++ )
		{
			// numbers at length bin
			T_Nlt( y, b )=value( T_Nat( y ) * P_al( b ));
		}
		//catch at length bin
		T_Clt( y ) = elem_prod( T_Nlt( y ), T_Ulength( y )); 
	}

	//survey obs error
	dvector surv_err( 1, nyt );
	surv_err.fill_randn( rnd3 );		//fill the vector with normal random numbers
    double sig = 0.4;
	
	int j,ii;
	for(j=1;j<=nyt;j++)
	{	
	   	ii=iyr(j);
	   	survB( j ) = value( sum( elem_prod( elem_prod(T_Nat((syr+ii-1)), va ), wa )));//survey biomass
	 	if(varon == 1) {survB( j )*= mfexp(surv_err( j )*sig-sig*sig/2);}
	 	survB( j ) /=1.0e+09;	// transform to 100s of tons	
	}


		// recruitment for age classess before start of collecting data
		T_parm(1,nage) = T_Nat( syr ); 
		T_Ro = Ro;
		T_kappa = reck;
		T_Ulength_fin(1,nlen) = T_Ulength(syr);


		for( int y = syr + 1; y <= eyr; y++ )
		{
			// recruitment for years in which data was collected: y-syr+nage indicates indexing from 1 to nyrs
			T_parm(y-syr+nage) = column( T_Nat, 1 )( y );
		}


	if( debug )	cout<<"L50year\t"<<T_L50year<<"\n Uyear\t"<<T_Uyear<<"\n Ulength\n"<<T_Ulength<<"\n Uage\n"<<T_Uage<<"\n Nat\n"<<T_Nat<<"\n Ct\n"<<T_Clt<<endl;

//_________________________________________
FUNCTION SRA
	dvar_vector eggs( syr, eyr ); // no of eggs produced
	dvar_matrix Upen( syr, eyr, 1, nlen ); 
	
	// INITIAL YEAR
	for( int a = 1; a <= nage; a++ ) 
	{
		Nat( syr, a ) = Ro * pow( Sa( a ), age( a ) - 1. );	// syr age-structure
	}				
	Nat( syr, nage ) /= 1. - Sa( nage );					//age plus group in syr
		
	eggs( syr ) = fec * Nat( syr );							// eggs from syr --- why is it indexed as 1? Just because it is, no special reason
	
	for( int b = 1; b <= nlen; b++ )
	{
		Nlt( syr, b ) = Nat( syr ) * P_al( b );												// length-structure in year-1
		Ulength( syr, b ) = Clt( syr, b ) / posfun( Nlt( syr, b ), Clt( syr, b ), fpen);	// exploitation by length, posfun is there to ensure that U is not negative and doesn't go above 1		
	}

	maxUy( syr ) = max( Ulength( syr ));		// maximum exploitation rate in the first year														
	
	for( int a = 1; a <= nage; a++ )
	{
		Uage( syr, a ) = Ulength( syr ) * P_la( a );	// exploitation by age: exploitation at lenght*prob of being at that length given you are of  given age - vector
	}

	// SUBSEQUENT YEARS
	for( int y = syr + 1; y <= eyr; y++ )
	{
		if( SR ==1 )
		{
			Nat( y, 1 ) = reca * eggs( y - 1 ) / ( 1. + recb * eggs( y - 1 )) * wt( y - 1 );	// recruitment
		}else{
			Nat( y, 1 ) = reca * eggs( y - 1 ) * mfexp( - recb * eggs( y - 1 )) * wt( y - 1 );
		}

		Nat( y )( 2, nage ) =++ posfun( elem_prod( elem_prod( Nat( y - 1 )( 1, Am1 ), Sa( 1, Am1 )), 1. - Uage( y - 1 )( 1, Am1 )), tiny, fpen );		// age-distribution post-recruitment
		Nat( y, nage ) += posfun( Nat( y - 1, nage ) * Sa( nage ) * ( 1. - Uage( y - 1, nage )), tiny, fpen);
		

		eggs( y ) = fec * Nat( y );																// eggs in year-y
		
		for( int b = 1; b <= nlen; b++ )
		{
			Nlt( y, b ) = Nat( y ) * P_al( b );												// length-distribution
			Ulength( y, b ) = Clt( y, b ) / posfun( Nlt( y, b ), Clt( y, b ), fpen);						// exploitation by length
		}
		
		for( int a = 1; a <= nage; a++ )
		{
			Uage( y, a ) = Ulength( y ) * P_la( a );										// exploitation by age
		}
			
		maxUy( y ) = max( Ulength( y ));
													// max exploitation across lengths
	}
	
  	for( int y = 1; y <= nyt; y++ )
	{
		psurvB( y ) = sum( elem_prod( elem_prod( Nat( y+syr-1 ), va ), wa ));
		psurvB( y ) /=1.0e+09;
	}
		
	muUl = colsum( Ulength ) / ( eyr - syr + 1 );										// mean exploitation rate across years
	
	for( int y = syr; y <= eyr; y++ )
	{
		Upen( y ) = pow( Ulength( y ) / posfun( maxUy( y ), tiny, fpen ) - muUl, 2. );					// penalty against dramatic changes in vulnerability
	}
	
	ssvul = sum( Upen );
	
	
FUNCTION observation_model
	

	zstat = log( elem_div( survB, psurvB ));												// z-statistic used for calculating MLE of q

	zstat -= mean(zstat);	        // posterior probability distribution of q
 																		// vulnerability penalty
	if(last_phase())
	{
		if( simeval )	
		{
			E_parm(1,nage) = Nat( syr );       // Estimated recruitment for age classess before start of collecting data 
			E_Ro = Ro;
			E_kappa = reck;
			E_Ulength_fin(1,nlen) = Ulength(syr);

			for( int y = syr + 1; y <= eyr; y++ )
			{
				E_parm(y-syr+nage) = column( Nat, 1 )( y ); //Estimated recruitment for years in which data was collected: y-syr+nage indicates indexing from 1 to nyrs
			}
		}
	}
	
//_________________________________________

FUNCTION objective_function 
	
	dvar_vector lvec(1,2);
	lvec.initialize();

	lvec(1)=dnorm(zstat,cv_it);
	lvec(2)=dnorm(wt,sigR);

	dvar_vector pvec(1,4);
	pvec.initialize();

	if(active(log_reck))
	{  
 		dvariable h=reck/(4.+reck);	
		pvec(1)=dbeta((h-0.2)/0.8,1.,1.);
   	}

   	if(last_phase())
	{
		pvec(2)=dnorm(log_wt,2.); 	// estimate recruitment deviations with dnorm function
	}
	else
	{
		pvec(3)=100.*norm2(log_wt); 	// assume deviations are at their normal mle analytical estimator  (obs-pred)^2
	}
	
	pvec(4) = ssvul/sigVul;
	nll = sum(lvec) + sum(pvec)*use_prior;



REPORT_SECTION	

	if( !simeval )
	{
		report<<"Ro\n"<<Ro<<endl;
		report<<"bo\n"<<Eo<<endl;
		report<<"kappa\n"<<reck<<endl;
		dvariable h = value(reck/(4.+reck));
		REPORT(h);
		REPORT(Linf);
		REPORT(k);
		REPORT(cvl);
 		REPORT(wt);
		REPORT(reca);
		REPORT(recb);
		REPORT(psurvB);
		REPORT(survB);
		report<<"yt\n"<<survB<<endl;
		report<<"bt\n"<<Nat*wa<<endl;
		dvar_vector wl = alw * pow( len, blw );
		REPORT(wl);
		REPORT(muUl);
		REPORT(maxUy);
		
		report<<"N\n"<<Nat<<endl;
		REPORT(Ulength);
		REPORT(Uage);
		REPORT(Clt);
		REPORT(Nlt);
 		dvector sbt=value(Nat.sub(syr,eyr)*fec);
 		REPORT(sbt);
		double depletion = sbt(eyr)/sbt(syr);
		REPORT(depletion);
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
		REPORT(va);
		REPORT(wa);
	}else{

		if(last_phase()){

		iter += 1;

		bias = elem_div( E_parm - T_parm, T_parm );
		biasRo = (E_Ro- T_Ro)/ T_Ro;
		biaskappa = (E_kappa- T_kappa)/ T_kappa;
		biasUlength_fin = elem_div( E_Ulength_fin - T_Ulength_fin, T_Ulength_fin );

		ofstream ofssra("SRA_bias.txt",ios::app);

		// outputs obj function value, minimum gradient, seed and bias values
		ofssra<<objective_function_value::pobjfun->gmax<<"\t"<<seed-12<<"\t"<<bias<<"\t"<< biasRo<<"\t"<< biaskappa<<"\t"<< biasUlength_fin <<endl;
		
		ofstream ofspar("SRA_param.txt",ios::app);

		// outputs obj function value, minimum gradient, seed and bias values
		ofspar<<objective_function_value::pobjfun->gmax<<"\t"<<seed-12<<"\t"<<E_parm<<"\t"<<T_parm<<"\t"<<E_Ro<<"\t"<<T_Ro<<"\t"<<E_kappa<<"\t"<<T_kappa<<"\t"<< E_Ulength_fin<<"\t"<< T_Ulength_fin <<endl;
		

		}
	}

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
	//#include <contrib.h>//IF you have ADMB-11
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
	if( !simeval )
	{
		cout<<"*******************************************"<<endl;
		cout<<"--Start time: "<<ctime(&start)<<endl;
		cout<<"--Finish time: "<<ctime(&finish)<<endl;
		cout<<"--Runtime: ";
		cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
		cout<<"*******************************************"<<endl;
	}
