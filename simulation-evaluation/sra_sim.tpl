//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Brett van Poorten
//Date:	June 21, 2013														 
//Purpose: length-based SRA / VPA based on Carl's spreadsheet
//Notes: 	basic code structure taken from Rob Ahrens - thanks for that			 
//							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	//int mod;										// if mod==1, do SRA; if mod==2, do VPA - given un the command line
	// given in directions.txt
	
	int phz_growth;		// phase for growth parameters
	int est;										// if est==1, estimate real data; if est==0, simulate-estimate
	int seed;										// seed to be used in the simulations
	int debug;										// if debug==1, only calculates the initial model run and spits out lots of numbers	
	int simeval									   //print true and predicted recruitment deviations
	number varon;								// if varon==1, includes variation in simulation parameters, if not, no variation
	int autooff;								//if autooff quits after calculating objective function			
	
	int tempmod;

	LOCAL_CALCS
		ifstream ifsin("directions.txt");       //input stream :  read from this file
		ifsin>>phz_growth>>est>>seed>>debug>>simeval>>varon>>autooff; // read in the following numbers
		//mod=1;          // set SRA mode on
		if( debug )	{ varon=0.;	autooff=1; 
			cout<<" phz_growth\t"<<phz_growth<<"\n est\t"<<est<<"\n seed\t"<<seed<<"\n debug\t"<<debug<<"\n simeval\t"<<simeval<<"\n varon\t"<<varon<<"\n autooff\t"<<autooff<<endl; // print all these variables to the screen
		}
		if( est == 0 ) 
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
		//int on,opt; 
		//if((on=option_match(ad_comm::argc,ad_comm::argv,"-mod",opt))>0) // if -mod is given
		//	mod=2; // set it to VPA mode
	END_CALCS
	

	//variables to be read in	
	init_int syr;								// start year of survey data
	init_int eyr;								// end year of survey data
	init_int nage;								// number of age-classes
	init_int nlen;								// number of length-bins (assume start at 0)
	init_int lstp;								// width of length-bins
	//init_int nityr;                           //numbers of years for which a survey is available
	init_int SR;								// stock-recruit relationship: 1==Beverton-Holt; 2==Ricker
	init_number m;								// natural mortality
	init_number alw;							// multiplier for length-weight relationship
	init_number blw;							// multiplier for length-weight relationship							// multiplier for length-weight relationship
	init_number wmat;							// initial weight at maturity
	

	init_number mat50;							// maturity parameter
	init_number matsd;							// maturity parameter
	init_number ahat;							// vulnerability parameter
	init_number ghat;							// vulnerability parameter
	init_vector va(1,nage);						// survey vulnerability
	init_int nyt
	init_vector survyrs(1,nyt);
	init_ivector iyr(1,nyt); 					// survey counter
	init_vector survB(1,nyt); 					// survey biomass

	init_matrix Clt(syr,eyr,1,nlen);			// catch at length and year	


	//init_vector Sa(1,nage);						// survival-at-age (assume constant across ages)
	
	//I'm including these on a separate file so no need to have them here.
	//init_number ilinf;
	//init_number ik
	//init_number ito
	//init_number icvl;
	//init_number ireck
	//init_number iRo
 	//init_vector iwt(syr,eyr-1);

	init_number cv_it;			// CV for cpue
	init_number sigR;			// sigma R
	init_number sigVul;
	init_int 	phz_reck;		// phase for reck
	init_int 	use_prior;		// add priors to the likehoods ? (1 or 0)


	init_int dend;							// eof
	
	LOCAL_CALCS	
		if( dend != 999 )
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}

		ad_exit(1);

	END_CALCS
	

	//simulation counters
	int iter;
	!! iter=0;
	int sim;                    // simulation switch

	vector age(1,nage);			// ages
	vector len(1,nlen);		    // initial length of each length bin
	int Am1;					// maximum age minus 1
	number tiny;				//small number used for penalize likelihoods	
	
	LOC_CALCS
		age.fill_seqadd( 1, 1 );  // fill in the age vector
		len.fill_seqadd( lstp, lstp ); //fill in length bins,  this case 1 by one
		Am1=nage-1;        // set AM1 to nage-1
		tiny=1.e-20;       // set tiny
	END_CALCS
	

	// import "true" parameter values
	!! ad_comm::change_datafile_name( "Base_pars.dat" ); // change file name to "Base_pars.dat"
	init_number lLinf;              //read in Linf
	init_number lk;				    // vonb k
	init_number to;					//vonb to
	init_number lcvl;				//cv for vonb curve
	init_number lreck;				// recruitment compensation ratio
	init_number lRo;				// average recruitment at virgin levels
	init_vector lwt(syr,eyr-1);		//recruitment deviations
	//init_vector tPplus(syr,eyr-1);  //plus group values
	//init_vector lNtermyear(1,nage); //number at terminal year, only needed for VPA
	init_number dend2;

	LOCAL_CALCS
		if( dend2 != 999 )
		{
			cout<<"Error reading base parameters.\n Fix it."<<endl;
			cout<<dend2<<endl;
			ad_exit(1);
		}
		if( !simeval )
		{
			cout<<"ESTIMATING WITH STOCK RECONSTRUCTION ANALYSIS (SRA)"<<endl;
		}
	END_CALCS

PARAMETER_SECTION
	init_number log_Linf(phz_growth); 		//log(linf) to be estimated in phase grow
	init_number log_k(phz_growth); 			// log(k), estimated in phase grow
	init_number to(phz_growth);				// to. to be estimated in phase grow
	init_number log_cvl(phz_growth);  		// cv at length
	init_number log_reck(phz_reck); 		// Goodyear compensation ratio only estimated in SRA
	init_number log_Ro;  					// average unfished recruitment - only estimated in SRA
	
	

	//init_number log_ptau(-1); // observation error, nor estimated in thei example why?
	init_bounded_dev_vector log_wt(syr,eyr-1,-10.,10.,2); // recruitment deviations. only estimated inSRA
	
	!! log_wt = log(lwt);

	//init_vector t_Pplus(syr,eyr-1,phz2);// vector of plus group at a year, proportion that survived within the group transformed to vary between +-inf
	
	objective_function_value nll;
	
	number fpen;
	number Linf;								// von Bertalanffy asymptotic length
	number k;									// von Bertalanffy metabolic parameter
	number cvl;									// coefficient of variation in length at age
	number reck;								// Goodyear recruitment compensation parameter
	number Ro;									// unfished recruitment
	number ssvul								//stores penalty for dramatic changes in vulnerability
	vector zstat(1,nyt);						// posterior probability distribution of q
	vector wt(syr,eyr-1);						// recruitment anomolies
	
	vector Sa(1,nage);		// survival-at-age (assume constant across ages)
 	vector vul(1,nage);		// age-specific vulnerabilities
	
	number Eo;									// unfished egg deposition
	number reca;								// alpha of stock-recruit relationship
	number recb;								// beta of stock-recruit relationship
	vector la(1,nage);							// length-at-age
	vector wa(1,nage);							// weight-at-age
	vector lxo(1,nage);							// unfished survivorship at age
	vector fec(1,nage);							// age-specific fecundity - used for stock-recruit function
	
	//pergunta: why were these removed?? might pertain to estimation
	vector T_parm(1,nage+eyr-syr);				// true initial age-structure and recruitment
	vector E_parm(1,nage+eyr-syr);				// estimated initial age-structure and recruitment
	vector bias(1,nage+eyr-syr);				// proportional bias in initial age-structure and recruitment estimates
	

	matrix P_la(1,nage,1,nlen);					// probability of being in a length bin given that you are at a given age age
	matrix P_al(1,nlen,1,nage);					// transpose of above
	
	
	//pergunta: what are all these?
	matrix Nat(syr,eyr,1,nage);				//numbers at age
	matrix Ulength(syr,eyr,1,nlen);			//harvest rate by length classes 
	matrix Uage(syr,eyr,1,nage);			// harvest rate by age
	matrix Nlt(syr,eyr,1,nlen);				
	vector maxUy(syr,eyr);
	vector muUl(1,nlen);
 	vector psurvB(syr,eyr);
	number q;


	



PRELIMINARY_CALCS_SECTION
    
	// if you are simulating-estimating data
	if( est==0 )
    {
        sim=1;					//turn the simulation mode on	
        gen_parms();			//generate parameters with or without variation	
        if( debug )				cout<<"parameters generated"<<endl;  
        trans_parms();			//transform parameters from logspace to normal space
        if( debug )				cout<<"parameters transformed"<<endl;
        incidence_functions();	// biological processes that will not be directly estimated in the assessment
        if( debug )				cout<<"incidence functions generated"<<endl;
        initialization();		// calculate recruitment parameters
        if( debug )				cout<<"recruitment parameters calculated"<<endl;
        sim_dynamics();			//simulate new data
        if( debug )				cout<<"dynamics simulated\n START ESTIMATING!!"<<endl;
        //exit(1);
    }

PROCEDURE_SECTION
	
	sim = 0;
	fpen = 0.;
	ssvul = 0.;
	trans_parms();
	if( phz_growth )	
	{	incidence_functions();	
	if( debug )	cout<<"incidence functions calculated"<<endl;	}
	initialization();
	if( debug )				cout<<"recruitment functions predicted"<<endl;
	objective_function();
	if( debug )			cout<<"objective function calculated"<<endl;		
	if( autooff )		exit(1); // stops all the procedure before the estimation starts
	
//_________________________________________
FUNCTION gen_parms

	// add variation to parameters
	random_number_generator rnd( seed );	//cout<<"rnd"<<endl;  create random number genetrators based on seed provided
	dvector r( 1, 4 );   					// creat a vector with 4 spaces	 to store variation in growth and recruitment parameters						
	r.fill_randn( rnd );					//fill the vector with normal random numbers
	dvector r2( syr, eyr-1 );					
	r2.fill_randn( rnd );
	dvector r3( 1, nage );						
	r3.fill_randn( rnd  );
	dvector r4( syr, eyr - 1 );	  			// vector to store recruitment deviations				
	r4.fill_randn( rnd  );        			//fill the vector with normally distributed deviations
	
	double cv1 = 0.1*varon,	cv2 = 0.2*varon, cv3 = 0.4*varon; 	// cv* variation switch, if 0 no variation
	log_Linf = lLinf * ( 1. + r( 3 ) * cv1 );  					// variation in  linfcomes from r(3)
	log_k = lk * ( 1. + r( 4 ) * cv1 );		   					// variation in k comes from r(4)
	log_cvl = lcvl;                            					// cv remains the same
	log_reck = lreck * ( 1. + r( 1 ) * cv1 );  					// recruitment compensation 
	log_Ro = lRo * ( 1. + r( 2 ) * cv2 );						// unexploited average recruitment
	if( varon )	lwt=0.;											// if variaton is on, initialize lwt
	log_wt = lwt + r4 * cv3;									// log_lwt= normally distributed error *cv
	//t_Pplus = tPplus + elem_prod( tPplus, r2 ) * cv2; 			//add error to the tPplus (varies between -inf and +inf)
	//log_Ntermyear = lNtermyear + elem_prod( lNtermyear, r3 ) * cv1;

	ofstream ofs3( "SRA.pin" ); //print real numbers to a pin file, initial guesses
	ofs3<<"#log_Linf\n"<<log_Linf<<endl;
	ofs3<<"#log_k\n"<<log_k<<endl;
	ofs3<<"#to\n"<<to<<endl;
	ofs3<<"#log_cvl\n"<<log_cvl<<endl;
	ofs3<<"#log_reck\n"<<log_reck<<endl;
	ofs3<<"#log_Ro\n"<<log_Ro<<endl;
	ofs3<<"#log_wt\n"<<lwt<<endl;
	ofs3<<"#log_ptau\n"<<0<<endl;
	//ofs3<<"#t_Pplus\n"<<t_Pplus<<endl;
	//ofs3<<"#log_Ntermyear\n"<<log_Ntermyear<<endl;
	ofs3<<"#dend\n"<<999<<endl;
	if( debug )	cout<<"parameters output"<<endl;

//_________________________________________
FUNCTION trans_parms
		//bbring the parameters back to normal space
	Linf = mfexp( log_Linf );
	k 	 = mfexp( log_k );
	cvl  = mfexp( log_cvl );
	reck = mfexp( log_reck ) + 1.;
	Ro 	 = mfexp( log_Ro );
	wt 	 = mfexp( log_wt );
	
	if( !varon ) wt = mfexp( lwt ); 		// if variation mode is not on, set wt=lwt (0?)
	//Pplus = 1. / ( 1. + mfexp( -t_Pplus ));	// transform t_Pplus to a scale between 0 and 1 
	//Ntermyear = mfexp( log_Ntermyear ); 	//numbers at age in the terminal year
	if( debug )	cout<<"reck "<<reck<<"\n Ro "<<Ro<<"\n wt "<<wt<<endl;


FUNCTION incidence_functions
	//  biological processes that will not be directly estimated in the assessment
	dvector z1( 1, nlen ); 
	dvector z2( 1, nlen );
	double zn;
	dvar_vector std( 1, nage );

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
		//if( fec( a ) < 0. )	fec( a ) = 0.; //set negative numbers to 0 - no need for this, plogis doesnt go negative
		
		z1 = (( len - lstp * 0.5 )-value( la( a )))/value( std( a ));//calculate for each length bin the probs tof being of specific ages
		// the "v(la)/std" is to bring the distribution back to the standard normal, and uses len as center of length bin, not it's lower limit
		z2 = (( len + lstp * 0.5 )-value( la( a )))/value( std( a ));

		for( int b=1; b<= nlen; b++ )
		{
			P_la( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b )); // calculates the proportion of being of a given age given your length
		}
	}
	
	P_al = trans( P_la ); //transpose matrix to length by age
	if( debug )	cout<<"la\t"<<la<<"\n wa\t"<<wa<<"\n lxo\t"<<lxo<<"\n fec\t"<<fec<<"\n len\t"<<len<<"\n P_la\n"<<P_la<<"\n P_al\n"<<P_al<<endl;

FUNCTION initialization
	Eo = Ro * sum( elem_prod( lxo, fec )); // calculate the number of eggs at a virgin level
	reca = reck * Ro / Eo; // parameter a of the BH recruitment function	
	
	if( SR )	recb = ( reck - 1. ) / Eo; // 
	else		recb = log( reck ) / Eo;
	
	if( debug )	cout<<"Eo\t"<<Eo<<"\n reca\t"<<reca<<"\n recb\t"<<recb<<endl;

//_________________________________________
FUNCTION sim_dynamics
	random_number_generator rnd3( seed );	//cout<<"rnd"<<endl; // overwrite rnd3 and 4 withdifferent seeds
	random_number_generator rnd4( seed );	//cout<<"rnd"<<endl;
	random_number_generator rnd5( seed );
	
	dvector r2( 1, 2 );		
	r2.fill_randn( rnd4 );
	double r = mfexp( log( 0.8 ) * ( 1. + r2( 1 ) * 0.1 * varon ));				//  growth rate of exploitation of 0.8 with cv=0.1
	double k = -log(value( mean( Sa )));												//carrying capacity of exploitation: M=F
	dvector Uyear( syr, eyr ); 													// vector of exploitation rates
	dvector L50year( syr, eyr ); 												// length at 50% vulnerability
	dvector vuln( 1, nage ); 													// vulnerability at age vector
	dvector eggs( syr, eyr ); 													//no of eggs produced, derived from spawning stock biomass* fecundity
	dmatrix Ulength( syr, eyr, 1, nlen ); 										// matrix of exploitation rate at length
	dmatrix Uage( syr, eyr, 1, nage );  										// matrix of exploitation rate at age
	dmatrix Nat( syr, eyr, 1, nage ); 											// matrix to store numbers at age
	dmatrix Nlt( syr, eyr, 1, nlen ); 											// matrix o store numbers at length
	Nat.initialize(); 															// make sure that Nat is empty
	
	// exploitation rate in the first year is set to very small number
	Uyear( syr ) = 0.001; 
	//logistic curve for development of fishery																		
	Uyear( syr + 1 ) = Uyear( syr ) + Uyear( syr ) * r * ( 1. - Uyear( syr ) / k ); 	
	Uyear( syr + 2 ) = Uyear( syr + 1 ) + Uyear( syr + 1) * r * ( 1. - mean( Uyear( syr, syr + 1)) / k );
	// L50 set to half of Linf
	L50year( syr ) = value( Linf ) * 0.5; 	
	// exploitation rate at length for the first year													
	Ulength( syr ) = Uyear( syr ) / ( 1. +mfexp( -1.7 * ( len - L50year( syr )) / 1.5 )); 
	//cumulative normal vulnerability with sd=1.5
	
	//logistic curve for development of fishery : use running 3 yr average to smooth out trajectory 
	for( int y = syr + 3; y <= eyr; y++ )
	{
		Uyear( y ) = Uyear( y - 1 ) + Uyear( y - 1 ) * r * ( 1. - mean( Uyear( y - 3, y - 1 )) / k ); 
	}
	
	// length at 50% vulnerability declines until that last 10 years according to a logistic curve with up side on the right
	//exploitation at length is cumulative normal approximation with sd 1.5
	for( int y = syr + 1; y <= eyr - 10; y++ )
	{
		L50year( y ) = L50year( y - 1 ) - L50year( y - 1 ) * 0.1 * ( 1. - L50year( y - 1 ) / value( Linf * 0.6 ));	
		Ulength( y ) = Uyear( y ) / ( 1. +mfexp( -1.7 * ( len - L50year( y )) / 1.5 ));
	}
	
	// length at 50% vulnerability increases quickly in the last 10 years;
	//exploitation at length is cumulative normal approximation with sd 1.5
	for( int y = eyr - 9; y <= eyr; y++ )
	{
		L50year( y ) = L50year( y - 1 ) + L50year( y - 1 ) * 0.5 * ( 1. - L50year( y - 1 ) / value( Linf * 0.6 ));	
		Ulength( y ) = Uyear( y ) / ( 1. + mfexp( -1.7 * ( len - L50year( y )) / 1.5 ));
	}
	
	// exploitation at age= exploitation at length *probability of being of a given length given you are from a age class.
	for( int y = syr; y <= eyr; y++ )
	{
		for( int a = 1; a <= nage; a++ )
		{
			Uage( y, a ) = Ulength( y ) * value( P_la( a ));				
		}
	}

	//numbers at age
	for( int a = 1; a <= nage; a++ ) 
	{
		Nat( syr, a ) = value( Ro  * pow( Sa( a ), age( a ) - 1 ));// in the first year 
	}

	Nat( syr, nage ) /= value(1. - Sa( nage ));// plus group in the first year
	eggs( syr ) = value( fec * Nat( syr ));//eggs or ssb in the first year
	
	for( int y = syr + 1; y <= eyr; y++ )
	{
		//recruitment
		Nat( y, 1 ) = value( reca *eggs( y - 1 )/( 1. + value( recb ) * eggs( y - 1 )) * value( wt( y - 1 )));
		//numbers at age
		Nat( y )( 2, nage ) =++ value(elem_prod( elem_prod( Nat( y - 1 )( 1, Am1 ), Sa( 1, Am1 )), 1. - Uage( y -  1 )( 1, Am1 )));
		// plus group
		Nat( y, nage ) += value(Nat( y - 1, nage ) * Sa( nage ) * ( 1. - Uage( y - 1, nage )));	
		eggs( y ) = value( fec * Nat( y ));
		
		for( int b = 1; b <= nlen; b++ )
		{
			// numbers at length bin
			Nlt( y, b )=value( Nat( y ) * P_al( b ));
		}
		//catch at length bin
		Clt( y ) = elem_prod( Nlt( y ), Ulength( y )); 
	}

	dvector r5( 1, nyt );
	r5.fill_randn( rnd5 );		//fill the vector with normal random numbers
    double sig = 0.4;
	
	// parei aqui nityr nao existe mais.
	int j,ii;
	for(j=1;j<=nyt;j++)
	{	
	   	ii=survyrs(j);
	   	survB( j ) = value( sum( elem_prod( elem_prod(Nat(ii), va ), wa )));//survey biomass
	 	survB( j )*= mfexp(r5( j )*sig-sig*sig/2);
		
	}
	
	
	if( simeval )	
	{
		// recruitment for age classess before start of collecting data
		T_parm(1,nage) = Nat( syr ); 
		for( int y = syr + 1; y <= eyr; y++ )
		{
			// recruitment for years in which data was collected: y-syr+nage indicates indexing from 1 to nyrs
			T_parm(y-syr+nage) = column( Nat, 1 )( y );
		}
	}
	if( debug )	cout<<"L50year\t"<<L50year<<"\n Uyear\t"<<Uyear<<"\n Ulength\n"<<Ulength<<"\n Uage\n"<<Uage<<"\n Nat\n"<<Nat<<"\n Ct\n"<<Clt<<endl;

//_________________________________________
FUNCTION SRA
	dvar_vector eggs( 1, nage ); // no of eggs produced by age class
	dvar_matrix Upen( syr, eyr, 1, nlen ); //penalty 
	
	// INITIAL YEAR
	for( int a = 1; a <= nage; a++ ) 
	{
		Nat( syr, a ) = Ro * pow( Sa( a ), age( a ) - 1. );	// syr age-structure
	}				
	Nat( syr, nage ) /= 1. - Sa( nage );					//age plus group in syr
		
	eggs( 1 ) = fec * Nat( syr );							// eggs from syr --- why is it indexed as 1? Just because it is, no special reason
	
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

		maxUy( y ) = max( Ulength( y ));													// max exploitation across lengths
	}
  
  /*	for( int y = syr; y <= eyr; y++ )
	{
		psurvB( y ) = sum( elem_prod( elem_prod( Nat( y ), va ), wa ));
		psurvB( y ) /=1.0e+09;
	}*/
		
	muUl = colsum( Ulength ) / ( eyr - syr + 1 );										// mean exploitation rate across years
	
	for( int y = syr; y <= eyr; y++ )
	{
		Upen( y ) = pow( Ulength( y ) / posfun( maxUy( y ), tiny, fpen ) - muUl, 2. );					// penalty against dramatic changes in vulnerability
	}
	
	psurvB = Nat * elem_prod(wa,va);
	ssvul = sum( Upen );
	//cout<<"muUl\t"<<muUl<<"\n maxUy\t"<<maxUy<<"\n Upen\n"<<Upen<<"\n Uage\n"<<Uage<<"\n Clt\n"<<Clt<<"\n Nat\n"<<Nat<<endl;	exit(1);
	//cout<<Nat<<endl;
	//cout<<survB<<endl;
	//int j,ii;
	//for(j=1;j<=nityr;j++)
	//{	
	//  ii=survyrs(j);
	//  psurvB( j ) = sum( elem_prod( elem_prod( Nat(ii)(1,nage), va ), wa ));
	//}
	//cout<<survB<<endl;
	//	cout<<psurvB<<endl;
	//cout<<zstat<<endl;
	
FUNCTION observation_model
	zstat = log( elem_div( survB, psurvB ));												// z-statistic used for calculating MLE of q
	
	zstat -= mean(zstat);	        // posterior probability distribution of q
 //	cout<<zstat<<endl;
																		// vulnerability penalty
	cout<<zstat<<endl;
	if( simeval )	
	{
		E_parm(1,nage) = Nat( syr );       // Estimated recruitment for age classess before start of collecting data 
		for( int y = syr + 1; y <= eyr; y++ )
		{
			E_parm(y-syr+nage) = column( Nat, 1 )( y ); //Estimated recruitment for years in which data was collected: y-syr+nage indicates indexing from 1 to nyrs
		}
	}
	if( debug )	cout<<"wt\n"<<wt<<"\n Ulength\n"<<Ulength<<"\n Uage\n"<<Uage<<"\n Nat\n"<<Nat<<"\n Nlt\n"<<Nlt<<"\n survB\t"<<survB<<"\n psurvB\t"<<psurvB<<"\n zstat\t"<<zstat<<endl;
	
//_________________________________________

FUNCTION objective_function 
	
	dvar_vector lvec(1,2);
	lvec.initialize();

	lvec(1)=dnorm(zstat,cv_it);
	lvec(2)=dnorm(wt,sigR);dvar_vector 

	pvec(1,4);
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


	
	//if( debug )	cout<<"sig "<<sig<<"\t pvec "<<pvec<<"\t var(wt) "<<var(wt)<<"\t fpen "<<fpen<<"\t ssvul "<<ssvul<<"\t nll "<<nll<<endl;
	//exit(1);
	

REPORT_SECTION
	
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
// 	REPORT(Nat);
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


	if( simeval )
	{
		bias = elem_div( E_parm - T_parm, T_parm );

			ofstream ofssra("SRA bias.txt",ios::app);
			ofssra<<objective_function_value::pobjfun->gmax<<"\t"<<seed-12<<"\t"<<bias<<endl;// outputs obj function value, minimum gradient, seed and bias values
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
