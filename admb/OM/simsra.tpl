//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Roberto Licandeo, Brett Van Poorten and Catarina Wor
//Date:	June 21, 2013; Update: 8 july 2014 
//Purpose: length-based SRA / VPA based on Carl's spreadsheet
//Notes: 	basic code structure taken from Rob Ahrens - thanks for that			 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	
	// Read new seed in
	int seed;

	LOC_CALCS
		ifstream ifs( "../seed.txt" ); // if this file is available
		ifs>>seed; //read in the seed
		seed += 3; // add 10 to the seed
		ofstream ofs( "../seed.txt" ); //put out to seed.txt
		ofs<<seed<<endl; //the new value of the seed


	END_CALCS

	//scenario name
	init_int scnNumber;

	//model dimensions
	init_int syr;
	init_int eyr;
	init_int rep_yr;
	init_int sage;
	init_int nage;
	init_int nlen;
	init_int slen;
	init_int lstp;
	init_int niyr;

	init_int proc_err;
	init_int obs_err;
	
	//true parameter values
	init_number sigR; 	    		// standard deviation for recruitment deviations
	init_number sigVul; 
	init_number tau;				// standard deviation for survey observation error
	init_number tau_length;			// standard deviation for observation error
	init_number Sa;					// natural survival
	
	init_number wmat;				// 50% maturity 
	//init_vector fec(sage,nage); 			// fecundity at age
	
	//Growth control parameters
	//init_int nlinfch;
	init_number Linf;				// vb Linfinity
	init_number k;					// vb k
	init_number to;					// vb to
	init_number cvl;				// coefficient of variantion for age at length curve
	init_number alw;				// parameter for length weight relationship
	init_number blw;				// parameter for length weight relationship
	init_number reck;				// recruitment compensarion ratio
	init_number Ro;					// Average unfished recruitment
	init_number q;					// catchability coefficient


	//Selectivity control parameters
	init_int nselch;			//number of times selectivitu changes in the time series
	init_matrix selControl(1,nselch,1,4);

	init_number Ulenmu;					// 50% selectivity at lenth
	init_number Ulensd;					// sd for selectivity at lenght


	init_vector va(sage,nage);        // vector of survey vulnerabilities at age
	init_vector ut(syr,eyr);		// exploitation rate pattern
	init_vector iyr(1,niyr);     	// years for which a survey is available

	init_int dend;					// end of file number

		LOCAL_CALCS
		if( dend != 999 )
		{
			cout<<"Error reading true values.\n Fix it."<<endl;
			cout<< niyr<<endl;
			cout<< iyr <<endl;
			cout<< dend <<endl;
			ad_exit(1);
		}
		
	END_CALCS


	// derived quantities
	number reca;					// BH parameter a
	number recb;					// BH parameter b
	number Eo;						// Average unfished egg production	
	number Am1;						// penultimate age
	number phie;					// Equilibrium fished fecundity
	number sbo;						// Unfished stock spawning biomass
	//number test;					// test to see if S-R function is calculated correctly, should return Ro 
	number tiny;						// very small number to be used in the fpen function
	number sampleSigR;

	//selectivity parameters
	vector selyr(1,nselch);
	vector selb(1,nselch);
	vector sela(1,nselch);
	vector selg(1,nselch);
	ivector indselyr(syr,eyr);

	matrix  sellen(1,nselch,1,nlen);


	LOC_CALCS
		selyr = column(selControl,1);
		selb  = column(selControl,2);
		sela  = column(selControl,3);
		selg  = column(selControl,4);

		int tmp;
       	tmp= 1;

       
       	for(int y=syr;y<=eyr;y++)
       	{
       		if(selyr(tmp)>=y){
       			indselyr(y)=tmp;
       		}else{
       			indselyr(y)=tmp+1;
       			tmp++;
       		}
       	}

       
	END_CALCS

	vector wt(syr,eyr);				// Recruitment deviations
	vector eps(syr,eyr); 			// Observation errors deviations
	vector age(sage,nage); 			// age vector
	vector vbt(syr,eyr); 			// vulnerable biomass (assuming a single gear)
	vector ct(syr,eyr);  			// catches
	//vector bt(syr,eyr);	  			// total biomass
	vector sbt(syr,eyr); 			// spawning biomass
	vector depl(syr,eyr);			// depletion based on spawning biomass
	vector maxUy(syr,eyr);   		// annual U based on vulberable biomass
	vector avgUy(syr,eyr);

	vector len(1,nlen); 			// length classes
	vector lxo(sage,nage); 			// Equilibrium unfished numbers ata age
	vector wa(sage,nage); 				// Weight at age
	vector fec(sage,nage); 			// fecundity at age
	vector la(sage,nage); 				// length at age
	vector std(sage,nage); 			// std for length at age curve

	vector umsy(rep_yr,eyr);
	vector msy(rep_yr,eyr);
	vector utarget(rep_yr,eyr);
	vector ytarget(rep_yr,eyr);
	

	matrix P_al(sage,nage,1,nlen); 	// matrix of proportion of age at length
	matrix P_la(1,nlen,sage,nage); 	// matrix of proportion of age at length (transpose)
	matrix Nat(syr,eyr,sage,nage); 	// Numbers at age
	matrix Ulength(syr,eyr,1,nlen); // Exploitation rate at length
	matrix ObsUlength(syr,eyr,1,nlen); // Exploitation rate at length based on vulnerable biomass
	matrix Uage(syr,eyr,sage,nage); // Exploitation rate at length
	matrix Nlt(syr,eyr,1,nlen); 	// Numbers at length
	matrix Clt(syr,eyr,1,nlen); 	// Catch at length
	matrix obsClt(syr,eyr,1,nlen); 	// Observed catch at length



	
	
	LOC_CALCS

		//generate random variables
		random_number_generator rng(seed);
		wt.fill_randn(rng);
		wt*=sigR;


		eps.fill_randn(rng);
		eps*=tau;
		age.fill_seqadd(sage,1);
		len.fill_seqadd(slen,lstp);
		Am1=nage-1;	
		tiny=1.e-40;


	END_CALCS
	
PARAMETER_SECTION

	init_number atemoia;

	number fpen;
	objective_function_value no_f; 

PRELIMINARY_CALCS_SECTION	  
 	
	
PROCEDURE_SECTION


	
	incidence_functions();
	propAgeAtLengh();
	initialYear();
	populationDynamics();

	calc_msy();
	calc_Fspr_target(.4);

	output_data();
	output_ctl();
	output_true();
	exit(1);


FUNCTION incidence_functions

	sampleSigR=sqrt(((wt-mean(wt))*(wt-mean(wt)))/size_count(wt));
	
			

 	la.initialize();
 	std.initialize();

	la = Linf*(1.-mfexp(-k*(age-to)));  //average length at age
	std = la*cvl;  		  //std for length at age


	lxo(sage) = 1.; //first age	
	for(int a = sage+1 ; a <= nage ; a++)
	{
		lxo(a) = lxo(a-1)*Sa; // proportion of individuals at age surviving M only
	}
	lxo(nage) /= (1.-Sa); // age plus group

	wa = alw * pow(la,blw); //weight at age


 	for(int w=sage; w<=nage;w++)
 	{
 		if(wa(w)>wmat)
 		{
 			fec(w)=wa(w)-wmat;
 		}
 	}
 	

	phie = lxo*fec; 

	
	reca = reck/phie; 
	recb = (reck - 1.)/(Ro*phie); 
	//cout<<"recb is"<<recb<<endl;
	sbo  = Ro*phie;

	fpen = 0.;

	calcSellen();

	//cout<<"OK after incidence_functions"<<endl;
	
FUNCTION propAgeAtLengh

	dvector z1(1,nlen); 				// intermediate steps for calculating proportion of age at length
	dvector z2(1,nlen); 				// intermediate steps for calculating proportion of age at length

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

	cout<<"OK after propAgeAtLengh"<<endl;


FUNCTION initialYear
	
	Nat( syr, sage )= Ro;
	for( int a = sage+1; a <= nage; a++ )
	{
		Nat( syr, a ) = Nat( syr, a - 1 ) * Sa ;	// initial age-structure
	}
	Nat( syr, nage ) /= (1.-Sa);
	

	Nlt(syr) = Nat(syr)*P_al;

	//Exploitation rate at length
	

	Ulength(syr) = ut(syr)*sellen(indselyr(syr));
	
	//Ulength(syr) = ut(syr)/(1.+mfexp(-1.7*(len-Ulenmu)/Ulensd)); 
	


	//Exploitation rate at age
	//Uage(syr) = Ulength(syr)*P_la;

	for( int au = sage; au <= nage; au++ ){
			Uage( syr )(au) = Ulength( syr ) * P_al(au);
		}

	
	Clt(syr) = elem_prod(Nlt(syr),Ulength(syr));

	ObsUlength(syr) = elem_div(Clt(syr),Nlt(syr));

	maxUy(syr)=max(ObsUlength(syr));
	avgUy(syr)=mean(ObsUlength(syr)(1,nlen));

	addErrorClt(syr);

	vbt(syr) = q * Nat(syr)*elem_prod(wa,va) * mfexp((eps(syr))*obs_err); // cpue -tau*tau/2.
	//vbt(syr) = q * Nat(syr)*elem_prod(wa,va) * mfexp((eps(syr))*obs_err); // cpue
		
	// Add process error to all ages in initial year
	//bt(syr) = Nat(syr)* wa * mfexp((eps(syr))*obs_err); //-tau*tau/2.				     // survey
	
	//spawning biomass
	sbt(syr) = fec * Nat(syr);

	//cout<<"Nat(syr)"<<endl<<Nat(syr)<<endl;
	//cout<<"fec"<<endl<<fec<<endl;

	//spawning biomass depletion
	depl(syr) = sbt(syr)/sbo;

	cout<<"OK after initialYear"<<endl;



FUNCTION populationDynamics	


	int i;
	for(i=syr;i<=eyr-1;i++)
	{	   
	    
	   
	    //recruitment
	    //Nat(i+1,sage) = (reca*sbt(i)/(1.+recb*sbt(i)))*mfexp((wt(i+1)-sigR*sigR/2.)*proc_err);
	   
	    Nat(i+1,sage) = (reca*(sbt(i))/(1.+recb*(sbt(i))))*mfexp((wt(i+1)-sigR*sigR/2.)*proc_err);//-sigR*sigR/2.
	    //ages 2 -nage
	    //Nat(i+1)(sage+1,nage) = ++elem_prod(Nat(i)(sage,nage-1)*Sa,1.-Uage(i)(sage,nage-1));
	    
	    for( int a = sage +1; a <= nage; a++ ){
			Nat( i+1 )( a ) =  Nat( i )( a - 1 )* Sa * (1. - Uage( i)( a -1 ));
				
		}
		
		//Nat( i+1, nage ) += Nat( i, nage ) *Sa*  (1. - Uage( i)( nage ));
		Nat( i+1, nage ) /= (1. - ( Sa* (1.-Uage(i)(nage) ) ) );
		
		//Proportion of individuals at length 
		//note admb matrix multiplication is yj = \sum_i xi * mij 
		Nlt(i+1) = Nat(i+1)*P_al;
		//Explitation rate at length
		//Ulength(i+1) = ut(i+1)/(1.+mfexp(-1.7*(len-4.)/0.1)); //4 is length of 50% mat hard coded in
		//calcUlength(i+1,indselyr(i+1));
		Ulength(i+1) = ut(i+1)*sellen(indselyr(i+1));
		//exploitation rate at age
		//Uage(i+1) = Ulength(i+1)*P_la;
		for( int au = sage; au <= nage; au++ ){
			Uage( i+1 )(au) = Ulength( i+1 ) * P_al(au);
		}
	
		Clt(i+1) = elem_prod(Nlt(i+1),Ulength(i+1));
		addErrorClt(i+1);
		ObsUlength(i+1) = elem_div(Clt(i+1),Nlt(i+1));
		maxUy(i+1)=max(ObsUlength(i+1));
		avgUy(i+1)=mean(ObsUlength(i+1)(1,nlen));
		// Vulnerable biomass for survey
		vbt(i+1) = q * Nat(i+1)*elem_prod(wa,va) * mfexp((eps(i+1))*obs_err); // cpue -tau*tau/2.
		//vbt(i+1) = q * Nat(i+1)*elem_prod(wa,va) * exp((eps(i+1))*obs_err); // cpue
		
		//Total biomass - what is this additional obs error representing? 
		//bt(i+1) = Nat(i+1)* wa * exp((eps(i+1))*obs_err); 
		
		//spawning biomass
		sbt(i+1) = fec * Nat(i+1) ;				    
		//spawning biomass depletion
		depl(i+1) = sbt(i+1)/sbo;
	}
	//cout<<"Nat"<<endl<<Nat<<endl;

	//cout<<"OK after populationDynamics"<<endl;

FUNCTION void addErrorClt(const int& ii)

       		
			obsClt(ii) = rmvlogistic(Clt(ii),tau_length,seed+ii)* sum(Clt(ii));
			//obsClt(ii) = Clt(ii);
			//cout<<"Clt(ii) is"<<endl<<Clt(ii)<<endl;
			//cout<<"obsClt(ii) is"<<endl<<obsClt(ii)<<endl;
			
      	

FUNCTION  calcSellen

	int b, si;
	

	for(int si=1;si<=nselch;si++){
		for(int b=1;b<=nlen;b++){
			
			sellen(si)(b) = (1/(1-selg(si)))*
							pow((1-selg(si))/selg(si),selg(si))*
							((exp(sela(si)*selg(si)*(selb(si)-len(b))))/
							(1+exp(sela(si)*(selb(si)-len(b)))));
			//sellen(si)(b) = 1.;///(1.+mfexp(-1.7*(len(b)-4.)/0.1));
		}
	}

	
	
	///(1.+mfexp(-1.7*(len-4.)/0.1));

FUNCTION calc_msy

	dvector utest(1,1001);
	utest.fill_seqadd(0,0.001);

 //This function calculates MSY in the lazy and slow way. 
 	int k, kk ;
	int NF=size_count(utest);
	
	dmatrix selage(rep_yr,eyr,sage,nage);

	

	for(int y=rep_yr+1;y<=eyr;y++){

		selage(y)= (Uage(y)/max(Uage(y)));

		//max((Uage(y)/mean(Uage(y))));


		dvector ye(1,NF);
		ye.initialize();
		
		for(k=1; k<=NF; k++)
		{
			dvector lz(sage,nage);
			dvariable phieq;
			dvariable phiz;
			dvariable req;
			
			lz.initialize();
			phieq.initialize();
			phiz.initialize();
			req.initialize();
			
			lz(sage) = 1.; //first age	
			for(int a = sage+1 ; a <= nage ; a++)
			{
				lz(a) = lz(a-1)*Sa*(1.-utest(k)*selage(y)(a-1)); // proportion of individuals at age surviving M only
			}
			lz(nage) /= (1.-Sa*(1.-utest(k)*selage(y)(nage))); // age plus group
			phiz= lz*fec;
			
			phieq = elem_prod(lz,selage(y))*wa;
			req = (Ro*(reck-phie/phiz)/(reck-1.));
			
			ye(k)= value(utest(k)*req*phieq);
		}

		
		msy(y)= max(ye);
		double mtest;	

		for(kk=1; kk<=NF; kk++)
		{
			mtest=ye(kk);
				
			if(mtest==msy(y)){
				umsy(y)=mean(utest(kk)*selage(y));
			} 
		}


	}


FUNCTION void calc_Fspr_target( double target )

  //This function calculates MSY in the lazy and slow way. 
 

	
	dvector ftest(1,1001);
	ftest.fill_seqadd(0,0.001);
 	int k, kk, a;
	int NF=size_count(ftest);

	dmatrix selage(rep_yr,eyr,sage,nage);




	for(int y=rep_yr+1;y<=eyr;y++){

		selage(y)= (Uage(y)/mean(Uage(y)));


		//selage= seltotal(ii);
		dvector tmp_phiz(1,NF);
		dvector tmp_target(1,NF);
		dvector ye(1,NF);
		
		tmp_phiz.initialize();
		tmp_target.initialize();
		ye.initialize();


		for(k=1; k<=NF; k++)
		{

	
			dvector lzt(sage,nage);
			dvariable phieq;
			dvariable req;
			
			lzt.initialize();
			phieq.initialize();
			req.initialize();

			lzt(sage)=1.0;

		
			for(int a = sage+1 ; a <= nage ; a++)
			{
				lzt(a) = lzt(a-1)*Sa*(mfexp(-ftest(k)*selage(y)(a-1))); // proportion of individuals at age surviving M only
			}
			
			lzt(nage) /= (1.-Sa*(mfexp(-ftest(k)*selage(y)(nage)))); // age plus group
		
			tmp_phiz(k) = lzt*fec;
		
			tmp_target(k) = fabs((tmp_phiz(k)/phie - target));

			phieq = elem_prod(lzt,(1.-mfexp(-ftest(k)*selage(y))))*wa;

			req = Ro*(reck-phie/tmp_phiz(k))/(reck-1.);
			
			ye(k)= value(req*phieq);
			
		}
		
		double ttest;

		ttest =  min(tmp_target);

		for(kk=1; kk<=NF; kk++)
		{
			if(tmp_target(kk)==ttest){
				utarget(y)=1.-mfexp(-ftest(kk));
				ytarget(y)=ye(kk);
			} 
		}
	}
				
		




FUNCTION output_ctl
	
	dvariable rb;
	rb=mean(column(Nat,sage));

	dvariable rini;
	rini=Nat(rep_yr,sage) / mfexp( wt(rep_yr)-sigR*sigR/2.); //

	dvector iwti(rep_yr+1,eyr);
	iwti.initialize();

	dvector iwtinit(rep_yr-(nage-sage),rep_yr);
	iwtinit.initialize();


	ofstream mfs("../SA/LSRA.ctl");
	mfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	mfs<<"## CONTROL FILE TEMPLATE                                                                ##"<< endl;
	mfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	mfs<<"## ------------------------------------------------------------------------------------ ##"<< endl;
	mfs<<"## CONTROLS FOR LEADING PARAMETERS                                                      ##"<< endl;
	mfs<<"##  Prior descriptions:                                                                 ##"<< endl;
	mfs<<"##                      -0 uniform      (0,0)                                           ##"<< endl;
	mfs<<"##                      -1 normal       (p1=mu,p2=sig)                                  ##"<< endl;
	mfs<<"##                      -2 lognormal    (p1=log(mu),p2=sig)                             ##"<< endl;
	mfs<<"##                      -3 beta         (p1=alpha,p2=beta)                              ##"<< endl;
	mfs<<"##                      -4 gamma        (p1=alpha,p2=beta)                              ##"<< endl;	
	mfs<<"##                      -5 no prior        pvec=0.0                              ##"<< endl;	
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"## npar"<<endl<< "3"<< endl;
	//mfs<<"## npar"<<endl<< "4"<< endl;
	mfs<<"## ival         		lb      	ub        phz     prior   p1      p2        #parameter            ##"<< endl;
	mfs<<"## ——————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<< log(Ro*0.8) <<"\t"<< 3.0 <<"\t"<< 7.0   <<"\t"<<  1  <<"\t"<< 0 <<"\t"<< 3.0	<<"\t"<< 7.0   	<<"\t"<<"#log_ro   	##"<<endl;
	//mfs<< log(rini*1.1)  	 <<"\t"<< 3.0 <<"\t"<< 7.0   <<"\t"<<  1  <<"\t"<< 0  <<"\t"<< 3.0 	<<"\t"<< 7.0   	<<"\t"<<"#log_rinit   	##"<<endl;
	//mfs<< log(Ro) <<"\t"<< 3.0 <<"\t"<< 7.0   <<"\t"<<  1  <<"\t"<< 5  <<"\t"<< 3.0	<<"\t"<< 7.0   	<<"\t"<<"#log_ro   	##"<<endl;
   	//mfs<< log(rini)  	 <<"\t"<< 3.0 <<"\t"<< 7.0   <<"\t"<<  1  <<"\t"<< 5  <<"\t"<< 3.0 	<<"\t"<< 7.0   	<<"\t"<<"#log_rinit   	##"<<endl;
   	//mfs<< log(reck*0.8) 	 <<"\t"<<  1.6 <<"\t"<< 4.0   <<"\t"<<  1  <<"\t"<< 0  <<"\t"<<  1.6 	<<"\t"<< 4.0  	<<"\t"<<"#log_reck  ##"<<endl;
   	mfs<< log(reck*0.8)  	  <<"\t"<<  1.6 <<"\t"<< 5.0   <<"\t"<<  1  <<"\t"<< 1  <<"\t"<<  log(reck)	<<"\t"<< 0.5	<<"\t"<<"#log_reck  ##"<<endl;
   	//mfs<< log(reck*0.8) 	  <<"\t"<<  1.6 <<"\t"<< 5.0   <<"\t"<<  1  <<"\t"<< 1  <<"\t"<<  log(reck)	<<"\t"<< 0.9  	<<"\t"<<"#log_reck  ##"<<endl;
   	//mfs<< log(8) 	 	 <<"\t"<<  1.0 <<"\t"<< 5.0   <<"\t"<<  1  <<"\t"<< 1  <<"\t"<<  log(reck)	<<"\t"<< 0.8 	<<"\t"<<"#log_reck  ##"<<endl;
   	//mfs<< log(Linf)   <<"\t"<< 1.3  <<"\t"<< 4.0   <<"\t"<<  -3  <<"\t"<< 0  <<"\t"<<  1.3 	<<"\t"<< 4.0 	<<"\t"<<"#log_Linf  ##"<<endl;
   	//mfs<< log(k)  <<"\t"<< -3.0 <<"\t"<< -0.2  <<"\t"<<  -3  <<"\t"<< 0  <<"\t"<< -3.0 	<<"\t"<< -0.2  	<<"\t"<<"#log_k  	##"<<endl;
   	//mfs<< to  	<<"\t"<< -2.0 <<"\t"<< 0.0   <<"\t"<<   -4  <<"\t"<< 0  <<"\t"<< -2.0 	<<"\t"<<  0.0  	<<"\t"<<"#to 	##"<<endl;
   	//mfs<< log(cvl)  <<"\t"<< -7.0 <<"\t"<< -0.1  <<"\t"<< 	-3  <<"\t"<< 0  <<"\t"<< -7.0 	<<"\t"<< -0.1	<<"\t"<<"#log_cvl   ##"<<endl;
   	mfs<< log(sigR) 	 <<"\t"<< -3.0 <<"\t"<< 8.0  <<"\t"<< 	-3  <<"\t"<< 5 <<"\t"<< -3.0 	<<"\t"<< 8.0	<<"\t"<<"#log_sigR   ##"<<endl;
   // mfs<< log(sigR) 	 <<"\t"<< -3.0 <<"\t"<< 8.0  <<"\t"<< 	-3  <<"\t"<< 5 <<"\t"<< -3.0 	<<"\t"<< 8.0	<<"\t"<<"#log_sigR_init   ##"<<endl;   
    //mfs<< log(sigVul) 	 <<"\t"<< -4.0 <<"\t"<< 15.0  <<"\t"<< 	3  <<"\t"<< 0  <<"\t"<< -4.0 	<<"\t"<< 15.0	<<"\t"<<"#log_sigR   ##"<<endl;
    //mfs<< log(tau)  <<"\t"<< -7.0 <<"\t"<< 8.0 <<"\t"<< 	2  <<"\t"<< 0  <<"\t"<< -7.0 	<<"\t"<< 8.0	<<"\t"<<"#log_cv_it   ##"<<endl;
    mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"##initial values for recruitment deviations ##"<< endl;
	//mfs<<"# wt "<< endl << mfexp(wt(rep_yr+1,eyr)*proc_err) <<endl;
	//mfs<<"# wt "<< endl << (wt(rep_yr-(nage-sage),eyr-(nage-sage+1)-1)*0.0) <<endl<<wt(eyr-(nage-sage+1),eyr)*0.0<< endl;
	mfs<<"# wt "<< endl << iwtinit << endl;
	mfs<<"# wt "<< endl << iwti << endl;
	//mfs<<"# wt "<< endl << wt(rep_yr-(nage-sage),eyr)<< endl;
	//mfs<<"# wt_init "<< endl << exp(wt(rep_yr-(nage-sage),rep_yr-1)*0) <<endl;
	mfs<<"##initial values for recruitment deviations in first year ##"<< endl;
	mfs<<"#log(q) prior - same codes as above##"<< endl;
	mfs<<"#prior   p1      p2  ##"<< endl;
  	mfs<< 5 <<"\t"<<	   0	 <<"\t"<<	 0.2 << endl;
	mfs<<"#closed loop ##"<<endl<<0<< endl;
	mfs<<"# eof " << endl << 999 <<endl;

	//cout<<"OK after otput_ctl"<<endl;
	  
FUNCTION output_data
	

	ofstream ofs("../SA/LSRA.dat");
	ofs<<"# syr " << endl << rep_yr <<endl;
	ofs<<"# eyr " << endl << eyr <<endl;
	ofs<<"# sage "<< endl << sage <<endl;
	ofs<<"# nage "<< endl << nage <<endl;
	ofs<<"# slen "<< endl << slen <<endl;
	ofs<<"# nlen "<< endl << nlen <<endl;
	ofs<<"# lstp "<< endl << lstp <<endl;
	ofs<<"# SR function " << endl << 1 <<endl;
	ofs<<"# m " << endl << -log(Sa) <<endl;
	ofs<<"# alw " << endl << alw <<endl;
	ofs<<"# blw "<< endl << blw <<endl;
	ofs<<"# waobs "<< endl << wa <<endl;
	ofs<<"# wasw "<< endl << 0 <<endl;

	//ofs<<"# mat50  "<< endl << feca <<endl;
	//ofs<<"# matsd " << endl << fecg <<endl;
	//ofs<<"# ahat " << endl << ahat <<endl;
	//ofs<<"# ghat "<< endl << ghat <<endl;
	ofs<<"# vul "<< endl << va <<endl;
	ofs<<"# fec "<< endl << fec <<endl;
	ofs<<"# nyt "<< endl << niyr <<endl;
	ofs<<"# iyr " << endl << iyr <<endl;
	ofs<<"# yt " << endl << vbt(rep_yr,eyr)  <<endl;
	ofs<<"# Clt"<< endl << obsClt.sub(rep_yr,eyr) <<endl;
	//ofs<<"# Clt"<< endl << Clt.sub(rep_yr,eyr) <<endl;
	ofs<<"# linf "<< endl << Linf <<endl;//     
	ofs<<"# k "<< endl << k <<endl;
	ofs<<"# to " << endl << to <<endl;
	ofs<<"# cvl " << endl << cvl <<endl;
	//ofs<<"# ireck "<< endl << reck <<endl;
	//ofs<<"# iRo "<< endl << Ro <<endl;

	ofs<<"# cv_it " << endl << tau <<endl;
	//ofs<<"# sigR " << endl << sigR <<endl;
	ofs<<"# sigVul " << endl << sigVul <<endl;
	//ofs<<"# phz_reck "<< endl << 2 <<endl;
	//ofs<<"# phz_growth  "<< endl << -4  <<endl;
	//ofs<<"# use_prior  "<< endl << 0 <<endl;
	ofs<<"# u_init " << endl << 0.0 <<endl; //ut(rep_yr-1)
	//ofs<<"# P_al " << endl << P_al <<endl;
	ofs<<"# eof " << endl << 999 <<endl;

	//cout<<"OK after output_dat"<<endl;
	
	
FUNCTION output_true
	  
	ofstream ofs("true_data_lsra.rep");

	double tRbar;

	for(int ni=rep_yr;ni<=eyr;ni++){

		
		tRbar += (Nat(ni)(sage))/mfexp((wt(ni))*proc_err);//-sigR*sigR/2.

	}
	tRbar /= (eyr-rep_yr+1);
	//	/mfexp(wt(rep_yr,eyr)*proc_err))/(eyr-rep_yr+1);

	ofs<<"scnNumber" << endl << scnNumber <<endl;
	ofs<<"ut" << endl << ut <<endl;
	ofs<<"Nat" << endl << Nat <<endl;
	ofs<<"Nlt" << endl << Nlt <<endl;
	ofs<<"Clt" << endl << Clt.sub(syr,eyr) <<endl;
	ofs<<"Ro" << endl << Ro <<endl;
	ofs<<"Rbar" << endl << tRbar <<endl;	
	ofs<<"Rinit" << endl << Nat(rep_yr)(sage) /mfexp((wt(rep_yr)-sigR*sigR/2.)*proc_err)<<endl;
	ofs<<"reck" << endl << reck <<endl;
	ofs<<"cvl" << endl << cvl <<endl;
	ofs<<"cv_it" << endl << tau <<endl;
	ofs<<"Linf" << endl << Linf <<endl;
	ofs<<"k" << endl << k <<endl;
	ofs<<"to" << endl << to <<endl;
	ofs<<"sbt" << endl << sbt <<endl;
	ofs<<"depl" << endl << depl <<endl;
	ofs<<"q" << endl << q <<endl;	
	ofs<<"sellen" << endl << sellen <<endl;	
	//ofs<<"selage" << endl << selage <<endl;
	ofs<<"syr" << endl << syr <<endl;
	ofs<<"eyr" << endl << eyr <<endl;
	ofs<<"rep_yr" << endl << rep_yr <<endl;	
	//ofs<<"wt" << endl << ((wt(rep_yr+1,eyr))*proc_err) <<endl;	//-sigR*sigR/2.
	ofs<<"wt" << endl << ((wt(rep_yr-(nage-sage),eyr))*proc_err) <<endl;
	//ofs<<"wt_init" << endl << wt(rep_yr-(nage-sage),rep_yr-1)*proc_err <<endl;	//-sigR*sigR/2.
	//ofs<<"expwt_init" << endl << exp(wt(rep_yr-(nage-sage),rep_yr-1)*proc_err) <<endl;
	ofs<<"lxo" << endl << lxo <<endl;
	ofs<<"fec" << endl << fec <<endl;	
	ofs<<"wa" << endl << wa <<endl;	
	ofs<<"la" << endl << la <<endl;	
	ofs<<"umsy" << endl << umsy<<endl;	
	ofs<<"msy" << endl << msy<<endl;
	ofs<<"utarget" << endl << utarget<<endl;
	ofs<<"ytarget" << endl << ytarget<<endl;
	ofs<<"Ulength" << endl << Ulength <<endl;	
	ofs<<"ObsUlength"<<endl << ObsUlength<< endl;
	ofs<<"maxUy"<<endl << maxUy<< endl;
	ofs<<"avgUy"<<endl << avgUy<< endl;
	ofs<<"reca"<<endl << reca<< endl;
	ofs<<"recb"<<endl << recb<< endl;
	ofs<<"phie"<<endl << phie << endl;
	ofs<<"it " << endl << vbt(rep_yr,eyr)  <<endl;
	ofs<<"P_al " << endl <<P_al  <<endl;
	ofs<<"sampleSigR " << endl <<sampleSigR  <<endl;



