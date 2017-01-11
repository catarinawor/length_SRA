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
		seed += 13; // add 10 to the seed
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
	//init_int slen;
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
	vector bt(syr,eyr);	  			// total biomass
	vector sbt(syr,eyr); 			// spawning biomass
	vector depl(syr,eyr);			// depletion based on spawning biomass
	vector maxUy(syr,eyr);   		// annual U based on vulberable biomass
	
	vector len(1,nlen); 			// length classes
	vector lxo(sage,nage); 			// Equilibrium unfished numbers ata age
	vector wa(sage,nage); 				// Weight at age
	vector fec(sage,nage); 			// fecundity at age
	vector la(sage,nage); 				// length at age
	vector std(sage,nage); 			// std for length at age curve

	vector umsy(syr,eyr);
	vector msy(syr,eyr);
	

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
		len.fill_seqadd(lstp,lstp);
		Am1=nage-1;	
		tiny=1.e-20;


	END_CALCS
	
PARAMETER_SECTION

	number fpen;
	objective_function_value no_f; 

PRELIMINARY_CALCS_SECTION	  
 	
	incidence_functions();
	propAgeAtLengh();
	initialYear();
	populationDynamics();

	calc_msy();

	output_data();
	output_ctl();
	output_true();
	exit(1);
	
PROCEDURE_SECTION

FUNCTION incidence_functions
	
			

 	la.initialize();
 	std.initialize();

	la = Linf*(1.-exp(-k*(age-to)));  //average length at age
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
	sbo  = Ro*phie;

	fpen = 0.;

	calcSellen();

	//cout<<"OK after incidence_functions"<<endl;
	
FUNCTION propAgeAtLengh

	dvector z1(1,nlen); 				// intermediate steps for calculating proportion of age at length
	dvector z2(1,nlen); 				// intermediate steps for calculating proportion of age at length
	

 	for( int a = sage; a <= nage; a++ )
	{
		z1 = (( len - lstp * 0.5 )-la( a ))/std( a );
		z2 = (( len + lstp * 0.5 )-la( a ))/std( a );
		for( int b=1; b<= nlen; b++ )
		{
			P_al( a, b )=cumd_norm( z2( b ))-cumd_norm( z1( b )); // calculates the proportion of a given age given your length
		}
	}
	
	P_la = trans(P_al); //transpose matrix to length by age

	cout<<"OK after propAgeAtLengh"<<endl;


FUNCTION initialYear
	
	Nat( syr, sage )= Ro;//*exp((wt(syr)-sigR*sigR/2)*proc_err);
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
	Uage(syr) = Ulength(syr)*P_la;
	
	Clt(syr) = elem_prod(Nlt(syr),Ulength(syr));

	ObsUlength(syr) = elem_div(Clt(syr),elem_prod(Nlt(syr),sellen(indselyr(syr))));

	maxUy(syr)=max(ObsUlength(syr));

	addErrorClt(syr);

	vbt(syr) = q * Nat(syr)*elem_prod(wa,va) * mfexp((eps(syr)-tau*tau/2.)*obs_err); // cpue
	//vbt(syr) = q * Nat(syr)*elem_prod(wa,va) * mfexp((eps(syr))*obs_err); // cpue
		
	// Add process error to all ages in initial year
	bt(syr) = Nat(syr)* wa * mfexp((eps(syr)-tau*tau/2.)*obs_err); 				     // survey
	
	//spawning biomass
	sbt(syr) = fec * Nat(syr);

	//cout<<"Nat(syr)"<<endl<<Nat(syr)<<endl;
	//cout<<"fec"<<endl<<fec<<endl;

	//spawning biomass depletion
	depl(syr) = sbt(syr)/sbt(syr);

	cout<<"OK after initialYear"<<endl;



FUNCTION populationDynamics	


	int i;
	for(i=syr;i<=eyr-1;i++)
	{	   
	    
	   
	    //recruitment
	    //Nat(i+1,sage) = (reca*sbt(i)/(1.+recb*sbt(i)))*mfexp((wt(i+1)-sigR*sigR/2.)*proc_err);
	   
	    Nat(i+1,sage) = (reca*sbt(i)/(1.+recb*sbt(i)))*mfexp((wt(i+1))*proc_err);
	    
	    //ages 2 -nage
	    Nat(i+1)(sage+1,nage) = ++elem_prod(Nat(i)(sage,nage-1)*Sa,1.-Uage(i)(sage,nage-1));
	    
		Nat( i+1, nage ) /= (1. - ( Sa* (1.-Uage(i)(nage) ) ) );
		
		//Proportion of individuals at length 
		//note admb matrix multiplication is yj = \sum_i xi * mij 
		Nlt(i+1) = Nat(i+1)*P_al;

		//Explitation rate at length
		//Ulength(i+1) = ut(i+1)/(1.+mfexp(-1.7*(len-4.)/0.1)); //4 is length of 50% mat hard coded in
		//calcUlength(i+1,indselyr(i+1));
		Ulength(i+1) = ut(i+1)*sellen(indselyr(i+1));

		//exploitation rate at age
		Uage(i+1) = Ulength(i+1)*P_la;
		Clt(i+1) = elem_prod(Nlt(i+1),Ulength(i+1));

		addErrorClt(i+1);

		ObsUlength(i+1) = elem_div(Clt(i+1),Nlt(i+1));

		maxUy(i+1)=max(ObsUlength(i+1));

		// Vulnerable biomass for survey
		vbt(i+1) = q * Nat(i+1)*elem_prod(wa,va) * mfexp((eps(i+1)-tau*tau/2.)*obs_err); // cpue
		//vbt(i+1) = q * Nat(i+1)*elem_prod(wa,va) * exp((eps(i+1))*obs_err); // cpue
		

		//Total biomass - what is this additional obs error representing? 
		bt(i+1) = Nat(i+1)* wa * exp((eps(i+1)-tau*tau/2.)*obs_err); 
		
		//spawning biomass
		sbt(i+1) = fec * Nat(i+1);				     // survey
		//spawning biomass depletion
		depl(i+1) = sbt(i+1)/sbt(syr);
	}
	//cout<<"Nat"<<endl<<Nat<<endl;

	cout<<"OK after populationDynamics"<<endl;

FUNCTION void addErrorClt(const int& ii)

       		//dvector ppl(1,nlen);
       		//ppl.initialize();
       		
			obsClt(ii) = rmvlogistic(Clt(ii),tau_length,seed+ii);
      	

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

	dvector utest(1,101);
	utest.fill_seqadd(0,0.01);

 //This function calculates MSY in the lazy and slow way. 
 	int k, kk ;
	int NF=size_count(utest);
	
	dmatrix selage(syr,eyr,sage,nage);

	
	
	

	for(int y=syr;y<=eyr;y++){

		selage(y)= Uage(y)/max(Uage(y));

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
				lz(a) = lz(a-1)*Sa*(1-utest(k)); // proportion of individuals at age surviving M only
			}
			lz(nage) /= (1.-Sa*(1-utest(k))); // age plus group

			phiz= lz*fec;
			
			phieq = elem_prod(lz,selage(y))*wa;
			req = Ro*(reck-phie/phiz)/(reck-1);
			
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


FUNCTION output_ctl
	

	ofstream mfs("../SA/length_sra.ctl");
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
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	//mfs<<"## npar"<<endl<< "7"<< endl;
	mfs<<"## npar"<<endl<< "6"<< endl;
	mfs<<"## ival         		lb      	ub        phz     prior   p1      p2        #parameter            ##"<< endl;
	mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<< log(1.50)  		 <<"\t"<< -4.0 <<"\t"<< 8.0   <<"\t"<<  1  <<"\t"<< 0  <<"\t"<< -4.0	<<"\t"<< 8.0   	<<"\t"<<"#log_ro   	##"<<endl;
	//mfs<< log(Ro)  		 <<"\t"<< -2 <<"\t"<< 8.0   <<"\t"<<  1  <<"\t"<< 0  <<"\t"<< -2	<<"\t"<< 8.0   	<<"\t"<<"#log_ro   	##"<<endl;
	//mfs<< 0.0  		 <<"\t"<< -4.0 <<"\t"<< 4.0   <<"\t"<<  1  <<"\t"<< 0  <<"\t"<< -4.0 	<<"\t"<< 4.0   	<<"\t"<<"#log_rbar   	##"<<endl;
   	//mfs<< 0.0  	 	 <<"\t"<< -4.0 <<"\t"<< 4.0   <<"\t"<<  1  <<"\t"<< 1  <<"\t"<< 0.0 	<<"\t"<< 0.5   	<<"\t"<<"#log_rinit   	##"<<endl;
   	mfs<< log(8) 	 <<"\t"<<  0.0 <<"\t"<< 4.0   <<"\t"<<  1  <<"\t"<< 0  <<"\t"<<  0.0 	<<"\t"<< 4.0  	<<"\t"<<"#log_reck  ##"<<endl;
   	//mfs<< log(reck) 	 <<"\t"<<  0.0 <<"\t"<< 4.0   <<"\t"<<  1  <<"\t"<< 0  <<"\t"<<  0.0 	<<"\t"<< 4.0  	<<"\t"<<"#log_reck  ##"<<endl;
   	mfs<< log(Linf)   <<"\t"<< 1.3  <<"\t"<< 4.0   <<"\t"<<  -3  <<"\t"<< 0  <<"\t"<<  1.3 	<<"\t"<< 4.0 	<<"\t"<<"#log_Linf  ##"<<endl;
   	mfs<< log(k)  <<"\t"<< -3.0 <<"\t"<< -0.2  <<"\t"<<  -3  <<"\t"<< 0  <<"\t"<< -3.0 	<<"\t"<< -0.2  	<<"\t"<<"#log_k  	##"<<endl;
   	mfs<< to  	<<"\t"<< -2.0 <<"\t"<< 0.0   <<"\t"<<   -4  <<"\t"<< 0  <<"\t"<< -2.0 	<<"\t"<<  0.0  	<<"\t"<<"#to 	##"<<endl;
   	mfs<< log(cvl)  <<"\t"<< -7.0 <<"\t"<< -0.1  <<"\t"<< 	-4  <<"\t"<< 0  <<"\t"<< -7.0 	<<"\t"<< -0.1	<<"\t"<<"#log_cvl   ##"<<endl;
    mfs<<"## ———————————————————————————————————————————————————————————————————————————————————— ##"<< endl;
	mfs<<"##initial values for recruitment deviations ##"<< endl;
	//mfs<<"# wt "<< endl << exp(wt(rep_yr+1,eyr)) <<endl;
	
	mfs<<"# wt "<< endl << exp(wt(rep_yr+1,eyr-(nage-sage+1)-1)) <<endl<<exp(wt(eyr-(nage-sage+1),eyr)*0.)<< endl;
	mfs<<"##initial values for recruitment deviations in first year ##"<< endl;
	//mfs<<"# wt_init "<< endl << exp(wt(rep_yr-(nage-sage),rep_yr-1)) <<endl;

	//cout<<"OK after otput_ctl"<<endl;
	  
FUNCTION output_data
	

	ofstream ofs("../SA/length_sra.dat");
	ofs<<"# syr " << endl << rep_yr <<endl;
	ofs<<"# eyr " << endl << eyr <<endl;
	ofs<<"# sage "<< endl << sage <<endl;
	ofs<<"# nage "<< endl << nage <<endl;
	ofs<<"# nlen "<< endl << nlen <<endl;
	ofs<<"# lstp "<< endl << lstp <<endl;
	ofs<<"# SR function " << endl << 1 <<endl;
	ofs<<"# m " << endl << -log(Sa) <<endl;
	ofs<<"# alw " << endl << alw <<endl;
	ofs<<"# blw "<< endl << blw <<endl;
	//ofs<<"# mat50  "<< endl << feca <<endl;
	//ofs<<"# matsd " << endl << fecg <<endl;
	//ofs<<"# ahat " << endl << ahat <<endl;
	//ofs<<"# ghat "<< endl << ghat <<endl;
	ofs<<"# vul "<< endl << va <<endl;
	ofs<<"# fec "<< endl << fec <<endl;
	ofs<<"# nyt "<< endl << niyr <<endl;
	ofs<<"# iyr " << endl << iyr <<endl;
	ofs<<"# yt " << endl << vbt(rep_yr,eyr)  <<endl;
	//ofs<<"# Clt"<< endl << obsClt.sub(rep_yr,eyr) <<endl;
	ofs<<"# Clt"<< endl << Clt.sub(rep_yr,eyr) <<endl;
	//ofs<<"# ilinf "<< endl << Linf <<endl;
	//ofs<<"# ik "<< endl << k <<endl;
	//ofs<<"# it0 " << endl << to <<endl;
	//ofs<<"# icvl " << endl << cvl <<endl;
	//ofs<<"# ireck "<< endl << reck <<endl;
	//ofs<<"# iRo "<< endl << Ro <<endl;

	ofs<<"# cv_it " << endl << tau <<endl;
	ofs<<"# sigR " << endl << sigR <<endl;
	ofs<<"# sigVul " << endl << sigVul <<endl;
	//ofs<<"# phz_reck "<< endl << 2 <<endl;
	//ofs<<"# phz_growth  "<< endl << -4  <<endl;
	//ofs<<"# use_prior  "<< endl << 0 <<endl;
	ofs<<"# u_init " << endl << 0.0 <<endl;
	ofs<<"# eof " << endl << 999 <<endl;

	cout<<"OK after otput_dat"<<endl;
	
	
FUNCTION output_true
	  
	ofstream ofs("true_data_lsra.rep");

	double tRbar;

	for(int ni=rep_yr;ni<=eyr;ni++){

		
		tRbar += (Nat(ni)(sage))/mfexp(wt(ni)*proc_err);

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
	ofs<<"Rinit" << endl << Nat(rep_yr)(sage)/mfexp(wt(rep_yr)*proc_err) <<endl;
	ofs<<"reck" << endl << reck <<endl;
	ofs<<"Linf" << endl << Linf <<endl;
	ofs<<"k" << endl << k <<endl;
	ofs<<"to" << endl << to <<endl;
	ofs<<"cvl" << endl << cvl <<endl;
	ofs<<"sbt" << endl << sbt <<endl;
	ofs<<"depl" << endl << depl <<endl;
	ofs<<"q" << endl << q <<endl;	
	ofs<<"sellen" << endl << sellen <<endl;	
	ofs<<"syr" << endl << syr <<endl;
	ofs<<"eyr" << endl << eyr <<endl;
	ofs<<"rep_yr" << endl << rep_yr <<endl;	
	ofs<<"wt" << endl << mfexp(wt(rep_yr+1,eyr)*proc_err) <<endl;	
	ofs<<"lxo" << endl << lxo <<endl;
	ofs<<"fec" << endl << fec <<endl;	
	ofs<<"wa" << endl << wa <<endl;	
	ofs<<"umsy" << endl << umsy<<endl;	
	ofs<<"msy" << endl << msy<<endl;
	ofs<<"Ulength" << endl << Ulength <<endl;	
	ofs<<"ObsUlength"<<endl << ObsUlength<< endl;
	ofs<<"maxUy"<<endl << maxUy<< endl;
	ofs<<"reca"<<endl << reca<< endl;
	ofs<<"recb"<<endl << recb<< endl;
	ofs<<"phie"<<endl << phie << endl;






	



	

