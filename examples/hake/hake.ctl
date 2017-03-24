## -------------------------------------------------------------------------- ##
## CONTROL FILE TEMPLATE                                                   ##
## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##
## CONTROLS FOR LEADING PARAMETERS                                         ##
##  Prior descriptions:                                                    ##
##                      -0 uniform      (0,0)                              ##
##                      -1 normal       (p1=mu,p2=sig)                     ##
##                      -2 lognormal    (p1=log(mu),p2=sig)                ##
##                      -3 beta         (p1=alpha,p2=beta)                 ##
##                      -4 gamma        (p1=alpha,p2=beta)                 ##
## ——————————————————————————————————————————————————————————————————————— ##
## npar
5
## ival    	lb      ub        phz     prior   p1      p2      #parameter   ##
## ——————————————————————————————————————————————————————————————————————— ##
7.0			5	10	 1	0	5	10	#log_ro   	##
6.2 		5	10	 1	0	5	10	#log_rinit   	##
#3.5 		1.3	4.5	 1	0	1.3	4.5	#log_reck  ##
3.0 		1.3	4.5	 2	1	3.3	0.5	#log_reck  ##
#0.3364722	-3	1	-3	0	-3	1	#log_SigR
-0.2231436	-3	1	-3	0	-3	1	#log_SigR
-6.907755	-10	0	3	0	-10	0	#log_can
## ——————————————————————————————————————————————————————————————————————— ##
##initial values for recruitment deviations ##
# wt
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1
#log(q) prior - same codes as above
#prior   p1      p2 
1		   0		 0.15
##initial values for recruitment deviations in first year ##
999