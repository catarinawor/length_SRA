##------------------------------------------------------------------------ ##
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
## ival		lb      ub      phz     prior 	 p1      p2        #parameter  ##
## ——————————————————————————————————————————————————————————————————————— ##
8.0			5.0		35	 	 1		0	5.0		35.0		#log_ro   	  ##
8.0 		5.0		35.0	 1		0	5.0		35.0	#log_rinit   ##
1.7 		1.0		6.0	 	 1		0	1.0		6.0		#log_reck    	  ##
#1.7 		1.0		6.0	 	 1		1	2.0		0.4		#log_reck    	  ##
-0.8209806	-3		1		-3		0	-3		1	#log_sigR   ##
#-0.4462871 	-3		1		-3		0	-3		1	#log_sigR   ##
-4.60517	-6	0	-3	0	-6	0	#log_can  ##
## ———————————————————————————————————————————————————————————————————————##
##initial values for recruitment deviations ##
# wt
1 1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1 1 
1 1 1 1 
#log(q) prior - same codes as above
#prior   p1      p2 
0		  -10	 20
##initial values for recruitment deviations in first year ##
999