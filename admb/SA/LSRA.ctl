## ------------------------------------------------------------------------------------ ##
## CONTROL FILE TEMPLATE                                                                ##
## ------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------ ##
## CONTROLS FOR LEADING PARAMETERS                                                      ##
##  Prior descriptions:                                                                 ##
##                      -0 uniform      (0,0)                                           ##
##                      -1 normal       (p1=mu,p2=sig)                                  ##
##                      -2 lognormal    (p1=log(mu),p2=sig)                             ##
##                      -3 beta         (p1=alpha,p2=beta)                              ##
##                      -4 gamma        (p1=alpha,p2=beta)                              ##
##                      -5 no prior        pvec=0.0                              ##
## ———————————————————————————————————————————————————————————————————————————————————— ##
## npar
3
## ival         		lb      	ub        phz     prior   p1      p2        #parameter            ##
## ——————————————————————————————————————————————————————————————————————————————————— ##
4.38203	3	7	1	0	3	7	#log_ro   	##
2.07944	1.6	4	1	1	2.30259	0.5	#log_reck  ##
-0.510826	-3	8	-3	5	-3	8	#log_sigR   ##
## ———————————————————————————————————————————————————————————————————————————————————— ##
##initial values for recruitment deviations ##
# wt 
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# wt 
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
##initial values for recruitment deviations in first year ##
#log(q) prior - same codes as above##
#prior   p1      p2  ##
1	0	0.5
#closed loop ##
0
# eof 
999
