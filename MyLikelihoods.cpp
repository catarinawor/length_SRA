#include<admodel.h>

dvariable dunif(const dvariable& x, const dvariable& pmin, const dvariable& pmax)
{
	dvariable p=log(1./(pmax-pmin));
	if(x<pmin||x>pmax)
	{
		p=-1e20;
	}
	return(p);
}

dvariable dnorm(dvector obs, dvar_vector mu, dvariable sig)
{
	int n=size_count(obs);
	return(n*log(sig)+norm2(obs-mu)/(2.0*sig*sig));

}

dvariable dnorm(dvar_vector res, dvariable mu, dvariable sig)
{
	int n=size_count(res);
	return(n*log(sig)+norm2(res-mu)/(2.0*sig*sig));
}

dvariable dnorm(const double& x, const dvariable& mu, const dvariable& sig)
{
	return(log(sig)+square(x-mu)/(2.0*sig*sig));
}

double dnorm(const double& x, const double& mu, const double& sig)
{
	return(log(sig)+square(x-mu)/(2.0*sig*sig));
}

dvariable dnorm(dvar_vector res, dvariable sig)
{
	int n=size_count(res);
	return(n*log(sig)+norm2(res)/(2.0*sig*sig));
}

dvar_vector dnorm2(const dvar_vector& x, const dvariable& mu,const dvariable& std)
{
	double pi=3.14159265358979323846;
	int n=size_count(x);
	return n/(sqrt(2.*pi)*std)*mfexp(-pow(x-mu,2.)/(2.*std*std));
}

dvariable full_dnorm(const dvariable& x, const dvariable& mu,const dvariable& std)
{
	double pi=3.14159265358979323846;
	//return log(1./(pow(2.*pi,0.5)*s)*exp(-square(x-m)/(2.*s*s)));
	return -(log(std)+0.5*log(2.*pi))-0.5*square(x-mu)/(std*std);
}

dvariable full_dnorm(const dvariable& res, const dvariable& s)
{
	double pi=3.14159265358979323846;
	return (-log(s)-0.5*log(2.*pi)-0.5*square(res)/(s*s));
}

dvariable full_dnorm(const dvar_vector& residual,const dvariable& std)
{
	double pi=3.14159265358979323846;
	int n=size_count(residual);
	return -n*(log(std)+0.5*log(2.*pi))-0.5*norm2(residual)/(std*std);
}

dvariable full_dnorm(const int& n, const dvar_vector& residual,const dvariable& std)
{
	double pi=3.14159265358979323846;
	return -n*(log(std)+0.5*log(2.*pi))-0.5*norm2(residual)/(std*std);
}

dvariable full_dnorm(const int& n, const dvar_vector& x, const dvar_vector mu,const dvariable& std)
{
	double pi=3.14159265358979323846;
	//int n=size_count(residual);
	return -(log(std)+0.5*log(2.*pi))-0.5*norm2(x-mu)/(std*std);
}

dvar_vector full_dnorm(const dvar_vector& x, const dvariable& mu,const dvariable& std)
{
	double pi=3.14159265358979323846;
	return -(log(std)+0.5*log(2.*pi))-0.5*pow(x-mu,2.)/(std*std);
}

dvariable full_dnorm(const dvariable& x, const double& mu, const double& std)
{
	double pi=3.14159265358979323846;
	return -(log(std)+0.5*log(2.*pi))-0.5*pow(x-mu,2.)/(std*std);
}

dvariable full_dnorm(const dvar_vector& x, const dvector& mu, const dvector& std)
{
	double pi=3.14159265358979323846;
	return sum(-(log(std)+0.5*log(2.*pi))-0.5*elem_div(pow(x-mu,2.),square(std)));
}

dvariable full_dnorm(const dvar_matrix& x, const dmatrix& mu, const dmatrix& std)
{
	double pi=3.14159265358979323846;
	return sum(-(log(std)+0.5*log(2.*pi))-0.5*elem_div(pow(x-mu,2.),square(std)));
}

dvariable full_dnorm(const dvar_matrix& residual,const dvariable& std)
{
	double pi=3.14159265358979323846;
	int n=size_count(residual);
	return n*(log(std)+0.5*log(2.*pi))+0.5*norm2(residual)/(std*std);
}

dvariable full_dlnorm(const dvariable& x, const dvariable& mu, const dvariable& std)
{
	double pi=3.14159265358979323846;
	return -0.5*log(2.*pi)-log(std)-log(x)-0.5*square(log(x)-mu)/(std*std);
}

dvariable dbeta(const dvariable& x, const dvariable& a, const dvariable& b)
{
	dvariable out;
	if(x<0. || x>1.)	out=-1e10;
	else	out=gammln(a+b)-(gammln(a)+gammln(b))+(a-1)*log(x)+(b-1)*log(1-x);
	return out;
}

dvariable dgamma(const dvar_vector& x, const dvariable& k, const dvariable& theta)
{
	int N=size_count(x);
	return((k-1.)*sum(log(x))-sum(x)/theta-N*k*log(theta)-N*gammln(k));
}

dvariable dgamma(const dvariable& x, const dvariable& k, const dvariable& theta)
{
	//return(k*log(theta)+(k-1)*log(x)-theta*x-gammln(k));
	return((k-1.)*log(x)-x/theta-k*log(theta)-gammln(k));
}

dvariable dlnorm(dvar_vector res, dvariable mu, dvariable sig)
{
	int n=size_count(res);
	return(n*log(sig)+norm2(log(res)-log(mu))/(2.0*sig*sig));
}

dvariable dlnorm(dvariable res, dvariable mu, dvariable sig)
{
	return(log(sig)+square(log(res)-log(mu))/(2.0*sig*sig));
}

dvar_vector dlnorm2(const dvar_vector& x, const dvariable& mu,const dvariable& std)
{
	double pi=3.14159265358979323846;
	int n=size_count(x);
	return elem_prod(n/(sqrt(2.*pi)*std*x),mfexp(-pow(log(x)-mu,2.)/(2.*std*std)));
}

//dvariable dlpois(const double& x,const dvariable& lambda)
//{
//	return -log(pow(lambda,x))+lambda+gammln(x+1);
//}

dvariable factorial(const double& x)
{
	double result=1;
	for(int k=x;k>1;k--)
	{
		result*=k;
	}
	return result;
}

dvector factorial(const dvector& x)
{
	int n=size_count(x);
	dvector result(1,n);
	result=1;
	for(int i=1;i<=n;i++)
	{
		for(int k=x(i);k>1;k--)
		{
			result(i)*=k;
		}
	}
	return result;
}

dvariable dpois(const double& x, const dvariable& lambda)
{
	if(lambda<0)	return(-lambda*100);
	else				return pow(lambda,x)*exp(-lambda)/factorial(x);
}

dvariable dlpois(const dvector& x, const dvar_vector& lambda)
{
	return sum(x*log(lambda)-lambda-gammln(x+1.));
}

dvariable dlpois(const double& x, const dvariable& lambda)
{
	return x*log(lambda)-lambda-gammln(x+1.);
}

dvariable dlpois(const dvar_vector& x, const dvar_vector& lambda)
{
	return(sum(x*log(lambda)-lambda-gammln(x+1.)));
}

dvariable carl_dpois(const double& x, const dvariable& lambda)
{
	return(-lambda+x*log(lambda));
}

dvariable carl_dpois(const dvar_vector& x, const dvar_vector& lambda)
{
	return(sum(-lambda+elem_prod(x,log(lambda))));
}


dvariable dbinom(const double& n, const double& k, const dvariable& p)
{
	if(p>1.)	return -1e10;
	else
		return log_comb(n,k)+k*log(p)+(n-k)*log(1-p);
}

dvariable dbinom(const dvariable& n, const double& k, const dvariable& p)
{
	if(p>1.)	return -1e10;
	else
		return log_comb(n,k)+k*log(p)+(n-k)*log(1-p);
}

dvar_vector dnbinom(const dvar_vector& x, const dvar_vector& mu, const dvariable tau)
{
	dvar_vector p(mu.indexmin(),mu.indexmax());
	dvar_vector temp(mu.indexmin(),mu.indexmax());
	temp=0.;
	for(int i=mu.indexmin();i<=mu.indexmax();i++)
	{	if(mu(i)<1.e-15)	temp(i)=1e-15;
		else						temp(i)=mu(i);}
	p=tau/(tau+temp);
	dvar_vector out=gammln(x+temp)-(gammln(temp)+gammln(x+1.))+elem_prod(mu,log(p))+elem_prod(x,log(1.-p));
	for(int i=mu.indexmin();i<=mu.indexmax();i++)
	{	if(x(i)==0.)	out(i)=0.;}
	return(out);
}

dvariable dmultinom(const dvar_vector& qi, const dvar_vector& qihat)
{
	dvar_vector pihat(qihat.indexmin(),qihat.indexmax());
	dvariable ans;
	dvariable tiny=1.e-100;
	if(sum(qi)>tiny)
	{
		pihat=qihat/sum(qihat);
		ans=sum(elem_prod(qi,log(pihat)));
	}
	else ans=0.;
	return(ans);
}

dvar_vector plogis(const dvector& x,dvariable& m,dvariable& s)
{
	return (1. / (1. + mfexp(-(x-m)/s)));

}

dvar_vector plogis(const dvector& x,const double& m,const double& s)
{
	return (1. / (1. + mfexp(-(x-m)/s)));
}

const double pi=3.141593;

dvariable get_F(const double& ct,const dvariable& M,dvariable& fpen,const dvar_vector& bt,const dvar_vector& va)
{				//this is the Baranov catch equation
	//t9 is the predicted catch
	//t16 is the derivative of the predicted catch
	double minsurv=0.01;
	dvariable ft;
	dvariable btmp=0.98*sum(elem_prod(bt,va));
	dvariable ctmp=ct;
	ft=ct/btmp;	//intial guess for ft
	if(1.-ft<minsurv){
		ft = 1.-posfun((1.-ft),minsurv,fpen);
		ctmp=ft*btmp;
	}

	for(int i = 1; i<= 7; i++)
	{
			//Barnov catch equation
			dvar_vector T3 = M + va * ft;
			dvar_vector T6 = mfexp(-T3);
			dvariable T9 = bt * elem_prod(elem_div(va * ft , T3) , (0.1e1 - T6));

			//derivative of the catch equation
			dvar_vector t2 = M + va * ft;
			dvar_vector t3 = 0.10e1 / t2;
			dvar_vector t5 = mfexp(-t2);
			dvar_vector t6 = 0.1e1 - t5;
			dvar_vector t8 = elem_prod(va , va);
			dvar_vector t9 = t8 * ft;
			dvar_vector t10 = elem_prod(t2 , t2);
			dvariable t16 =bt * (elem_prod(elem_prod(va , t3) , t6)
										- elem_prod(elem_div(t9 , t10) , t6)
										+ elem_prod(elem_prod(t9 , t3) , t5));

			//newton step
			ft -= (T9-ctmp)/t16;
	}
	return (ft);
}

dvariable negfun(const dvariable&x,const double eps,dvariable& pen)
{
	if (x<=eps) {
		return x;
	} else {
		pen+=.01*square(x-eps);
		return eps+(eps-eps/(2.-eps/x))/100.;
	}
}

dvariable posfun2(const dvariable&x,const double eps,dvariable& pen)
{
	if (x>=eps) {
		return x;
	} else {
		pen+=.01*square(x-eps);
		return eps/(1-x/eps);
	}
}

dvariable posfun3(const dvariable &x, const double eps, dvariable& pen)	//This allows eps=0
{
	if (x>=eps)	return x;
	else {
		pen+=0.01*square(x-eps);
		return x*x+eps;
	}
}
