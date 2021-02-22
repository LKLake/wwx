#include "matvec.h"
void gmres(SparseMtx<double> a,Vcr<double>& x, const Vcr<double>& b, double& eps, int& iter, int pn, int m)
{
	 int ret =  a.GMRES(x, b,eps,iter,pn,m);
	// return ret;

}