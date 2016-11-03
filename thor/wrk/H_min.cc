#define NO 5

#include "thor_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


extern double       b2_max;
extern tps          K, g;
extern ss_vect<tps> Map, A0, A1, Map_res;

const int  max_ind = 10, n_prm_max = 20;

bool         h[max_ind][max_ind][max_ind][max_ind][max_ind], fit_chrm;
int          n_b3, b3s[mpole_max], n_bn, *bns_fam, *bns_n, n_cell;
int          check_range, adj_tune, adj_chrom, n_steps;
long int     beta_loc1, beta_loc2, beta_loc3;
double       Jx, Jy, delta, ksi1[2], nu0_x, nu0_y, eps_nu, chi2;
double       beta1[2], beta2[2], beta3[2];
double       bnL_max[mpole_max], eps_ksi, *bns;
double       scl_dnu, scl_ksi_nl, scl_dnuddelta, step, scl_dnudJ;
double       nu_x_min, nu_x_max, nu_y_min, nu_y_max;
ss_vect<tps> Id_scl;

const double max_Ax = 5e-3, max_Ay = 5e-3, max_delta = 3e-2;

const double scl_ksi1[] = { 1e3, 1e3 };
const bool   mirror_sym = true;


int ncom;
double *pcom,*xicom,(*nrfunc)(double []);
void (*nrdfun)(double [], double []);


double amotry(double **p, double y[], double psum[], int ndim,
	      double (*funk)(double []), int ihi, double fac)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry=dvector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_dvector(ptry,1,ndim);
	return ytry;
}


#define TINY 1.0e-10
#define NMAX 5000
#define GET_PSUM \
  for (j=1;j<=ndim;j++) {				\
    for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];	\
    psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void amoeba(double **p, double y[], int ndim, double ftol,
	    double (*funk)(double []), int *nfunk)
{
  int i,ihi,ilo,inhi,j,mpts=ndim+1;
  double rtol,sum,swap,ysave,ytry,*psum;

  psum=dvector(1,ndim);
  *nfunk=0;
  GET_PSUM
    for (;;) {
      ilo=1;
      ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
      for (i=1;i<=mpts;i++) {
	if (y[i] <= y[ilo]) ilo=i;
	if (y[i] > y[ihi]) {
	  inhi=ihi;
	  ihi=i;
	} else if (y[i] > y[inhi] && i != ihi) inhi=i;
      }
      rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
//      printf("%12.5e %12.5e %12.5e %12.5e\n", y[ihi], y[ilo], rtol, ftol);
      if (rtol < ftol) {
	SWAP(y[1],y[ilo])
	  for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
	    break;
      }
      if (*nfunk >= NMAX) {
	nrerror("NMAX exceeded");
	return;
      }
      *nfunk += 2;
      ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
      if (ytry <= y[ilo])
	ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
      else if (ytry >= y[inhi]) {
	ysave=y[ihi];
	ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
	if (ytry >= ysave) {
	  for (i=1;i<=mpts;i++) {
	    if (i != ilo) {
	      for (j=1;j<=ndim;j++)
		p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
	      y[i]=(*funk)(psum);
	    }
	  }
	  *nfunk += ndim;
	  GET_PSUM
	    }
      } else --(*nfunk);
    }
  free_dvector(psum,1,ndim);
}
#undef TINY
#undef GET_PSUM
#undef NMAX
#undef SWAP


double f1dim(double x)
{
  int j;
  double f,*xt;

  xt=dvector(1,ncom);
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free_dvector(xt,1,ncom);
  return f;
}


double df1dim(double x)
{
  int j;
  double df1=0.0;
  double *xt,*df;

  xt=dvector(1,ncom);
  df=dvector(1,ncom);
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  (*nrdfun)(xt,df);
  for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
  free_dvector(df,1,ncom);
  free_dvector(xt,1,ncom);
  return df1;
}


#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
	    double *fc,	double (*func)(double))
{
  double ulim,u,r,q,fu,dum;

  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
      }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	  SHFT(*fb,*fc,fu,(*func)(u))
	  }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
      }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT


#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	    double *xmin)
{
  int iter;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }
  }
  nrerror("Too many iterations in brent");
  *xmin=x;
  return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT


#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

double dbrent(double ax, double bx, double cx, double (*f)(double),
	      double (*df)(double), double tol, double *xmin)
{
  int iter,ok1,ok2;
  double a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
  double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  dw=dv=dx=(*df)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol1=tol*fabs(x)+ZEPS;
    tol2=2.0*tol1;
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      d1=2.0*(b-a);
      d2=d1;
      if (dw != dx) d1=(w-x)*dx/(dx-dw);
      if (dv != dx) d2=(v-x)*dx/(dx-dv);
      u1=x+d1;
      u2=x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde=e;
      e=d;
      if (ok1 || ok2) {
	if (ok1 && ok2)
	  d=(fabs(d1) < fabs(d2) ? d1 : d2);
	else if (ok1)
	  d=d1;
	else
	  d=d2;
	if (fabs(d) <= fabs(0.5*olde)) {
	  u=x+d;
	  if (u-a < tol2 || b-u < tol2)
	    d=SIGN(tol1,xm-x);
	} else {
	  d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
      } else {
	d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if (fabs(d) >= tol1) {
      u=x+d;
      fu=(*f)(u);
    } else {
      u=x+SIGN(tol1,d);
      fu=(*f)(u);
      if (fu > fx) {
	*xmin=x;
	return fx;
      }
    }
    du=(*df)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      MOV3(v,fv,dv, w,fw,dw)
	MOV3(w,fw,dw, x,fx,dx)
	MOV3(x,fx,dx, u,fu,du)
	} else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	MOV3(v,fv,dv, w,fw,dw)
	  MOV3(w,fw,dw, u,fu,du)
	  } else if (fu < fv || v == x || v == w) {
	MOV3(v,fv,dv, u,fu,du)
	  }
    }
  }
  nrerror("Too many iterations in routine dbrent");
  return 0.0;
}
#undef ITMAX
#undef ZEPS
#undef MOV3


#define TOL 2.0e-4

void linmin(double p[], double xi[], int n, double *fret,
	    double (*func)(double []))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=dvector(1,n);
	xicom=dvector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_dvector(xicom,1,n);
	free_dvector(pcom,1,n);
}
#undef TOL


#define TOL 2.0e-4

void dlinmin(double p[], double xi[], int n, double *fret,
	     double (*func)(double []), void (*dfunc)(double [], double []))
{
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;

  ncom=n;
  pcom=dvector(1,n);
  xicom=dvector(1,n);
  nrfunc=func;
  nrdfun=dfunc;
  for (j=1;j<=n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin);
  for (j=1;j<=n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free_dvector(xicom,1,n);
  free_dvector(pcom,1,n);
}
#undef TOL


#define TINY 1.0e-25
#define ITMAX 200

void powell(double p[], double **xi, int n, double ftol, int *iter,
	    double *fret, double (*func)(double []))
{
  int i,ibig,j;
  double del,fp,fptt,t,*pt,*ptt,*xit;

  pt=dvector(1,n);
  ptt=dvector(1,n);
  xit=dvector(1,n);
  *fret=(*func)(p);
  for (j=1;j<=n;j++) pt[j]=p[j];
  for (*iter=1;;++(*iter)) {
    fp=(*fret);
    ibig=0;
    del=0.0;
    for (i=1;i<=n;i++) {
      for (j=1;j<=n;j++) xit[j]=xi[j][i];
      fptt=(*fret);
      linmin(p,xit,n,fret,func);
      if (fptt-(*fret) > del) {
	del=fptt-(*fret);
	ibig=i;
      }
    }
    if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY) {
      free_dvector(xit,1,n);
      free_dvector(ptt,1,n);
      free_dvector(pt,1,n);
      return;
    }
    if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
    for (j=1;j<=n;j++) {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(ptt);
    if (fptt < fp) {
      t=2.0*(fp-2.0*(*fret)+fptt)*DSQR(fp-(*fret)-del)-del*DSQR(fp-fptt);
      if (t < 0.0) {
	linmin(p,xit,n,fret,func);
	for (j=1;j<=n;j++) {
	  xi[j][ibig]=xi[j][n];
	  xi[j][n]=xit[j];
	}
      }
    }
  }
}
#undef TINY
#undef ITMAX


#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	    double (*func)(double []), void (*dfunc)(double [], double []))
{
  int j,its;
  double gg,gam,fp,dgg;
  double *g,*h,*xi;

  g=dvector(1,n);
  h=dvector(1,n);
  xi=dvector(1,n);
  fp=(*func)(p);
  (*dfunc)(p,xi);
  for (j=1;j<=n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    linmin(p,xi,n,fret,func);
//    dlinmin(p,xi,n,fret,func,dfunc);
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
      FREEALL
	return;
    }
    fp= *fret;
    (*dfunc)(p,xi);
    dgg=gg=0.0;
    for (j=1;j<=n;j++) {
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
    if (gg == 0.0) {
      FREEALL
	return;
    }
    gam=dgg/gg;
    for (j=1;j<=n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }
  nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL


#define ALF 1.0e-4
#define TOLX 1.0e-7

void lnsrch(int n, double xold[], double fold, double g[], double p[],
	    double x[], double *f, double stpmax, int *check,
	    double (*func)(double []))
{
  int i;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0,rhs1,rhs2,slope,sum,temp,
    test,tmplam;

  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
  test=0.0;
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*func)(x);
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return;
    else {
      if (alam == 1.0)
	tmplam = -slope/(2.0*(*f-fold-slope));
      else {
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc < 0.0) tmplam=0.5*alam;
	  else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
	  else tmplam=-slope/(b+sqrt(disc));
	}
	if (tmplam > 0.5*alam)
	  tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    alam=FMAX(tmplam,0.1*alam);
  }
}
#undef ALF
#undef TOLX


#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_dvector(xi,1,n);free_dvector(pnew,1,n);		      \
  free_dmatrix(hessin,1,n,1,n);free_dvector(hdg,1,n);free_dvector(g,1,n);     \
  free_dvector(dg,1,n);

void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	    double(*func)(double []), void (*dfunc)(double [], double []))
{
  int check,i,its,j;
  double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
  double *dg,*g,*hdg,**hessin,*pnew,*xi;

  dg=dvector(1,n);
  g=dvector(1,n);
  hdg=dvector(1,n);
  hessin=dmatrix(1,n,1,n);
  pnew=dvector(1,n);
  xi=dvector(1,n);
  fp=(*func)(p);
  (*dfunc)(p,g);
  for (i=1;i<=n;i++) {
    for (j=1;j<=n;j++) hessin[i][j]=0.0;
    hessin[i][i]=1.0;
    xi[i] = -g[i];
    sum += p[i]*p[i];
  }
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
    fp = *fret;
    for (i=1;i<=n;i++) {
      xi[i]=pnew[i]-p[i];
      p[i]=pnew[i];
    }
    test=0.0;
    for (i=1;i<=n;i++) {
      temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) {
      FREEALL
      cout << scientific << setprecision(3)
	   << "dfpmin (TOLX): " << test << " (" << TOLX << ")" << endl;
      return;
    }
    for (i=1;i<=n;i++) dg[i]=g[i];
    (*dfunc)(p,g);
    test=0.0;
    den=FMAX(*fret,1.0);
    for (i=1;i<=n;i++) {
      temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
      if (temp > test) test=temp;
    }
    if (test < gtol) {
      FREEALL
      cout << scientific << setprecision(3)
	   << "dfpmin (gtol): " << test << " (" << gtol << ")" << endl;
      return;
    }
    for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
    for (i=1;i<=n;i++) {
      hdg[i]=0.0;
      for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
    }
    fac=fae=sumdg=sumxi=0.0;
    for (i=1;i<=n;i++) {
      fac += dg[i]*xi[i];
      fae += dg[i]*hdg[i];
      sumdg += SQR(dg[i]);
      sumxi += SQR(xi[i]);
    }
    if (fac > sqrt(EPS*sumdg*sumxi)) {
      fac=1.0/fac;
      fad=1.0/fae;
      for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
      for (i=1;i<=n;i++) {
	for (j=i;j<=n;j++) {
	  hessin[i][j] += fac*xi[i]*xi[j]
	    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
	  hessin[j][i]=hessin[i][j];
	}
      }
    }
    for (i=1;i<=n;i++) {
      xi[i]=0.0;
      for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
    }
  }
  nrerror("too many iterations in dfpmin");
  FREEALL
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL


void prt_b3s(double bns[])
{
  int  i, j;
  FILE *outf;

  const string file_name = "b3s.out";

  outf = file_write(file_name.c_str());
  
  for (i = 1; i <= ncom; i++)
    for (j = 1; j <= get_n_Kids(abs(bns_fam[i])); j++) {
      if (bns_fam[i] > 0) 
	fprintf(outf, "%9s(%2d) = %10.5f %1d\n",
		get_Name(abs(bns_fam[i])), j, get_bnL(bns_fam[i], 1, bns_n[i]),
		bns_n[i]);
   }

  fclose(outf);
}

void select_h(void)
{
  int  i, j, k, l, m;

  for (i = 0; i <= no_tps-1; i++)
    for (j = 0; j <= no_tps-1; j++)
      for (k = 0; k <= no_tps-1; k++)
	for (l = 0; l <= no_tps-1; l++)
	  for (m = 0; m <= no_tps-1; m++)
	    // enable all non-chromatic driving terms
	    h[i][j][k][l][m] = ((i+j+k+l <= 4) && (m == 0))? true : false;

  // linear chromaticity
  h[1][1][0][0][1] = true; h[0][0][1][1][1] = true;

  // 2nd order amplitude dependent tune shift
  h[2][2][0][0][0] = true; h[1][1][1][1][0] = true; h[0][0][2][2][0] = true;

  // 4th order amplitude dependent tune shift
  // balance nonlinear terms
  h[3][3][0][0][0] = true; h[2][2][1][1][0] = true;
  h[0][0][3][3][0] = true; h[1][1][2][2][0] = true;

  // 6th order amplitude dependent tune shift
  h[4][4][0][0][0] = true; h[0][0][4][4][0] = true;
  h[3][3][1][1][0] = true; h[1][1][3][3][0] = true;
  h[2][2][2][2][0] = true;

  if (true) {
    // nonlinear chromaticity
    for (m = 2; m <= no_tps-3; m++) {
      h[1][1][0][0][m] = true; h[0][0][1][1][m] = true;
    }

    // balance nonlinear terms
    //    h[1][1][0][0][3] = false; h[0][0][1][1][3] = false;
    //    h[1][1][0][0][5] = false; h[0][0][1][1][5] = false;
  }

  if (true) {
    // amplitude dependent chromaticity

    for (m = 1; m <= no_tps-5; m++) {
      h[2][2][0][0][m] = true; h[0][0][2][2][m] = true;
      h[1][1][1][1][m] = true;
    }

    for (m = 1; m <= no_tps-7; m++) {
      h[3][3][0][0][m] = true; h[0][0][3][3][m] = true;
      h[2][2][1][1][m] = true; h[1][1][2][2][m] = true;
    }
  }

  // exclude tune
  h[1][1][0][0][0] = false; h[0][0][1][1][0] = false;

  // exclude delta dependence of pathlength
  for (m = 2; m <= no_tps-1; m++)
    h[0][0][0][0][m] = false;

  if (true) {
    // higher order dispersion
    for (m = 2; m <= no_tps-2; m++) {
      h[1][0][0][0][m] = true; h[0][1][0][0][m] = true;
    }

    // delta dependance of the beta functions
    for (m = 1; m <= no_tps-3; m++) {
      h[2][0][0][0][m] = true; h[0][2][0][0][m] = true;
      h[0][0][2][0][m] = true; h[0][0][0][2][m] = true;
    }
  }
}


double H_fun(double bns[])
{
  const int  m_max = 500;  // max no of constraints

  char    hs[m_max][max_str];
  int     i, j, k, l, m, m1;
  double  h2, b[m_max];
  tps     r, K_re, K_im, g_re, g_im;

  const int  n_prt = 10;

  prt_b3s(bns);

  for (i = 1; i <= ncom; i++)
    set_bn(bns_fam[i], bns_n[i], bns[i]);

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  CtoR(K*Id_scl, K_re, K_im); CtoR(g*Id_scl, g_re, g_im);

  // mirror symmetric cell => g_re = 0
  m1 = 0; h2 = 0.0;
  for (i = 0; i <= no_tps-1; i++)
    for (j = 0; j <= no_tps-1; j++)
      for (k = 0; k <= no_tps-1; k++)
	for (l = 0; l <= no_tps-1; l++)
	  for (m = 0; m <= no_tps-1; m++) {
	    if ((0 < i+j+k+l+m) && (i+j+k+l+m <= no_tps-1) &&
		h[i][j][k][l][m] &&
		((fabs(h_ijklm(g_im, i, j, k, l, m)) > 0.0) ||
		 (fabs(h_ijklm(K_re, i, j, k, l, m)) > 0.0))) {
	      sprintf(hs[++m1-1], "h_%d%d%d%d%d", i, j, k, l, m);

	      if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m)) {
		// horizontal linear chromaticity
		b[m1-1] = scl_ksi1[X_]*(ksi1[X_]*M_PI*2.0*Jx*delta
					+h_ijklm(K_re, i, j, k, l, m));
		h2 += sqr(b[m1-1]);
		  
	      } else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m)) {
		// vertical linear chromaticity
		b[m1-1] = scl_ksi1[Y_]*(ksi1[Y_]*M_PI*2.0*Jy*delta
					+h_ijklm(K_re, i, j, k, l, m));
		h2 += sqr(b[m1-1]);
	      } else if (is_h_ijklm(2, 2, 0, 0, 0, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 2, 2, 0, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 1, 1, 0, i, j, k, l, m) ||
			 is_h_ijklm(3, 3, 0, 0, 0, i, j, k, l, m) ||
			 is_h_ijklm(2, 2, 1, 1, 0, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 2, 2, 0, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 3, 3, 0, i, j, k, l, m) ||
			 is_h_ijklm(4, 4, 0, 0, 0, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 4, 4, 0, i, j, k, l, m) ||
			 is_h_ijklm(3, 3, 1, 1, 0, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 3, 3, 0, i, j, k, l, m) ||
			 is_h_ijklm(2, 2, 2, 2, 0, i, j, k, l, m)) {
		// amplitude dependent tune shift
		b[m1-1] = scl_dnu*h_ijklm(K_re, i, j, k, l, m);
		h2 += sqr(b[m1-1]);
	      } else if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 0, 0, 4, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 1, 1, 4, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 0, 0, 5, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 1, 1, 5, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 0, 0, 6, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 1, 1, 6, i, j, k, l, m)) {
		// nonlinear chromaticity
		b[m1-1] = scl_ksi_nl*h_ijklm(K_re, i, j, k, l, m);
		h2 += sqr(b[m1-1]);
	      } else if (is_h_ijklm(2, 2, 0, 0, 1, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 2, 2, 1, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 1, 1, 1, i, j, k, l, m) ||
			 is_h_ijklm(2, 2, 0, 0, 2, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 2, 2, 2, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 1, 1, 2, i, j, k, l, m) ||
			 is_h_ijklm(2, 2, 0, 0, 3, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 2, 2, 3, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 1, 1, 3, i, j, k, l, m) ||
			 is_h_ijklm(2, 2, 0, 0, 4, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 2, 2, 4, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 1, 1, 4, i, j, k, l, m) ||
			 is_h_ijklm(3, 3, 0, 0, 1, i, j, k, l, m) ||
			 is_h_ijklm(2, 2, 1, 1, 1, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 2, 2, 1, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 3, 3, 1, i, j, k, l, m) ||
			 is_h_ijklm(3, 3, 0, 0, 2, i, j, k, l, m) ||
			 is_h_ijklm(2, 2, 1, 1, 2, i, j, k, l, m) ||
			 is_h_ijklm(1, 1, 2, 2, 2, i, j, k, l, m) ||
			 is_h_ijklm(0, 0, 3, 3, 2, i, j, k, l, m)) {
		// amplitude dependent chromaticity
		b[m1-1] = scl_ksi_nl*h_ijklm(K_re, i, j, k, l, m);
		h2 += sqr(b[m1-1]);
	      } else {
		sprintf(hs[m1-1], "i_%d%d%d%d%d", i, j, k, l, m);
		b[m1-1] = h_ijklm(g_im, i, j, k, l, m);
		h2 += sqr(b[m1-1]);
		if (!mirror_sym) {
		  sprintf(hs[++m1-1], "r_%d%d%d%d%d", i, j, k, l, m);
		  b[m1-1] = h_ijklm(g_re, i, j, k, l, m);
		  h2 += sqr(b[m1-1]);
		}
	      }
	    }
	  }

  if (false || (h2 < chi2)) {
    i = 0;
    do {
      cout << endl;
      for (j = 1; j <= n_prt; j++) {
	i++;
	if (i <= m1) cout << setw(10) << hs[i-1];
      }
      cout << endl;

      i -= n_prt;
      for (j = 1; j <= n_prt; j++) {
	i++;
	if (i <= m1)
	  cout << scientific << setprecision(2) << setw(10) << b[i-1];
      }
      cout << endl;
    } while (i < m1);

    cout << endl;
    cout << scientific << setprecision(3)
	 << "chi2:" << setw(10) << h2 << endl;

    for (i = 1; i <= ncom; i++)
      cout << scientific << setprecision(3)
	   << setw(11) << get_bnL(bns_fam[i], 1, bns_n[i]);
    cout << endl;

    chi2 = h2;
  }

  return h2;
}


double get_dh(const tps &h,
	      const int i, const int j, const int k, const int l, const int m)
{

  return 2.0*h_ijklm(h, i, j, k, l, m)*h_ijklm_p(h, i, j, k, l, m, 7);
}


void H_dfun(double bns[], double dh2[])
{
  int      i1, i, j, k, l, m;
  tps      r, K_re, K_im, g_re, g_im;
  ofstream outf;

  for (i = 1; i <= ncom; i++)
    set_bn(bns_fam[i], bns_n[i], bns[i]);

  // mirror symmetric cell => g_re = 0
  for (i1 = 1; i1 <= ncom; i1++) {
    set_bn_par(bns_fam[i1], bns_n[i1], 7);

    danot_(no_tps-1);
    get_Map();
    danot_(no_tps);
    K = MapNorm(Map, g, A1, A0, Map_res, 1);
    CtoR(K*Id_scl, K_re, K_im); CtoR(g*Id_scl, g_re, g_im);

    dh2[i1] = 0.0;
    for (i = 0; i <= no_tps-1; i++)
      for (j = 0; j <= no_tps-1; j++)
	for (k = 0; k <= no_tps-1; k++)
	  for (l = 0; l <= no_tps-1; l++)
	    for (m = 0; m <= no_tps-1; m++) {
	      if ((0 < i+j+k+l+m) && (i+j+k+l+m <= no_tps-1) &&
		  h[i][j][k][l][m] &&
		  ((fabs(h_ijklm(g_im, i, j, k, l, m)) > 0.0) ||
		   (fabs(h_ijklm(K_re, i, j, k, l, m)) > 0.0))) {
		if (is_h_ijklm(1, 1, 0, 0, 1, i, j, k, l, m)) {
		  // horizontal linear chromaticity
		  dh2[i1] +=
		    sqr(scl_ksi1[X_])
		    *2.0*(M_PI*2.0*Jx*delta+h_ijklm(K_re, i, j, k, l, m))
		    *h_ijklm_p(K_re, i, j, k, l, m, 7);
		} else if (is_h_ijklm(0, 0, 1, 1, 1, i, j, k, l, m)) {
		  // vertical linear chromaticity
		  dh2[i1] += sqr(scl_ksi1[Y_])*get_dh(K_re, i, j, k, l, m);
		} else if (is_h_ijklm(2, 2, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 0, i, j, k, l, m) ||
			   is_h_ijklm(4, 4, 0, 0, 0, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 4, 4, 0, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 1, 1, 0, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 3, 3, 0, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 2, 2, 0, i, j, k, l, m)) {
		  // amplitude dependent tune shift
		  dh2[i1] +=
		    sqr(scl_ksi1[Y_])
		    *2.0*(M_PI*2.0*Jy*delta+h_ijklm(K_re, i, j, k, l, m))
		    *h_ijklm_p(K_re, i, j, k, l, m, 7);
		} else if (is_h_ijklm(1, 1, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 3, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 3, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 4, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 4, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 5, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 5, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 0, 0, 6, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 1, 1, 6, i, j, k, l, m)) {
		  // nonlinear chromaticity
		  dh2[i1] += sqr(scl_ksi_nl)*get_dh(K_re, i, j, k, l, m);
		} else if (is_h_ijklm(2, 2, 0, 0, 1, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 1, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 1, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 3, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 3, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 3, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 0, 0, 4, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 2, 2, 4, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 1, 1, 4, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 1, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 1, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 1, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 1, i, j, k, l, m) ||
			   is_h_ijklm(3, 3, 0, 0, 2, i, j, k, l, m) ||
			   is_h_ijklm(2, 2, 1, 1, 2, i, j, k, l, m) ||
			   is_h_ijklm(1, 1, 2, 2, 2, i, j, k, l, m) ||
			   is_h_ijklm(0, 0, 3, 3, 2, i, j, k, l, m)) {
		  // amplitude dependent chromaticity
		  dh2[i1] += sqr(scl_ksi_nl)*get_dh(K_re, i, j, k, l, m);
		} else {
		  dh2[i1] += get_dh(g_im, i, j, k, l, m);
		  if (!mirror_sym)
		    dh2[i1] += get_dh(g_re, i, j, k, l, m);
		}
	      }
	    }

    clr_bn_par(bns_fam[i1], bns_n[i1]);
  }

  if (false) {
    cout << endl;
    cout << "dh2" << endl;
    for (i = 1; i <= ncom; i++)
      cout << scientific << setprecision(3) << setw(11) << dh2[i];
    cout << endl;
  }
}


void H_min(const int k)
{
  int    i, j, iter;
  double **xi, *y, **p, fret;

  const double d_bn = 1e-1;

  select_h();

  ncom = n_bn; nrfunc = H_fun; nrdfun = H_dfun;

  chi2 = 1e10;

  switch (k) {
  case 1:\
    printf("\nDownhill simplex\n");

    y = dvector(1, n_bn); p = dmatrix(1, n_bn+1, 1, n_bn);

    for (i = 1; i <= n_bn+1; i++) {
      for (j = 1; j <= n_bn; j++) {
	p[i][j] = bns[j];
	if (i == j) p[i][j] += d_bn;
      }
      y[i] = H_fun(p[i]);
    }

    amoeba(p, y, n_bn, 1e-10, H_fun, &iter);

    free_dvector(y, 1, n_bn); free_dmatrix(p, 1, n_bn+1, 1, n_bn);
    break;
  case 2:
    printf("\nPowell's method\n");

    xi = dmatrix(1, n_bn, 1, n_bn);

    for (i = 1; i <= n_bn; i++)
      for (j = 1; j <= n_bn; j++)
	xi[i][j] = (i == j)? 1.0 : 0.0; 

    powell(bns, xi, n_bn, 1e-10, &iter, &fret, H_fun);

    free_dmatrix(xi, 1, n_bn, 1, n_bn);
    break;
  case 3:
    printf("\nConjugated gradient\n");

    frprmn(bns, n_bn, 1e-15, &iter, &fret, H_fun, H_dfun);
    break;
  case 4:
    printf("\nMetric method\n");

    dfpmin(bns, n_bn, 1e-30, &iter, &fret, H_fun, H_dfun);
    break;
  default:
    printf("\n*** H_min: undefined method %d (1-4)\n", k);
    exit(1);
  }
}


void chk_lat(double nu[], double ksi[])
{
  double        alpha1[2];
  ss_vect<tps>  nus;

  // get_Map();
  get_COD(10, 1e-10, 0.0, true);
  K = MapNorm(Map, g, A1, A0, Map_res, 1);
  nus = dHdJ(K); get_nu_ksi(nus, nu, ksi); get_ab(alpha1, beta1, 0);

  printf("\nalpha_x  = %10.3f, alpha_y = %10.3f"
	 ", beta_x =  %10.3f, beta_x =  %10.3f\n",
	 alpha1[X_], alpha1[Y_], beta1[X_] , beta1[Y_]);
  prt_nu(nus);
}


void get_locs()
{
  double  alpha1[2], alpha2[2], alpha3[2];

  beta_loc1 = get_loc(get_Fnum("mp"), 1); get_ab(alpha1, beta1, beta_loc1);
  beta_loc2 = get_loc(get_Fnum("ss"), 1); get_ab(alpha2, beta2, beta_loc2);
//  beta2[X_] = 1.0; beta2[Y_] = 1.0;
  beta_loc3 = get_loc(get_Fnum("ls"), 1); get_ab(alpha3, beta3, beta_loc3);
//  beta3[X_] = 15.0; beta3[Y_] = 3.0;

  printf("\nalpha1_x = %6.3f, alpha1_y = %6.3f"
	 ", beta1_x = %6.3f, beta1_y  = %6.3f",
	 alpha1[X_], alpha1[Y_], beta1[X_], beta1[Y_]);
  printf("\nalpha2_x = %6.3f, alpha2_y = %6.3f"
	 ", beta2_x = %6.3f, beta2_y = %6.3f",
	 alpha2[X_], alpha2[Y_], beta2[X_], beta2[Y_]);
  printf("\nalpha3_x = %6.3f, alpha3_y = %6.3f"
	 ", beta3_x = %6.3f, beta3_y = %6.3f",
	 alpha3[X_], alpha3[Y_], beta3[X_], beta3[Y_]);
}


void no_mpoles(void)
{
  int j, k;

  cout << endl;
  cout << "zeroing multipoles" << endl;
  cout << endl;
  for (j = 0; j < n_elem; j++)
    if (elem[j].kind == Mpole)
      for (k = Sext; k < mpole_max; k++) {
//	cout << "zeroing " << elem[j].Name << endl;
	set_bn(elem[j].Fnum, elem[j].Knum, k, 0.0);
      }
}


void get_prm(const char *file_name)
{
  char      line[max_str];      
  ifstream  prm_in;

  file_rd(prm_in, file_name);

  do
    prm_in.getline(line, max_str);
  while (strstr(line, "#") != NULL);

  sscanf(line, "%*s %d %lf %lf %lf", &adj_tune, &nu0_x, &nu0_y, &eps_nu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d %lf %lf %lf",
	 &adj_chrom, &ksi1[X_], &ksi1[Y_], &eps_ksi);
  fit_chrm = eps_ksi < 0.0; eps_ksi = fabs(eps_ksi);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d", &check_range);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_x_min);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_x_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_y_min);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &nu_y_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d", &n_steps);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %d", &n_cell);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &ds_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &b2_max);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Sext]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Oct]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &bnL_max[Dec]);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnu);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_ksi_nl);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnuddelta);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &scl_dnudJ);

  prm_in.getline(line, max_str);
  sscanf(line, "%*s %lf", &step);

  cout << endl;
  cout << fixed << setprecision(6)
       << "fit_tune      = " << adj_tune
       << ", nu0_x  = " << nu0_x << ", nu0_y  = " << nu0_y
       << scientific << setprecision(1) << ", eps_nu = " << eps_nu << endl;
  cout << fixed << setprecision(6)
       << "fit_chrom     = " << adj_chrom
       << ", ksi0_x = " << ksi1[X_] << ", ksi0_y = " << ksi1[Y_]
       << scientific << setprecision(1) << ", eps_ksi = " << eps_ksi << endl;
  cout << "check_range   = " << check_range << endl;
  cout << endl;
  cout << fixed << setprecision(6)
       << "n_steps   = " << n_steps
       << ", nu_x_min = " << nu_x_min << ", nu_x_max = " << nu_x_max
       << ", nu_y_min = " << nu_y_min << ", nu_y_max = " << nu_y_max << endl;
  cout << endl;
  cout << "n_cell        = " << n_cell << endl;
  cout << fixed << setprecision(2)
       << "ds_max        = " << ds_max << endl;
  cout << fixed << setprecision(1)
       << "b2_max        = " << b2_max << endl;
  cout << fixed << setprecision(1)
       << "b3L_max       = " << bnL_max[Sext] << endl;
  cout << fixed << setprecision(1)
       << "b4L_max       = " << bnL_max[Oct] << endl;
  cout << fixed << setprecision(1)
       << "b5L_max       = " << bnL_max[Dec] << endl;
  cout << fixed << setprecision(1)
       << "scl_dnu       = " << scl_dnu << endl;
  cout << fixed << setprecision(1)
       << "scl_ksi_nl    = " << scl_ksi_nl << endl;
  cout << scientific << setprecision(1)
       << "scl_dnuddelta = " << scl_dnuddelta << endl;
  cout << scientific << setprecision(1)
       << "scl_dnuddJ    = " << scl_dnudJ << endl;
  cout << fixed << setprecision(2)
       << "step          = " << step << endl;
}


void get_b2s(int &n_b2, int b2_Fams[])
{

  n_b2 = 0;

  b2_Fams[n_b2++] = get_Fnum("ql1");
  b2_Fams[n_b2++] = get_Fnum("ql2");
  b2_Fams[n_b2++] = get_Fnum("ql3");

  b2_Fams[n_b2++] = get_Fnum("qh1");
  b2_Fams[n_b2++] = get_Fnum("qh2");
  b2_Fams[n_b2++] = get_Fnum("qh3");
}


void get_bns()
{
  int  k;

  bns_fam = ivector(1, n_prm_max); bns_n = ivector(1, n_prm_max);
  bns = dvector(1, n_prm_max);

  n_bn = 0;

  bns_fam[++n_bn] = get_Fnum("sf");  bns_n[n_bn] = Sext;
  bns_fam[++n_bn] = get_Fnum("sd");  bns_n[n_bn] = Sext;

  bns_fam[++n_bn] = get_Fnum("o1");  bns_n[n_bn] = Oct;
  bns_fam[++n_bn] = get_Fnum("o2");  bns_n[n_bn] = Oct;
  bns_fam[++n_bn] = get_Fnum("o3");  bns_n[n_bn] = Oct;
  bns_fam[++n_bn] = get_Fnum("o4");  bns_n[n_bn] = Oct;

  if (n_bn > n_prm_max) {
    cout << "get_bns: n_prm_max exceeded " << n_bn << "(" << n_prm_max
	 << ")" << endl;
    exit(0);
  }

  for (k = 1; k <= n_bn; k++)
    bns[k] = get_bn(bns_fam[k], 1, bns_n[k]);

  cout << endl;
  cout << "get_bns: no of multipole families " << n_bn << endl;
}


int main(int argc, char *argv[])
{
  int              b2_Fams[n_prm_max], n_b2 = 0;
  double           alpha[2], beta[2], nu[2], ksi[2];
  tps              Hnl, H2, gn, h, h_re, h_im, H_num, dH, H, H_re, H_im;
  tps              K_re, K_im;
  tps              g_re, g_im;
  ss_vect<double>  x;
  ss_vect<tps>     A_inv, ps, ps_re, ps_im, nus, R, R_inv, Map2;
  ss_vect<tps>     Map_Fl, Mn, Mk, J, eta, M1, M2, M3, M4, Mp;
  ofstream         outf, K_out, nus_out, A1_out, J_out;
  ifstream         inf;


  danot_(no_tps-1);

  rad_on    = false; H_exact        = false; totpath_on   = false;
  cavity_on = false; quad_fringe_on = false; emittance_on = false;
  IBS_on    = false;


  rd_mfile(argv[1], elem); rd_mfile(argv[1], elem_tps);

  // initialize the symplectic integrator after energy has been defined
  ini_si();

  // disable from TPSALib- and LieLib log messages
  idprset(-1);

  danot_(3);

  if (true) chk_lat(nu, ksi);

  get_ab(alpha, beta, 0);
  cout << endl;
  cout << fixed << setprecision(3)
       << "alpha_x  = " << setw(6) << alpha[X_]
       << ", alpha_y = " << setw(6) << alpha[Y_]
       << ", beta_x = " << setw(6) << beta[X_]
       << ", beta_y  = " << setw(6) << beta[Y_] << endl;

  if (true) prt_alphac();

  Jx = sqr(max_Ax)/(2.0*beta1[X_]); Jy = sqr(max_Ay)/(2.0*beta1[Y_]);
  delta = max_delta;

  Id_scl.identity();
  Id_scl[x_] *= sqrt(2.0*Jx); Id_scl[px_] *= sqrt(2.0*Jx);
  Id_scl[y_] *= sqrt(2.0*Jy); Id_scl[py_] *= sqrt(2.0*Jy);
  Id_scl[delta_] *= delta;

  get_prm("tune_scan.prm");

  if (adj_tune) {
    get_locs();

    get_b2s(n_b2, b2_Fams);

    fit_tune(nu0_x, nu0_y,
	     beta1[X_], beta1[Y_], beta_loc1,
	     beta2[X_], beta2[Y_], beta_loc2,
	     beta3[X_], beta3[Y_], beta_loc3,
	     n_b2, b2_Fams, eps_nu, true);
  }

  if (adj_chrom) {
    n_b3 = 0;
    b3s[n_b3++] = get_Fnum("sf"); b3s[n_b3++] = get_Fnum("sd");

    if (fit_chrm) {
      danot_(3);
      no_mpoles();
      cavity_on = false;
      fit_chrom(ksi1[X_], ksi1[Y_], n_b3, b3s, true);
    }

    get_bns(); H_min(atoi(argv[2]));
  }

  danot_(no_tps-1);
  get_Map();
  danot_(no_tps);
  K = MapNorm(Map, g, A1, A0, Map_res, no_tps);
  CtoR(K*Id_scl, K_re, K_im); CtoR(g*Id_scl, g_re, g_im);
  CtoR(get_h()*Id_scl, h_re, h_im);
  CtoR(get_H()*Id_scl, H_re, H_im);
  nus = dHdJ(K);

  if (false) Id_scl.identity();

  file_wr(outf, "map.dat"); outf << Map; outf.close();
  file_wr(outf, "K.dat"); outf << K_re*Id_scl; outf.close();
  file_wr(outf, "g.dat"); outf << g_im*Id_scl; outf.close();
  file_wr(outf, "nus.dat"); outf << nus[3] << nus[4]; outf.close();
  file_wr(outf, "h.dat"); outf << h_re*Id_scl; outf.close();
  file_wr(outf, "H.dat"); outf << H_re*Id_scl; outf.close();
}
