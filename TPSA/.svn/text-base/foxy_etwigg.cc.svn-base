#include <math.h>
#include "tracy.h"


/* Grammar:

     //DA { "B"egin, "D"efinition, "E"nd }
       { "D"efinition, "V"ariable, "F"unction, "SUB"routine }
       { "RI", "RE"al, "IN"teger, "DA" }
       { "INTernal, "EXT"ernal, "COM", "FUN"ction }
       { variable name } { order } { no of variables } { array size }
*/

void get_Axy(double z, double BoBrho, double kx, double ky, double kz,
	     DAvect *x, DAvect *AxoBrho, DAvect *AyoBrho, long radia)
{

//DA  B D;
//DA  D V DA EXT x 6 6 6 ;
//DA  D V DA EXT AxoBrho 6 6 4 ; D V DA EXT AyoBrho 6 6 4 ;
//DA  D V DA INT cx 6 6 ;        D V DA INT sx 6 6 ;
//DA  D V RE INT kx ;            D V RE INT ky ; D V RE INT kz ;
//DA  D V DA INT chy 6 6 ;       D V DA INT shy 6 6 ;
//DA  D V RE INT sz ;            D V RE INT cz ;
//DA  D V IN INT i ;
//DA  D V RE INT BoBrho ;
//DA  E D ;
//DA  0 ;

  int    i;
  double cz;

//--------------------------------------------------------------------------
//DA  cx = cos(kx*x[1]); sx = sin(kx*x[1]);                             //DA
  for (i = 0; i <= 3; ++i) {
//DA  AxoBrho[i] = 0.0; AyoBrho[i] = 0.0;                               //DA
  }
  for (i = 1; i <= 1; ++i) {
    /* the sum over harmonics assumes kx zero */
//DA  chy = cosh(i*ky*x[3]); shy = sinh(i*ky*x[3]);                     //DA
    sz = sin(i*kz*z);
//DA  AxoBrho[0] = AxoBrho[0] + BoBrho/kz*cx*chy*sz;                    //DA
//DA  AyoBrho[0] = AyoBrho[0] + BoBrho*kx/(ky*kz)*sx*shy*sz;            //DA
    /* derivatives with respect to x */                                 //DA
//DA  AxoBrho[1] = AxoBrho[1] - BoBrho*kx/kz*sx*chy*sz;                 //DA
//DA  AyoBrho[1] = AyoBrho[1] + BoBrho*(kx*kx)/(ky*kz)*cx*shy*sz;       //DA
    /* derivatives with respect to y */
//DA  AxoBrho[2] = AxoBrho[2] + BoBrho*ky/kz*cx*shy*sz;                 //DA
//DA  AyoBrho[2] = AyoBrho[2] + BoBrho*kx/kz*sx*chy*sz;                 //DA

    if (globval.radiation) {
      cz = cos(kz*z);
      /* derivatives with respect to z */
//DA    AxoBrho[3] = AxoBrho[3] + BoBrho*cx*chy*cz;                     //DA
//DA    AyoBrho[3] = AyoBrho[3] + BoBrho*kx/ky*sx*shy*cz;               //DA
    }
  }

//DA  dadal
}


void etwigg(long nstep, double len, double lambda, double BoBrho,
	    double kx, DAvect *map, double crad, long pthlen, long radia)
{
  int    i;
  double hodp, b[3], h;
  double z, b2, d1, d2, x2;
  double a11, a12, a21, a22, c11, c12, c21, c22, dp, xf[6];
  double ky, kz, xp, yp;
  double det, AxoBrho[4], AyoBrho[4];

  const double zero = 0.0;

//DA  B D;
//DA  D V DA EXT map 6 6 6 ;
//DA  D V RE INT x_ ;            D V RE INT px_ ;
//DA  D V RE INT y_ ;            D V RE INT py_ ;
//DA  D V RE INT delta_ ;        D V RE INT ct_ ;
//DA  D V DA INT AxoBrho 6 6 4 ; D V DA INT AyoBrho 6 6 4 ;
//DA  D V DA INT dp 6 6 ;        D V DA INT hodp 6 6 ; D V RE INT h ;
//DA  D V DA INT a11 6 6 ;       D V DA INT a12 6 6 ;
//DA  D V DA INT a21 6 6 ;       D V DA INT a22 6 6 ;
//DA  D V DA INT det 6 6 ;
//DA  D V DA INT d1 6 6 ;        D V DA INT d2 6 6 ;
//DA  D V DA INT c11 6 6 ;       D V DA INT c12 6 6 ;
//DA  D V DA INT c21 6 6 ;       D V DA INT c22 6 6 ;
//DA  D V DA INT x2 6 6 ;
//DA  D V DA INT xp 6 6 ;        D V DA INT yp 6 6 ;
//DA  D V DA INT xf 6 6 6 ;
//DA  D V DA INT B 6 6 3 ;       D V DA INT B2 6 6 ;
//DA  D F RI get_Axy 9 ;
//DA  D F DA B2perp 5 ;
//DA  D F RI pow 2 ;
//DA  D V RE INT crad ;
//DA  E D ;
//DA  0 ;

  /* first order symplectic integrator for wiggler using expanded Hamiltonian*/

  kz = M_PI*2.0/lambda; ky = sqrt(pow(kz, 2) + pow(kx, 2));
  h = len/nstep; z = 0.0;
  for (i = 1; i <= nstep; ++i) {
    get_Axy(z, BoBrho, kx, ky, kz, &map[1], AxoBrho, AyoBrho, radia);

//--------------------------------------------------------------------------
//DA  dp = map[5] + 1.0; hodp = h/dp;                                   //DA
//DA  a11 = hodp*AxoBrho[1]; a12 = hodp*AyoBrho[1];                     //DA
//DA  a21 = hodp*AxoBrho[2]; a22 = hodp*AyoBrho[2];                     //DA
//DA  det = 1.0 - a11 - a22 + a11*a22 - a12*a21;                        //DA
//DA  d1 = hodp*AxoBrho[0]*AxoBrho[1]; d2 = hodp*AxoBrho[0]*AxoBrho[2]; //DA
//DA  c11 = (1.0-a22)/det; c12 = a12/det;                               //DA
//DA  c21 = a21/det; c22 = (1.0-a11)/det;                               //DA
//DA  x2 = c11*(map[2]-d1) + c12*(map[4]-d2);                           //DA

//DA  map[4] = c21*(map[2]-d1) + c22*(map[4]-d2); map[2] = x2;          //DA
//DA  map[1] = map[1] + hodp*(map[2]-AxoBrho[0]);                       //DA
//DA  map[3] = map[3] + hodp*map[4];                                    //DA
//DA  map[6] = map[6] + h*(((map[2]-AxoBrho[0])/dp)**2
//DA           +((map[4]-AyoBrho[0])/dp)**2)/2.0;                       //DA

    if (pthlen == 1)
//DA    map[6] = map[6] + h;                                            //DA

    if (globval.radiation) {
//DA    xp = x[px_]/dp; yp = x[py_]/dp;                                 //DA
//DA    B[0] = -AyoBrho[3]; B[1] = AxoBrho[3];                          //DA
//DA    B[2] = AyoBrho[1] - AxoBrho[2];                                 //DA
//DA    B2 = B2perp(0.0, B[0], x[x_], xp, yp);                          //DA
//DA    xf[py_] = -(crad)*pow(dp, 2)*B2                                 //DA
//DA              *((pow(xp, 2)+pow(yp, 2))/2.0+1.0);                   //DA
      /* good in large machine: conservation of the dx/dz and dy/dz */
//DA    xf[x_] = xp*xf[py_]; xf[y_] = yp*xf[py_];                       //DA
//DA    x[px_] = x[px_] + h*xf[x_]; x[py_] = x[py_] + h*xf[y_];         //DA
//DA    x[delta_] = x[delta_] + h*xf[py_];                              //DA
    }

    z += h;
  }

//DA  dadal
}
