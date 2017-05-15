#define NO 7

#include "thor_lib.h"

int  no_tps   = NO,
     ndpt_tps = 5;


int main()
{

  cout << endl;
  cout << "allocate a, b, and c" << endl;
  tps  a, b, c;
//  tps  a = tps(0.0, 1), b = tps(0.0, 2), c;

  cout << endl;
  cout << "initialize a and b" << endl;
  a = tps(1.0, 1); b = tps(2.0, 2);
  cout << a << b;

  cout << endl;
  cout << "a+b" << endl;
  c = a + b;
  cout << c;

  // test copy constructor
  cout << endl;
  cout << "a=1.111" << endl;
  a = 1.111;
  cout << a;

  // test copy constructor
  cout << endl;
  cout << "a=a" << endl;
  a = a;
  cout << a;

  cout << endl;
  cout << "a=a+1" << endl;
  a = a + 1;
  cout << a;

  cout << endl;
  cout << "allocate A, B, C" << endl;
  ss_vect<double>  A, B, C;

  cout << endl;
  cout << "allocate U, V, and W" << endl;
  ss_vect<tps>     U, V, W;

  cout << endl;
  cout << "U.zero()" << endl;
  U.zero();
  cout << U;

  cout << endl;
  cout << "V.identity()" << endl;
  V.identity();
  cout << V;

  cout << endl;
  cout << "V += V" << endl;
  V += V;
  cout << V;

  cout << endl;
  cout << "V = V + V" << endl;
  V = V + V;
  cout << V;

  cout << endl;
  cout << "V = V" << endl;
  V = V;
  cout << V;
}
