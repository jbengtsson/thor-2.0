ps = 0; eps = 0; contour = 1; N = 15;

res_1 = 1; res_2 = 1;

# working point
nu_x = 33.47; nu_y = 15.66;
set label "{/Symbol \264}" at nu_x, nu_y centre textcolor lt 3;
set label "{/Symbol \264}" at nu_x, nu_y centre textcolor lt 3;
nu_x = 33.15; nu_y = 15.69;
set label "{/Symbol \264}" at nu_x, nu_y centre textcolor lt 3;
set label "{/Symbol \264}" at nu_x, nu_y centre textcolor lt 3;
nu_x = 33.15; nu_y = 16.30;
set label "{/Symbol \264}" at nu_x, nu_y centre textcolor lt 3;
set label "{/Symbol \264}" at nu_x, nu_y centre textcolor lt 3;

if (!ps) set terminal x11;
#if (!ps) set terminal x11 dashed;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lw 1 lc rgb "red";
set style line 2 lw 1 lc rgb "cyan";
set style line 3 lw 1 lc rgb "blue";

if (contour) \
  set nosurface; \
  set noztics;

set clabel "%4.1f"; set key left;
#unset clabel;
set noztics; unset colorbox;
#set cbrange [0:1];

#set cntrparam level 18;
set cntrparam level 25;

# x <-> horizontal, y <-> vertical, z <-> perpendicular to screen
# rot_x, rot_z, scale, scale_z
if (!contour) set view 65, 15, 1, 1;
if (contour) set view 0, 0, 1, 1;

set palette rgbformulae 22, 13, -31 negative;

if (ps) set output "tune_scan.ps"

set title "Normalized Dynamic Aperture";
set xlabel "{/Symbol n}_x"; set ylabel "{/Symbol n}_y";

xmin = 2.15; xmax = 2.25; ymin = 1.01; ymax = 1.1;
#xmin = 2.15; xmax = 2.3; ymin = 1.01; ymax = 1.1; \
if (N == 1) \
  scl = 1.0/15.0;
if (N != 1) \
  scl = 1.0; \
  xmin = N*xmin; xmax = N*xmax; ymin = N*ymin; ymax = N*ymax;

nx = int(xmin/scl); ny = int(ymin/scl);

set xrange [xmin:xmax]; set yrange [ymin:ymax];

r10(n) = scl*n; r01(n) = scl*n;

r11(x, n) = scl*n - x; r1m1(x, n) = x - scl*n;
r11_inv(y, n) = scl*n - y; r1m1_inv(y, n) = scl*n + y;

r20(n) = scl*n/2.0; r02(n) = scl*n/2.0;

r30(n) = scl*n/3.0;
r12(x, n) = (scl*n-x)/2.0; r1m2(x, n) = (x-scl*n)/2.0;
r12_inv(y, n) = scl*n - 2.0*y; r1m2_inv(y, n) = scl*n + 2.0*y;

r40(n) = scl*n/4.0; r04(n) = scl*n/4.0;
r22(x, n) = scl*n/2.0 - x; r2m2(x, n) = x - scl*n/2.0;
r22_inv(y, n) = scl*n/2.0 - y; r2m2_inv(y, n) = scl*n/2.0 + y;

#  set arrow from r10(nx+2), ymin to r10(nx+2), ymax nohead ls 1; \
#  set arrow from r20(2*nx+3), ymin to r20(2*nx+3), ymax nohead ls 1; \
#  set arrow from r30(3*nx+7), ymin to r30(3*nx+7), ymax nohead ls 1; \
#  set arrow from r40(4*nx+9), ymin to r40(4*nx+9), ymax nohead ls 3; \

if (res_1) \
  set arrow from r10(nx+1), ymin to r10(nx+1), ymax nohead ls 1; \
  set arrow from r20(2*nx+1), ymin to r20(2*nx+1), ymax nohead ls 1; \
  set arrow from r30(3*nx+1), ymin to r30(3*nx+1), ymax nohead ls 1; \
  set arrow from r30(3*nx+2), ymin to r30(3*nx+2), ymax nohead ls 1; \
  set arrow from r30(3*nx+4), ymin to r30(3*nx+4), ymax nohead ls 1; \
  set arrow from r30(3*nx+5), ymin to r30(3*nx+5), ymax nohead ls 1; \
\
  set arrow from xmin, r01(ny+1) to xmax, r01(ny+1) nohead ls 1; \
  set arrow from xmin, r02(2*ny+1) to xmax, r02(2*ny+1) nohead ls 1; \
  set arrow from xmin, r02(2*ny+3) to xmax, r02(2*ny+3) nohead ls 1; \
\
  set arrow from xmin, r11(xmin, nx+ny+1), ymax to \
    r11_inv(ymin, nx+ny+1), ymin nohead ls 2; \
  set arrow from r11_inv(ymax, nx+ny+2), ymax to \
    xmax, r11(xmax, nx+ny+2), ymin nohead ls 2; \
  set arrow from r11_inv(ymax, nx+ny+3), ymax to \
    xmax, r11(xmax, nx+ny+3) nohead ls 2; \
  set arrow from xmin, r1m1(xmin, nx-ny-1) to \
    r1m1_inv(ymax, nx-ny-1), ymax nohead ls 2; \
  set arrow from xmin, r1m1(xmin, nx-ny) to \
    r1m1_inv(ymax, nx-ny), ymax nohead ls 2; \
  set arrow from r1m1_inv(ymin, nx-ny+1), ymin to \
    xmax, r1m1(xmax, nx-ny+1), ymax nohead ls 2; \
\
  set arrow from xmin, r12(xmin, nx+2*ny+1) to \
     r12_inv(ymin, nx+2*ny+1), ymin nohead ls 1; \
  set arrow from xmin, r12(xmin, nx+2*ny+2) to \
     r12_inv(ymin, nx+2*ny+2), ymin nohead ls 1; \
  set arrow from xmin, r12(xmin, nx+2*ny+3) to \
    xmax, r12(xmax, nx+2*ny+3) nohead ls 1; \
  set arrow from r12_inv(ymax, nx+2*ny+4), ymax to \
    xmax, r12(xmax, nx+2*ny+4) nohead ls 1; \
  set arrow from xmin, r1m2(xmin, nx-2*ny-2) to \
    r1m2_inv(ymax, nx-2*ny-2), ymax nohead ls 1; \
  set arrow from xmin, r1m2(xmin, nx-2*ny-1) to \
    xmax, r1m2(xmax, nx-2*ny-1) nohead ls 1; \
  set arrow from r1m2_inv(ymin, nx-2*ny), ymin to \
    xmax, r1m2(xmax, nx-2*ny) nohead ls 1; \
  set arrow from r1m2_inv(ymin, nx-2*ny+1), ymin to \
    xmax, r1m2(xmax, nx-2*ny+1) nohead ls 1;

# redundant
#set arrow from xmin, r22(xmin, 98) to r22_inv(ymin, 98), ymin nohead ls 3;
#set arrow from r22_inv(ymax, 100), ymax to xmax, r22(xmax, 100) nohead ls 3;
#set arrow from xmin, r2m2(xmin, 34) to r2m2_inv(ymax, 34), ymax nohead ls 3;
#set arrow from r2m2_inv(ymin, 36), ymin to xmax, r2m2(xmax, 36) nohead ls 3;
#  set arrow from r22_inv(ymax, 2*nx+2*ny+7), ymax to \
#    xmax, r22(xmax, 2*nx+2*ny+7) nohead ls 3; \

if (res_2) \
  set arrow from r40(4*nx+1), ymin to r40(4*nx+1), ymax nohead ls 3; \
  set arrow from r40(4*nx+3), ymin to r40(4*nx+3), ymax nohead ls 3; \
  set arrow from r40(4*nx+5), ymin to r40(4*nx+5), ymax nohead ls 3; \
  set arrow from r40(4*nx+6), ymin to r40(4*nx+6), ymax nohead ls 3; \
  set arrow from r40(4*nx+7), ymin to r40(4*nx+7), ymax nohead ls 3; \
  set arrow from xmin, r04(4*ny+1) to xmax, r04(4*ny+1) nohead ls 1; \
  set arrow from xmin, r04(4*ny+3) to xmax, r04(4*ny+3) nohead ls 1; \
  set arrow from xmin, r04(4*ny+5) to xmax, r04(4*ny+5) nohead ls 1; \
  set arrow from xmin, r22(xmin, 2*nx+2*ny+3) to \
    r22_inv(ymin, 2*nx+2*ny+3), ymin nohead ls 3; \
  set arrow from r22_inv(ymax, 2*nx+2*ny+5), ymax to \
    xmax, r22(xmax, 2*nx+2*ny+5), ymin nohead ls 3; \
  set arrow from xmin, r2m2(xmin, 2*nx-2*ny-1) to \
    r2m2_inv(ymax, 2*nx-2*ny-1), ymax nohead ls 3; \
  set arrow from r2m2_inv(ymin, 2*nx-2*ny+1), ymin to \
    xmax, r2m2(xmax, 2*nx-2*ny+1) nohead ls 3; \
  set arrow from r2m2_inv(ymin, 2*nx-2*ny+3), ymin to \
    xmax, r2m2(xmax, 2*nx-2*ny+3) nohead ls 3;

splot "tune_scan.dat" using ($5*N):($7*N):9 notitle with lines lt palette z;

if (!ps) pause(-1);
