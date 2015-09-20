ps = 0; eps = 0;

N = 15;

if (!ps) set terminal x11;
#if (!ps) set terminal x11 dashed;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

file_name = "`echo $TRACY_LIB`/gnuplot/jet.dat";
# Load 64-color palette for Jet
set palette model RGB file file_name using ($1/255):($2/255):($3/255);

set cntrparam level 25;
set view map;
set pm3d;

if (ps) set output "nu2_1.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "{/Symbol n}_x";
set xlabel "2J_x"; set ylabel "2J_y";
splot "get_dnu2_1.out" using (1e6*$1):(1e6*$2):(log(abs(N*$3))) notitle \
      w l lt palette z;

set origin 0.0, 0.0;
set title "{/Symbol n}_y";
set xlabel "2J_x"; set ylabel "2J_y";
splot "get_dnu2_1.out" using (1e6*$1):(1e6*$2):(log(abs(N*$4))) notitle \
      w l lt palette z;

unset multiplot;
if (!ps) pause(-1);

if (ps) set output "nu2_1.ps"

set multiplot;

set size 1.0, 0.5; set origin 0.0, 0.5;
set title "{/Symbol n}_x";
set xlabel "{/Symbol d}"; set ylabel "2J_x";
splot "get_dnu2_2.out" using (1e2*$1):(1e6*$2):(log(abs(N*$3))) notitle \
      w l lt palette z;

set origin 0.0, 0.0;
set title "{/Symbol n}_y";
set xlabel "{/Symbol d}"; set ylabel "2J_x";
splot "get_dnu2_2.out" using (1e2*$1):(1e6*$2):(log(abs(N*$4))) notitle \
      w l lt palette z;

unset multiplot;
if (!ps) pause(-1);
