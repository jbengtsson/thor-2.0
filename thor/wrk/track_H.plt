ps = 0; eps = 0;

if (!ps) set terminal x11;
if (ps && !eps) \
  set terminal postscript enhanced color solid lw 2 "Times-Roman" 20;
if (ps && eps) \
  set terminal postscript eps enhanced color solid lw 2 "Times-Roman" 20;

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

min(a, b) = (a < b)? a : b;

if (ps) set output "track_H_1.ps"
set title "J_x";
set xlabel "{/Symbol f}_x";
set ylabel "[mm{/Symbol \327}mrad]";
set yrange [0:];
plot "track_H_1.dat" using 7:(min($6, 50.0)) notitle with points ls 1, \
     "track_H_2.dat" using 7:(min($6, 50.0)) notitle with points ls 1, \
     "track_H_3.dat" using 7:(min($6, 50.0)) notitle with points ls 1, \
     "track_H_4.dat" using 7:(min($6, 50.0)) notitle with points ls 1, \
     "track_H_5.dat" using 7:(min($6, 50.0)) notitle with points ls 1, \
     "track_H_6.dat" using 7:(min($6, 50.0)) notitle with points ls 1, \
     "track_H_7.dat" using 7:(min($6, 50.0)) notitle with points ls 1, \
     "track_H_1.dat" using 11:(min($10, 50.0)) notitle with points ls 2, \
     "track_H_2.dat" using 11:(min($10, 50.0)) notitle with points ls 2, \
     "track_H_3.dat" using 11:(min($10, 50.0)) notitle with points ls 2, \
     "track_H_4.dat" using 11:(min($10, 50.0)) notitle with points ls 2, \
     "track_H_5.dat" using 11:(min($10, 50.0)) notitle with points ls 2, \
     "track_H_6.dat" using 11:(min($10, 50.0)) notitle with points ls 2, \
     "track_H_7.dat" using 11:(min($10, 50.0)) notitle with points ls 2;
if (!ps) pause -1;

if (ps) set output "track_H_2.ps"
set title "J_y";
set xlabel "{/Symbol f}_y";
set ylabel "[mm{/Symbol \327}mrad]";
set yrange [0:];
plot "track_H_1.dat" using 9:(min($8, 20.0)) notitle with points ls 3, \
     "track_H_2.dat" using 9:(min($8, 20.0)) notitle with points ls 3, \
     "track_H_3.dat" using 9:(min($8, 20.0)) notitle with points ls 3, \
     "track_H_4.dat" using 9:(min($8, 20.0)) notitle with points ls 3, \
     "track_H_5.dat" using 9:(min($8, 20.0)) notitle with points ls 3, \
     "track_H_6.dat" using 9:(min($8, 20.0)) notitle with points ls 3, \
     "track_H_7.dat" using 9:(min($8, 20.0)) notitle with points ls 3, \
     "track_H_1.dat" using 13:(min($12, 20.0)) notitle with points ls 1, \
     "track_H_2.dat" using 13:(min($12, 20.0)) notitle with points ls 1, \
     "track_H_3.dat" using 13:(min($12, 20.0)) notitle with points ls 1, \
     "track_H_4.dat" using 13:(min($12, 20.0)) notitle with points ls 1, \
     "track_H_5.dat" using 13:(min($12, 20.0)) notitle with points ls 1, \
     "track_H_6.dat" using 13:(min($12, 20.0)) notitle with points ls 1, \
     "track_H_7.dat" using 13:(min($12, 20.0)) notitle with points ls 1;
if (!ps) pause -1;

if (ps) set output "track_H_3.ps"
set title "H";
set xlabel "Turn no";
set ylabel "1e-6";
set yrange [0.01:];
set logscale y;
plot "track_H_1.dat" using 1:(min($14, 100.0)) notitle with lines ls 1, \
     "track_H_2.dat" using 1:(min($14, 100.0)) notitle with lines ls 1, \
     "track_H_3.dat" using 1:(min($14, 100.0)) notitle with lines ls 1, \
     "track_H_4.dat" using 1:(min($14, 100.0)) notitle with lines ls 1, \
     "track_H_5.dat" using 1:(min($14, 100.0)) notitle with lines ls 1, \
     "track_H_6.dat" using 1:(min($14, 100.0)) notitle with lines ls 1, \
     "track_H_7.dat" using 1:(min($14, 100.0)) notitle with lines ls 1;
if (!ps) pause -1;
