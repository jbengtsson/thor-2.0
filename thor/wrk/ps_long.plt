ps = 0;

f_s = 14; l_w = 2;
if (ps == 0) \
  set terminal x11; \
else if (ps == 1) \
  set terminal postscript enhanced color solid lw l_w "Times-Roman" f_s; \
  ext = "ps"; \
else if (ps == 2) \
  set terminal postscript eps enhanced color solid lw l_w "Times-Roman" f_s; \
  ext = "eps"; \
else if (ps == 3) \
  set terminal pdf enhanced color solid lw l_w font "Times-Roman f_s"; \
  ext = "pdf"; \
else if (ps == 4) \
  set term pngcairo enhanced color solid lw l_w font "Times-Roman f_s"; \
  ext = "png";

set grid;

set style line 1 lt 1 lw 1 lc rgb "blue";
set style line 2 lt 1 lw 1 lc rgb "green";
set style line 3 lt 1 lw 1 lc rgb "red";

if (ps) set output "ps_long.".(ext)

set multiplot;

set size 0.5, 0.5; set origin 0.0, 0.5;
set title "Hor Phase-Space";
set xlabel "x [m]"; set ylabel "p_x [rad]";
#set format x "%.1e"; set format y "%.1e";
#set xtics 5e-6; set ytics 2.5e-6;
#set xtics 1e-5; set ytics 5e-6;
plot "ps_long.out" using 2:3 notitle with points ls 1;

set origin 0.5, 0.5;
set title "Ver Phase-Space";
set xlabel "y [m]"; set ylabel "p_y [rad]";
#set format x "% g"; set format y "% g";
#set xtics 1e-5; set ytics 1e-6;
#set xtics 5e-5; set ytics 2e-6;
plot "ps_long.out" using 4:5 notitle with points ls 3;

set origin 0.5, 0.0;
set title "Long Phase-Space";
set xlabel "ct [m]"; set ylabel "{/Symbol d}";
#set format x "%.1e"; set format y "%.1e";
#set xtics 1e-10; set ytics 2e-4;
#set xtics 5e-9; set ytics 2e-4;
#set xrange [-1e-10:1e-10];
#set xrange [-5e-9:5e-9];
#set xtics autofreq;
plot "ps_long.out" using 7:6 notitle with points ls 2;

unset multiplot;
if (!ps) pause mouse "click on graph to cont.\n";
