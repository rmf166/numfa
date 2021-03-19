#!/bin/sh
rm final.pdf
sols="DD SC LD LC"
for sol in $sols
do
  for sw in 1 2 3 4
  do
    for prb in 1 2 3 4 5 6
    do
      rm plot.p
      echo 'set autoscale' >> plot.p
      echo 'unset logscale; set logscale x' >> plot.p
      echo 'unset label' >> plot.p
      echo 'set xtic auto' >> plot.p
      echo 'set ytic auto' >> plot.p
      if [ ${sw} == 1 ]
      then
        echo 'set title "'${sol}', p = '${prb}'"' >> plot.p
      else
        echo 'set title "'${sol}', p = '${prb}', s = '${sw}'"' >> plot.p
      fi
      echo 'set xlabel "{/Symbol s_T}h" enhanced' >> plot.p
      echo 'set ylabel "{/Symbol r}" enhanced' >> plot.p
      echo 'set yr [0:1]' >> plot.p
      echo 'set xr [0.01:100]' >> plot.p
      echo 'plot    "numres-p'${prb}'-s'${sw}'-'${sol}'.dat" using 1:2 notitle with points pointtype 1 ps 1 lc rgb "red", \' >> plot.p
      echo '        "numres-p'${prb}'-s'${sw}'-'${sol}'.dat" using 1:3 notitle with points pointtype 1 ps 1 lc rgb "green", \' >> plot.p
      echo '        "numres-p'${prb}'-s'${sw}'-'${sol}'.dat" using 1:4 notitle with points pointtype 1 ps 1 lc rgb "blue", \' >> plot.p
      echo '        "numres-p'${prb}'-s'${sw}'-'${sol}'.dat" using 1:5 notitle with points pointtype 1 ps 1 lc rgb "violet", \' >> plot.p
      echo '        "result-p'${prb}'-s'${sw}'-'${sol}'.dat" using 1:2 title "c=0.8"  with lines linetype 1 lc rgb "red", \' >> plot.p
      echo '        "result-p'${prb}'-s'${sw}'-'${sol}'.dat" using 1:3 title "c=0.9"  with lines linetype 1 lc rgb "green", \' >> plot.p
      echo '        "result-p'${prb}'-s'${sw}'-'${sol}'.dat" using 1:4 title "c=0.99" with lines linetype 1 lc rgb "blue", \' >> plot.p
      echo '        "result-p'${prb}'-s'${sw}'-'${sol}'.dat" using 1:5 title "c=0.9999" with lines linetype 1 lc rgb "violet"' >> plot.p
      echo 'set key left top' >> plot.p
      echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
      echo 'set output "plot-'${prb}'-'${sw}'-'${sol}'.pdf"' >> plot.p
      echo 'replot' >> plot.p
      echo 'set terminal x11' >> plot.p
      gnuplot plot.p
    done
  done
done
for sol in $sols
do
  for sw in 1 2 3 4
  do
    for prb in 1 2 3 4 5 6
    do
      rm plot.p
      echo 'set autoscale' >> plot.p
      echo 'unset logscale; set logscale x' >> plot.p
      echo 'unset label' >> plot.p
      echo 'set xtic auto' >> plot.p
      echo 'set ytic auto' >> plot.p
      if [ ${sw} == 1 ]
      then
        echo 'set title "'${sol}'(lp), p = '${prb}'"' >> plot.p
      else
        echo 'set title "'${sol}'(lp), p = '${prb}', s = '${sw}'"' >> plot.p
      fi
      echo 'set xlabel "{/Symbol s_T}h" enhanced' >> plot.p
      echo 'set ylabel "{/Symbol r}" enhanced' >> plot.p
      echo 'set yr [0:1]' >> plot.p
      echo 'set xr [0.01:100]' >> plot.p
      echo 'plot    "numres-lp'${prb}'-s'${sw}'-'${sol}'.dat" using 1:2 notitle with points pointtype 1 ps 1 lc rgb "red", \' >> plot.p
      echo '        "numres-lp'${prb}'-s'${sw}'-'${sol}'.dat" using 1:3 notitle with points pointtype 1 ps 1 lc rgb "green", \' >> plot.p
      echo '        "numres-lp'${prb}'-s'${sw}'-'${sol}'.dat" using 1:4 notitle with points pointtype 1 ps 1 lc rgb "blue", \' >> plot.p
      echo '        "numres-lp'${prb}'-s'${sw}'-'${sol}'.dat" using 1:5 notitle with points pointtype 1 ps 1 lc rgb "violet", \' >> plot.p
      echo '        "result-lp'${prb}'-s'${sw}'-'${sol}'.dat" using 1:2 title "c=0.8"  with lines linetype 1 lc rgb "red", \' >> plot.p
      echo '        "result-lp'${prb}'-s'${sw}'-'${sol}'.dat" using 1:3 title "c=0.9"  with lines linetype 1 lc rgb "green", \' >> plot.p
      echo '        "result-lp'${prb}'-s'${sw}'-'${sol}'.dat" using 1:4 title "c=0.99" with lines linetype 1 lc rgb "blue", \' >> plot.p
      echo '        "result-lp'${prb}'-s'${sw}'-'${sol}'.dat" using 1:5 title "c=0.9999" with lines linetype 1 lc rgb "violet"' >> plot.p
      echo 'set key left top' >> plot.p
      echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
      echo 'set output "plot-l'${prb}'-'${sw}'-'${sol}'.pdf"' >> plot.p
      echo 'replot' >> plot.p
      echo 'set terminal x11' >> plot.p
      gnuplot plot.p
    done
  done
done
pdfunite plot*DD.pdf plot*SC.pdf plot*LD.pdf plot*LC.pdf final.pdf
