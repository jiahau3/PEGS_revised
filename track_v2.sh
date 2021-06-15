#!/bin/bash
# script : track.sh tmp/*
[[ -d track ]] || mkdir track;
awk -v r=8100 -v ymax=820 -v tmax=706 -v OFS='\t' 'function dsq(x1,y1,x2,y2){return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)} \
    NR==FNR{x[NR]=int($1);y[NR]=int($2);nb=NR; fn=sprintf("track/%03d",nb); print x[NR],y[NR],$3,0 > fn} \
    NR>FNR{ min=1e9; \
       for(i=1;i<=nb;i++){ if(int($2)>ymax) continue;\
           d2=dsq(x[i],y[i],$1,$2); \
           if(d2<min){imin=i;min=d2;}\
#if(FILENAME=="tmp/026")print "track",i,imin,min,x[i],y[i],$1,$2,$3,FILENAME; \
       } \
       if(min<r && $3<tmax) {fnn=sprintf("track/%03d",imin); print int($1),int($2),$3,min > fnn; x[imin]=$1; y[imin]=$2; } \
         else if(int($2)<ymax && $3<tmax) {nb++;x[nb]=int($1);y[nb]=int($2); fnnn=sprintf("track/%03d",nb); print x[nb],y[nb],$3,min > fnnn; } \
#       if(nb==9&&FILENAME=="tmp/026")print "*****************",imin,min,x[nb],y[nb],$1,$2,$3; \
       }' \
$*
