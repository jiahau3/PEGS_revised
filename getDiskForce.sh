#!/bin/bash
# get each disk information from .csv through matching the xy value with tracked disk file and obtaining the same time.
fn=$1;
[[ -d disk ]] || mkdir disk 
for file in `ls track`;
do [[ -d disk/disk$file ]] || mkdir disk/disk$file;
done
for m in `ls track`; do for a in `awk '{print $1,$2,$3}' track/$m | tr ' ' 'x'`; do t=$((${a##*x})); xy=${a%x*}; awk 'NR>1' $fn | tr ',' ' ' | cat  | awk '{printf("%sx%s ",int($4),int($5));print}' | grep $xy | awk -v t=$t '$2==t'; done | awk '{for(i=2; i<NF; i++) printf "%s ",$i; print $NF}' > disk/disk$m/disk$m; done
