for m in `ls disk`; do
awk -v OFS='\t' '{t=$1; x[t]=$2; y[t]=$3; f=$4; fx[t]+=$5; fy[t]+=$6; torq[t]+=$8}END{for (i in fx) print i, x[i], y[i], fx[i], fy[i], torq[i]}' disk/$m/txyF | sort -n > disk/$m/totF_torq
done
