for m in `ls disk`; do
awk -v OFS='\t' '{$2=$3=$7=$8=$13=""; $4=int($4); $5=int($5); $6=sprintf("%.3f", $6); $9=sprintf("%.3f", $9); $10=sprintf("%.3f", $10); $11=int($11); $12=int($12); f=$6; fx=f*cos($9+$10); fy=f*sin($9+$10); fn=f*cos($9); ft=f*sin($9); fx=sprintf("%.3f", fx); fy=sprintf("%.3f", fy); fn=sprintf("%.3f", fn); ft=sprintf("%.3f", ft); frame=110;Dm=8; Dpx=98; mmpPx=Dm/Dpx; alpha=$9; beta=$10; cx=int($4+$11); cy=int($5+$12); t=$1; x=int($4); y=int($5); print t,x,y,f,fx,fy,fn,ft,alpha,beta,cx,cy}' disk/$m/$m > disk/$m/txyF
done
