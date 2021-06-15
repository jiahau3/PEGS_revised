#!/bin/bash
# script : data_analysis.sh for analyzing calibrataion data
Dir=`basename $0 .sh`; [[ -d $Dir ]] || mkdir $Dir;
cd $Dir;
for d in raw center m; do [[ -d $d ]] || mkdir $d; done;

function cp2raw { cp ../$imgdir/???.png raw; }
function cp2m { cp ../../scripts/*{.m,.py} m; }
function fndcntr { cd m; python3 fndcntr.py $r $R $CS ../raw/ $png ../$cntdir/; cd ..;}
function getcontact { cd m; matlab $mopt "getcontact('../raw/$png', '../$cntctInD/', '../$cntctOutD/',$Dm, $Dpx, $g2guess, $FS, $DT, $conR, $cG2Thrsd, $ctrstL, $ctrstH, $shift4calibration); quit"; cd ..; }
function calForce { cd m; matlab $mopt "calForce('../$calFInD/$mat', '../$calFOutD/', $rMask, $tF1, $tF2, $tF3, $tF4, $g2_cal_a, $g2_cal_b, $exitX, $exitY, $D_arch2exit, $given_force, $calibrate, $optimization); quit"; cd ..; }
function ForceAdjMat { cd m; matlab $mopt "ForceAdjMat('../$FadjInD/$mat','../raw/','../$FadjOutD/', $E_starS, $E_starL, $Fc); quit"; cd ..; }

function action { return; }

## default setting ##
mopt="-nodesktop -nosplash -r"; 
imgdir=4b_0043_30s_fps20; 
r=48; R=52; CS=13; png=001.png; cntdir=center; 
cntctInD=center/txt; cntctOutD=s1; Dm=0.008; Dpx=100; g2guess=100; FS=6.25e4*2e-3; DT=5; conR=10; cG2Thrsd=0.0001; ctrstL=0.005; ctrstH=0.995; shift4calibration=0;
calFInD=s1/mat; calFOutD=s2; mat=001.mat; rMask=0.45; tF1=1.0; tF2=2.0; tF3=0.0; tF4=0.0; g2_cal_a=90/0.88; g2_cal_b=0; exitX=427; exitY=475; D_arch2exit=260; given_force=0; calibrate=0; optimization=1;
FadjInD=s2/maskR${rMask}_mat; FadjOutD=s3;
E_starS=3.76e6;
E_starL=2.28e6;
Fc=0.20;

while true; do case $1 in
      -imgdir)   imgdir=$2;       shift 2;;  # raw image directory
         -png)   png=$2;          shift 2;;  # .png filename
      -cntdir)   cntdir=$2;       shift 2;;  # create center directory containing txt/, png/
    -cntctInD)   cntctInD=$2;     shift 2;;  # center txt directory
   -cntctOutD)   cntctOutD=$2;    shift 2;;  # create directory for disk infomration (stage 1) (containing contact points, initial guess force)
          -Dm)   Dm=$2            shift 2;;  # Diameter of disk in meter
         -Dpx)   Dpx=$2           shift 2;;  # Diameter size in pixel
     -g2guess)   g2guess=$2;      shift 2;;  # give G2 value for estimating initial force per disk
          -FS)   FS=$2            shift 2;;  # photoelastic stress coefficient
          -DT)   DT=$2            shift 2;;  # distance tolerance(pixel) for considering as contacts
        -conR)   conR=$2          shift 2;;  # setting radius of contact circle 
    -cG2Thrsd)   cG2Thrsd=$2      shift 2;;  # G2 threshold for determing if a contact is valid
      -ctrstL)   ctrstL=$2        shift 2;;
      -ctrstH)   ctrstH=$2        shift 2;;
-shift4calibration) shift4calibration=$2 shift 2;;
     -calFInD)   calFInD=$2;      shift 2;;  # .mat file directory in stage 1 directory
         -mat)   mat=$2;          shift 2;;  # .mat filename
       -rMask)   rMask=$2;        shift 2;;  # ratio of mask radius to whole disk
    -tF)   shift; tF1=$1; tF2=$2; tF3=$3; tF4=$4;     shift 4;;  # adjust initial force for contact points > 2
    -g2_cal_a)   g2_cal_a=$2      shift 2;;
    -g2_cal_b)   g2_cal_b=$2      shift 2;;
       -exitX)   exitX=$2         shift 2;;
       -exitY)   exitY=$2         shift 2;;
 -D_arch2exit)   D_arch2exit=$2   shift 2;;
 -given_force)   given_force=$2   shift 2;;
   -calibrate)   calibrate=$2;    shift 2;;
-optimization)   optimization=$2  shift 2;;
    -calFOutD)   calFOutD=$2;     shift 2;;  # create directory for disk information (stage 2) (containing force images for each disk)
     -FadjInD)   FadjInD=$2;      shift 2;;  # .mat file directory in stage 2 directory
    -FadjOutD)   FadjOutD=$2;     shift 2;;  # create directory for disk information (stage 3) 
     -E_starS)   E_starS=$2;      shift 2;;
     -E_starL)   E_starL=$2;      shift 2;;
          -Fc)   Fc=$2;           shift 2;;
           -r)   r=$2;            shift 2;;  # inner radius
           -R)   R=$2;            shift 2;;  # outer radius
          -CS)   CS=$2;           shift 2;;
            *)   break;;
  esac
done

#echo imgdir=$imgdir r=$r R=$R f=$f 

while true; do case $1 in 
            cp2m)  cp2m;            shift  ;;
          cp2raw)  cp2raw;          shift  ;;
         fndcntr)  fndcntr;         shift  ;;
      findcenter)  findcenter       shift  ;;
      getcontact)  getcontact;      shift  ;;
        calForce)  calForce;        shift  ;;
     ForceAdjMat)  ForceAdjMat;     shift  ;;
          action)  action;          shift  ;;
               *) break;;
   esac
done
