#!/bin/bash
#SBATCH -J peel
#SBATCH --no-requeue
#SBATCH -c 32
#SBATCH -t 00:40:00
#SBATCH -p skylake

snr=10
selfcal_rnd=2
obsname="QSO-J1429-0104"

#SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_DIR="/rds/user/pmk46/hpc-work/19A-056/pipeline"
ROOT="${SCRIPT_DIR}/.."
SELFCALDIR="${ROOT}/imaging/${obsname}/selfcal/"
msorig="${ROOT}/ms/${obsname}.ms"
imfits="${ROOT}/imaging/${obsname}/clean/${obsname}-selfcal-${selfcal_rnd}.im-MFS-image.fits"
srcs="${ROOT}/imaging/${obsname}/clean/${obsname}-selfcal-${selfcal_rnd}.im-sources-all.txt"
first="${ROOT}/regions/FIRST/${obsname}-regions.txt"
solint="${SELFCALDIR}/solint.pickle"

pythonpath="${SCRIPT_DIR}/casa/bin/python3"
container="apptainer exec --pid /rds/user/pmk46/hpc-work/19A-056/pipeline/work/singularity/pmkeller-imaging-tools.img"
cnt=0

name=$(printf "SRC%02d" $cnt)
#rm -rf $name
#mkdir $name
cd $name

msres="${obsname}-res.ms"
ms="${obsname}.ms"

#cp -rf $msorig $ms
#
#$pythonpath ../preppeel.py --ms $ms --imfits $imfits --srcs $srcs --first $first --selfcal $SELFCALDIR --rnds $(($selfcal_rnd + 1))
#
## predict sources
#for i in $(seq 0 8);
#do
#    $container DP3 msin="$ms" msout="$ms" msin.band=$i msin.datacolumn="MODEL_DATA" msout.datacolumn="MODEL_DATA" predict.operation="replace" steps=[predict] predict.sourcedb="srcsout.txt"
#done 

readarray -t region_array < "region.txt"
coord="${region_array[0]}"
size="${region_array[1]}"
size=$((size*8))

echo "${coord}, ${size}"

## subtract source model from corrected visibilities
#$pythonpath -c "import casatasks; 
#casatasks.uvsub('$ms'); 
#casatasks.split('$ms', '$msres', datacolumn='corrected'); 
#casatasks.phaseshift(vis='$msres', outputvis='shifted.ms', phasecenter='J2000 $coord')
#casatasks.uvsub('$ms', reverse=True);"
#
## image with WSClean
#$container wsclean -name "${name,,}" \
#  -size $size $size \
#  -niter 1000 \
#  -threshold 0.00015 \
#  -mgain 0.8 \
#  -gain 0.1 \
#  -weight briggs -1.0 \
#  -scale 0.2asec \
#  -taper-gaussian 1.0asec \
#  -use-wgridder \
#  -multiscale \
#  -join-channels \
#  -channels-out 9 \
#  -deconvolution-channels 9 \
#  -fit-spectral-pol 3 \
#  -no-update-model-required \
#  -save-source-list \
#  -abs-mem 6 \
#  "shifted.ms"
#
#rm -r shifted.ms
#rm -r *-0*-*.fits
#
#$pythonpath -c "import casatasks; casatasks.setjy('$msres', fluxdensity=[0, 0, 0, 0], usescratch=True)"
#
## predict sources
#for i in $(seq 0 8);
#do
#    $container DP3 msin="$msres" msout="$msres" msin.band=$i msin.datacolumn="MODEL_DATA" msout.datacolumn="MODEL_DATA" predict.operation="replace" steps=[predict] predict.sourcedb="${name,,}-sources.txt"
#done 
#
#solintmin=$(< "tself.txt")
#$pythonpath ../optsolint.py --ms $msres --solint_min $solintmin
#$pythonpath ../selfcal.py --ms $msres --solintfile $solint --calmode "ap"
#$pythonpath ../applycal.py --ms $msres --dir "" --nself 1
#
#msself="${obsname}-res-selfcal-0.ms"
#$pythonpath -c "import casatasks; casatasks.phaseshift(vis='$msself', outputvis='shifted.ms', phasecenter='J2000 $coord')"
#$container wsclean -name "${name,,}-selfcal" \
#  -size $size $size \
#  -niter 1000 \
#  -threshold 0.00015 \
#  -mgain 0.8 \
#  -gain 0.1 \
#  -weight briggs -1.0 \
#  -scale 0.2asec \
#  -taper-gaussian 1.0asec \
#  -use-wgridder \
#  -multiscale \
#  -join-channels \
#  -channels-out 9 \
#  -deconvolution-channels 9 \
#  -fit-spectral-pol 3 \
#  -no-update-model-required \
#  -save-source-list \
#  -abs-mem 6 \
#  "shifted.ms"
#
#rm -r shifted.ms
#rm -r *-0*-*.fits
#
#$pythonpath -c "import casatasks; casatasks.setjy('$msself', fluxdensity=[0, 0, 0, 0], usescratch=True)"
#
## predict sources
#for i in $(seq 0 8);
#do
#    $container DP3 msin="$msself" msout="$msself" msin.band=$i msin.datacolumn="MODEL_DATA" msout.datacolumn="MODEL_DATA" predict.operation="replace" steps=[predict] predict.sourcedb="${name,,}-selfcal-sources.txt"
#done 

mscp="${obsname}-cp.ms"
msself="${obsname}-res-selfcal-0.ms"
#cp -rf $msres $mscp
#$pythonpath ../putmod.py --take $msres --add $msself
#$pythonpath ../peel.py --ms $mscp

$pythonpath -c "import casatasks; casatasks.phaseshift(vis='$msself', outputvis='shifted.ms', phasecenter='J2000 $coord'); casatasks.uvsub('shifted.ms')"
#$pythonpath ../peel.py --ms "shifted.ms"

$container wsclean -name "${name,,}-peeled" \
  -size 128 128 \
  -niter 1000 \
  -auto-mask 3 \
  -auto-threshold 0.3 \
  -mgain 0.8 \
  -gain 0.1 \
  -weight briggs -1.0 \
  -scale 0.3asec \
  -taper-gaussian 1.5asec \
  -use-wgridder \
  -multiscale \
  -join-channels \
  -channels-out 9 \
  -deconvolution-channels 9 \
  -fit-spectral-pol 3 \
  -no-update-model-required \
  -save-source-list \
  -abs-mem 6 \
  "shifted.ms"

rm -r shifted*

#
# make final image
#$container wsclean -name "${imname}" \
#  -size 12288 12288 \
#  -niter 100000 \
#  -auto-mask 5 \
#  -auto-threshold 1 \
#  -mgain 0.8 \
#  -gain 0.1 \
#  -weight briggs -0.5 \
#  -scale 0.3arcsec \
#  -taper-gaussian 1.5asec \
#  -multiscale \
#  -use-wgridder \
#  -join-channels \
#  -channels-out 9 \
#  -deconvolution-channels 9 \
#  -fit-spectral-pol 3 \
#  -parallel-deconvolution 2048 \
#  -deconvolution-threads 32 \
#  -parallel-gridding 32 \
#  -parallel-reordering 32 \
#  -no-update-model-required \
#  -save-source-list \
#  $ms