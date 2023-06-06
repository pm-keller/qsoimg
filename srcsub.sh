msres=$1
msmod=$2
automask=$4
autothresh=$5
ncpu=$6

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
pythonpath=/opt/conda/envs/minicasa/bin/python3

# read source coordinates
readarray srcs < "$3"

# image and subtract sources
cnt=0

for src in "${srcs[@]}"
do
    ((cnt++))

    # format string
    src=${src%$'\n'}
    src=${src%$'\r'}

    # separate RA and DEC
    ra=$(echo $src | cut -f1 -d ' ')
    dec=$(echo $src | cut -f2 -d ' ')
    size=$(echo $src | cut -f3 -d ' ')

    # make source name "SRC??"
    name=$(printf "SRC%02d" $cnt)
    rm -rf $name
    mkdir $name
    cd $name

    echo "$ra $dec" 

    $pythonpath -c "import casatasks; casatasks.phaseshift(vis='../$msres', outputvis='shifted.ms', phasecenter='J2000 $ra $dec')"
    
    # image with WSClean
    wsclean -name "${name,,}" \
    -size $size $size \
    -niter 1000 \
    -auto-mask $automask \
    -auto-threshold $autothresh \
    -mgain 0.8 \
    -gain 0.1 \
    -weight briggs -0.5 \
    -scale 0.3asec \
    -taper-gaussian 1.5asec \
    -use-wgridder \
    -join-channels \
    -channels-out 9 \
    -deconvolution-channels 9 \
    -fit-spectral-pol 3 \
    -deconvolution-threads $ncpu \
    -parallel-gridding $ncpu \
    -parallel-reordering $ncpu \
    -no-update-model-required \
    -save-source-list \
    "shifted.ms"

    rm -r shifted.ms
    rm -r *-0*-*.fits

    # predict model and add to model column
    for i in $(seq 0 8);
    do
        DP3 msin="../$msres" msout="../$msres" msin.band=$i msin.datacolumn="MODEL_DATA" msout.datacolumn="MODEL_DATA" predict.operation="replace" steps=[predict] predict.sourcedb="${name,,}-sources.txt" numthreads=$ncpu
    done 

    cd ../

    $pythonpath "${SCRIPT_DIR}/addmod.py" --take $msres --add $msmod

    # subtract source model from corrected visibilities
    $pythonpath -c "import casatasks; casatasks.uvsub('$msres')"

done