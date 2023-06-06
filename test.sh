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

    # image with WSClean
    wsclean -name "${name,,}" -size 128 128 -niter 1000000 -auto-mask 10.0 -auto-threshold 0.3 -mgain 0.9 -gain 0.2 -weight briggs 0.0 -scale 0.2asec -taper-gaussian 1.5asec -join-channels -channels-out 9 -fit-spectral-pol 4 -deconvolution-channels 9 -deconvolution-threads $ncpu -parallel-gridding $ncpu -parallel-reordering $ncpu -shift "$ra" "$dec" -no-negative -save-source-list "../$msres"

    rm -rf *-0*-*.fits

    cd ../
done