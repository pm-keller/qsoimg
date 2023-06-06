import casatasks
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--ms', type=str, help='Path to measurement set.')
parser.add_argument('--mscp', type=str, help='Path to store copy of measurement set.')
parser.add_argument('--imname', type=str, help='Image Name.')
parser.add_argument('--outlierfile', type=str, help='Outlier File.')
args = parser.parse_args()

# make a copy of the measurement set
shutil.copytree(args.ms, args.mscp)

# fill model column
casatasks.tclean(args.mscp, 
    imagename=f"{args.imname}",
    imsize=8192, 
    cell="0.3arcsec", 
    pblimit=-0.1,
    niter=0, 
    gridder="wproject",
    weighting="briggs",
    robust=0.0,
    #aterm=False,
    #psterm=True,
    #cfcache=f"{args.imname}.cf",
    uvtaper="140klambda",
    #uvrange=["3~200 klambda"],
    wprojplanes=128,
    specmode="mfs",
    deconvolver="mtmfs",
    nterms=2,
    outlierfile=args.outlierfile,
    savemodel="modelcolumn",
    calcres=False, 
    calcpsf=False
)