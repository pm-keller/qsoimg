params.root = projectDir.resolve('../')
params.msdir = params.root.resolve('ms')
params.aostrat = projectDir.resolve('aoflagger_strategy.lua')
params.srcs = projectDir.resolve('')
params.pythonpath = [container: "/opt/conda/envs/minicasa/bin/python3", projdir: "$projectDir/casa/bin/python3"]
params.imsize = 12288 // pixels
params.cell = 0.3 // cell size
params.imreg = params.imsize * params.cell / (2 * 60) 
// radius in arcmin of region to be imaged

process mkdir {
  cache = 'lenient'
  errorStrategy 'ignore'
  
  input:
  path ms

  script:
  obsname = ms.getSimpleName()
  imagingDir = params.root.resolve("imaging/")
  cleanDir = imagingDir.resolve("${obsname}/clean/")
  flaggingDir = imagingDir.resolve("${obsname}/flagging/")
  selfcalDir = imagingDir.resolve("${obsname}/selfcal/")

  """
  mkdir -p $imagingDir
  mkdir -p $cleanDir
  mkdir -p $flaggingDir
  mkdir -p $selfcalDir
  """
}

process aoflagger {
  cache = 'lenient'
  errorStrategy 'ignore'
  
  input:
  path ms
  path aostrat

  output:
  path ms

  when:
  obsname = ms.getSimpleName()
  finished = file(projectDir.resolve('finished.txt'))
  obsnames = finished.readLines()
  !(obsnames.contains(obsname))

  script:
  flaggingDir = params.root.resolve("imaging/${obsname}/flagging/")

  """
  # aoflagger removes flags on completely flagged columns.
  # Save the current flags and merge them later with the flags from aoflagger.
  if [ -d "$params.msdir/${obsname}.flagversions/flags.pre-ao"]
  then 
    $params.pythonpath.container -c "import casatasks; casatasks.flagmanager('$ms', mode='restore', versionname='before_ao_flags')"
  else
    $params.pythonpath.container $projectDir/saveflags.py --ms $ms --name "ao" --when "before"
  fi

  # run aoflagger
  aoflagger -v -skip-flagged -strategy $aostrat $ms 

  # merge flags
  $params.pythonpath.container -c "import casatasks; casatasks.flagmanager('$ms', mode='restore', versionname='before_ao_flags', merge='or')"

  # save flags
  $params.pythonpath.container $projectDir/saveflags.py --ms $ms --name "ao" --when "after"

  # copy results to field flagging directory
  cp -rf *npy $flaggingDir
  cp -rf *png $flaggingDir
  """
}

process wscleanInspect {
  publishDir {params.root.resolve("inspect/${obsname}")}, mode: 'copy'
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms

  output:
  path imfile

  when:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  imname = "${obsname}.inspect"
  imfile = file(params.root.resolve("inspect/${obsname}/${imname}-MFS-image.fits"))
  !(imfile.exists())  

  """
  wsclean -name '${imname}' \
  -size 8192 8192 \
  -join-channels \
  -channels-out 3 \
  -scale 1asec \
  -taper-gaussian 1asec \
  -use-wgridder \
  -parallel-gridding 4 \
  -weight briggs 0.0 \
  $ms

  rm -rf *.im-000*-*.fits
  rm -rf *.im*dirty.fits
  rm -rf *.im*psf.fits
  """
}


process wscleanDirty {
  cpus {"${ncpu}"}
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms
  val cell
  val size
  val briggs
  val ncpu

  output:
  path ms, emit: ms
  path '*image.fits', emit: dirty

  when:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  finished = file(projectDir.resolve('finished.txt'))
  obsnames = finished.readLines()
  !(obsnames.contains(obsname))

  script:
  imname = "${ms.getSimpleName()}.im"
  cleanDir = params.root.resolve("imaging/${obsname}/clean/")

  """
  wsclean -name "${imname}" \
  -size $size $size \
  -niter 0 \
  -weight briggs $briggs \
  -scale "${cell}arcsec" \
  -taper-gaussian 1.5asec \
  -use-wgridder \
  -join-channels \
  -channels-out 3 \
  -parallel-gridding $ncpu \
  -parallel-reordering $ncpu \
  -no-update-model-required \
  $ms

  rm -rf *.im-000*-*.fits
  rm -rf *.im*dirty.fits
  rm -rf *.im*psf.fits
  cp -rf *.fits $cleanDir
  """
}


process getRegions {
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms
  path imfile
  val rout // radius in arcmin above which sources should be treated as outliers

  output:
  path ms, emit: ms
  path "*-regions-mask.txt", emit: regions
  path "*-regions-outliers.txt", emit: outliers
  path "*-regions-mask.im", emit: mask

  when:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  regions = file(params.root.resolve("inspect/${obsname}/${obsname}-regions"))
  first = file(params.root.resolve("regions/FIRST/${obsname}-regions.txt"))
  regions.exists()

  """
  $params.pythonpath.projdir $projectDir/regions.py --regions $regions --first $first --ms $ms --radius $rout --image $imfile
  cp -rf ${regions}* .

  """
}

process wscleanStokesV {
  publishDir {params.root.resolve("imaging/${obsname}/clean/")}, mode: 'copy'
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms 
  val cell
  val size
  val briggs
  val opt

  output:
  path '*MFS-image.fits', emit: image

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  imname = "${ms.getSimpleName()}-Stokes-V.im"

  """
  wsclean -name "${imname}" \
  -size $size $size \
  -niter 0 \
  -weight briggs $briggs \
  -scale "${cell}arcsec" \
  -taper-gaussian 1.5asec \
  -use-wgridder \
  -join-channels \
  -channels-out 9 \
  -no-update-model-required \
  -pol V \
  $opt \
  $ms

  # remove images to free up space
  rm -rf *.im-000*-*.fits
  rm -rf *.im-*dirty.fits
  rm -rf *.im-*psf.fits
  rm -rf *.im-*model.fits
  rm -rf *.im-*residual.fits
  """
}

process wsclean {
  cpus {"${ncpu}"}
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms 
  val cell
  val size
  val mask
  val automask
  val autothresh
  val briggs
  val ncpu
  val opt

  output:
  path ms, emit: ms
  path '*MFS-image.fits', emit: image
  path '*MFS-image-pb.fits', emit: image_pb, optional: true
  path '*MFS-residual.fits', emit: res
  path '*sources.txt', emit: srcs

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  imname = "${ms.getSimpleName()}.im"
  cleanDir = params.root.resolve("imaging/${obsname}/clean/")

  if(mask == true) {
    mask = file(params.root.resolve("inspect/${obsname}/${obsname}-regions-mask.im"))
    opt = opt + " -casa-mask $mask"
  }

  """
  # remove copies of intermediate measurement sets to free up space
  rm -rf $params.root/imaging/$obsname/*/*.ms

  wsclean -name "${imname}" \
  -size $size $size \
  -niter 100000 \
  -auto-mask $automask \
  -auto-threshold $autothresh \
  -mgain 0.8 \
  -gain 0.1 \
  -weight briggs $briggs \
  -scale "${cell}arcsec" \
  -taper-gaussian 1.5asec \
  -use-wgridder \
  -join-channels \
  -channels-out 9 \
  -deconvolution-channels 9 \
  -fit-spectral-pol 3 \
  -parallel-deconvolution 2048 \
  -deconvolution-threads $ncpu \
  -parallel-gridding $ncpu \
  -parallel-reordering $ncpu \
  -no-update-model-required \
  -save-source-list \
  $opt \
  $ms

  # remove output channel images and dirty images to free up space
  rm -rf *.im-000*-*.fits
  rm -rf *.im-*dirty.fits
  rm -rf *.im-*psf.fits
  cp -rf *.fits $cleanDir
  cp -rf *.txt $cleanDir

  """
}

process predict {
  cpus {"${ncpu}"}
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms
  path srcs
  val ncpu

  output:
  path msmod, emit: ms

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  imname = "${ms.getSimpleName()}.im"
  msmod = "${obsname}-mod.ms"

  """ 
  # make copy of measurement set to store residual
  cp -rL $ms $msmod

  # predict sources from main field
  $params.pythonpath.container -c "import casatasks; casatasks.setjy('$msmod', fluxdensity=[0, 0, 0, 0], usescratch=True)"

  for i in \$(seq 0 8);
  do
      DP3 msin="$msmod" msout="$msmod" msin.band=\$i msin.datacolumn="MODEL_DATA" msout.datacolumn="MODEL_DATA" predict.operation="replace" steps=[predict] predict.sourcedb=$srcs numthreads=$ncpu
  done 
  """
}

process imageOutliers {
  cpus {"${ncpu}"}
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms
  path srcs
  val automask
  val autothresh
  val ncpu

  output:
  path 'SRC*'
  path "${imname}-sources-all.txt", emit: srcs
  path "${imname}-outliers.txt", emit: outliers
  path msmod, emit: ms

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  imname = "${ms.getSimpleName()}.im"
  msres = "${obsname}-res.ms"
  msmod = "${obsname}-mod.ms"
  outliers = file(params.root.resolve("inspect/${obsname}/${obsname}-regions-outliers.txt"))
  cleanDir = params.root.resolve("imaging/${obsname}/clean/")

  """ 
  # make copy of measurement set to store residual
  cp -rL $ms $msres

  # predict sources from main field
  $params.pythonpath.container -c "import casatasks; casatasks.setjy('$msres', fluxdensity=[0, 0, 0, 0], usescratch=True)"

  for i in \$(seq 0 8);
  do
      DP3 msin="$msres" msout="$msres" msin.band=\$i msin.datacolumn="MODEL_DATA" msout.datacolumn="MODEL_DATA" predict.operation="replace" steps=[predict] predict.sourcedb=$srcs numthreads=$ncpu
  done 

  # make copy of measurement set to store residual
  cp -rL $msres $msmod

  # compute residuals and image outliers
  $params.pythonpath.container -c "import casatasks; casatasks.uvsub('$msres')"
  bash $projectDir/srcsub.sh $msres $msmod $outliers $automask $autothresh $ncpu

  # merge source list
  bash $projectDir/mergesourcedb.sh $imname

  # remove residual measurement set
  rm -rf $msres

  cp -rf SRC* $cleanDir
  cp -rf *.txt $cleanDir
  """
}

process subtractOutliers {
  publishDir {params.root.resolve("imaging/${obsname}/clean/")}
  cpus {"${ncpu}"}
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms
  val ncpu

  output:
  path msres, emit: ms

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  msres = "${obsname}-res.ms"
  outliers = file(params.root.resolve("imaging/${obsname}/clean/${obsname}-selfcal-2.im-outliers.txt"))

  """ 
  # make copy of measurement set to store residual
  cp -rL $ms $msres

  # predict sources from main field
  $params.pythonpath.container -c "import casatasks; casatasks.setjy('$msres', fluxdensity=[0, 0, 0, 0], usescratch=True)"

  for i in \$(seq 0 8);
  do
      DP3 msin="$msres" msout="$msres" msin.band=\$i msin.datacolumn="MODEL_DATA" msout.datacolumn="MODEL_DATA" predict.operation="replace" steps=[predict] predict.sourcedb=$outliers numthreads=$ncpu
  done 

  $params.pythonpath.container -c "import casatasks; casatasks.uvsub('$msres')"
  """
}

process getInitSolInt {
  publishDir {params.root.resolve("imaging/${obsname}/selfcal/")}, mode: 'copy'
  errorStrategy 'ignore'

  input:
  path ms
  path imfile
  val snr

  output:
  path "tself.txt"

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  stat = params.root.resolve("imaging/${obsname}/clean/stat.pickle")

  """
  #! $params.pythonpath.projdir

  import pickle
  import numpy as np
  import casatools
  import casatasks

  # get observation time
  tb = casatools.table()
  tb.open("$ms")
  time_array = tb.getcol("TIME")
  tb.close()
  ttot = max(time_array) - min(time_array)

  # get image statistics
  stat = casatasks.imstat("$imfile")

  # per-antenna RMS. 5~sqrt(N_antennas - 3)
  rms_ant = stat['medabsdevmed'] * 5.0

  # self-calibration solution interval
  tself = ttot / (stat['max'] / rms_ant / $snr)**2 

  # double up tself, because it will be used to calibrate on polarisations separately
  tself = max(2, 2 * tself[0])

  print(tself)

  with open('tself.txt', 'w') as f:
    f.write('%d' % tself)
  """
}

process getSolInt {
  errorStrategy 'ignore'

  input:
  path ms
  val calmode

  output:
  path ms, emit: ms
  path "solint.pickle", emit: solint

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  selfcalDir = params.root.resolve("imaging/${obsname}/selfcal/")
  tself = selfcalDir.resolve("tself.txt")

  """
  tmin=\$(< $tself)

  $params.pythonpath.projdir $projectDir/optsolint.py --ms $ms --solint_min \$tmin --calmode $calmode

  # copy solution interval to self-calibration directory
  cp -rf solint.pickle selfcalDir
  """
}


process selfCal {
  cache = 'lenient'
  errorStrategy 'ignore'

  input:
  path ms
  val rnd
  val calmode
  val snr

  output:
  path "${obsname}-selfcal-${rnd}.ms", emit: ms
  path "selfcal-${rnd}*.tb"

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  selfcal = "${obsname}-selfcal-${rnd}.ms"
  selfcalDir = params.root.resolve("imaging/${obsname}/selfcal")
  solint = selfcalDir.resolve("solint.pickle")

  """
  #! $params.pythonpath.container

  import os
  import shutil
  import casatasks
  import casatools
  import pickle
  
  shutil.copytree("$ms", "${selfcal}.tmp")

  # make solution intervals that are a fraction of the total observation time
  with open("$solint", "rb") as file:
    solint_dict = pickle.load(file)
  
  # get observation time
  tb = casatools.table()
  tb.open("$ms")
  time_array = tb.getcol("TIME")
  tb.close()
  ttot = max(time_array) - min(time_array)

  print(solint_dict, ttot)

  # G: gaincal without combining data, T: combine polarisations, S: combine spectral windows
  # The solution interval needs to be adjusted accordingly and is doubled for calmode = "ap"
  N = 1
  if "${calmode}" == "ap":
    N = 2

  solintG = N * max(solint_dict["G"], 36 * solint_dict["ST"])
  solintT = N * max(solint_dict["ST"], 18 * solint_dict["ST"])
  solintS = N * max(solint_dict["S"], 2 * solint_dict["ST"])
  solintST = min(N * solint_dict["ST"], ttot)

  gaintables, spwmaps = [], []

  if (ttot // solintG):
      solintG = (ttot + 2) / (ttot // solintG)

      # SPW self-calibrationm combine polarisations
      casatasks.gaincal \
      ( \
        "${selfcal}.tmp", \
        caltable="selfcal-${rnd}-G.tb", \
        solint=f"{solintG}s", \
        refant="ea19", \
        calmode="${calmode}", \
        solnorm=True, \
        solmode="R", \
        gaintype="G", \
        minsnr=$snr \
      )
      gaintables.append("selfcal-${rnd}-G.tb")
      spwmaps.append([0,1,2,3,4,5,6,7,8])

  if (ttot // solintT):
      solintT = (ttot + 2) / (ttot // solintT)
      
      # SPW self-calibrationm combine polarisations
      casatasks.gaincal \
      ( \
        "${selfcal}.tmp", \
        caltable="selfcal-${rnd}-T.tb", \
        solint=f"{solintT}s", \
        refant="ea19", \
        calmode="${calmode}", \
        solnorm=True, \
        solmode="R", \
        gaintype="T", \
        minsnr=$snr, \
        gaintable=gaintables,
        spwmap=spwmaps \
      )
      gaintables.append("selfcal-${rnd}-T.tb")
      spwmaps.append([0,1,2,3,4,5,6,7,8])

  if (ttot // solintS):
      solintS = (ttot + 2) / (ttot // solintS)

      # SPW self-calibrationm combine polarisations
      casatasks.gaincal \
      ( \
        "${selfcal}.tmp", \
        caltable="selfcal-${rnd}-S.tb", \
        solint=f"{solintS}s", \
        refant="ea19", \
        calmode="${calmode}", \
        solnorm=True, \
        solmode="R", \
        gaintype="G", \
        combine="spw", \
        minsnr=$snr, \
        gaintable=gaintables,
        spwmap=spwmaps \
      )
      gaintables.append("selfcal-${rnd}-S.tb")
      spwmaps.append([0,0,0,0,0,0,0,0,0])

  if (ttot // solintST):
      # SPW self-calibrationm combine polarisations
      casatasks.gaincal \
      ( \
        "${selfcal}.tmp", \
        caltable="selfcal-${rnd}-ST.tb", \
        solint=f"{solintST}s", \
        refant="ea19", \
        calmode="${calmode}", \
        solnorm=True, \
        solmode="R", \
        gaintype="T", \
        combine="spw", \
        minsnr=$snr, \
        gaintable=gaintables,
        spwmap=spwmaps \
      )
      gaintables.append("selfcal-${rnd}-ST.tb")
      spwmaps.append([0,0,0,0,0,0,0,0,0])

  print(f"applying calibration tables {gaintables}")
  casatasks.applycal \
  ( \
    "${selfcal}.tmp", \
    gaintable = gaintables, \
    spwmap = spwmaps \
  )

  # freeze self-calibrated data
  casatasks.split("${selfcal}.tmp", "$selfcal", datacolumn="corrected")
  casatasks.setjy("${selfcal}", fluxdensity=[0, 0, 0, 0], usescratch=True)
  
  # copy model data column
  tb = casatools.table()
  tb.open("${selfcal}.tmp")
  model = tb.getcol("MODEL_DATA")
  tb.close()
  tb.open("${selfcal}", nomodify=False)
  tb.putcol("MODEL_DATA", model)
  tb.close()

  shutil.rmtree("${selfcal}.tmp")

  # copy self-calibration tables to self-calibration directory
  cp -rf *.tb $selfcalDir
  """
}

process flagRFI {
  cache = 'lenient'
  errorStrategy 'ignore'
  
  input:
  path ms
  path aostrat

  output:
  path ms, emit: ms
  path "*.png"
  path "*.npy"

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]
  flaggingDir = params.root.resolve("imaging/${obsname}/flagging/")

  """
  $params.pythonpath.container $projectDir/flagresidual.py --ms $ms  

  # copy results to field flagging directory
  cp -rf *npy $flaggingDir
  cp -rf *png $flaggingDir
  """
}

process removeFiles {
  errorStrategy 'ignore'

  input:
  path ms

  script:
  obsname = (ms =~ /QSO-J[0-9]+[+-]+[0-9]+/)[0]

  """
  echo $obsname >> $projectDir/finished.txt

  sleep 60

  mv $params.root/imaging/$obsname/clean/$obsname-res.im-MFS-image.fits $params.root/final/images/
  mv $params.root/imaging/$obsname/clean/$obsname-res.im-MFS-image-pb.fits $params.root/final/images-pb/
  mv $params.root/imaging/$obsname/clean/$obsname-res.im-MFS-residual.fits $params.root/final/residuals/
  mv $params.root/imaging/$obsname/clean/$obsname-res-Stokes-V.im-MFS-image.fits $params.root/final/stokes-V/

  rm -rf $projectDir/work/*/*/$obsname*.ms
  rm -rf $projectDir/work/*/*/$obsname*.fits
  rm -rf $projectDir/work/*/*/$obsname*tmp*
  rm -rf $projectDir/work/*/*/$obsname*psf*
  rm -rf $projectDir/work/*/*/$obsname*beam*
  rm -rf $projectDir/work/*/*/$obsname*selfcal*
  rm -rf $params.root/imaging/$obsname/*/*.ms
  rm -rf $params.root/inspect/$obsname/*mask.im

  """
}