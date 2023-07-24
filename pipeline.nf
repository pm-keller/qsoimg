params.root = projectDir.resolve('../')
params.msdir = params.root.resolve('ms')
params.aostrat = projectDir.resolve('aoflagger_strategy.lua')
params.srcs = projectDir.resolve('')
params.pythonpath = [container: "/opt/conda/envs/minicasa/bin/python3", projdir: "$projectDir/casa/bin/python3"]
params.imsize = 12288 // pixels
params.cell = 0.3 // cell size
params.imreg = params.imsize * params.cell / (2 * 60) 
// radius in arcmin of region to be imaged

// PROCESSES FOR FIRST ROUND OF IMAGING  
include { mkdir; aoflagger; wscleanDirty; wscleanInspect; getRegions; wsclean as wsclean1; imageOutliers as imageOutliers1; getInitSolInt; getSolInt; selfCal as selfCal1} from './modules.nf'
// PROCESSES FOR SECOND ROUND OF IMAGING  
include { wsclean as wsclean2; imageOutliers as imageOutliers2; selfCal as selfCal2} from './modules.nf'
// PROCESSES FOR THIRD ROUND OF IMAGING  
include { wsclean as wsclean3; imageOutliers as imageOutliers3; selfCal as selfCal3; flagRFI} from './modules.nf'
// PROCESSES FOR FOURTH ROUND OF IMAGING  
include { wsclean as wsclean4; imageOutliers as imageOutliers4; selfCal as selfCal4} from './modules.nf'
// PROCESSES FOR FINAL IMAGES  
include { wsclean as wscleanFinal1; wsclean as wscleanFinal2; subtractOutliers; wscleanStokesV; removeFiles} from './modules.nf'

workflow imagingPipeline {
  take:
  ms 
  
  main:

  // FIRST ROUND
  //mkdir(ms)
  //wscleanInspect(ms, 1, 8192, -0.5, 4)
  aoflagger(ms, params.aostrat)
  wscleanDirty(ms, params.cell, params.imsize, -0.5, 32)
  getRegions(wscleanDirty.out.ms, params.imreg.round(2))
  wsclean1(getRegions.out.ms, params.cell, params.imsize, true, 10, 1, -0.5, 32, "-no-negative -local-rms")
  imageOutliers1(wsclean1.out.ms, wsclean1.out.srcs, 10, 1, 32)
  getInitSolInt(wsclean1.out.ms, wsclean1.out.image, 10.0)
  getSolInt(imageOutliers1.out.ms, "p")
  selfCal1(getSolInt.out.ms, 0, "p", 5)
  
  // SECOND ROUND
  wsclean2(selfCal1.out.ms, params.cell, params.imsize, true, 5, 0.5, -0.5, 32, "-no-negative -local-rms")
  imageOutliers2(wsclean2.out.ms, wsclean2.out.srcs, 5, 0.5, 32)
  selfCal2(imageOutliers2.out.ms, 1, "p", 5)

  // THIRD ROUND
  wsclean3(selfCal2.out.ms, params.cell, params.imsize, true, 5, 0.5, -0.5, 32, "-no-negative -local-rms")
  imageOutliers3(wsclean3.out.ms, wsclean3.out.srcs, 5, 0.5, 32)
  selfCal3(imageOutliers3.out.ms, 2, "ap", 5.5)
  flagRFI(selfCal3.out.ms, params.aostrat)

  // FOURTH ROUND
  wsclean4(flagRFI.out.ms, params.cell, params.imsize, false, 5, 0.5, -0.5, 32, "-local-rms")
  imageOutliers4(wsclean4.out.ms, wsclean4.out.srcs, 5, 0.5, 32)
  selfCal4(imageOutliers4.out.ms, 3, "ap", 5)

  // FINAL IMAGING
  subtractOutliers(selfCal4.out.ms, 32)
  wscleanFinal2(subtractOutliers.out.ms, params.cell, params.imsize, "final", 4, 0.3, -0.5, 32, "-apply-primary-beam -multiscale") 
  wscleanStokesV(wscleanFinal2.out.ms, params.cell, params.imsize, -0.5, 1, "") 
  removeFiles(wscleanStokesV.out.ms)

}

workflow imagingPipelineNoSelfCal {
  take:
  ms 
  
  main:

  // FIRST ROUND
  //mkdir(ms)
  //wscleanInspect(ms, 1, 8192, -0.5, 4)
  wscleanDirty(ms, params.cell, params.imsize, -0.5, 32)
  getRegions(wscleanDirty.out.ms, params.imreg.round(2))
  wsclean1(getRegions.out.ms, params.cell, params.imsize, true, 10, 1, -0.5, 32, "-no-negative -local-rms")
  imageOutliers1(wsclean1.out.ms, wsclean1.out.srcs, 10, 1, 32)
  subtractOutliers(imageOutliers1.out.ms, 32)
  flagRFI(subtractOutliers.out.ms, params.aostrat)
  

  // FINAL IMAGING
  wscleanFinal2(flagRFI.out.ms, params.cell, params.imsize, "final", 4, 0.3, -0.5, 32, "-apply-primary-beam -multiscale") 
  wscleanStokesV(wscleanFinal2.out.ms, params.cell, params.imsize, -0.5, 1, "") 
  removeFiles(wscleanStokesV.out.ms)

}

workflow {
  mslist = [
    params.msdir.resolve("QSO-J0818+1722.ms"), 
    params.msdir.resolve("QSO-J0828+2633.ms"),
    params.msdir.resolve("QSO-J0844-0052.ms")
  ]

  msChannel = Channel.fromPath(mslist, type: 'dir')

  imagingPipeline(msChannel)

}
