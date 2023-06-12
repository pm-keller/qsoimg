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
include { mkdir; aoflagger; wscleanDirty; getRegions; wsclean as wsclean1; imageOutliers as imageOutliers1; getInitSolInt; getSolInt; selfCal as selfCal1} from './modules.nf'
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
  mkdir(ms)
  aoflagger(ms, params.aostrat)
  wscleanDirty(ms, params.cell, params.imsize, -0.5, 32)
  getRegions(wscleanDirty.out.ms, wscleanDirty.out.dirty, params.imreg.round(2))
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
  wscleanFinal2(subtractOutliers.out.ms, params.cell, params.imsize, false, 4, 0.3, -0.5, 32, "-apply-primary-beam") 
  wscleanStokesV(wscleanFinal2.out.ms, params.cell, params.imsize, -0.5, "") 
  removeFiles(wscleanFinal2.out.ms)

}

workflow {
  mslist = [
    //params.msdir.resolve("QSO-J0001+0000.ms"),
    //params.msdir.resolve("QSO-J0001+0006.ms"),
    //params.msdir.resolve("QSO-J0004-0049.ms"),
    //params.msdir.resolve("QSO-J0009+3252.ms"),
    //params.msdir.resolve("QSO-J0024+3913.ms"),
    //params.msdir.resolve("QSO-J0028+0457.ms"),
    //params.msdir.resolve("QSO-J0033-0125.ms"),
    //params.msdir.resolve("QSO-J0038-1025.ms"),
    //params.msdir.resolve("QSO-J0045+0901.ms"),
    //params.msdir.resolve("QSO-J0050+3445.ms"),
    //params.msdir.resolve("QSO-J0055+0146.ms"),
    //params.msdir.resolve("QSO-J0100+2802.ms"),
    //params.msdir.resolve("QSO-J0106-0030.ms"),
    //params.msdir.resolve("QSO-J0109-3047.ms"),
    //params.msdir.resolve("QSO-J0112+0110.ms"),
    //params.msdir.resolve("QSO-J0136+0226.ms"),
    //params.msdir.resolve("QSO-J0142-3327.ms"),
    //params.msdir.resolve("QSO-J0159-3633.ms"),
    //params.msdir.resolve("QSO-J0202-0251.ms"),
    //params.msdir.resolve("QSO-J0206-0255.ms"),
    //params.msdir.resolve("QSO-J0210-0456.ms"),
    //params.msdir.resolve("QSO-J0216-0455.ms"),
    //params.msdir.resolve("QSO-J0217-0208.ms"),
    //params.msdir.resolve("QSO-J0221-0802.ms"),
    //params.msdir.resolve("QSO-J0226+0302.ms"),
    //params.msdir.resolve("QSO-J0227-0605.ms"),
    //params.msdir.resolve("QSO-J0231-2850.ms"),
    //params.msdir.resolve("QSO-J0235-0532.ms"),
    // AUTOMASK 4.0, STOCKESV
    //params.msdir.resolve("QSO-J0303-0019.ms"),
    //params.msdir.resolve("QSO-J0305-3150.ms"),
    //params.msdir.resolve("QSO-J0330-4025.ms"),
    //params.msdir.resolve("QSO-J0353+0104.ms"),
    //params.msdir.resolve("QSO-J0402+2451.ms"),
    //params.msdir.resolve("QSO-J0421-2657.ms"),
    //params.msdir.resolve("QSO-J0422-1927.ms"),
    //params.msdir.resolve("QSO-J0559-1535.ms"),
    //params.msdir.resolve("QSO-J0818+1722.ms"),
    //params.msdir.resolve("QSO-J0828+2633.ms"),
    //params.msdir.resolve("QSO-J0842+1218.ms"),
    //params.msdir.resolve("QSO-J0844-0052.ms"),
    //params.msdir.resolve("QSO-J0844+0226.ms"),
    //params.msdir.resolve("QSO-J0853+0139.ms"),
    //params.msdir.resolve("QSO-J0857+0056.ms"),
    //params.msdir.resolve("QSO-J0859+0022.ms"),
    //params.msdir.resolve("QSO-J0902+0155.ms"),
    //params.msdir.resolve("QSO-J0905+0300.ms"),
    //params.msdir.resolve("QSO-J0911+0152.ms"),
    //params.msdir.resolve("QSO-J0935-0110.ms"),
    //params.msdir.resolve("QSO-J1030+0524.ms"),
    //params.msdir.resolve("QSO-J1034-1425.ms"),
    //params.msdir.resolve("QSO-J1036-0232.ms"),
    //params.msdir.resolve("QSO-J1048-0109.ms"),
    //params.msdir.resolve("QSO-J1048+4637.ms"),
    //params.msdir.resolve("QSO-J1110-1329.ms"),
    //params.msdir.resolve("QSO-J1120+0641.ms"),
    //params.msdir.resolve("QSO-J1137+0045.ms"),
    //params.msdir.resolve("QSO-J1137+3549.ms"),
    //params.msdir.resolve("QSO-J1148+0702.ms"),
    //params.msdir.resolve("QSO-J1148+5251.ms"),
    //params.msdir.resolve("QSO-J1152+0055.ms"),
    //params.msdir.resolve("QSO-J1201+0133.ms"),
    //params.msdir.resolve("QSO-J1205-0000.ms"),
    //params.msdir.resolve("QSO-J1207-0005.ms"),
    //params.msdir.resolve("QSO-J1207+0630.ms"),
    //params.msdir.resolve("QSO-J1208-0200.ms"),
    //params.msdir.resolve("QSO-J1212+0505.ms"),
    //params.msdir.resolve("QSO-J1217+0131.ms"),
    //params.msdir.resolve("QSO-J1231+0052.ms"),
    //params.msdir.resolve("QSO-J1250+3130.ms"),
    //params.msdir.resolve("QSO-J1254-0014.ms"),
    //params.msdir.resolve("QSO-J1257+6349.ms"),
    //params.msdir.resolve("QSO-J1306+0356.ms"),
    //params.msdir.resolve("QSO-J1319+0950.ms"),
    //params.msdir.resolve("QSO-J1342+0928.ms"),
    //params.msdir.resolve("QSO-J1344+0128.ms"),
    //params.msdir.resolve("QSO-J1347-0157.ms"),
    //params.msdir.resolve("QSO-J1350-0027.ms"),
    //params.msdir.resolve("QSO-J1401+2749.ms"),
    //params.msdir.resolve("QSO-J1402+4024.ms"),
    //params.msdir.resolve("QSO-J1416+0015.ms"),
    //params.msdir.resolve("QSO-J1417+0117.ms"),
    //params.msdir.resolve("QSO-J1423-0018.ms"),
    //params.msdir.resolve("QSO-J1425-0015.ms"),
    //params.msdir.resolve("QSO-J1426-0128.ms"),
    //params.msdir.resolve("QSO-J1427+3312.ms"),
    //params.msdir.resolve("QSO-J1428-1602.ms"),
    //params.msdir.resolve("QSO-J1429-0002.ms"),
    //params.msdir.resolve("QSO-J1429-0104.ms"),
    //params.msdir.resolve("QSO-J1429+5447.ms"),
    //params.msdir.resolve("QSO-J1431-0724.ms"),
    //params.msdir.resolve("QSO-J1440+0107.ms"),
    //params.msdir.resolve("QSO-J1448+4333.ms"),
    //params.msdir.resolve("QSO-J1509-1749.ms"),
    //params.msdir.resolve("QSO-J1512+4422.ms"),
    //params.msdir.resolve("QSO-J1516+4228.ms"),
    //params.msdir.resolve("QSO-J1525+4303.ms"),
    //params.msdir.resolve("QSO-J1526-2050.ms"),
    //params.msdir.resolve("QSO-J1558-0724.ms"),
    //params.msdir.resolve("QSO-J1559+2212.ms"),
    //params.msdir.resolve("QSO-J1602+4228.ms"),
    //params.msdir.resolve("QSO-J1603+5510.ms"),
    //params.msdir.resolve("QSO-J1609+3041.ms"),
    //params.msdir.resolve("QSO-J1609+5328.ms"),
    //params.msdir.resolve("QSO-J1612+5559.ms"),
    //params.msdir.resolve("QSO-J1623+3112.ms"),
    //params.msdir.resolve("QSO-J1629+2407.ms"),
    //params.msdir.resolve("QSO-J1630+4012.ms"),
    //params.msdir.resolve("QSO-J1641+3755.ms"),
    //params.msdir.resolve("QSO-J1724+1901.ms"),
    //params.msdir.resolve("QSO-J1932+7139.ms"),
    //params.msdir.resolve("QSO-J2032-2114.ms"),
    //params.msdir.resolve("QSO-J2054-0005.ms"),
    //params.msdir.resolve("QSO-J2100-1715.ms"),
    //params.msdir.resolve("QSO-J2132+1217.ms"),
    //params.msdir.resolve("QSO-J2201+0155.ms"),
    //params.msdir.resolve("QSO-J2211-3206.ms"),
    //params.msdir.resolve("QSO-J2215+2606.ms"),
    //params.msdir.resolve("QSO-J2216-0016.ms"),
    params.msdir.resolve("QSO-J2219+0102.ms"),
    params.msdir.resolve("QSO-J2228+0128.ms"),
    params.msdir.resolve("QSO-J2228+0152.ms"),
    params.msdir.resolve("QSO-J2229+1457.ms"),
    params.msdir.resolve("QSO-J2232+0012.ms"),
    params.msdir.resolve("QSO-J2232+2930.ms"),
    params.msdir.resolve("QSO-J2236+0032.ms"),
    params.msdir.resolve("QSO-J2239+0207.ms"),
    params.msdir.resolve("QSO-J2240-1839.ms"),
    params.msdir.resolve("QSO-J2255+0503.ms"),
    params.msdir.resolve("QSO-J2310+1855.ms"),
    params.msdir.resolve("QSO-J2315-0023.ms"),
    params.msdir.resolve("QSO-J2318-0246.ms"),
    params.msdir.resolve("QSO-J2318-3029.ms"),
    //params.msdir.resolve("QSO-J2318-3113.ms"),
    //params.msdir.resolve("QSO-J2329-0301.ms"),
    params.msdir.resolve("QSO-J2348-3054.ms"),
    params.msdir.resolve("QSO-J2356+0017.ms"),
    //params.msdir.resolve("QSO-J2356+0023.ms"),
    //params.msdir.resolve("QSO-J2356-0622.ms"),
]

  msChannel = Channel.fromPath(mslist, type: 'dir')

  imagingPipeline(msChannel)
}
