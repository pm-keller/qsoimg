projectDir=$PWD
rootDir=$PWD/../

while read obsname; do
  rm -rf $projectDir/work/*/*/$obsname*.ms
  rm -rf $projectDir/work/*/*/$obsname*.fits
  rm -rf $projectDir/work/*/*/$obsname*tmp*
  rm -rf $projectDir/work/*/*/$obsname*psf*
  rm -rf $projectDir/work/*/*/$obsname*beam*
  rm -rf $projectDir/work/*/*/$obsname*selfcal*
  rm -rf $rootDir/imaging/$obsname/*/*.ms
  rm -rf $rootDir/inspect/$obsname/*mask.im
done <finished.txt
