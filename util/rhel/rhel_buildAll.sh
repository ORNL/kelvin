#!/bin/sh

if [[ $# -eq 1 ]]; then
   echo $1
   export kelvin_baseDir=$(echo $1 | tr -d '\r')
   export kelvin_installDir=$(echo $kelvin_baseDir/install | tr -d '\r')
   export kelvin_depBuildDir=$(echo $kelvin_baseDir/deps-build | tr -d '\r')
   mkdir $kelvin_installDir
   mkdir $kelvin_depBuildDir
   mkdir $kelvin_baseDir/build
else
   echo "Please specify the base directory for the Kelvin source code."
   exit
fi

cd $kelvin_baseDir/util/rhel
sh ./rhel_install_deps.sh 
sh ./rhel_build_parsers.sh $kelvin_depBuildDir $kelvin_installDir 
sh ./rhel_build_mfem.sh $kelvin_depBuildDir $kelvin_installDir 
sh ./rhel_build_kelvin.sh $kelvin_baseDir/build $kelvin_installDir 

