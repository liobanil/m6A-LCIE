#!/bin/bash --login

#SBATCH -J consensusPeakset
#SBATCH -N 4
#SBATCH -e conPeak.%J.err
#SBATCH -o conPeak.%J.out
#SBATCH --mem=32G
#SBATCH --time=36:00:00

helpFunction()
{
   echo ""
   echo "Usage: sbatch $0 -r" 
   echo -e "\t-r path to peaks files (.bed or any variant)"
   exit 1 # Exit script after printing help
}

echo 'Getting args...'
while getopts "r:" opt
do
   case "$opt" in
      r ) peaksPath="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done
echo 'Done reading args...'

#Print helpFunction in case parameters are empty
if [ -z "$peaksPath" ] 
then
   helpFunction
fi

conda activate basictools-env

#Experiments
exp=( "control" "t6" "t24" )

#bedpos - operation "everything" cats and sorts beds followed by merge
for c in ${exp[@]}; do 
echo "[bedops: sorting and merging beds to get all replicate peak set]"
sort-bed  $peaksPath/*${c}*1*.narrowPeak $peaksPath/*${c}*2*.narrowPeak $peaksPath/*${c}*3*.narrowPeak > $peaksPath/${c}_allreps.bed
bedops --merge $peaksPath/${c}_allreps.bed > $peaksPath/${c}_allreps_mergePeaks.bed
done

echo "[Making global peak set...]"
bedops --merge $peaksPath/control_allreps_mergePeaks.bed $peaksPath/t6_allreps_mergePeaks.bed $peaksPath/t24_allreps_mergePeaks.bed > $peaksPath/globalPeaks.bed
