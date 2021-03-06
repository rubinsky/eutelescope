#!/bin/sh

first="17"
last="17"

RUNLIST="runlist.csv"


modus="straight"
#modus="daf"
#modus="gbl"

DRY="--dry-run"

for i in `seq $first $last`; do

 jobsub.py $DRY -c config.cfg -csv $RUNLIST converter  $i
 jobsub.py $DRY -c config.cfg -csv $RUNLIST clustering $i
##### jobsub.py  $DRY -c config.cfg -csv $RUNLIST filter $i

if [[ $modus == "straight" ]]; then
 jobsub.py $DRY -c config.cfg -csv $RUNLIST hitmaker   $i
# alignment using straight line assumption
 jobsub.py $DRY -c config.cfg -csv $RUNLIST alignStraightLinePR      $i
# fitter using broken line implementation by F.Zarnezki
 jobsub.py $DRY -c config.cfg -csv $RUNLIST trackTestFitter $i

elif [[ $modus == "daf" ]]; then
 jobsub.py $DRY -c config.cfg -csv $RUNLIST hitmaker   $i
 jobsub.py $DRY -c config.cfg -csv $RUNLIST alignDafFitter   $i
 jobsub.py $DRY -c config.cfg -csv $RUNLIST trackDafFitter   $i

elif [[ $modus == "gbl" ]]; then
# FOR NEW gbl ONE NEEDS TO GET HITS IN LOCAL COORDINATE SYSTEM:
 jobsub.py  $DRY -c config.cfg -csv $RUNLIST hitlocal $i
# Exhautsive and Helix track search results are not identical - to be investigated (perhaps trivially explained)

#### jobsub.py  $DRY -c config.cfg -csv $RUNLIST tracksearchExh $i
 jobsub.py  $DRY -o MaxMissingHitsPerTrack="2"  -c config.cfg -csv $RUNLIST tracksearchHelix $i
 jobsub.py  $DRY                                -c config.cfg -csv $RUNLIST gbltraj          $i
 jobsub.py  $DRY                                -c config.cfg -csv $RUNLIST gbltrajmille     $i
# echo jobsub.py  $DRY -c config.cfg -csv $RUNLIST aligngbl $i
#
# bash alignm26.sh -r ${i} -l ${RUNLIST} -c 100
# bash align1fei4.sh -r ${i} -l ${RUNLIST} -c 1000

#### jobsub.py  $DRY -c config.cfg -csv $RUNLIST aligngbl $i
#
# 10 iterations to aligne 6 planes (6D)
#   jobsub.py  $DRY -c config.cfg -o GearAlignedFile="gear-${i}-20-21.xml" -csv $RUNLIST trackgbl $i
# jobsub.py  $DRY -c config.cfg -o GearAlignedFile="gear-${i}-10.xml" -o Chi2Cut="5000" -csv $RUNLIST trackgbl $i

fi
done

#
