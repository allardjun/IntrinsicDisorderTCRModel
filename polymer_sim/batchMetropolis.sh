# Script to run single parameter set for all phosphostates

#####################################################
# Choose how much stiffening per side of binding site
STIFFENRANGE=5 # -1 means don't stiffen

# Name output files
OUTPUTFILENAME="../data/TCRZeta.Stiffenrange.$STIFFENRANGE"

# Find total phosphostates
TOTALITERATIONS=`wc -l < OccupiediSitesZeta.txt`

echo "Length of file is $TOTALITERATIONS"

# run simulation for each phosphostate
for ((IT=1; IT<=$TOTALITERATIONS; IT++))
do
      # read specified line of text files
      OCCUPIEDSITES="`awk 'NR==iter' iter=$IT OccupiediSitesZeta.txt`"

      OCCUPIEDSITESNOSPACE="`awk 'NR==iter' iter=$IT OccupiediSitesZetaNoSpace.txt`"

      # print to screen the line read
      echo "Line $IT of file is $OCCUPIEDISITES"

      # run program with specified parameters
      ./metropolis.out parameters.txt $OUTPUTFILENAME.$IT $OCCUPIEDSITES $OCCUPIEDSITESNOSPACE $STIFFENRANGE &

      # print to screen the process ID and the name of the run
      echo "PID of simulation $IT is $!"
done

echo "Done calling metropolis."

# wait for all background processes to finish before concatenating files
wait

echo "Done waiting for processes to finish."

# loop through all files, concatenate them into one file
for ((IT=1; IT<=$TOTALITERATIONS; IT++))
do

cat $OUTPUTFILENAME.$IT >> $OUTPUTFILENAME.cat

done
