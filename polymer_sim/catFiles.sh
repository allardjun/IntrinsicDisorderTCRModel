##############################################
# Concatenate output files
# Tar individual files
# 08 02 2019

##############################################
# Set up variables
#ITERArray=( 4 8 12 16 20 )
ITERArray=( 20 )
#ITER2Array=( 1 2 3 5 9 10 )
ITER2Array=( 10 )

I2MAX=1024

FILENAME=TCRMembrane1FilSweep
FOLDERNAME=NITAM
FOLDERNAME2=NFIL

##############################################
# Make folders for concatenated files

if [ ! -d CatFiles ]
then
    mkdir CatFiles

    mkdir CatFiles/MissingFiles
fi

##############################################
# loop through all files, concatenate them into one file

for ITER in "${ITERArray[@]}"
do
    cd $FOLDERNAME$ITER

    for ITER2 in "${ITER2Array[@]}"
    do

      cd $FOLDERNAME2$ITER2

      echo "We are here: $PWD"

      # if cat file does not exist
      if [[ ! -e "$FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.cat" ]];
      then
          for ((IT=1;IT<=$I2MAX;IT++))
          do
              # if individual file exists, add data to cat file. If not, report as missing.
              if [ -e "$FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.$IT" ];
              then
                 # if there is more than one line in the file, just take first line
                  if [[ "$(wc -l < $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.$IT)" -gt "1" ]];
                  then
                      echo -e "Too many lines in output file $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.$IT.\n" >> ../CatFiles/CatFiles.err

                      awk 'NR==1' $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.$IT >> $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.cat
                  else
                      cat $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.$IT >> $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.cat
                  fi

              else
                  echo -e "$FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.$IT\n" >> MissingFiles.$ITER.$FOLDERNAME2.$ITER2
              fi


          done

          # if more lines in cat file than expected
          if [[ "$(wc -l < $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.cat)" -gt "$I2MAX" ]];
          then
              echo -e "Too many lines in cat file $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.cat.\n" >> ../../CatFiles/CatFiles.err
          fi

          # if less lines in cat file than expected
          if [[ "$(wc -l < $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.cat)" -lt "$I2MAX" ]];
          then
              echo -e "Missing lines in cat file $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.cat.\n" >> ../../CatFiles/CatFiles.err
          fi

          # copy cat file to folder
          cp $FILENAME$FOLDERNAME.$ITER.$FOLDERNAME2.$ITER2.cat ../../CatFiles/

          # copy missing files file to folder
          if [[ -e MissingFiles.$ITER.$FOLDERNAME2.$ITER2 ]];
          then
              cp MissingFiles.$ITER.$FOLDERNAME2.$ITER2 ../../CatFiles/MissingFiles/
          fi

          echo "Finished concatenating output files $FOLDERNAME$ITER.$FOLDERNAME2$ITER2."

      fi

      cd ..
    done

    cd ..

done

echo "Done concatenating output files."

echo "We are here: $PWD"

##############################################
# copy concatenated files to independent CatFiles folder

cp -R CatFiles/ ../CatFiles/

##############################################
