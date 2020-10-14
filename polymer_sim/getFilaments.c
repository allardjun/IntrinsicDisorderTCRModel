/*** Allard Group jun.allard@uci.edu                    ***/

void getFilaments();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void getFilaments()
{
    /********* INITIALIZE Filaments *******************/
    
    switch (filamentInputMethod)
    {
            
        case 0:  // use identical filaments, number and length set from parameters.txt file

            for(nf=0;nf<NFil;nf++)
            {
                N[nf]=Ntemp;
                if (TALKATIVE) printf("This is number of rods in filament %ld: %ld\n",nf, N[nf]);
            }
            
            break;
            
        case 1: //filaments from file
            
            filList = fopen(filamentFilename, "r");
            char line[200];
            nf=0;
            
            while (fgets(line, sizeof(line), filList))
            {
                N[nf]=atoi(line);
                nf++;
            }
            
            fclose(filList);
            
            // count number of filaments
            NFil=nf;

            break;

    }
    
    //for debugging - prints a list of the filament lengths
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            printf("Filament: %ld\n", nf);
            fflush(stdout);
            
            printf("N: %ld\n", N[nf]);
            fflush(stdout);
        }
    }
    
}

/********************************************************************************************************/

