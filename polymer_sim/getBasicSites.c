/*** Allard Group jun.allard@uci.edu                    ***/

void getBasicSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void getBasicSites()
{
    /********* INITIALIZE BASIC SITES *******************/
    

    basicSiteList = fopen(basicSiteFilename, "r");
    char line[200];
    nf=0;
    
    while (fgets(line, sizeof(line), basicSiteList))
    {
        iy=0;
        
        // if line has something on it, set iSites equal to parts of line
        if(line[0] != '\n')
        {
            char * linepart;
            linepart = strtok(line," ,");
            while(linepart != NULL)
            {
                
                if(atoi(linepart)!=-1)
                {
                    basicSite[nf][iy] = atoi(linepart);
                    linepart = strtok(NULL, " ,");
                    iy++;
                }
                else
                {
                    linepart = strtok(NULL, " ,");
                }
            }
        }
        
        // set basicSiteTotal and update filament
        basicSiteTotal[nf]=iy;
        if(basicSiteTotal[nf]==0)
        {
            printf("Filament %ld has no basic sites.\n",nf);
            fflush(stdout);
        }
        nf++;
        
    }
    
    if(nf!=NFil)
    {
        printf("Error! Number of filaments mismatch between filaments and  basicSites!\n");
        fflush(stdout);
        exit(0);
    }
    
    fclose(basicSiteList);

    
    // warning for basic sites located past filament end
    for(nf=0;nf<NFil;nf++)
    {
        for (iBasic=0; iBasic<basicSiteTotal[nf];iBasic++)
        {
            if (basicSite[nf][iBasic] >= N[nf])
            {
                printf("Warning! Basic site is located past end of polymer number %ld !",nf);
                fflush(stdout);
            }
        }
    }
    
    //initializes basic residues to 0
    for(nf=0;nf<NFil;nf++)
    {
        for(i=0;i<N[nf];i++)
        {
            BasicSitesYN[nf][i] = 0;
        }
    }
    
    //Full residue list of basic yes or no
    for(nf=0;nf<NFil;nf++)
    {
        for (i=0;i<N[nf];i++)
        {
            for(iBasic=0;iBasic<basicSiteTotal[nf];iBasic++)
            {
                if(i==basicSite[nf][iBasic])
                {

                    BasicSitesYN[nf][i]=1; //set that joint to "basic"

                }
            }
        }
    }
    
    //for debugging - prints a list of the iSites
    
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            printf("Filament: %ld\n", nf);
            fflush(stdout);
            
            for (iBasic=0;iBasic<basicSiteTotal[nf];iBasic++)
            {
                printf("basicSite: %ld\n", basicSite[nf][iBasic]);
                fflush(stdout);
            }
            
            printf("basicSiteTotal: %ld\n", basicSiteTotal[nf]);
            fflush(stdout);
        }
    }


}

/********************************************************************************************************/

