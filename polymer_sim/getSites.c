/*** Allard Group jun.allard@uci.edu                    ***/

void getSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/
long siteCounter,bSiteCounter;

/********************************************************************************************************/
void getSites()
{
    /********* INITIALIZE ISITES *******************/


    switch (iSiteInputMethod) //switch in Batch Script to specify which set of iSites you want to use
    {

        case 0:  // iSites initialized for human CD3Zeta-Chain

            for(nf=0;nf<NFil;nf++)
            {
                iSiteTotal[nf] = 6;

                for (iy=0;iy<iSiteTotal[nf];iy++) //initializes iSite array
                {
                    iSite[nf][iy]=0;
                }
            }

            for(nf=0;nf<NFil;nf++)
            {
                //specify iSite locations by rod number (located at n-1) (i.e. if you have a polymer of N=50, and want an iSite at the 36th rod, input iSite[]=35)

                iSite[nf][0]=20;
                iSite[nf][1]=31;
                iSite[nf][2]=59;
                iSite[nf][3]=71;
                iSite[nf][4]=90;
                iSite[nf][5]=101;
            }
            break;

        case 1: //do nothing, use command line input

            break;

        case 2: //input iSites from file
            // use -1 on a line to denote no iSites for that filament
            // a blank line will also denote no iSites, but doesn't work if final filament is one without iSites

            iSiteList = fopen(iSiteFilename, "r");
            char line[200];
            nf=0;

            while (fgets(line, sizeof(line), iSiteList))
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
                            iSite[nf][iy] = atoi(linepart);
                            linepart = strtok(NULL, " ,");
                            iy++;
                        }
                        else
                        {
                            linepart = strtok(NULL, " ,");
                        }
                    }
                }

                // set iSiteTotal and update filament
                iSiteTotal[nf]=iy;
                if(iSiteTotal[nf]==0)
                {
                    printf("Filament %ld has no iSites.\n",nf);
                    fflush(stdout);
                }
                nf++;

            }

            if(nf!=NFil)
            {
                printf("Error! Number of filaments mismatch between filaments and iSites!\n");
                fflush(stdout);
                exit(0);
            }

            fclose(iSiteList);

            break;

        case 3: // use last site as only iSite

            for(nf=0;nf<NFil;nf++)
            {
                iSiteTotal[nf] = 1;
                for(iy=0;iy<iSiteTotal[nf];iy++)
                {
                    iSite[nf][iy]=0;
                }
            }

            for(nf=0;nf<NFil;nf++)
            {
                iSite[nf][0] = N[nf]-1;
            }

            break;

    }

    // Determine total number of iSites for system
    NumberiSites = 0;
    for(nf=0;nf<NFil;nf++)
        NumberiSites += iSiteTotal[nf];

    //Warning for possible user error
    for(nf=0;nf<NFil;nf++)
    {
        for (iy=0; iy<iSiteTotal[nf];iy++)
        {
            if (iSite[nf][iy] >= N[nf])
            {
                printf("Warning! Site is located past end of polymer in filament %ld!\n",nf);
                fflush(stdout);
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

            for (iy=0;iy<iSiteTotal[nf];iy++)
            {
                printf("iSite: %ld\n", iSite[nf][iy]);
                fflush(stdout);
            }

            printf("iSiteTotal: %ld\n", iSiteTotal[nf]);
            fflush(stdout);
        }

        printf("Number of iSites in system: %ld\n",NumberiSites);
        fflush(stdout);
    }

    /*****************************************************/
    /********* INITIALIZE BOUND ISITES *******************/
    /*****************************************************/

    if (MULTIPLE) //if looking at multiple binding (i.e. MULTIPLE set to 1 in driveM)
    {
        switch (bSiteInputMethod)
        {
            case 0:
                for(nf=0;nf<NFil;nf++)
                {
                    bSiteTotal[nf] = 4; //total number of iSites bound

                    for (ib=0;ib<bSiteTotal[nf];ib++) //initializes iSite array
                    {
                        bSite[nf][ib]=0;
                    }
                }
                for(nf=0;nf<NFil;nf++)
                {
                    bSite[nf][0]=10;
                    bSite[nf][1]=10;
                    bSite[nf][2]=40;
                    bSite[nf][3]=40;
                }
                break;

            case 1: //do nothing, use command line input

                break;

            case 2: //bSites for multiple binding of ZAP-70

                siteCounter = 0;
                for(nf=0;nf<NFil;nf++)
                {
                    bSiteCounter=0;
                    for (iy=0;iy<iSiteTotal[nf];iy++)
                    {
                        if (occupied[siteCounter]==1)
                        {
                            bSite[nf][bSiteCounter]=iSite[nf][iy];
                            bSiteCounter++;
                        }
                        siteCounter++;
                    }
                    bSiteTotal[nf]=bSiteCounter;
                }

                break;

            case 3: // read bound sites from file
                // use -1 on a line to denote no bSites for that filament
                // a blank line will also denote no bSites, but doesn't work if final filament is one without bSites

                bSiteList = fopen(bSiteFilename, "r");
                char line[200];
                nf=0;

                while (fgets(line, sizeof(line), bSiteList))
                {
                    ib=0;

                    // if line has something on it, set bSites equal to parts of line
                    if(line[0] != '\n')
                    {
                        char * linepart;
                        linepart = strtok(line," ,");
                        while(linepart != NULL)
                        {

                            if(atoi(linepart)!=-1)
                            {
                                bSite[nf][ib] = atoi(linepart);
                                linepart = strtok(NULL, " ,");
                                ib++;
                            }
                            else
                            {
                                linepart = strtok(NULL, " ,");
                            }
                        }
                    }

                    // set iSiteTotal and update filament
                    bSiteTotal[nf]=ib;
                    if(bSiteTotal[nf]==0)
                    {
                        printf("Filament %ld has no bound sites.\n",nf);
                        fflush(stdout);
                    }
                    nf++;

                }

                if(nf!=NFil)
                {
                    printf("Error! Number of filaments mismatch between filaments and bSites!\n");
                    fflush(stdout);
                    exit(0);
                }

                fclose(bSiteList);

                break;

            case 4: // use last site as only bSite

                for(nf=0;nf<NFil;nf++)
                {
                    bSiteTotal[nf] = 1;
                    for(ib=0;ib<bSiteTotal[nf];ib++)
                    {
                        bSite[nf][ib]=0;
                    }
                }
                for(nf=0;nf<NFil;nf++)
                {
                    bSite[nf][0] = N[nf]-1;
                }

                break;

        }

        // Determine total number of iSites for system
        NumberbSites = 0;
        for(nf=0;nf<NFil;nf++)
            NumberbSites += bSiteTotal[nf];

        //Warning for possible user error
        for(nf=0;nf<NFil;nf++)
        {
            for (iy=0;iy<bSiteTotal[nf];iy++)
            {
                if (bSite[nf][iy] >= N[nf])
                {
                    printf("Warning! Bound site is located past end of filament %ld!\n",nf);
                    fflush(stdout);
                }
            }
        }

        //for debugging - prints a list of the bSites
        if (TALKATIVE)
        {
            for(nf=0;nf<NFil;nf++)
            {
                printf("Filament: %ld \n", nf);
                fflush(stdout);

                for (ib=0; ib<bSiteTotal[nf]; ib++)
                {
                    printf("bSite[%ld][%ld] is %ld \n",nf, ib, bSite[nf][ib]);
                    fflush(stdout);
                }

                printf("bSiteTotal = %ld \n", bSiteTotal[nf]);
                fflush(stdout);
            }

            printf("Number of bound sites in system = %ld \n", NumberbSites);
            fflush(stdout);
        }

    }





}

/********************************************************************************************************/
