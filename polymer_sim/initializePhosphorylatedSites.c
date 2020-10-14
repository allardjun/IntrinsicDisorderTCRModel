/*** Allard Group jun.allard@uci.edu                    ***/

void initializePhosphorylatedSites();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

int totalPhosphorylated[NFILMAX];
/*******************************************************************************/
void initializePhosphorylatedSites()
{

    //print which of the iSites are phosphorylated
    if (TALKATIVE)
    {
        printf("Phosphorylated sites:\n");
        fflush(stdout);
        for (i=0; i<NumberiSites; i++)
        {
            printf("occupied[ %ld ] =  %f\n",i, occupied[i]);
            fflush(stdout);
        }
    }

    //initializes phosphorylated residues to 0 (none phosphorylated)
    for(nf=0;nf<NFil;nf++)
    {
        for(i=0;i<N[nf];i++)
        {
            PhosphorylatedSites[nf][i] = 0;
        }
    }

    //Phosphorylate sites
    // Note this is ONLY correct if PhosElectroRange does not extend past end points of polymer
    int siteCounter = 0;
    for(nf=0;nf<NFil;nf++)
    {
        for(ty=0;ty<iSiteTotal[nf];ty++)
        {
            if(occupied[siteCounter]==1) //might want to check the truth value on this - equals for double?
            {
                // include exit in case phosphorylation goes too far (can fix this to work appropriately if needed in the future)
                if( (iSite[nf][ty]-PhosElectroRange < 0) || (iSite[nf][ty]+PhosElectroRange+1 > N[nf]) )
                {
                    printf("Error! Phosphorylation extends beyond polymer %ld !",nf);
                    fflush(stdout);
                    exit(0);
                }
                // assign which segments are phosphorylated
                for(i=(iSite[nf][ty]-PhosElectroRange);i<(iSite[nf][ty]+PhosElectroRange+1);i++)
                {
                    PhosphorylatedSites[nf][i]=1; //set that joint to "phosphorylated"
                }
            }
            siteCounter++;
        }
    }

    // for debugging, count and print total segments "phosphorylated"
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            totalPhosphorylated[nf] = 0;
            for (i=0;i<N[nf];i++)
            {
                if (PhosphorylatedSites[nf][i]==1)
                {
                    totalPhosphorylated[nf]++;
                }
            }
            printf("Total Phosphorylated on filament %ld : %d\n", nf,totalPhosphorylated[nf]);
            fflush(stdout);
        }
    }
}
/***********************************************************************************/

