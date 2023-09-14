#include <vector>
#include <math.h>
#include "countsAndMeans.h" 
#include "flowTypeFilter.h" 
#include <iostream>
#include <R_ext/Print.h>

using namespace Rcpp ;

/*****************************************************************************
 * Knuth's ancient algorithm for calculating combinations 
 ****************************************************************************/

unsigned long combinations(unsigned int n, unsigned int k)
{
     if (k > n)
         return 0;
     unsigned long r = 1;
     for (unsigned int d = 1; d <= k; ++d)
     {
         r *= n--;
         r /= d;
     }
     return r;
}


/*****************************************************************************
 * Calculate number of total populations given max channels per pop 
 * and total channels 
 ****************************************************************************/


unsigned int calcNumPop (int maxPopSize, int numChannels)
{
  int total = 0;
  for(int i = 1; i <=maxPopSize; i++)
  {
    total = total + combinations(numChannels,i) * (1 << i);
  }
  return total;
}

/*****************************************************************************
 * Convert thresholds from an R list of double vectors to a C++ vector<vector<double> >
 ****************************************************************************/

vector<vector<double> > convertThresholds (List RThresholds)
{
    vector<vector<double> > thresholds = vector<vector<double> >();


    for (int i = 0; i < RThresholds.length(); i++)
    {
        //Make a new vector, then iterate over the converted R vector pushing values in
        thresholds.push_back(vector<double>());
        NumericVector chanThresh = RThresholds[i];

        for (int j = 0; j < chanThresh.length(); j++)
            thresholds[i].push_back(chanThresh[j]);
    }
    

    return thresholds;
}

/*****************************************************************************
 * Main function called from R
 * Translates R data types into CPP equivalents
 * Calls down to deeper, pure CPP function to calculate counts and means
 ****************************************************************************/
    SEXP countCells (   //Inputs:
                        SEXP RPartitionsPerMarker, //no of partitions per channel
                        SEXP RThresholds,   //threshold points per channel; list of vectors
                        SEXP RmaxPopSize,
                        SEXP RPropMarkers, // labels of threshold markers
                        SEXP RMfiData, // expr data from the flowFrame, only including mfi channels
                        SEXP RFrameExprData, //expr data from the flowFrame, only including threshold channels
                        SEXP RNumPops, //total number of populations
			SEXP Rverbose //talk?
                        )
    {
      //dump_pid_status();
        //Perform the easier conversions using built-in functions:
        int maxPopSize = as<int>(RmaxPopSize);
        List RcppThresholds(RThresholds);
        NumericVector propMarkers(RPropMarkers);
        unsigned int numPops = as<int>(RNumPops); 
        unsigned int verbose = as<int>(Rverbose); 

        IntegerVector partitionsPerMarker (RPartitionsPerMarker);
        NumericMatrix mfiData(RMfiData);
        NumericMatrix frameExpr(RFrameExprData);
        int nRowFrame = frameExpr.nrow();
     
        //Convert thresholds using a custom function:
        vector<vector<double> > thresholds = convertThresholds(RcppThresholds);

        if (verbose)
	  Rprintf("\nCalculating phenotype counts and MFIs.");
        //cout << "\nCalculating counts for " << numPops << " phenotypes made up of " << partitionsPerMarker[1] << " partitions of 1-" <<  maxPopSize << " markers each, from a total of " << propMarkers[1] << "markers.\n";
       

        //Initialise R space for returns:
        IntegerVector counts(numPops);
        NumericMatrix mfis(numPops, mfiData.ncol());
        CharacterVector popCodes(numPops); //should be a "character matrix", but can use array indexing later on
     
        for (int i=0; i<partitionsPerMarker.length(); i++)
            partitionsPerMarker[i]+=1; //"base" required by calculateCountsAndMeans is partitions+1

        //Call down:
    
        calculateCountsAndMeans(
                propMarkers.length(), 
                partitionsPerMarker.begin(), 
                nRowFrame, 
                maxPopSize,
                thresholds, 
                frameExpr.begin(), 
                mfiData,
		verbose,
                counts, 
                mfis, 
                popCodes);
        //Return as R List
        return List::create(Named("counts") = counts,
                     Named("mfis") = mfis,
                     Named("codes") = popCodes);

    }




