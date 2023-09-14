#ifndef _COUNTSANDMEANS_
#define _COUNTSANDMEANS_
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

//void dump_pid_status();

void calculateCountsAndMeans(
                        const int marker_count, //total number of markers to threshold
                        const int *base, //number of thresholds for each channel 
                        const int data_count, //number of rows in list mode flow data 
                        const int max_level, //only generate phenotypes of up to this many markers
                        vector<vector<double> > &thresholds, //values of thresholds defining partitions, given per channel 
                        double *data, //"list mode" flow data (matrix of expression values for cells by channel)
                        NumericMatrix mfi_data, //"list mode" flow data, only MFI channels
			const int verbose,

                        //Outputs:
                        IntegerVector &counts, // 
                        NumericMatrix &mfis, // 
                        StringVector &pheno_codes //byte codes representing the immunophenotypes

                         );




#endif

