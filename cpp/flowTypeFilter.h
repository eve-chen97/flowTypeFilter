#ifndef _flowType_FLOWTYPE_
#define _flowType_FLOWTYPE_

#include <Rcpp.h>

RcppExport  SEXP countCells (   //Inputs:
                        SEXP RPartitionsPerChannel, //no of partitions per channel
                        SEXP RThresholds,   //threshold points per channel; list of vectors
                        SEXP RmaxPopSize,
                        SEXP RnumChannels, // needs to be not just numChannels, but labels of channels
                        SEXP RFrameExprData, //expr data from the flowFrame
                        SEXP RNrowFrameData, //Number of rows in flowFrame
                        SEXP RNumPops, //total number of populations
                        SEXP Rverbose //total number of populations
                        );


#endif
