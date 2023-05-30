/** @file GseaRcpp.cc
 * @brief GseaRcpp implemenation file */

#include <Rcpp.h>
#include "gsea.hh"
using namespace Rcpp;

/**
 * @class GseaRcpp
 * @brief Translates Rcpp data structures to C++ data structures
 */
class GseaRcpp
{
private:
    Gsea *gsea;

public:
    /**
    * @brief Gsea creator function used when GSEA is using Gsea$runChunked()
    * @param sampleIds sample ids as a character vector of the expression matrix
    * @param geneIds gene ids as a character vector of the expression matrix
    * @param geneSets gene ids as a character vector of the expression matrix
    * @param nThreads number of threads used to run GSEA, 0 if all CPU threads want to be used
    * @post Gene sets, sample ids and gene ids are initialised
    */
    GseaRcpp(CharacterVector sampleIdsRcpp,
             CharacterVector geneIdsRcp,
             List geneSetsRcpp,
             uint nThreads);

    /**
    * @brief Gsea creator function when GSEA is run using Gsea$run()
    * @param geneSets gene sets as list
    * @param expressionMatrix expression matrix as numeric matrix, it has no null rownames and colnames
    * @param nThreads number of threads used to run GSEA, 0 if all CPU threads want to be used
    * @post Gene sets, sample ids, gene ids and expression matrix are initialised
    */
    GseaRcpp(NumericMatrix expressionMatrix,
             List geneSets,
             uint nThreads);

    /**
    * @brief Runs GSEA for the given expression matrix using the gene sets initialised in the Gsea creator function
    * @param expressionMatrix numeric matrix containing counts in the cells, samples in the rows and genes in the columns
    * @post The correspinding chunk file contains the ES score for each sample and gene set
    */
    void runChunked(const NumericMatrix &countMatrixRcpp);

    /**
    * @brief Filter the chunked GSEA results by selecting the nFilteredGeneSets gene sets with more variance across the samples
    * @param nFilteredGeneSets number of gene sets selected to be written in the filtered results file
    * @pre runChunked has been run at least one time
    * @post filtered-results.csv contains the filtered results
    */

    void filterResults(uint nFilteredGeneSets, string chunksPath, string outFilename);

    /**
    * @brief Runs GSEA for the given expression matrix and gene sets in the Gsea creator function
    * @param outFileName name of the output file
    * @param ioutput number of gene sets between status output
    * @post ES results are written into outFileName file
    */
    void run(string outFileName, uint ioutput);

    /**
    * @brief Normalize the expression matrix using rpm and centering the samples
    * @post The expression matrix is normalized
    */
    void normalizeExprMatrix();

    ~GseaRcpp();
};
