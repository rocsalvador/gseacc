/** @file gsea.hh
 * @brief Gsea header file */

#ifndef GSEA_HH
#define GSEA_HH

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <thread>
#include <chrono>
#include <filesystem>
#include <cassert>

using namespace std;
using namespace chrono;

/** @struct GseaSample
 * @brief Gene count with a pointer to the position of the gene id in Gsea::geneIds */
struct GeneSample
{
    uint32_t geneId;
    float count;
};

/** @struct GseaSet
 * @brief Gene set struct containing the gene set id and its genes in a set */
struct GeneSet
{
    string geneSetId;
    unordered_set<string> geneSet;
};

/** @struct GseaSetPtr
 * @brief Gene set variance value with a pointer to the position of the gene set in Gsea::geneSets */
struct GeneSetPtr
{
    uint geneSetPtr;
    float value;
};

/** @class Gsea
 * @brief Runs Gsea indepentdently of Rcpp  */
class Gsea
{
private:
    // gsea.config variables
    string expressionMatrixFilename;
    string geneSetsFilename;
    string outputFilename;
    char expressionMatrixSep;
    char geneSetsSep;
    char outputSep;
    uint ioutput;
    uint batchSize;
    uint chunk;
    uint nThreads;
    bool normalizedData;
    bool scRna;

    /// Thread in charge of printing the status
    uint logThread;

    /// Variable to keep track of the current sample while running runChunked()
    uint currentSample;
    /// Path to the folder where chunks are saved
    filesystem::path chunksPath;

    /// Array containing the gene sets
    vector<GeneSet> geneSets;

    /// Matrix containing the gene counts
    vector<vector<GeneSample>> expressionMatrix;
    /// Array containing sample ids
    vector<string> sampleIds;
    /// Array containing gene ids
    vector<string> geneIds;

    /// Matrix containing GSEA results
    vector<vector<float>> results;

    /// Number of genes in the expression matrix
    uint nGenes;
    /// Number of samples in the expression matrix
    uint nSamples;
    /// Number of gene sets in the gene sets
    uint nGeneSets;

    /// Time point when GSEA was started
    system_clock::time_point startGSEATime;

    /**
    * @brief GeneSample comparator function to sort gene samples in decreasing order
    * @param g1 first GeneSample
    * @param g2 second GeneSample
    * @return True if g1.count is bigger than g2.count, false otherwise
    */
    static bool geneSampleComp(const GeneSample &g1, const GeneSample &g2);

    /**
    * @brief GeneSetPtr comparator function to sort gene samples in decreasing order
    * @param g1 first GeneSetPtr
    * @param g2 second GeneSetPtr
    * @return True if g1.value is bigger than g2.value, false otherwise
    */
    static bool geneSetPtrComp(const GeneSetPtr &g1, const GeneSetPtr &g2);

    void readRna();

    void readScRna();

    void runScRna();

    void runRna();

    /**
    * @brief Reads the gsea.config file
    * @post All configuration variables are set up
    */
    void readConfig();

    /**
    * @brief Rpm the expression matrix
    * @post Each expression matrix row sums 1 million
    */
    void rpm();

    /**
    * @brief Mean center the expression matrix
    * @post The expression matrix is mean centered rowwise
    */
    void meanCenter();

    void sortColumnsJob(uint columnStart, uint columnEnd);

    /**
    * @brief Runs the gsea for all the expression matrix, dividing it in nThreads
    * @pre expressionMatrix rows contain genes, expressionMatrix columns contain samples
    * @post The samples in the results matrix contain the ES
    */
    void enrichmentScore();

    /**
    * @brief Runs the gsea from lineStart to lineEnd samples, assuming genes in the rows and samples in the columns
    * @param startSample start sample
    * @param endSample end sample
    * @pre expressionMatrix rows contain genes, expressionMatrix columns contain samples
    * @post The samples startSample to endSample in the results matrix contain the ES
    */
    void enrichmentScoreJob(uint sampleStart, uint sampleEnd);

    /**
    * @brief Runs the gsea from startSample to endSample samples, assuming samples in the rows and genes in the columns
    * @param startSample start sample
    * @param endSample end sample
    * @pre expressionMatrix rows contain samples, expressionMatrix columns contain genes
    * @post The samples startSample to endSample in the results matrix contain the ES
    */
    void scEnrichmentScoreJob(uint sampleStart, uint sampleEnd);

    /**
    * @brief Writes the results into outputFilename
    * @post Results are written into outputFilenName
    */
    void writeResults();

public:
    /**
    * @brief Gsea creator function to use the class without R, it reads the configuration from
    * gsea.config, and all input data from files
    * @post Gene sets and the expression matix are initialised
    */
    Gsea();

    /**
    * @brief Gsea creator function used when GSEA is using Gsea::runChunked()
    * @param sampleIds sample ids of the expression matrix
    * @param geneIds gene ids of the expression matrix
    * @param geneSets array of gene sets
    * @param nThreads number of threads used to run GSEA, 0 if all CPU threads want to be used
    * @post Gene sets, sample ids and gene ids are initialised
    */
    Gsea(vector<string> &sampleIds,
         vector<string> &geneIds,
         vector<GeneSet> &geneSets,
         uint nThreads);
    /**
    * @brief Gsea creator function when GSEA is run using Gsea::run()
    * @param geneSets array of gene sets
    * @param expressionMatrix expression matrix
    * @param geneIds gene ids of the expression matrix
    * @param sampleIds sample ids of the expression matrix
    * @param nThreads number of threads used to run GSEA, 0 if all CPU threads want to be used
    * @param scRna true if rows contain samples and columns genes, false otherwise
    * @post Gene sets, sample ids, gene ids and expression matrix are initialised
    */
    Gsea(vector<GeneSet> &geneSets,
         vector<vector<GeneSample>> &expressionMatrix,
         vector<string> &geneIds,
         vector<string> &sampleIds,
         uint threads,
         bool scRna);

    /**
    * @brief Runs GSEA for the given expression matrix and gene sets in the Gsea creator function
    * @param outFileName name of the output file
    * @param ioutput how many gene sets between status output
    * @post ES results are written into outFileName file
    */
    void run(string outFileName = "", uint ioutput = 10);

    /**
    * @brief Runs GSEA for the given expression matrix using the gene sets initialised in the Gsea creator function
    * @param expressionMatrix matrix containing counts in the cells, samples in the rows and genes in the columns
    * @post The correspinding chunk file contains the ES score for each sample and gene set
    */
    void runChunked(vector<vector<GeneSample>> &expressionMatrix);


    /**
    * @brief Filter the chunked GSEA results by selecting the nFilteredGeneSets gene sets with more variance across the samples
    * @param nFilteredGeneSets number of gene sets selected to be written in the filtered results file
    * @param chunksPath path where the chunks are stored, "" if chunks where generated in the same session with runChunked()
    * @param outFileName name of the filtered results file
    * @pre runChunked has been run at least one time
    * @post outFileName contains the filtered results
    */
    void filterResults(uint nFilteredGeneSets, string chunksPath, string outFilenName);

    /**
    * @brief Normalize the expression matrix using rpm and centering the samples
    * @post The expression matrix is normalized
    */
    void normalizeExprMatrix();

    ~Gsea();
};

#endif
