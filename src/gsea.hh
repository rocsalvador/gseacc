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

struct GeneSample
{
    uint32_t geneId;
    float count;
};

struct GeneSet
{
    string geneSetId;
    unordered_set<string> geneSet;
};

struct GeneSetPtr
{
    uint geneSetPtr;
    float value;
};

/** @class Gsea
 * @brief  */
class Gsea
{
private:
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
    uint logThread;

    bool normalizedData;
    bool scRna;

    uint currentSample;

    vector<GeneSet> geneSets;

    vector<vector<GeneSample>> expressionMatrix;
    vector<string> sampleIds;
    vector<string> geneIds;

    vector<vector<float>> results;

    uint nGenes;
    uint nSamples;
    uint nGeneSets;

    system_clock::time_point startGSEATime;

    static bool geneSampleComp(const GeneSample &g1, const GeneSample &g2);

    static bool geneSetPtrComp(const GeneSetPtr &g1, const GeneSetPtr &g2);

    void readRna();

    void readScRna();

    void runScRna();

    void runRna();

    void readConfig();

    void rpm();

    void meanCenter();

    void sortColumnsJob(uint columnStart, uint columnEnd);

    void enrichmentScore();

    /**
    * @brief Runs the gsea from lineStart to lineEnd samples, assuming genes in the rows and samples in the columns
    * @param startSample start sample
    * @param endSample end sample
    * @pre expressionMatrix rows contain genes, expressionMatrix columns contain samples
    * @post The samples startSample to endSample in the results matrix contain the ES
    */
    //' @name Double$new
    //' @title Constructs a new Double object
    //' @param v A value to encapsulate
    void enrichmentScoreJob(uint sampleStart, uint sampleEnd);

    /**
    * @brief Runs the gsea from startSample to endSample samples, assuming samples in the rows and genes in the columns
    * @param startSample start sample
    * @param endSample end sample
    * @pre expressionMatrix rows contain samples, expressionMatrix columns contain genes
    * @post The samples startSample to endSample in the results matrix contain the ES
    */
    void scEnrichmentScoreJob(uint sampleStart, uint sampleEnd);

    void writeResults();

public:
    /**
    * @brief Gsea creator function to use the class without R, it reads the configuration from
    * gsea.config, and all input data from files
    * @post Gene sets and the expression matix (only in Rna-seq case) are initialised
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
    * @pre runChunked has been run at least one time
    * @post Filtered-results.csv contains the filtered results
    */
    void filterResults(uint nFilteredGeneSets);

    /**
    * @brief Normalize the expression matrix using rpm and centering the samples
    * @post The expression matrix is normalized
    */
    void normalizeExprMatrix();

    ~Gsea();
};

#endif
