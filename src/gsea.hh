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

    void enrichmentScoreJob(uint sampleStart, uint sampleEnd);

    void scEnrichmentScoreJob(uint sampleStart, uint sampleEnd);

    void writeResults();

public:
    Gsea();

    Gsea(vector<string> &sampleIdsRcpp,
         vector<string> &geneIdsRcp,
         vector<GeneSet> &geneSets,
         uint nThreads);

    Gsea(vector<GeneSet> &geneSets,
         vector<vector<GeneSample>> &expressionMatrix,
         vector<string> &geneIds,
         vector<string> &sampleIds,
         uint threads,
         bool scRna);

    void run(string outFileName, uint ioutput);

    void runChunked(vector<vector<GeneSample>> &expressionMatrix);

    void filterResults(uint nFilteredGeneSets);

    void normalizeExprMatrix();

    ~Gsea();
};

#endif
