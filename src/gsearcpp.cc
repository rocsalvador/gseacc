/** @file gsearcpp.cc
 * @brief GseaRcpp implementation file */

#include "gsearcpp.hh"


GseaRcpp::GseaRcpp(CharacterVector sampleIdsRcpp,
                   CharacterVector geneIdsRcpp,
                   List geneSetsRcpp,
                   uint nThreads)
{
    vector<string> sampleIds = as<vector<string>> (sampleIdsRcpp);
    vector<string> geneIds = as<vector<string>> (geneIdsRcpp);

    CharacterVector geneSetsIdsRcpp = geneSetsRcpp.names();
    vector<GeneSet> geneSets = vector<GeneSet> (geneSetsRcpp.size());
    vector<string> geneSetsIds = as<vector<string>>(geneSetsIdsRcpp);
    for (uint i = 0; i < geneSetsRcpp.length(); ++i)
    {
        CharacterVector geneSetRcpp = geneSetsRcpp[i];
        vector<string> geneVector = as<vector<string>>(geneSetRcpp);
        unordered_set<string> geneSet = unordered_set<string>(geneVector.begin(), geneVector.end());
        geneSets[i] = {geneSetsIds[i], geneSet};
    }

    gsea = new Gsea(sampleIds, geneIds, geneSets, nThreads);
}

GseaRcpp::GseaRcpp(NumericMatrix expressionMatrixRcpp,
                   List geneSetsRcpp,
                   uint threads)
{
    vector<string> sampleIds, geneIds;
    sampleIds = as<vector<string>> (colnames(expressionMatrixRcpp));
    geneIds = as<vector<string>> (rownames(expressionMatrixRcpp));
    double nGenes = expressionMatrixRcpp.ncol();
    double nSamples = expressionMatrixRcpp.nrow();
    vector<vector<GeneSample>> expressionMatrix(nSamples, vector<GeneSample>(nGenes));
    for (uint i = 0; i < nSamples; ++i)
    {
        for (uint j = 0; j < nGenes; ++j)
        {
            expressionMatrix[i][j] = {i, float(expressionMatrixRcpp(i, j))};
        }
    }

    CharacterVector geneSetsIdsRcpp = geneSetsRcpp.names();
    vector<GeneSet> geneSets = vector<GeneSet> (geneSetsRcpp.size());
    vector<string> geneSetsIds = as<vector<string>>(geneSetsIdsRcpp);
    for (uint i = 0; i < geneSetsRcpp.length(); ++i)
    {
        CharacterVector geneSetRcpp = geneSetsRcpp[i];
        vector<string> geneVector = as<vector<string>>(geneSetRcpp);
        unordered_set<string> geneSet = unordered_set<string>(geneVector.begin(), geneVector.end());
        geneSets[i] = {geneSetsIds[i], geneSet};
    }

    gsea = new Gsea(geneSets, expressionMatrix, geneIds, sampleIds, threads, false);
}

void GseaRcpp::runChunked(const NumericMatrix &countMatrixRcpp)
{
    double nGenes = countMatrixRcpp.ncol();
    double nSamples = countMatrixRcpp.nrow();
    vector<vector<GeneSample>> expressionMatrix(nSamples, vector<GeneSample>(nGenes));
    for (uint i = 0; i < nSamples; ++i)
    {
        for (uint j = 0; j < nGenes; ++j)
        {
            expressionMatrix[i][j] = {i, float(countMatrixRcpp(i, j))};
        }
    }

    gsea->runChunked(expressionMatrix);
}

void GseaRcpp::filterResults(uint nFilteredGeneSets, string chunksPath, string outFileName)
{
    gsea->filterResults(nFilteredGeneSets, chunksPath, outFileName);
}

void GseaRcpp::run(string outFileName, uint ioutput)
{
    gsea->run(outFileName, ioutput);
}

void GseaRcpp::normalizeExprMatrix()
{
    gsea->normalizeExprMatrix();
}

GseaRcpp::~GseaRcpp()
{
    delete gsea;
}

