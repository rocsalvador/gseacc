#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include "Rcpp/vector/instantiation.h"
#include "gsea.hh"
#include "gsearcpp.hh"
using namespace Rcpp;
using namespace std;

RCPP_MODULE(GseaModule) {
    class_<GseaRcpp>("GseaRcpp")
    .constructor<CharacterVector, CharacterVector, List, uint>()
    .constructor<NumericMatrix, List, uint>()
    .method("runChunked", &GseaRcpp::runChunked)
    .method("filterResults", &GseaRcpp::filterResults)
    .method("run", &GseaRcpp::run)
    .method("normalizeExprMatrix", &GseaRcpp::normalizeExprMatrix)
    ;
}

// [[Rcpp::export]]
void writeGeneSets(List geneSetsRcpp, String fileName)
{
    ofstream file(fileName);
    unordered_map<string, unordered_set<string>> geneSets;
    CharacterVector geneSetsIdsRcpp = geneSetsRcpp.names();
    vector<string> geneSetsIds = as<vector<string>>(geneSetsIdsRcpp);
    for (uint i = 0; i < geneSetsRcpp.length(); ++i)
    {
        file << geneSetsIds[i] << ",";
        CharacterVector geneSetRcpp = geneSetsRcpp[i];
        vector<string> geneVector = as<vector<string>>(geneSetRcpp);
        for (uint j = 0; j < geneVector.size(); ++j)
        {
            if (j != 0)
                file << ",";
            file << geneVector[j];
        }
        file << endl;
    }
}

// [[Rcpp::export]]
List readGeneSets(std::string fileName)
{
    ifstream file(fileName);
    std::string line;
    getline(file, line);
    stringstream ssLine(line);
    List geneSets;
    CharacterVector geneSetIds;
    while (getline(file, line))
    {
        list<string> geneSet;
        stringstream ssLine(line);
        string rowName;
        getline(ssLine, rowName, ',');
        geneSetIds.push_back(rowName);

        string valueStr;
        while (getline(ssLine, valueStr, ','))
        {
            geneSet.push_back(valueStr);
        }
        geneSets.push_back(geneSet);
    }
    geneSets.names() = geneSetIds;
    return geneSets;
}

// [[Rcpp::export]]
NumericMatrix readCsv(String fileName, char sep = ',', bool hasRowNames = true, bool hasColNames = true) {
    ifstream file(fileName);
    string line;


    uint nCols = 0;
    uint nRows = 0;
    bool first = true;
    while (getline(file, line))
    {
        if (first) {
            stringstream ssLine(line);
            string valueStr;
            while (getline(ssLine, valueStr, sep))
                ++nCols;
            first = false;
        }
        ++nRows;
    }

    if (hasRowNames) --nRows;

    file.clear();
    file.seekg(0, ios::beg);


    NumericMatrix matrix(nRows, nCols);
    CharacterVector rowNames(nRows);
    CharacterVector colNames(nCols);
    uint i = 0;
    if (hasColNames) {
        getline(file, line);
        stringstream ssLine(line);
        string valueStr;
        while (getline(ssLine, valueStr, sep))
        {
            colNames(i) = valueStr;
            ++i;
        }
    }

    i = 0;
    while (getline(file, line))
    {
        stringstream ssLine(line);

        if (hasRowNames) {
            string rowName;
            getline(ssLine, rowName, sep);
            rowNames(i) = rowName;
        }

        string valueStr;
        uint j = 0;
        while (getline(ssLine, valueStr, sep))
        {
            matrix(i, j) = stof(valueStr);
            ++j;
        }
        ++i;
    }

    if (hasRowNames)
        rownames(matrix) = rowNames;
    if (hasColNames)
        colnames(matrix) = colNames;

    return matrix;
}
