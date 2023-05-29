/** @file rcpp.cc
 * @brief Rcpp exports file */

#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include "gsea.hh"
#include "gsearcpp.hh"
using namespace Rcpp;
using namespace std;

RCPP_EXPOSED_CLASS(GseaRcpp)
RCPP_MODULE(GseaModule) {
    class_<GseaRcpp>("Gsea")
    .constructor<CharacterVector, CharacterVector, List, uint>()
    .constructor<NumericMatrix, List, uint>()
    .method("runChunked", &GseaRcpp::runChunked, "Run GSEA for the expression matrix chunk")
    .method("filterResults", &GseaRcpp::filterResults)
    .method("run", &GseaRcpp::run)
    .method("normalizeExprMatrix", &GseaRcpp::normalizeExprMatrix)
    ;
}


/** @brief Write the gene sets list into fileName file, each line has as a first element the gene set id and then gene ids of the gene set
* @param geneSetsRcpp list that has a gene set in every elemnt
* @param fileName output file name
* @pre fileName is not a null string
* @post Gene sets written into fileName file
*/
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

/**
 * @brief Read the gene sets list from the fileName file, each line has as a first element the gene set id and then gene ids of the gene set
 * @param fileName input file name
 * @pre fileName is a valid gene sets file
 * @post Return the gene sets as a list of lists
 */
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

/**
 * @brief Read a csv file as a numeric matrix
 * @param fileName input file name
 * @param sep Csv element separator, default: sep = ','
 * @param hasRowNames true if the csv has row names, false otherwise
 * @param hasColNames true if the csv has column names, false otherwise
 * @pre fileName is a valid csv file
 * @post Return the csv file as a numeric matrix
 */
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

    if (hasColNames) --nRows;
    if (hasRowNames) --nCols;

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
