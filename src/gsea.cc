/** @file gsea.cc
 * @brief Gsea implementation file */

#include "gsea.hh"
#include <chrono>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void printTime(system_clock::time_point timePoint)
{
    char timeString[9];
    time_t timePointC = system_clock::to_time_t(timePoint);
    struct tm *tm = localtime(&timePointC);
    strftime(timeString, sizeof(timeString), "%H:%M:%S", tm);
    cout << "[" << timeString << "]";
}

Gsea::Gsea()
{
    system_clock::time_point startIOTime = system_clock::now();

    readConfig();

    if (!scRna)
        readRna();
    else
        readScRna();

    system_clock::time_point endIOTime = system_clock::now();
    cout << "IO elapsed time: " << duration_cast<milliseconds>(endIOTime - startIOTime).count() / 1000.0 << " s" << endl;
}

Gsea::Gsea(vector<string> &sampleIds,
           vector<string> &geneIds,
           vector<GeneSet> &geneSets,
           uint nThreads)
{
    currentSample = 0;
    chunk = 0;

    if (nThreads == 0)
        this->nThreads = thread::hardware_concurrency();
    else
        this->nThreads = nThreads;
    this->geneIds = geneIds;
    this->sampleIds = sampleIds;

    nGenes = geneIds.size();
    nSamples = sampleIds.size();
    nGeneSets = geneSets.size();

    for (uint k = 0; k < nGeneSets; ++k)
    {
        bool hit = false;
        for (uint i = 0; i < nGenes and not hit; ++i)
        {
            if (geneSets[k].geneSet.find(geneIds[i]) != geneSets[k].geneSet.end())
                hit = true;
        }
        if (hit)
            this->geneSets.push_back(geneSets[k]);
    }

    nGeneSets = this->geneSets.size();

    cout << "[GSEA input size]" << endl;
    cout << "Sampled genes: " << nGenes << endl;
    cout << "Samples:       " << nSamples << endl;
    cout << "Gene sets:     " << geneSets.size() << endl;
    cout << endl;
}

Gsea::Gsea(vector<GeneSet> &geneSets,
           vector<vector<GeneSample>> &expressionMatrix,
           vector<string> &geneIds,
           vector<string> &sampleIds,
           uint threads,
           bool scRna)
{
    this->geneSets = geneSets;
    this->expressionMatrix = expressionMatrix;
    nGenes = geneIds.size();
    nSamples = sampleIds.size();
    nGeneSets = geneSets.size();
    this->sampleIds = sampleIds;
    this->geneIds = geneIds;
    if (nThreads == 0)
        this->nThreads = thread::hardware_concurrency();
    else
        this->nThreads = threads;
    this->scRna = scRna;
    this->outputSep = ',';
    results = vector<vector<float>>(geneSets.size(), vector<float>(nSamples));
}

void Gsea::readConfig()
{
    ifstream file("./gsea.config");
    if (!file.is_open())
    {
        file.close();
        ofstream outFile("./gsea.config");
        string aux;
        outFile << "expression-matrix-file:     expression-matrix.csv" << endl;
        outFile << "sep:                        ," << endl;
        outFile << "gene-sets-file:             gene-sets.csv" << endl;
        outFile << "sep:                        ," << endl;
        outFile << "output-file:                results.csv" << endl;
        outFile << "sep:                        ," << endl;
        outFile << "threads-used:               0" << endl;
        outFile << "normalized-data:            0" << endl;
        outFile << "ioutput:                    100" << endl;
        outFile << "scrna:                      0" << endl;
        outFile << "batch-size:                 50" << endl;

        expressionMatrixFilename = "expression-matrix.csv";
        expressionMatrixSep = ',';
        geneSetsFilename = "gene-sets.csv";
        geneSetsSep = ',';
        outputFilename = "results.csv";
        outputSep = ',';
        nThreads = 0;
        normalizedData = false;
        ioutput = 0;
        scRna = 0;
        batchSize = 50;
        outFile.close();
    }
    else
    {
        string aux;
        file >> aux >> expressionMatrixFilename;
        file >> aux >> expressionMatrixSep;
        if (expressionMatrixSep == 't')
            expressionMatrixSep = '\t';
        file >> aux >> geneSetsFilename;
        file >> aux >> geneSetsSep;
        file >> aux >> outputFilename;
        file >> aux >> outputSep;
        file >> aux >> nThreads;
        file >> aux >> normalizedData;
        file >> aux >> ioutput;
        file >> aux >> scRna;
        file >> aux >> batchSize;
    }

    if (nThreads == 0)
        nThreads = thread::hardware_concurrency();

    cout << "[GSEA config]" << endl;
    cout << "expression-matrix-file: " << expressionMatrixFilename << endl;
    cout << "sep:                    " << expressionMatrixSep << endl;
    cout << "gene-sets-file:         " << geneSetsFilename << endl;
    cout << "sep:                    " << geneSetsSep << endl;
    cout << "output-file:            " << outputFilename << endl;
    cout << "sep:                    " << outputSep << endl;
    cout << "threads-used:           " << nThreads << endl;
    cout << "normalized-data:        " << normalizedData << endl;
    cout << "ioutput:                " << ioutput << endl;
    cout << "scrna:                  " << scRna << endl;
    cout << "batch-size:             " << batchSize << endl;
    cout << endl;

    file.close();
}

void Gsea::readRna()
{
    ifstream file(expressionMatrixFilename);
    string line;

    // Read first row (sample ids)
    getline(file, line);
    stringstream ssLine(line);
    string colName;
    uint i = 0;
    while (getline(ssLine, colName, expressionMatrixSep))
    {
        sampleIds.push_back(colName);
        ++i;
    }

    i = 0;
    while (getline(file, line))
    {
        stringstream ssLine(line);

        // Read first column (gene id)
        string rowName;
        getline(ssLine, rowName, expressionMatrixSep);
        geneIds.push_back(rowName);

        string valueStr;
        bool first = true;
        bool nullRow = true;
        uint j = 0;
        while (getline(ssLine, valueStr, expressionMatrixSep))
        {
            if (first)
            {
                expressionMatrix.push_back(vector<GeneSample>(sampleIds.size()));
                first = false;
            }
            float count = stof(valueStr);
            if (nullRow and count != 0)
                nullRow = false;
            expressionMatrix[i][j] = {i, count};
            ++j;
        }
        if (nullRow)
        {
            expressionMatrix.pop_back();
            geneIds.pop_back();
        }
        else
            ++i;
    }

    nGenes = expressionMatrix.size();
    if (nGenes > 0)
        nSamples = expressionMatrix[0].size();
    file.close();

    file = ifstream(geneSetsFilename);
    while (getline(file, line))
    {
        stringstream ssLine(line);

        string rowName;
        getline(ssLine, rowName, geneSetsSep);

        unordered_set<string> genes;
        string valueStr;
        while (getline(ssLine, valueStr, geneSetsSep))
        {
            genes.insert(valueStr);
        }
        geneSets.push_back({rowName, genes});
    }
    file.close();

    nGeneSets = geneSets.size();

    results = vector<vector<float>>(geneSets.size(), vector<float>(nSamples));
}

void Gsea::readScRna()
{
    // Read gene sets file
    ifstream file = ifstream(geneSetsFilename);
    string line;
    uint i = 0;
    while (getline(file, line))
    {
        stringstream ssLine(line);

        string rowName;
        getline(ssLine, rowName, geneSetsSep);

        unordered_set<string> genes;
        string valueStr;
        while (getline(ssLine, valueStr, geneSetsSep))
        {
            genes.insert(valueStr);
        }
        geneSets.push_back({rowName, genes});
        ++i;
    }
    nGeneSets = geneSets.size();
    file.close();
}

void Gsea::runScRna()
{
    ifstream file = ifstream(expressionMatrixFilename);
    ofstream oFile = ofstream(outputFilename);
    string line;

    // Ignore first row, already read
    getline(file, line);

    bool first = true;
    for (auto geneSet : geneSets)
    {
        if (first)
            first = false;
        else
            oFile << ",";
        oFile << geneSet.geneSetId;
    }
    oFile << endl;

    uint nLines = batchSize;
    uint totalLines = nThreads * nLines;
    expressionMatrix = vector<vector<GeneSample>>(totalLines);
    vector<string> sampleNames = vector<string>(totalLines);
    results = vector<vector<float>>(totalLines, vector<float>(geneSets.size()));

    uint i = 0;
    first = true;
    while (getline(file, line))
    {
        stringstream ssLine(line);

        // Read first column (sample id)
        string rowName;
        getline(ssLine, rowName, expressionMatrixSep);
        sampleNames[i % totalLines] = rowName;

        string valueStr;
        uint j = 0;
        while (getline(ssLine, valueStr, expressionMatrixSep))
        {
            if (first)
                expressionMatrix[i % totalLines].push_back({j, stof(valueStr)});
            else
                expressionMatrix[i % totalLines][j] = {j, stof(valueStr)};
            ++j;
        }

        if (i != 0 and i % (totalLines - 1) == 0)
        {
            first = false;
            vector<thread> threads = vector<thread>(nThreads);
            for (uint t = 0; t < nThreads; ++t)
            {
                uint startLine = nLines * t;
                uint endLine = startLine + nLines;
                threads[t] = thread(&Gsea::scEnrichmentScoreJob, this, startLine, endLine);
            }

            for (thread &t : threads)
                t.join();

            for (uint t = 0; t < totalLines; ++t)
            {
                oFile << sampleNames[t];
                for (uint l = 0; l < geneSets.size(); ++l)
                {
                    oFile << results[t][l] << ",";
                }
                oFile << endl;
            }

            system_clock::time_point now = system_clock::now();
            printTime(now);
            cout << " Sample " << i;
            uint ETA = (nSamples - i) * duration_cast<milliseconds>(now - startGSEATime).count() / (i * 60 * 1000);
            cout << " ETA: " << ETA << " min" << endl;
        }
        ++i;
    }

    uint offset = i % totalLines;
    uint linesPerThread = offset / nThreads;
    uint offsetLines = linesPerThread % nThreads;
    vector<thread> threads = vector<thread>(nThreads);
    for (uint t = 0; t < nThreads; ++t)
    {
        uint startLine = linesPerThread * t;
        uint endLine = startLine + linesPerThread;
        if (t == nThreads - 1)
            endLine += offsetLines;
        threads[t] = thread(&Gsea::scEnrichmentScoreJob, this, startLine, endLine);
    }

    for (thread &t : threads)
        t.join();

    for (uint t = 0; t < offset; ++t)
    {
        oFile << sampleNames[t];
        for (uint l = 0; l < geneSets.size(); ++l)
        {
            oFile << results[t][l] << ",";
        }
        oFile << endl;
    }

    file.close();
}

void Gsea::rpm()
{
    for (uint j = 0; j < nSamples; ++j)
    {
        float sum = 0;

        for (uint i = 0; i < nGenes; ++i)
            sum += expressionMatrix[i][j].count;

        float multFactor = 1000000 / sum;
        for (uint i = 0; i < nGenes; ++i)
            expressionMatrix[i][j].count *= multFactor;
    }
}

void Gsea::meanCenter()
{
    for (uint i = 0; i < nGenes; ++i)
    {
        float mean = 0;
        for (uint j = 0; j < nSamples; ++j)
            mean += expressionMatrix[i][j].count;
        mean /= nSamples;

        for (uint j = 0; j < nSamples; ++j)
            expressionMatrix[i][j].count -= mean;
    }
}

bool Gsea::geneSampleComp(const GeneSample &g1, const GeneSample &g2)
{
    return g1.count > g2.count;
}

bool Gsea::geneSetPtrComp(const GeneSetPtr &g1, const GeneSetPtr &g2)
{
    return g1.value > g2.value;
}

void Gsea::sortColumnsJob(uint startSample, uint endSample)
{
    for (uint j = startSample; j < endSample; ++j)
    {
        vector<GeneSample> column = vector<GeneSample>(nGenes);
        for (uint i = 0; i < nGenes; ++i)
            column[i] = expressionMatrix[i][j];
        sort(column.begin(), column.end(), &Gsea::geneSampleComp);
        for (uint i = 0; i < nGenes; ++i)
            expressionMatrix[i][j] = column[i];
    }
}

void Gsea::enrichmentScore()
{
    uint samplesPerThread = nSamples / nThreads;
    uint offset = nSamples % samplesPerThread;
    vector<thread> threads = vector<thread>(nThreads);
    for (uint i = 0; i < nThreads; ++i)
    {
        uint startSample = i * samplesPerThread;
        uint endSample = startSample + samplesPerThread;
        if (i == nThreads - 1)
            endSample += offset;
        threads[i] = thread(&Gsea::enrichmentScoreJob, this, startSample, endSample);
        if (i == nThreads - 1)
        {
            auto threadId = threads[i].get_id();
            logThread = *static_cast<unsigned int *>(static_cast<void *>(&threadId));
        }
    }

    for (thread &t : threads)
        t.join();
}

void Gsea::enrichmentScoreJob(uint startSample, uint endSample)
{
    assert(endSample <= nSamples);

    sortColumnsJob(startSample, endSample);

    auto threadId = this_thread::get_id();
    uint id = *static_cast<unsigned int *>(static_cast<void *>(&threadId));

    for (uint k = 0; k < nGeneSets; ++k)
    {
        float geneSetSize = geneSets[k].geneSet.size();
        float posScore = sqrt((nGenes - geneSetSize) / geneSetSize);
        float negScore = -sqrt((geneSetSize / (nGenes - geneSetSize)));
        for (uint j = startSample; j < endSample; ++j)
        {
            float currentValue = 0;
            float maxValue = 0;
            for (uint i = 0; i < nGenes; ++i)
            {
                if (geneSets[k].geneSet.find(geneIds[expressionMatrix[i][j].geneId]) != geneSets[k].geneSet.end())
                {
                    currentValue += posScore;
                }
                else
                    currentValue += negScore;

                if (i == 0)
                    maxValue = currentValue;
                maxValue = max(currentValue, maxValue);
            }
            results[k][j] = maxValue;
        }
        if (k != 0 and id == logThread and k % ioutput == 0)
        {
            system_clock::time_point now = system_clock::now();
            printTime(now);
            cout << " Gene set " << k;

            ulong ETA = (geneSets.size() - k) * duration_cast<milliseconds>(now - startGSEATime).count() / (k * 60 * 1000);
            cout << " ETA: " << ETA << " min" << endl;
        }
    }
}

void Gsea::scEnrichmentScoreJob(uint startSample, uint endSample)
{
    assert(endSample <= nSamples);

    for (uint i = startSample; i < endSample; ++i)
    {
        sort(expressionMatrix[i].begin(), expressionMatrix[i].end(), &Gsea::geneSampleComp);
    }

    for (uint i = startSample; i < endSample; ++i)
    {
        for (uint k = 0; k < nGeneSets; ++k)
        {
            float geneSetSize = geneSets[k].geneSet.size();
            float posScore = sqrt((nGenes - geneSetSize) / geneSetSize);
            float negScore = -sqrt((geneSetSize / (nGenes - geneSetSize)));
            float currentValue = 0;
            float maxValue = 0;
            for (uint j = 0; j < nGenes; ++j)
            {
                if (expressionMatrix[i][j].count == 0)
                    break;

                if (geneSets[k].geneSet.find(geneIds[expressionMatrix[i][j].geneId]) != geneSets[k].geneSet.end())
                {
                    currentValue += posScore;
                }
                else
                    currentValue += negScore;

                if (j == 0)
                    maxValue = currentValue;
                maxValue = max(currentValue, maxValue);
            }
            results[i][k] = maxValue;
        }
    }
}

void Gsea::writeResults()
{
    ofstream file(outputFilename);
    bool first = true;
    for (string &sampleId : sampleIds)
    {
        if (first)
            first = false;
        else
            file << outputSep;
        file << sampleId;
    }
    file << endl;

    vector<string> geneSetNames = vector<string>(geneSets.size());
    uint i = 0;
    for (auto &it : geneSets)
    {
        geneSetNames[i] = it.geneSetId;
        ++i;
    }

    i = 0;
    for (vector<float> &row : results)
    {
        file << geneSetNames[i];
        for (float value : row)
        {
            file << outputSep << value;
        }
        file << endl;
        ++i;
    }
    file.close();
}

void Gsea::runRna()
{
    enrichmentScore();

    writeResults();
}

void Gsea::run(string outFileName, uint ioutput)
{
    if (outFileName != "")
        this->outputFilename = outFileName;
    
    if (ioutput != 10)
        this->ioutput = ioutput;

    cout << "[GSEA input size]" << endl;
    cout << "Sampled genes: " << nGenes << endl;
    cout << "Samples:       " << nSamples << endl;
    cout << "Gene sets:     " << geneSets.size() << endl;
    cout << endl;

    startGSEATime = system_clock::now();

    printTime(startGSEATime);
    cout << " Started GSEA" << endl;

    if (scRna)
        runScRna();
    else
        runRna();

    cout << endl
         << "Elapsed time: " << duration_cast<minutes>(system_clock::now() - startGSEATime).count() << " min" << endl;
    cout << "Results written in " << outputFilename << endl;
}

void Gsea::runChunked(vector<vector<GeneSample>> &expressionMatrix)
{
    if (currentSample == 0)
        startGSEATime = system_clock::now();

    uint chunkSamples = expressionMatrix.size();
    if (chunkSamples > 0)
        nGenes = expressionMatrix[0].size();
    this->expressionMatrix = expressionMatrix;

    if (chunk == 0)
    {
        filesystem::path tmpPath = filesystem::temp_directory_path();
        ulong id = duration_cast<seconds>(startGSEATime.time_since_epoch()).count();
        filesystem::path chunksFolder = filesystem::path("chunks" + to_string(id));
        chunksPath = tmpPath / chunksFolder;
        if (not filesystem::exists(chunksPath))
            filesystem::create_directory(chunksPath);
        cout << "Chunks path: " << chunksPath << endl
             << endl;

        printTime(system_clock::now());
        cout << " Started GSEA" << endl;
    }

    vector<thread> threads = vector<thread>(nThreads);
    uint samplesPerThread = chunkSamples / nThreads;
    uint offset = chunkSamples % nThreads;

    results = vector<vector<float>>(chunkSamples, vector<float>(nGeneSets));

    for (uint t = 0; t < nThreads; ++t)
    {
        uint threadStartSample = samplesPerThread * t;
        uint threadEndSample = threadStartSample + samplesPerThread;
        if (t == nThreads - 1)
        {
            threadEndSample += offset;
            auto threadId = threads[t].get_id();
            logThread = *static_cast<unsigned int *>(static_cast<void *>(&threadId));
        }
        threads[t] = thread(&Gsea::scEnrichmentScoreJob, this, threadStartSample, threadEndSample);
    }

    for (thread &t : threads)
        t.join();

    filesystem::path chunkFile = filesystem::path(to_string(chunk));
    filesystem::path chunkPath = chunksPath / chunkFile;
    ofstream resultsFile(chunkPath);
    for (uint k = 0; k < nGeneSets; ++k)
    {
        for (uint i = 0; i < chunkSamples; ++i)
        {
            if (i != 0)
                resultsFile << ",";
            resultsFile << results[i][k];
        }
        resultsFile << endl;
    }
    resultsFile.close();

    system_clock::time_point now = system_clock::now();
    printTime(now);
    currentSample += chunkSamples;
    ++chunk;
    ulong ETA = (nSamples - currentSample) * duration_cast<milliseconds>(now - startGSEATime).count() / (currentSample * 60 * 1000);
    cout << " Sample: " << currentSample << " ETA: " << ETA << " min" << endl;
}

void Gsea::filterResults(uint nFilteredGeneSets, string chunksPathStr, string outFileName)
{
    assert(nFilteredGeneSets < nGeneSets);
    vector<GeneSetPtr> geneSetsVar = vector<GeneSetPtr>(nGeneSets);
    string line;

    if (chunksPathStr != "")
        chunksPath = filesystem::path(chunksPathStr);

    uint nChunks = 0;
    if (filesystem::exists(chunksPath))
        nChunks = (std::size_t)std::distance(std::filesystem::directory_iterator{chunksPath}, std::filesystem::directory_iterator{});

    if (nChunks == 0)
    {
        cerr << "No chunk files found in chunks/ directory" << endl;
        return;
    }

    vector<ifstream> chunkFiles = vector<ifstream>(nChunks);
    for (uint i = 0; i < nChunks; ++i)
    {
        filesystem::path chunkFile = filesystem::path(to_string(i));
        filesystem::path chunkPath = chunksPath / chunkFile;
        chunkFiles[i].open(chunkPath);
    }

    for (uint i = 0; i < nGeneSets; ++i)
    {
        float mean = 0;
        uint k = 0;
        vector<float> rowValues(nSamples);
        for (uint j = 0; j < nChunks; ++j)
        {
            getline(chunkFiles[j], line);
            stringstream ssLine(line);
            string valueStr;
            while (getline(ssLine, valueStr, ','))
            {
                float value = stof(valueStr);
                rowValues[k] = value;
                mean += rowValues[k];
                ++k;
            }
        }
        mean /= nSamples;
        float variance = 0;
        for (uint j = 0; j < nSamples; ++j)
        {
            variance += pow((rowValues[j] - mean), 2);
        }
        variance /= nSamples;

        geneSetsVar[i] = {i, variance};
    }

    sort(geneSetsVar.begin(), geneSetsVar.end(), &Gsea::geneSetPtrComp);
    ofstream variance("var");
    for (auto x : geneSetsVar)
        variance << geneSets[x.geneSetPtr].geneSetId << " " << x.value << endl;

    unordered_set<string> filteredSets;
    for (uint i = 0; i < nFilteredGeneSets; ++i)
        filteredSets.insert(geneSets[geneSetsVar[i].geneSetPtr].geneSetId);

    ofstream filteredResultsFile(outFileName);
    for (uint i = 0; i < nSamples; ++i)
    {
        if (i != 0)
            filteredResultsFile << ",";
        filteredResultsFile << sampleIds[i];
    }
    filteredResultsFile << endl;

    for (uint i = 0; i < nChunks; ++i)
    {
        chunkFiles[i].clear();
        chunkFiles[i].seekg(0, ios::beg);
    }
    for (uint i = 0; i < nGeneSets; ++i)
    {
        bool filtered = filteredSets.find(geneSets[i].geneSetId) != filteredSets.end();
        if (filtered)
            filteredResultsFile << geneSets[i].geneSetId;
        for (uint j = 0; j < nChunks; ++j)
        {
            getline(chunkFiles[j], line);
            if (filtered)
                filteredResultsFile << "," << line;
        }
        if (filtered)
            filteredResultsFile << endl;
    }
}

void Gsea::normalizeExprMatrix()
{
    rpm();

    meanCenter();
}

Gsea::~Gsea() {}
