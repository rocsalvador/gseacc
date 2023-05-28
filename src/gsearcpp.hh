#include <Rcpp.h>
#include "gsea.hh"
using namespace Rcpp;

class GseaRcpp
{
private:
    Gsea *gsea;

public:
    GseaRcpp(CharacterVector sampleIdsRcpp,
             CharacterVector geneIdsRcp,
             List geneSetsRcpp,
             uint nThreads);

    GseaRcpp(NumericMatrix expressionMatrix,
             List geneSets,
             uint threads = 1);

    void runChunked(const NumericMatrix &countMatrixRcpp);

    void filterResults(uint nFilteredGeneSets);

    void run(string outFileName, uint ioutput);

    void normalizeExprMatrix();

    ~GseaRcpp();
};
