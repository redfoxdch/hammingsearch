/**
    @function: kmeans based on neiborhood
    @author: chenghao deng
    @version: 2.0
    @date: 2016-2-27
    @institute: Xiamen University


**/
#ifndef NNGKMeans_H
#define NNGKMeans_H

#include "abstractkmeans.h"
#include <set>
#include <unordered_set>

class NNGKMeans: public AbstractKMeans
{
    const static unsigned NBIter, TreeNum;
    const static unsigned NIter;
    const static unsigned RangeFirst;
    const static unsigned RangeSecond;
    const static float Err0;

    bool outNN = false;

    unsigned long longloc;///long compute;

    double *Es, *arrayD;
    int     *Ns;
    double *Ds;
    unsigned *topList;
    float *topDist;

    unordered_set<unsigned> *visitList;

    float *radius;


    const static unsigned long Treesize;
    const static unsigned long TopK;
    const static unsigned long NNTop;

    int adjustNN();
    int optzI2(bool verbose, const unsigned NNTop,  const unsigned niter0);
    int NNKM(bool verbose, const unsigned NNTop, const unsigned iter0);
    int initNNGraph();
    int buildNNGraph();
    int updateNNGraph();
    double appRecall();

public :
    NNGKMeans();
    virtual ~NNGKMeans();
    bool   init(const char *srcfn);
    bool   init(unsigned char *mat, const int row, const int dim);
    bool   init(float *mat, const int row, const int dim);

    int recall();
    int NNI2(bool verbose, const unsigned niter0);
    int buildNN(const char *srcfn, const char *dstfn);
    int buildCluster(unsigned* NN, float* data, unsigned row, unsigned ndim, unsigned clustnumb);
    int buildCluster(const char *srcfn, const char *nnfn, const char *dstfn, unsigned clustnumb);

    int ITNNGraph();

    int copy2TopList(unsigned *top);

    int copyFromTopList(unsigned *top, unsigned topNum);

    bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);

    int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);
    bool   rndSeeds(const int k, int rseeds[], const int bound);
    bool   kppSeeds(const int k, int rseeds[], const int bound);
    void   saveCenters(const char *dstfn, bool append);
    void   saveKNNGraph(const char *dstfn, unsigned int k0);
    void   saveKNNGraph(const char *dstfn, const unsigned int n, const unsigned int k0, vector<unsigned int> &ids);
    int    fetchCenters(float *centers);

    static void test();

};





#endif
