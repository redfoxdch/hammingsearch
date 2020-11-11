/**
    @function: kmeans based on neiborhood
    @author: chenghao deng
    @version: 2.0
    @date: 2017-2-27
    @institute: Xiamen University


**/
#ifndef NNG_H
#define NNG_H
#include "util/distance.hpp"


class NNG
{
    bool EXTData;
    float *data;
    unsigned count, ndim;
    unsigned TreeNum, TreeSize, HubNum, PoolSize, NNIter, TopK, MaxSize, Check;

    const DISTANCE::Distance<float>* distance;

    unsigned *topList;
    float *topDist;

    int *labels;
    float *radius;
    unsigned clnumb;

    int adjustNN();
    int initMemory();
    int getTopDistance();
    int buildTreeGraph();
    int initNNGraph();
    int treeUp();
    int NNExtent();
    int NNDescent();
    int Extent();
    int NNHub();

    int label2clusters(vector<unsigned> *);

public :
    NNG(const DISTANCE::Distance<float>* d);
    virtual ~NNG();
    bool loadData(const char *srcfn);
    bool loadData(unsigned char *mat, const int row, const int dim);
    bool loadData(float *mat, const int row, const int dim);
    int showRecall();

    int reRank();
    int initNNGraph(unsigned *TopList);
    int getTopList(unsigned *&top, unsigned &topNum);
    int NNBuild(const char *srcfn, const char* dstfn, unsigned TreeNum, unsigned TreeSize, unsigned HubNum, unsigned PoolSize, unsigned NNIter, unsigned TopK, unsigned MaxSize);
    int NNBuild(unsigned type = 0);
    int saveNNGraph(const char *dstfn);

    static void test();
    static void CHTest();


};





#endif
