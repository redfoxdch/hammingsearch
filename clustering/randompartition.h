/**
    @function: clustering based on random parition
    @author: chenghao deng
    @version:1.0
    @date:2016-11-27
    @institute: XiaMen University


**/
#ifndef RANDOMPARTITION_H
#define RANDOMPARTITION_H

#include "abstractkmeans.h"
#include "util/distance.hpp"
#include "nnitem.h"

#include <set>
#include <map>
#include <unordered_set>

class RandomPartition: public AbstractKMeans
{
    const DISTANCE::Distance<float>* distance;
    const static unsigned NTRAILS;
    const static unsigned NIter;
    const static float Err0;
    float *D1, *D2, *C1, *C2, *tmpC1, *tmpC2;
    double *arrayD;
    int     *Ns;

    vector<unsigned> *clusters;

public :
    RandomPartition();
    RandomPartition(const DISTANCE::Distance<float>* );
    virtual ~RandomPartition();
    bool   init(const char *srcfn);
    bool   init(unsigned char *mat, const int row, const int dim);
    bool   init(float *mat, const int row, const int dim);

    int I2P(unsigned label, unsigned nlabel);
    int RP(unsigned label, unsigned nlabel);
    int KD(unsigned label, unsigned nlabel);
    int var(unsigned label, float *v);

    int random(unsigned label, unsigned nlabel);

    bool   config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose);

    int    clust(const unsigned int clust_num, const char *dstfn, const int verbose);
    bool   rndSeeds(const int k, int rseeds[], const int bound);
    bool   kppSeeds(const int k, int rseeds[], const int bound);
    void   saveCenters(const char *dstfn, bool append);
    int    fetchCenters(float *centers);

    static void test();


};





#endif
