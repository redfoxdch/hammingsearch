#ifndef BRUTEFORCENN_H
#define BRUTEFORCENN_H

#include "util/distance.hpp"
#include <map>

//#define NTHREADS 4

using namespace std;

class BruteForceNN
{
    unsigned NTHREADS ;
    const DISTANCE::Distance<float>* distance;
    const float *centers;
    unsigned crow = 20000;

public:

    BruteForceNN(const float *centers,unsigned crow,const DISTANCE::Distance<float>* d,unsigned NTHREADS = 0);
    BruteForceNN();
    virtual ~BruteForceNN();

    pair<int, double> * buildNNG(float *data, long row, int col, int topK);
    unsigned *BOWQuantization(const float *sift, long row, int col,const  float *centers, int crow)const ;

    void BOWQuantization(const char *sift_name, const char *word_name) const;


    static void test();
};



#endif
