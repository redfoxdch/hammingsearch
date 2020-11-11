#ifndef HE_H
#define HE_H

#include "util/distance.hpp"
#include <vector>
#include <thread>


#define LEN 64
#define COL 128
#define WORD_LEN 20000
//#define NTHREADS 4

using namespace std;

class HE
{
    unsigned NTHREADS = 40;

    const float *projection;///64*128
    const double *median;///64*20000 ts

    const DISTANCE::Distance<float> * distance;
    static double *transpose(const double *median, int row, int col);

    static void miniCoding(const float *sift, const unsigned *words, long row, unsigned *heCodes,const DISTANCE::Distance<float> * distance,const float *projection,const double *median);
    static void miniAdd(const float *sift, const unsigned *words, long row, vector<double>  *median,const DISTANCE::Distance<float> * distance,const float *projection);

public :
    HE() {}
    HE(const float *projection, const DISTANCE::Distance<float> *distance = new DISTANCE::L2DistanceAVXr4<float>, unsigned NTHREADS = 0)
    {
        this->projection = projection;
        this->median = NULL;
        this->distance = distance;
        this->NTHREADS = NTHREADS;
        if(NTHREADS == 0)
            this->NTHREADS = thread::hardware_concurrency();
        if(this->NTHREADS < 1)
            this->NTHREADS = 1;
    }
    HE(const float *projection, const double *median, const DISTANCE::Distance<float> *distance, unsigned NTHREADS = 0)
    {
        this->projection = projection;
        //this->median = transpose(median,LEN,WORD_LEN);
        this->median = median;
        this->distance = distance;
        this->NTHREADS = NTHREADS;
        if(NTHREADS == 0)
            this->NTHREADS = thread::hardware_concurrency();
        if(this->NTHREADS < 1)
            this->NTHREADS = 1;
    }
    void genMedian(const char *sift_file_name, const char *words_file, const char *dst) const;
    void enCoding(const char *sift_file_name, const char *words_file, const char *dst) const;
    unsigned *enCoding(const float *sift, const unsigned *words, long row) const;

    float static heScore(const unsigned *a, const unsigned *b);
    float static heScore(const pair<unsigned,unsigned> a, const pair<unsigned, unsigned> b);


    void static test();

    virtual ~HE()
    {
        if(median != NULL)
        {
           // delete []median;
            median = NULL;
        }
    }
};





#endif
