#ifndef GEOHENN_H
#define GEOHENN_H

#include "he.h"
#include "geo.h"
#include "info.h"
#include <string>
#include <thread>

#define WORD_LEN 20000
#define COL 128
#define SPARSE_THR 10

class GEOHENN
{
    unsigned NTHREADS = 40;
    float *projection;
    double *median;
    float *centers;
    const DISTANCE::Distance<float>* distance;
    unsigned *word;
    unsigned *he;
    unsigned row;



    static void miniTOPK( long st,const  unsigned *word, const  unsigned *he,const  float *geo,const   unsigned *wmap,  unsigned row,
                          const  unsigned *words, const  unsigned *hes,const   float *geos, const unsigned * sz,  long num, float *score);

public :

    INFO *storeInfo;///keep in memory

    GEOHENN();


    /////used in python

    GEOHENN(string projectionfile, string medianfile, string centersfile, string pic_info, string szfile,
            string wordsfile, string hesfile, string geosfile);
    GEOHENN(float *projection, double *median, float *centers, const DISTANCE::Distance<float> *distance,unsigned NTHREADS = 0);

    GEOHENN(float *projection, double *median, float *centers, const DISTANCE::Distance<float> *distance,const char *szfile,
            const char * wordsfile, const char * hesfile, const char * geosfile,unsigned NTHREADS = 0);

    unsigned *getWord()
    {
        return this->word;
    }
    unsigned *getHE()
    {
        return this->he;
    }


///
    pair<long, float> *findTop(const float *sift, const float *geo, unsigned row,  int &TOPK, float &value)
    {
        pair<long, float> * top = findTop(sift,geo,row,this->storeInfo,TOPK, value);
        return top;
    }
    pair<long, float> *findTop(const float *sift, const float *geo, unsigned row, INFO *storedInfo, int &TOPK, float &value);

    pair<long, float> *findTop(const float *sift,const float *geo,unsigned row, const unsigned *words,
                               const unsigned *hes, const float *geos, const unsigned *sz, long num, int &TOPK, float &value);

    static pair<long, float> *findTop(const  unsigned *word, const  unsigned *he, const float *geo, unsigned row, const  unsigned *words,
                                      const  unsigned *hes, const  float *geos,  const unsigned * sz,  long num,  int &TOPK,float &value, unsigned NTHREADS = 0);
    static float PairScore(const unsigned *word, const unsigned *he,const float *geo, const unsigned *wmap,const  int row, const unsigned *words,
                           const unsigned *hes, const float *geos, const int srow);
    static pair<long, float> *findTop(const  unsigned *word, const  unsigned *he, const float *geo, unsigned row, INFO *storedInfo, int &TOPK, float &value);

    virtual ~GEOHENN();

};



#endif
