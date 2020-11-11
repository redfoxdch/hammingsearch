#ifndef MISSIONAGENT_H
#define MISSIONAGENT_H

/**
* in charge of distributing the tasks
*
*
*@author:  Cheng-Hao Deng
*@date:    2018.7.19
*@version: 1.0
*
**/

#include <cstring>
#include <string>
#include <map>
#include "info.h"
#include "nn/geohenn.h"

using namespace std;

class MissionAgent
{

public:


/////used in python
    static string * matchANDupdate(unsigned uid, unsigned entry_id, string url,float *sift,float *geo, unsigned row, GEOHENN *mgnn, int &TOPK);
    MissionAgent() {}
    virtual ~MissionAgent() {}
    static bool genCenters(const char *sift_npy, const char * centers, int clustNum);
    static bool HEQuantization(const char *sift_npy,const char *projectionfile,const char *medianfile, const char *wordfile,const char *HEWords);

    static bool findGEOTop(const char *siftfile,const char *geofile, const char *shapefile, const char * projectionfile,
                           const char * medianfile, const char * centersfile, const char * wordsfile, const char * hesfile, const char *geosfile, const char *szfile, int TOPK);
    static bool nn(const char *sift_npy, const char *centers, const char *ALLwords, const char *ALLHEWords);
    static bool BruteForceSearch(const char *name, const char * centers, const char *words);


    static bool *findGEOTop(const char * projectionfile, const char * medianfile, const char * centersfile,
                            const char * wordsfile, const char * hesfile, const char *geosfile, const char *szfile, int TOPK);


};

#endif
