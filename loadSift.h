#ifndef LOADSIFT_H
#define LOADSIFT_H

#include <vector>
#include <string>

/***

@function: load SIFT feature from npy file
@author :Cheng-Hao Deng
@version : 1.0
@date: 2018-7-19
@institute: akulaku
/***/

using namespace std;

class loadSift
{


public :
    loadSift();
    static  int scanFiles(string inputDirectory, vector<string> &fileList);
    static  vector<string> loadFileName(const char *fileName);
    static  float *loadFromNameList(const char * filename, long &row, int &col,int MAX_LEN = RAND_MAX);
    static  float *loadFromFile(const char * file, long &row, int &col);
    static  float *loadFromDir(const char * dir, long &row, int &col);
    static  float *loadGeo(const char * dir, long &row, int &col);
    static void saveGeo(const char *file,float *geo, unsigned &row);

    static  float *loadGeoDir(const char * dir, long &row, int &col);
    static  float *loadGeoList(const char * dir, long &row, int &col);
    static  float *loadProjection(const char *file, int &row,int &col);
    static  float *loadCenters(const char *file, int &row,int &col);
    static  double *loadMedian(const char *file, int &row,int &col);
    static  unsigned *loadWords(const char *file, long &row);
    static  vector<unsigned>loadShape(const char *file, long &row);

    void static test();
    virtual ~loadSift();
};








#endif // LOADSIFT_H
