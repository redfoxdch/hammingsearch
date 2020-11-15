#ifndef IOAGENT_H
#define IOAGENT_H

#include "sparsematrix.h"

#include <string>
#include <vector>
#include <map>
#include <set>

using std::set;
using std::map;
using std::vector;
using std::string;

class IOAgent
{
public:
    IOAgent() {}
    virtual ~IOAgent() {}
    static float *load_fvecs(const char *fvecfn, unsigned int &r, unsigned int &d);
    static float *load_fvecs(vector<string> &filenames, unsigned int &r, unsigned int &d);
    static float *load_bvecs(const char *bvecfn, unsigned int &r, unsigned int &d);
    static float *load_bvecs(const char *bvecfn, unsigned int &r, unsigned int &d, unsigned row);

    static float *loadMatrix(const char *fn,     unsigned int &row, unsigned int &col);
    static float *loadMatrix(vector<string> &filenames,     unsigned int &row, unsigned int &col);
    static float *loadDat(const char *fn,     unsigned int &row, unsigned int &col);
    static float *loadDat(vector<string> &filenames,     unsigned int &row, unsigned int &col);
    static unsigned char *loadPQ(const char *srcfn, unsigned int &row, unsigned int&col);
    static float* loadPQInfo(const char *srcfn, unsigned int &pqn, unsigned int &pql, unsigned &pqm);
    static unsigned *loadNN(const char *nnfn, unsigned &row, unsigned col);

    static float *loadItms(const char *fn,  const char *idxKey,  unsigned int &row, unsigned int &col);
    static float *loadItms(const char *fn,  const char *idxKey,  const unsigned int line, unsigned int &row, unsigned int &col);
    static SparseMatrix loadBOVW(const char *fn, unsigned int &row, unsigned int &col);
    static SparseMatrix loadSparse(const char *fn, unsigned int &row, unsigned int &col);

    static unsigned int getMatrixInfo(const char *matfn, unsigned int &col, unsigned int &row);

    static void saveDat(const char *fn, const unsigned int &row, const unsigned int &col,float *Dat);
    static void saveDat(const char* fn, const unsigned int& row, const unsigned int &col, unsigned char *Dat);
    static void saveSparse(const char* fn, const unsigned int& row, const unsigned int &col, const SparseMatrix &sdata);
    static float *loadTruncMat(const char *srcFn, unsigned int &row, unsigned int &dim, vector<unsigned int> &ids);
    static void saveTruncMat(const char *dstFn, const unsigned int dim, float *mat, set<unsigned int> &ids);
    bool static readFileName(const char* srcdir, vector<string >& filename);
    static void saveEx(vector<string> &filenames, float *e, float *p, const char *result);
    static void addRslt(const char *filename, float *e, float *p, const char *result);
    static void saveAsFullMat(const float *dat, unsigned int &row, unsigned int &col, const char *file);
    static void saveknnGraph(const unsigned *knnGraph, unsigned int row, unsigned int dim,
                             map<unsigned int, unsigned int> &idMap, const char *dstFn);

    static void loadClust(const char *fn, map<unsigned int, set<int> *> &cluster);

    static void loadmClust(const char *srcFn, map<unsigned int, set<int> *> &clusts);
    static void saveClust(const char* fn, map<unsigned int, set<int> *>&cluster);

    static void test();
};

#endif // IOAGENT_H
