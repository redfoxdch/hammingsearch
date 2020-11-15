#ifndef ABSTRACTKMEANS_H
#define ABSTRACTKMEANS_H

#include "util/sparsematrix.h"
#include "util/clsinfo.h"
#include "util/nnitem.h"

#include <unordered_set>
#include <vector>
#include <map>

struct CLSIndex
{
public:
    CLSIndex()
    {
        index = 0;
        val   = 0.0;
    }
    CLSIndex(const unsigned int idx, const float val0)
    {
        index = idx;
        val   = val0;
    }
    unsigned int index;  /// number of members
    float val;           /// sum of intra-sim

    static bool LtComp (const CLSIndex *a,const CLSIndex *b)
    {
        return (a->val < b->val);
    }
    static bool LgComp (const CLSIndex *a,const CLSIndex *b)
    {
        return (a->val < b->val);
    }
};


enum KMN_OPT {_tkmn_ = 0, _kppkmn_ = 1, _xbkmn_ = 2, _xtkmn_ = 3, _mnkmn_ = 4, _olkmn_ = 5, _rbkmn_, _rrkmn_, _lvq_, _hek_, _imk_, _hec_, _prc_, _pqc_};
enum OPTZ {_I1_ = 0, _I2_ = 1, _I3_, _I4_, _E1_, _E2_,_T1_};
enum SEED {_rnd_ = 0, _kpp_ = 1, _non_};

class AbstractKMeans
{
protected:
    AbstractKMeans();

protected:
    static const unsigned int paraBound;
    unsigned long count, ndim, clnumb;
    unsigned int nITER;
    unsigned int nThrd;
    vector<NNItem*>  sorted_stack;
    static const float smallVal0;
    double  *lens, allEg;
    bool    _INIT_,  _REFER_, _REFINE_;
    char    srcmatfn[1024];
    char    mthStr[8];
    unsigned long mvs[2];
    double  dstort;
    SEED    seed;
    float   *data;

    SparseMatrix sdata;

    CLSInfo *infoMap;
    int     *labels;
    KMN_OPT kmMtd;
    bool dataType;
    bool SHUFFLE;
    OPTZ  myoptz;
public:
    unsigned int buildcluster(const char *srcfn, const char *dstfn,  const char *_seed_, const char *lg_first,
                              const char *crtrn, const int num, bool _refine_);
    unsigned int buildcluster(float *mat, const int row, const int dim, const char *dstfn,
                              const char *_seed_, const char  *lg_first, const char *crtrn, const int num, bool _refine_);

    bool  initMemry(const unsigned int dNum, const unsigned int clustNum);
    virtual bool  init(const char *srcfn) = 0;
    virtual bool  init(float *data, const int row, const int dim) = 0;
    virtual bool  config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose) = 0;
    virtual int   clust(const unsigned int clust_num, const char *dstfn, const int verbose) = 0;
    double  getI2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns);

    double  calAVGDist(const double *Ds, const int *Ns, const unsigned int clustNum, CLSInfo *infos);
    double  calAVGDist(const double *Ds, const int *Ns, const unsigned int clustNum);
    float  *getCluster(const unsigned int clabel0, unsigned int &row, unsigned int &dim);

    virtual void saveCenters(const char *dstfn, bool append) = 0;
    void    saveCenters(const char *dstfn, double *arrayD, bool append, bool norm);
    virtual int  fetchCenters(float *centers) = 0;

    int    idx2clabel(const int i);

    void   printClusters(const char *dstdir);
    void   printCluster(const char  *dstfn);
    bool   save_clust(const char    *dstfn);
    bool label2cluster(vector<unsigned> *clusters);
    bool getLabel(int *labels);

    static double normVects(float *vects, double *lens, const unsigned int d0, const unsigned int n0);
    static void normVects(const SparseMatrix &svects, double *lens, const unsigned int d0, const unsigned int n0);
    static void normVects(float *vects, const unsigned int d0, const unsigned int n0, double *lens);
    static void normVects(SparseMatrix &vects, const unsigned int d0, const unsigned int n0, double *lens);

    void resetlabels(unsigned int &clustnum);

    virtual ~AbstractKMeans();

};



#endif
