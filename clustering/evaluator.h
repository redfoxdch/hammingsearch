#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <vector>
#include <map>
#include <set>

using namespace std;


/**
        @author:chenghao deng
        @data:2016-11-18
        @version:1.0
**/

class Evaluator
{
public:
    Evaluator() {}
    virtual ~Evaluator() {}
    /** clustering evaluator **/

    static void loadClust(const char *srcFn, map<unsigned int, set<int> *> &clusts);//clust from zero;
    static float entropy(map<unsigned int, set<int> *> &clusts, map<unsigned int, set<int> *> &grdTrth);
    static float purity(map<unsigned int, set<int> *> &clusts,  map<unsigned int, set<int> *> &grdTrth);
    static float entropy(set<int> &clust, map<unsigned int, set<int> *> &grdTrth);
    static float purity(set<int> &clust, map<unsigned int, set<int> *> &grdTrth);
    static float entropy(const char *srcFn, const char *grdFn);
    static float purity(const char *srcFn, const char *grdFn);
    static int   interSect(const set<int> &s1, const set<int> &set2);
    static float funScore(const char *clustFn, const char *datFn, const char *fun);
    static void getSEP(const char *clustFn, const char *grdFn, const char *datFn, float &score, float &e, float &p, const char *fun);
    static void jaccard(map<unsigned int, set<int> *> &clusts, map<unsigned int, set<int> *> &grdTrh, const char *fn);
    static void jaccard(const char *setfn, const char* grdfn, const char *dstfn);
    static double distortion(const char* srcfn, const char* clusterfn);



    /** quantization evaluator **/
    static double recall(const char *nnfn, const char *grdfn);
    static double recall(const char *srcfn, const char *nnfn, unsigned top);

    /** image evaluator **/

    static void test();
};

#endif // EVALUATOR_H
