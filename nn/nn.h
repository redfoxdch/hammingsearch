#ifndef NN_H
#define NN_H

#include<vector>
using namespace std;

static inline unsigned InsertIntoKnn(vector<pair<float, unsigned> > &topmap, unsigned L, unsigned id, float dist, float &radius)
{
    if(topmap.size() < L && dist > radius)
    {
        topmap.push_back(pair<float, unsigned>(dist, id));
        radius = dist;
        return 0;
    }
    if(dist > radius)
    {
        return 0;
    }

    // find the location to insert
    unsigned j;
    unsigned K = topmap.size()-1;

    unsigned i = K;
    while (i > 0)
    {
        j = i - 1;
        if (topmap[j].first <= dist) break;
        i = j;
    }
    // check for equal ID
    unsigned l = i + 1;
    while (l > 0)
    {
        j = l - 1;
        if (topmap[j].first < dist) break;
        if (topmap[j].second == id) return K + 1;
        l = j;
    }
    // i <= K-1
    if(topmap.size() < L)
    {
        topmap.push_back(pair<float, unsigned>(topmap[K].first, topmap[K].second));
    }

    j = K;
    while (j > i)
    {
        topmap[j].first = topmap[j-1].first;
        topmap[j].second = topmap[j-1].second;
        --j;
    }
    topmap[i].first = dist;
    topmap[i].second = id;
    radius = topmap[topmap.size()-1].first;
    return 0;
}

static inline unsigned InsertIntoKnn (vector<unsigned> &topmap, unsigned id)
{

    unsigned i, j;
    unsigned K = topmap.size() - 1;
    unsigned l;

    if(topmap.size() == 0 || topmap[K] < id)
    {
        topmap.push_back(id);
        return 0;
    }
    if(topmap[K] == id)
        return 1;

    i = K;

    while (i > 0)
    {
        j = i - 1;
        if(topmap[j] <= id) break;
        i = j;
    }

    // check for equal ID
    l = i;
    while (l > 0)
    {
        j = l - 1;
        if(topmap[j] == id) return 1;
        l = j;
    }

    topmap.push_back(topmap[K]);

    j = K;
    while (j > i)
    {
        topmap[j] = topmap[j-1];
        --j;
    }
    topmap[i] = id;

    return 0;
}

static inline unsigned InsertIntoKnn (vector<pair<float, unsigned> > &topmap, unsigned id, float dist, float &radius)
{
// find the location to insert
    unsigned j;
    unsigned K = topmap.size()-1;

    unsigned i = K;
    while (i > 0)
    {
        j = i - 1;
        if (topmap[j].first <= dist) break;
        i = j;
    }
    // check for equal ID
    unsigned l = i;
    while (l > 0)
    {
        j = l - 1;
        if (topmap[j].first < dist) break;
        if (topmap[j].second == id) return K + 1;
        l = j;
    }
    // i <= K-1

    j = K;
    while (j > i)
    {
        topmap[j].first = topmap[j-1].first;
        topmap[j].second = topmap[j-1].second;
        --j;
    }
    topmap[i].first = dist;
    topmap[i].second = id;
    radius = topmap[K].first;
    return 0;
}

///recommand!
static inline unsigned InsertIntoKnn (unsigned *NNG, float *topDist, unsigned topK, float dist, unsigned id, float &radius)
{
    // find the location to insert
    unsigned j;
    unsigned K = topK;
    unsigned i = K;
    while (i > 0)
    {
        j = i - 1;
        if (NNG[j] == id) return K + 1;
        i = j;
    }
    // check for equal ID
    i = topK;
    while (i > 0)
    {
        j = i - 1;
        if (topDist[j] == radius) break;
        i = j;
    }
    // i <= K-1

    topDist[j] = dist;
    NNG[j] = id;

    i = K;
    radius = 0;
    while (i > 0)
    {
        j = i - 1;
        if (topDist[j] > radius)
            radius = topDist[j];
        i = j;
    }

    return i;
}


static inline unsigned InsertIntoKnn (unsigned *NNG, float *topDist, unsigned topK, unsigned id, float dist, float &radius)
{
    // find the location to insert
    unsigned j;
    unsigned K = topK-1;
    unsigned i = K;
    while (i > 0)
    {
        j = i - 1;
        if (topDist[j] <= dist) break;
        i = j;
    }
    // check for equal ID
    unsigned l = i;
    while (l > 0)
    {
        j = l - 1;
        if (topDist[j] < dist) break;
        if (NNG[j] == id) return K + 1;
        l = j;
    }
    // i <= K-1
    j = K;
    while (j > i)
    {
        topDist[j] = topDist[j-1];
        NNG[j] = NNG[j-1];
        --j;
    }
    topDist[i] = dist;
    NNG[i] = id;
    radius = topDist[K-1];
    return i;
}

class NN
{
    unsigned k;
    float *data;
    unsigned row;
    unsigned col;
    unsigned *topList;
    double recall;

public:
    NN();
    virtual ~NN();
    virtual void init() = 0;
    virtual void buildNN(const char *srcfn, const char *dstfn, unsigned &topk) = 0;
    virtual void buildNN(float *, unsigned &row, unsigned &col, unsigned &topk) = 0;

    ///NN file : N*K
    ///top1 top2 ...
    ///top1 top2 ...
    ///....
    bool saveNN(const char *dstfn);
    double getRecall(const char *grdfn);


    static double getRecall(const char *grdfn, const char *nnfn, unsigned topk);
    void getTopList(unsigned *topList);
    static void test();
};


#endif
