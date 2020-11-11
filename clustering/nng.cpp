#include "randompartition.h"
#include "util/distance.hpp"
#include "ioagent.h"
#include "vstring.h"
#include "timer.h"
#include "nng.h"
#include "nn/nn.h"

#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <map>


#define SIFT

static unsigned RangeFirst = 100;
static unsigned RangeSecond = 100;
static unsigned topRange = 50;


bool static inline fuLs(const pair<float, unsigned> &a, const pair<float, unsigned> &b)
{
    if(a.first == b.first)
        return a.second < b.second;
    return a.first < b.first;
}

NNG::NNG(const DISTANCE::Distance<float>* d):distance(d)
{
    this->EXTData = false;
    this->data = NULL;
    this->topList = NULL;
    this->topDist = NULL;
    this->count = 0;
    this->ndim = 0;
    this->NNIter = 6;
    this->Check = 0;
    this->HubNum = 30;
    this->PoolSize = 50;
    this->MaxSize = 30;
    this->TopK = 50;
    this->TreeNum = 8;
    this->TreeSize = 50;

    this->labels = NULL;
    this->radius = NULL;
}

NNG::~NNG()
{
    if(!this->EXTData)
    {
        delete [] this->data;
    }
    delete [] this->topList;
    delete [] this->topDist;
    if(this->labels != NULL)
    {
        delete [] labels;
        this->labels = NULL;
    }
}

int NNG::initMemory()
{
    this->topList = new unsigned[this->count*TopK];
    memset(this->topList, 0, sizeof(unsigned)*this->count*TopK);
    this->topDist = new float[this->count*TopK];
    memset(this->topDist, 0, sizeof(float)*this->count*TopK);
    this->radius = new float[this->count];
}

int NNG::getTopDistance()
{
    if(this->topDist = NULL)
    {
        topDist = new float[this->count];
    }
}

int NNG::adjustNN()
{
    unsigned i, j;
    vector<pair<float, unsigned> > topmap;
    topmap.resize(TopK);
    for(i = 0; i < this->count; i++)
    {
        for(j = 0; j < this->TopK; j++)
        {
            topmap[j].second = topList[i*TopK+j] ;
            topmap[j].first = topDist[i*TopK+j];
        }
        sort(topmap.begin(), topmap.end(), fuLs);
        for(j = 0; j < this->TopK; j++)
        {
            topList[i*TopK+j] = topmap[j].second;
            topDist[i*TopK+j] = topmap[j].first;
        }
    }

}

bool NNG::loadData(const char *srcfn)
{
    cout<<"Loading matrix ................... ";
    assert(srcfn);

    if(VString::endWith(srcfn, ".txt"))
    {
        this->data = IOAgent::loadMatrix(srcfn, this->count, this->ndim);
    }
    else if(VString::endWith(srcfn, ".fvecs"))
    {
        this->data = IOAgent::load_fvecs(srcfn, this->count, this->ndim);
    }
    else if(VString::endWith(srcfn, ".mat"))
    {
        this->data = IOAgent::loadDat(srcfn, this->count, this->ndim);
    }
    else
    {
        this->data = IOAgent::loadItms(srcfn, "fsvtab", 10000000, this->count, this->ndim);
    }

    if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }

    cout<<this->count<<"x"<<this->ndim<<endl;

    initMemory();

    cout<<"TopK\t"<<this->TopK<<endl;
    return true;
}

bool NNG::loadData(unsigned char *mat, const int row, const int dim)
{

}

bool NNG::loadData(float *mat, const int row, const int dim)
{

}

int NNG::buildTreeGraph()
{
    if(this->topDist == NULL)
    {
        this->getTopDistance();
    }

    unsigned i;
    this->clnumb = this->count/TreeSize;
    this->clnumb = pow(2, floor(log(this->clnumb)/log(2)));
    //cout<<this->clnumb<<endl;

    for(int iter = 0; iter < this->TreeNum; iter++)
    {
        cout<<"initialize cluster begin .......\n";
        if(this->labels == NULL)
        {
            this->labels = new int[this->count];
        }
        RandomPartition *xbk = new RandomPartition(this->distance);
        xbk->buildcluster(this->data, this->count, this->ndim, "tmpfile.txt",  "non", "large", "i2", this->clnumb, 1);
        xbk->getLabel(this->labels);
        delete xbk;
        cout<<"Initialize cluster finished\n";
        Timer *mt = new Timer;
        mt->start();
        this->treeUp();
        mt->end(true);
    }
    if(this->labels != NULL)
    {
        delete [] labels;
        labels = NULL;
    }
}

int NNG::initNNGraph()
{
    /**/
    unsigned i, j, k, id, rid, r1 = 0, r2 = 0;
    float dst = 0, maxdst;
    unsigned *rl = new unsigned[this->count];
    for(i = 0; i < this->count; i++)
    {
        rl[i] = i;
    }
    random_shuffle(rl, rl+this->count);

    vector<pair<float, unsigned> > topmap;
    topmap.resize(TopK);

    for(i = 0; i < this->count; i++)
    {
        rid = rand()%this->count;
        maxdst = 0;
        for(j = 0; j < TopK; j++)
        {
            id = rl[(j+rid)%this->count];
            if(id == i)
            {
                id = rl[(rid+TopK+1)%this->count];
            }

            topList[i*TopK+j] = id;
            r1  = i*this->ndim;
            r2  = id*this->ndim;
            dst = this->distance->compare(data+r1, data+r2, ndim);
            topDist[i*TopK+j] = dst;
            if(dst > maxdst)
                maxdst = dst;
        }
        radius[i] = maxdst;
    }
    this->adjustNN();
    /**/

    /**
    unsigned i, j, k, id, rid, r1 = 0, r2 = 0;
    float dst = 0;
    unsigned *rl = new unsigned[this->count];
    for(i = 0; i < this->count; i++)
    {
        rl[i] = i;
    }
    random_shuffle(rl, rl+this->count);


    for(i = 0; i < this->count; i++)
    {
        rid = rand()%this->count;
        for(j = 0; j < TopK; j++)
        {
            id = rl[(j+rid)%this->count];
            if(id == i)
            {
                id = rl[(rid+TopK+1)%this->count];
            }
            topList[i*TopK+j] = id;
            topDist[i*TopK+j] = RAND_MAX;
        }
    }
    /**/
}

int NNG::label2clusters(vector<unsigned>* clusters)
{
    int i;
    for(i = 0; i < this->clnumb; i++)
    {
        clusters[i].clear();
    }
    for(i = 0; i < this->count; i++)
    {
        clusters[this->labels[i]].push_back(i);
    }

    return 0;
}

int NNG::treeUp()
{
    unsigned sz = 0;
    vector<unsigned> *clusters = new vector<unsigned>[this->clnumb];
    unsigned i = 0, j = 0, k = 0, id, r, col, r1 = 0, r2 = 0;
    float dst = 0;
    this->label2clusters(clusters);
    cout<<"update tree\n";

    unsigned mxsize = 0;
    for(i = 0; i < this->clnumb; i++)
    {
        if(clusters[i].size() > mxsize)
            mxsize = clusters[i].size();
    }

    float *cdst = new float[mxsize*mxsize];

    for(i = 0; i < this->clnumb; i++)
    {
        sz     = clusters[i].size();
        for(j = 0; j < sz; j++)
        {
            id = clusters[i][j];
            for(k = j + 1; k < sz; k++)
            {
                r = clusters[i][k];
                dst = 0;
                r1  = r*ndim;
                r2  = id*ndim;
                dst = this->distance->compare(data+r1, data+r2, ndim);
                cdst[j*sz+k] = dst;
            }
        }
        for(j = 0; j < sz; j++)
        {
            id = clusters[i][j];
            for(k = j + 1; k < sz; k++)
            {
                r = clusters[i][k];
                dst = cdst[j*sz+k];
                if(dst <= radius[id])
                    InsertIntoKnn(topList+id*TopK, topDist+id*TopK, TopK, dst, r, radius[id]);
                if(dst <= radius[r])
                    InsertIntoKnn(topList+r*TopK, topDist+r*TopK, TopK, dst, id, radius[r]);
            }
        }

        clusters[i].clear();
    }
    this->adjustNN();

    delete [] cdst;
    cdst = NULL;
    delete [] clusters;
    clusters = NULL;
    cout<<"update tree finished\n";
    return 0;
}

int NNG::NNExtent()
{
    if(this->topDist == NULL)
    {
        this->getTopDistance();
    }
}

int NNG::NNDescent()
{
    if(this->topDist == NULL)
    {
        this->getTopDistance();
    }

    unsigned i, j, k, id, r, all, d;
    unsigned nwn;
    unsigned odn;
    unsigned dstSize = RangeSecond*(RangeSecond+RangeFirst);
    vector<unsigned> *nw = new vector<unsigned>[this->count];
    vector<unsigned> *od = new vector<unsigned>[this->count];
    float tmpdst, rdst;
    float *ndst = new float[dstSize];
    float *radius = new float[this->count];
    float *rradius = new float[this->count];
    bool *state = new bool[this->count];
    memset(state, 0, sizeof(bool)*this->count);
    memset(rradius, 0, sizeof(float)*this->count);

    vector<unsigned> *rIDMap = new vector<unsigned>[this->count];


    //set<unsigned> tmpset;
    unordered_set<unsigned> tmpset, rtmpset;
    vector<unsigned> tmpseq, oldseq, nseq, rseq;

    if(Check > TopK)
        d = TopK;
    else
        d = Check;
    if(MaxSize == 0)
        d = 0;
    //d = 15;
    for(i = 0; i < this->count; i++)
    {
        if(d > 0)
        {
            nwn = d;
            nw[i].resize(d);
            for(j = 0; j < d; j++)
            {
                nw[i][j] = topList[i*TopK+j];
            }
            rradius[i] = topDist[i*TopK+d-1];
        }
        radius[i] = topDist[i*TopK+TopK-1];
    }

    for(i = 0; i < this->count; i++)
    {
        k = 0;
        for(j = 0; j < TopK; j++)
        {
            if(k == HubNum)
                break;
            id = topList[i*TopK+j];
            nwn = nw[id].size();
            tmpdst = topDist[i*TopK+j];
            rdst = rradius[i];
            if(tmpdst > rdst)
            {
                rIDMap[id].push_back(i);
                k++;
            }
        }
    }

    /**/
    ///hub points
    for(i = 0; i < this->count; i++)
    {
        tmpseq.clear();
        for(auto it = rIDMap[i].begin(); it != rIDMap[i].end(); it++)
        {
            tmpseq.push_back(*it);
        }
        rIDMap[i].clear();
        random_shuffle(tmpseq.begin(), tmpseq.end());
        if(tmpseq.size() > PoolSize)
            tmpseq.resize(PoolSize);
        rtmpset.clear();
        nwn = nw[i].size();
        for(j = 0; j < nwn; j++)
        {
            id = nw[i][j];
            state[id] = true;
        }
        for(auto it = tmpseq.begin(); it != tmpseq.end(); it++)
        {
            if(!state[id])
            {
                nw[i].push_back(*it);
            }
        }
        for(j = 0; j < nwn; j++)
        {
            id = nw[i][j];
            state[id] = false;
        }
    }
    /**/

    for(unsigned iter = 0; iter < NNIter; iter++)
    {
        cout<<"iter bgin\n";
        Timer *mytm = new Timer;
        mytm->start();
        for(i = 0; i < this->count; i++)
        {
            ///initialize
            nwn = nw[i].size();
            odn = od[i].size();
            all = (nwn + odn);
            ///caculate distance
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                for(k = j+1; k < nwn; k++)
                {
                    r = nw[i][k];
                    tmpdst = this->distance->compare(data+r*ndim, data+id*ndim, ndim);
                    ndst[j*all+k] = tmpdst;
                }
            }
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                for(k = 0; k < odn; k++)
                {
                    r = od[i][k];
                    tmpdst = this->distance->compare(data+r*ndim, data+id*ndim, ndim);
                    ndst[j*all+nwn+k] = tmpdst;
                }
            }
            /**/
            /**upnn**/
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                for(k = j+1; k < nwn; k++)
                {
                    r = nw[i][k];
                    tmpdst = ndst[j*all+k];
                    if(tmpdst < radius[id])
                    {
                        InsertIntoKnn(topList+id*TopK, topDist+id*TopK, TopK, tmpdst, r, radius[id]);
                    }
                    if(tmpdst < radius[r])
                    {
                        InsertIntoKnn(topList+r*TopK, topDist+r*TopK, TopK, tmpdst, id, radius[r]);
                    }
                }
            }
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                for(k = 0; k < odn; k++)
                {
                    r = od[i][k];
                    tmpdst = ndst[j*all+nwn+k];
                    if(tmpdst < radius[id])
                    {
                        InsertIntoKnn(topList+id*TopK, topDist+id*TopK, TopK, tmpdst, r, radius[id]);
                    }
                    if(tmpdst < radius[r])
                    {
                        InsertIntoKnn(topList+r*TopK, topDist+r*TopK, TopK, tmpdst, id, radius[r]);
                    }
                }
            }
        }
        this->adjustNN();
        ///hubpoints
        ///danymic check TopK num
        /**/
        for(i = 0; i < this->count; i++)
        {
            nwn = nw[i].size();
            odn = od[i].size();
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                state[id] = true;
            }
            for(j = 0; j < odn; j++)
            {
                id = od[i][j];
                state[id] = true;
            }
            r = 0;
            k = 0;
            for(j = 0; j < TopK; j++)
            {
                if(r >= Check || r+k >= MaxSize)
                    break;
                id = topList[i*TopK+j];
                if(!state[id])
                {
                    r++;
                }
                else
                {
                    k++;
                }
            }
            if(r+k > 0)
                rradius[i] = topDist[r+k-1];
            else
                rradius[i] = 0;
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                state[id] = false;
            }
            for(j = 0; j < odn; j++)
            {
                id = od[i][j];
                state[id] = false;
            }
        }
        for(i = 0; i < this->count; i++)
        {
            k = 0;
            for(j = 0; j < TopK; j++)
            {
                if(k == HubNum)
                    break;
                id = topList[i*TopK+j];
                tmpdst = topList[i*TopK+j];
                rdst = rradius[id];
                if(tmpdst > rdst)
                {
                    rIDMap[id].push_back(i);
                    k++;
                }
            }
        }
        /**/
        ///update nw od;
        /**/
        for(i = 0; i < this->count; i++)
        {
            nwn = nw[i].size();
            odn = od[i].size();
            tmpseq.clear();
            tmpseq.reserve(PoolSize);
            for(auto it = rIDMap[i].begin(); it != rIDMap[i].end(); it++)
            {
                tmpseq.push_back(*it);
            }
            rIDMap[i].clear();
            random_shuffle(tmpseq.begin(), tmpseq.end());
            if(tmpseq.size() > PoolSize)
                tmpseq.resize(PoolSize);
            oldseq.clear();
            oldseq.reserve(nwn+odn);
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                oldseq.push_back(id);
                state[id] = true;
            }
            for(j = 0; j < odn; j++)
            {
                id = od[i][j];
                oldseq.push_back(id);
                state[id] = true;
            }
            nw[i].clear();
            od[i].clear();
            r = 0;
            k = 0;
            nseq.clear();
            nseq.reserve(RangeFirst);
            for(j = 0; j < TopK; j++)
            {
                if(r >= Check || r+k >= MaxSize)
                    break;
                id = topList[i*TopK+j];
                if(!state[id])
                {
                    nw[i].push_back(id);
                    r++;
                }
                else
                {
                    od[i].push_back(id);
                    k++;
                }
                nseq.push_back(id);
            }
            for(auto it = oldseq.begin(); it != oldseq.end(); it++)
            {
                id = *it;
                state[id] = false;
            }
            for(auto it = nseq.begin(); it != nseq.end(); it++)
            {
                id = *it;
                state[id] = true;
            }
            rseq.clear();
            rseq.reserve(tmpseq.size());
            for(auto it = tmpseq.begin(); it != tmpseq.end(); it++)
            {
                id = *it;
                if(!state[id])
                {
                    rseq.push_back(id);
                }
            }
            for(auto it = nseq.begin(); it != nseq.end(); it++)
            {
                id = *it;
                state[id] = false;
            }
            for(auto it = oldseq.begin(); it != oldseq.end(); it++)
            {
                id = *it;
                state[id] = true;
            }
            for(auto it = rseq.begin(); it != rseq.end(); it++)
            {
                id = *it;
                if(!state[id])
                {
                    nw[i].push_back(id);
                    r++;
                }
                else
                {
                    od[i].push_back(id);
                    k++;
                }
            }
            for(auto it = oldseq.begin(); it != oldseq.end(); it++)
            {
                id = *it;
                state[id] = false;
            }
        }

        /**/
        cout<<"iter "<<iter<<endl;
        mytm->end(true);
        /**/
    }

    cout<<"end\n";
///clear top K;
    for(i = 0; i < this->count; i++)
    {
        nw[i].clear();
        od[i].clear();
    }
    delete [] rIDMap;
    delete [] ndst;
    delete [] radius;
    delete [] rradius;
    delete [] state;
    delete [] nw;
    delete [] od;

    return 0;
}

int NNG::NNHub()
{
    if(this->topDist == NULL)
    {
        this->getTopDistance();
    }

    unsigned i, j, k, id, r, all, d;
    unsigned nwn;
    unsigned odn;
    unsigned dstSize = RangeSecond*(RangeSecond+RangeFirst);
    vector<unsigned> *nw = new vector<unsigned>[this->count];
    vector<unsigned> *od = new vector<unsigned>[this->count];
    vector<pair<float, unsigned> >::iterator it;
    float tmpdst, rdst;
    float *ndst = new float[dstSize];
    bool *state = new bool[this->count];
    memset(state, 0, sizeof(bool)*this->count);


    vector<unsigned> *rIDMap = new vector<unsigned>[this->count];
    // unsigned *oldNN = new unsigned[this->count*this->HubNum];

    unordered_set<unsigned> tmpset, rtmpset;
    vector<unsigned> tmpseq, oldseq;


    for(i = 0; i < this->count; i++)
    {
        k = 0;
        for(j = 0; j < TopK; j++)
        {
            if(k == Check)
                break;
            id = topList[i*TopK+j];
            rIDMap[id].push_back(i);
            k++;
        }
    }

    /**/
    ///hub points
    for(i = 0; i < this->count; i++)
    {
        tmpseq.clear();
        for(auto it = rIDMap[i].begin(); it != rIDMap[i].end(); it++)
        {
            tmpseq.push_back(*it);
        }
        rIDMap[i].clear();
        random_shuffle(tmpseq.begin(), tmpseq.end());
        if(tmpseq.size() > PoolSize)
            tmpseq.resize(PoolSize);
        rtmpset.clear();
        for(auto it = tmpseq.begin(); it != tmpseq.end(); it++)
        {
            nw[i].push_back(*it);
        }
    }
    /**/


    for(unsigned iter = 0; iter < NNIter; iter++)
    {
        cout<<"iter bgin\n";
        Timer *mytm = new Timer;
        mytm->start();

        for(i = 0; i < this->count; i++)
        {
            ///initialize
            nwn = nw[i].size();
            odn = od[i].size();
            all = (nwn + odn);

///caculate distance
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                for(k = j+1; k < nwn; k++)
                {
                    r = nw[i][k];
                    tmpdst = this->distance->compare(data+r*ndim, data+id*ndim, ndim);
                    ndst[j*all+k] = tmpdst;
                }
            }
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                for(k = 0; k < odn; k++)
                {
                    r = od[i][k];
                    tmpdst = this->distance->compare(data+r*ndim, data+id*ndim, ndim);
                    ndst[j*all+nwn+k] = tmpdst;
                }
            }
            /**/

            /**upnn**/
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                for(k = j+1; k < nwn; k++)
                {
                    r = nw[i][k];

                    tmpdst = ndst[j*all+k];
                    if(tmpdst < radius[id])
                    {
                        InsertIntoKnn(topList+id*TopK, topDist+id*TopK, TopK, tmpdst, r, radius[id]);
                    }

                    if(tmpdst < radius[r])
                    {
                        InsertIntoKnn(topList+r*TopK, topDist+r*TopK, TopK, tmpdst, id,  radius[r]);
                    }

                }
            }
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                for(k = 0; k < odn; k++)
                {
                    r = od[i][k];
                    tmpdst = ndst[j*all+nwn+k];
                    if(tmpdst < radius[id])
                    {
                        InsertIntoKnn(topList+id*TopK, topDist+id*TopK, TopK, tmpdst, r,  radius[id]);
                    }
                    if(tmpdst < radius[r])
                    {
                        InsertIntoKnn(topList+r*TopK, topDist+r*TopK, TopK, tmpdst, id,  radius[r]);
                    }
                }
            }
        }

        ///hubpoints
        if(this->Check < TopK)
        {
            this->adjustNN();
        }

        for(i = 0; i < this->count; i++)
        {
            k = 0;
            for(j = 0; j < TopK; j++)
            {
                if(k == Check)
                    break;
                id = topList[i*TopK+j];
                rIDMap[id].push_back(i);
                k++;
            }
        }
        Check += MaxSize;
        if(Check > HubNum)
            Check = HubNum;

        /**/

        /**/
        ///update nw od;
        /**/
        for(i = 0; i < this->count; i++)
        {
            nwn = nw[i].size();
            odn = od[i].size();
            all = nwn+odn;
            tmpseq.clear();
            tmpseq.reserve(PoolSize);
            for(auto it = rIDMap[i].begin(); it != rIDMap[i].end(); it++)
            {
                tmpseq.push_back(*it);
            }
            rIDMap[i].clear();
            random_shuffle(tmpseq.begin(), tmpseq.end());

            if(tmpseq.size() > PoolSize)
                tmpseq.resize(PoolSize);

            oldseq.clear();
            oldseq.reserve(nwn+odn);
            for(j = 0; j < nwn; j++)
            {
                id = nw[i][j];
                oldseq.push_back(id);
                state[id] = true;
            }
            for(j = 0; j < odn; j++)
            {
                id = od[i][j];
                oldseq.push_back(id);
                state[id] = true;
            }
            nw[i].clear();
            od[i].clear();

            r = 0;
            k = 0;
            for(auto it = tmpseq.begin(); it != tmpseq.end(); it++)
            {
                id = *it;
                if(!state[id])
                {
                    nw[i].push_back(id);
                    r++;
                }
                else
                {
                    od[i].push_back(id);
                    k++;
                }
            }
            for(auto it = oldseq.begin(); it != oldseq.end(); it++)
            {
                id = *it;
                state[id] = false;
            }
        }
        /**/
        cout<<"iter "<<iter<<endl;
        mytm->end(true);
        /**/
    }

    cout<<"end\n";
///clear top K;
    for(i = 0; i < this->count; i++)
    {
        nw[i].clear();
        od[i].clear();
    }
    this->adjustNN();

    delete [] rIDMap;
    delete [] ndst;

    delete [] state;
    delete [] nw;
    delete [] od;
    //delete [] oldNN;

    return 0;

}

int NNG::Extent()
{
    unsigned i, j , k, id, d, num, sz;
    double dst;
    multimap<double, unsigned> *dmap = new multimap<double, unsigned>[this->count];
    set<unsigned> tops;
    set<unsigned> *visitList = new set<unsigned>[this->count];

    for(i = 0; i < this->count; i++)
    {
        for(j = 0; j < this->TopK; j++)
        {
            id = this->topList[i*TopK+j];
            visitList[i].insert(id);
        }
    }

    for(i = 0; i < this->count; i++)
    {
        tops.clear();
        num = 0;
        for(j = 0; j < MaxSize; j++)
        {
            for(k = 0; k < Check; k++)
            {
                id = topList[topList[i*TopK+j]*TopK+k];
                if(id == i)
                    continue;
                sz = visitList[i].size();
                visitList[i].insert(id);
                if(sz == visitList[i].size())
                    continue;
                tops.insert(id);
                num++;
                if(num >= MaxSize)
                    break;
            }
            if(num >= MaxSize)
                break;
        }

        for(auto it = tops.begin(); it != tops.end(); it++)
        {
            id = *it;

            dst = 0;
            for(d = 0; d < ndim; d++)
            {
                dst += (data[i*ndim+d] - data[id*ndim+d])*(data[i*ndim+d] - data[id*ndim+d]);
            }
            dmap[i].insert(pair<double, unsigned>(dst, id));
            if(i < id && 0)
            {
                sz = visitList[id].size();
                visitList[id].insert(i);
                if(sz != visitList[id].size())
                {
                    dmap[id].insert(pair<double, unsigned>(dst, i));
                }
            }
        }

        for(j = 0; j < TopK; j++)
        {
            dst = topDist[i*TopK+j];
            id = topList[i*TopK+j];
            dmap[i].insert(pair<double, unsigned>(dst, id));
        }
        j = 0;
        for(auto it = dmap[i].begin(); it != dmap[i].end(); it++)
        {
            topDist[i*TopK+j] = it->first;
            topList[i*TopK+j] = it->second;
            j++;
            if(j == TopK)
                break;
        }
        dmap[i].clear();
    }

    for(i = 0; i < this->count; i++)
    {
        visitList[i].clear();
    }

    delete [] visitList;
    delete [] dmap;
    return 0;
}

int NNG::initNNGraph(unsigned *topList)
{
    memcpy(this->topList, topList, sizeof(unsigned)*this->count*TopK);
    this->getTopDistance();
}

int NNG::NNBuild(unsigned type)
{
    Timer *alt = new Timer;
    alt->start();
    Timer *tt = new Timer;
    tt->start();
    cout<<"init topList\n";
    this->initNNGraph();
    cout<<"init finished\n";
    cout<<"build tree\n";
    this->buildTreeGraph();
    cout<<"build tree finished\nbuild nnhub\n";
    tt->end(true);
    //this->showRecall();
    Timer *nt = new Timer;
    nt->start();
    if(type == 0)
        this->NNHub();
    else if(type == 1)
        this->NNDescent();
    nt->end(true);
    cout<<"all cost time:\n";
    alt->end(true);
    // this->showRecall();
    return 0;
}

int NNG::NNBuild(const char *srcfn, const char* dstfn, unsigned TreeNum, unsigned TreeSize, unsigned HubNum, unsigned PoolSize, unsigned NNIter, unsigned TopK, unsigned MaxSize)
{
    this->TreeNum = TreeNum;
    this->TreeSize = TreeSize;
    this->HubNum = HubNum;
    this->PoolSize = PoolSize;
    this->NNIter = NNIter;
    this->TopK = TopK;
    this->MaxSize = MaxSize;

    this->loadData(srcfn);
    this->NNBuild(0);
    this->saveNNGraph(dstfn);
    this->showRecall();
}

int NNG::showRecall()
{
    int i = 0, j = 0, id = 0;
    /**/
    unsigned tmpid[200];
    unsigned *rtop = new unsigned[200*this->count];
    ifstream is;
    if(this->count == 1000000 && this->ndim == 128)
        is.open("sift_base200nn.txt");
    else if(this->count == 100000 && this->ndim == 128)
        is.open("nn.txt");
    else if(this->count >= 5000 && this->count <= 6000 && this->ndim == 2048)
        is.open("ox5k200nn.txt");
    else if(this->count == 50000 && this->ndim == 2048)
        is.open("ox50k200nn.txt");
    else if(this->count >= 100000 && this->count <= 100000+1000 && this->ndim == 2048)
        is.open("ox100k200nn.txt");
    else if(this->count == 500000 && this->ndim == 960)
        is.open("gist_learn200nn.txt");
    else if(this->count == 1000000 && this->ndim == 960)
        is.open("gist_200nn.txt");
    else
    {
        cout<<"erro no ground truth for count "<<this->count<<endl;
        return -1;
    }

    if(!is.is_open())
    {
        cout<<"cant not open groundturth \n";
        return -1;
    }

    for(i = 0; i < this->count; i++)
    {
        for(j = 0; j < 200; j++)
        {
            is>>id;
            tmpid[j] = id;
        }
        for(j = 0; j < 200; j++)
            rtop[i*200+j] = tmpid[j];
    }
    is.close();
    unsigned r[4] = {0};
    unordered_set<unsigned> tmpset;
    unsigned b[4] = {1,3,5,10};
    unsigned a = 0;
    for(i = 0; i < this->count; i++)
    {
        for(id = 0; id < 4; id++)
        {
            a = b[id];
            tmpset.clear();
            for(j = 0; j < TopK; j++)
            {
                if(j >=a)
                    break;
                tmpset.insert(topList[i*TopK+j]);
            }
            for(j = 0; j < 1; j++)
            {
                tmpset.insert(rtop[i*200+j]);
            }
            if(tmpset.size() < 1+a)
            {
                r[id]++;
            }
            if(a == 0)
                cout<<a<<endl;
        }
    }
    unsigned cr[5] = {0};
    unsigned c[5] = {5,10,25,50,100};
    for(i = 0; i < this->count; i++)
    {
        for(id = 0; id < 5; id++)
        {
            a = c[id];
            tmpset.clear();
            for(j = 0; j < TopK; j++)
            {
                if(j >=a)
                    break;
                tmpset.insert(topList[i*TopK+j]);
            }
            for(j = 0; j < a; j++)
            {
                tmpset.insert(rtop[i*200+j]);
            }
            cr[id] += 2*a - tmpset.size();
        }
    }
    double ctop5 = (cr[0]+0.0)/this->count/5;
    double ctop10 = (cr[1]+0.0)/this->count/10;//10;
    double ctop15 = (cr[2]+0.0)/this->count/25;//50;
    double ctop20 = (cr[3]+0.0)/this->count/50;//100;
    double ctop25 = (cr[4]+0.0)/this->count/100;//100;

    double top1 = (r[0]+0.0)/this->count;
    double top10 = (r[1]+0.0)/this->count;//10;
    double top50 = (r[2]+0.0)/this->count;//50;
    double top100 = (r[3]+0.0)/this->count;//100;
    cout<<"recall 1 "<<top1<<" 3 "<<top10<<" 5 "<<top50<<" 10 "<<top100<<endl;
    cout<<"cover 5 "<<ctop5<<" 10 "<<ctop10<<" 25 "<<ctop15<<" 50 "<<ctop20<<" 100 "<<ctop25<<endl;

    /**/

    delete [] rtop;
    return 0;
}

int NNG::saveNNGraph(const char* dstfn)
{
    ofstream os(dstfn);
    if(!os.is_open())
    {
        cout<<"erro "<<"can not open file "<<dstfn<<" for save NNGraph\n";
        return -1;
    }
    os<<this->count<<" "<<this->TopK<<endl;
    for(int i = 0; i < this->count; i++)
    {
        os<<i<<" ";
        for(int j = 0; j < this->TopK; j++)
        {
            os<<this->topList[i*TopK+j]<<" ";
        }
        os<<"\n";
    }
    os.close();
    return 0;
}

int NNG::getTopList(unsigned *&top, unsigned &topNum)
{
    top = new unsigned[this->count*TopK];
    topNum = this->TopK;
    memcpy(top, this->topList, sizeof(unsigned)*TopK*this->count);
    return 0;
}

void NNG::CHTest()
{
    const char *srcMatFn1 = "/home/chenghaod/datasets/bignn/mat/sift1m/sift_base.txt";
    const char *srcMatFn2 = "/home/chenghaod/datasets/glove.fvecs";
    const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.itm";
    const char *dstfn = "/home/chenghaod/sift_learn_test.out";
    //const char *srcfn = "/home/chenghaod/dataset/bow/ox5k_nrsift_65k_hesaff_bow.txt";
    //const char *srcfn = "/home/chenghaod/src/finished/minhash/etc/data/holidays_nrsift_20k_hess_bow.mat";
    //const char *srcfn = "/home/chenghaod/data/sift/raw/sift_base.fvecs";
    // const char *srcfn = "sift_base.fvecs";

    //const char *srcfn = "/home/chenghaod/data/sift/raw/sift_learn.txt";
    const char *srcfn = "1mMat.txt";
    ///const char *srcfn = "/home/chenghaod/datasets/gist/gist_base.fvecs";
    //const char *srcfn = "/home/chenghaod/datasets/gist/gist_learn.fvecs";

    //paraTest("/home/chenghaod/data/sift/raw/sift_learn.txt", "tmpfile.txt");
    //exit(0);

    ///const char *nnfn = "vlad10m50.txt";
    const char *nnfn = "sift_grd.txt";


    NNG *mykm = new NNG(new DISTANCE::FastL2Distance<float>);
    mykm->loadData(srcfn);
    mykm->NNBuild(1);
    mykm->saveNNGraph(nnfn);
    delete mykm;
    cout<<"time\n";
}

void NNG::test()
{

}



