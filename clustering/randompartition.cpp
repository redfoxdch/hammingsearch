#include "randompartition.h"

#include "print2scrn.h"
#include "ioagent.h"
#include "cleaner.h"
#include "vstring.h"
#include "nnitem.h"
#include "pqmath.h"
#include "timer.h"

#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <limits>
#include <cassert>
#include <cstdlib>
#include <bitset>
#include <vector>
#include <cmath>
#include <ctime>

using namespace std;

const unsigned int RandomPartition::NTRAILS= 1;
const unsigned int RandomPartition::NIter = 20;
const float RandomPartition::Err0 = 0;
static unsigned ch = 10;
static unsigned t1 = 0, t2 = 0, t3 = 0;

inline bool mapLs(const pair<float, unsigned> &a, const pair<float, unsigned> &b)
{
    return a.first < b.first;
}

RandomPartition::RandomPartition()
{
    this->D1      = NULL;
    this->arrayD  = NULL;
    this->Ns = NULL;
    this->_INIT_  = false;
    this->_REFER_ = false;
    this->distance = new DISTANCE::FastL2Distance<float>;
    kmMtd   = _xtkmn_;
    strcpy(mthStr, "_xtk_");
    cout<<"Method ........................... random partition clustering\n";
}

RandomPartition::RandomPartition(const DISTANCE::Distance<float>* d):distance(d)
{
    this->D1      = NULL;
    this->arrayD  = NULL;
    this->_INIT_  = false;
    this->_REFER_ = false;
    kmMtd   = _xtkmn_;
    strcpy(mthStr, "_xtk_");
    cout<<"Method ........................... random partition clustering\n";
}

bool RandomPartition::init(const char *srcfn)
{
    cout<<"Loading matrix ................... ";
    assert(srcfn);
    strcpy(srcmatfn, srcfn);

    if(VString::endWith(srcfn, ".txt"))
    {
        unsigned ct, nd;
        this->data = IOAgent::loadMatrix(srcfn, ct, nd);
        this->count = ct;
        this->ndim = nd;
        this->dataType = 0;
    }
    else if(VString::endWith(srcfn, ".fvecs"))
    {
        unsigned ct, nd;
        this->data = IOAgent::load_fvecs(srcfn, ct, nd);
        this->count = ct;
        this->ndim = nd;
        this->dataType = 0;
    }
    else if(VString::endWith(srcfn, ".mat"))
    {
        unsigned ct, nd;
        this->data = IOAgent::loadDat(srcfn, ct, nd);

        this->count = ct;
        this->ndim = nd;
        /** this->sdata    = IOAgent::loadSparse(srcfn, this->count, this->ndim);
        this->dataType = 1;/**/
    }
    else if(VString::endWith(srcfn, ".tvecs"))
    {
        unsigned ct, nd;
        vector<unsigned int> ids;
        this->data = IOAgent::loadTruncMat(srcfn, ct, nd, ids);
        this->count = ct;
        this->ndim = nd;
        ids.clear();
        /** this->sdata    = IOAgent::loadSparse(srcfn, this->count, this->ndim);
        this->dataType = 1;/**/
    }
    else
    {
        unsigned ct, nd;
        this->data = IOAgent::loadItms(srcfn, "fsvtab", 10000000, ct, nd);
        this->count = ct;
        this->ndim = nd;
        //cout<<"Unrecognizable input file format!!!\n";
        //this->data = NULL;
        //exit(0);
    }

    if(this->data == NULL)
    {
        cout<<"Exceptions ....................... Loading matrix failed!\n";
        exit(0);
    }


    cout<<this->count<<"x"<<this->ndim<<endl;
    this->_INIT_  = true;
    this->_REFER_ = false;
    this->Ns      = (int*)calloc(this->clnumb, sizeof(int));
    this->arrayD = (double*)calloc(this->clnumb, this->ndim*sizeof(double));


    cout<<"clust num\t"<<this->clnumb<<endl;
    return true;
}

bool RandomPartition::init(float *mat, const int row, const int dim)
{
    assert(mat);

    this->data   = mat;
    this->count  = row;
    this->ndim   = dim;

    this->_INIT_  = true;
    this->_REFER_ = true;
    this->Ns      = (int*)calloc(this->clnumb, sizeof(int));
    this->arrayD = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
    cout<<"clust num\t"<<this->clnumb<<endl;
    return true;
}


bool RandomPartition::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
{
    if(verbose)
        cout<<"Distance function ................ l2\n";

    if(verbose)
        cout<<"Seeds ............................ ";
    if(!strcmp(_seed_, "rnd"))
    {
        if(verbose)
            cout<<"rand\n";
        seed = _rnd_;
    }
    else if(!strcmp(_seed_, "kpp"))
    {

        if(verbose)
            cout<<"kpp\n";
        seed = _kpp_;
    }
    else
    {
        if(verbose)
            cout<<"non\n";
        seed = _non_;
    }

    if(verbose)
        cout<<"Optimization function ............ ";
    if(!strcmp(crtrn, "i1"))
    {
        myoptz     = _I1_;
        if(verbose)
            cout<<"I1\n";
    }
    else if(!strcmp(crtrn, "i2"))
    {
        myoptz  = _I2_;
        if(verbose)
            cout<<"I2\n";
    }
    else if(!strcmp(crtrn, "i3"))
    {
        myoptz  = _I3_;
        if(verbose)
            cout<<"I3\n";
    }
    else if(!strcmp(crtrn, "i4") )
    {
        myoptz  = _I4_;
        if(verbose)
            cout<<"e1\n";
    }
    else if(!strcmp(crtrn, "e1") )
    {
        myoptz  = _E1_;
        if(verbose)
            cout<<"e1\n";
    }
    else if(!strcmp(crtrn, "e2") )
    {
        myoptz  = _E2_;
        if(verbose)
            cout<<"e2\n";
    }
    else if(!strcmp(crtrn, "t1"))
    {
        myoptz = _T1_;
        if(verbose)
            cout<<"t1\n";
    }
    else if(!strcmp(crtrn, "i4"))
    {
        myoptz = _I4_;
        if(verbose)
            cout<<"i4\n";
    }
    else
    {
        cout<<"Unkown optimize option '"<<crtrn<<"'!\n";
        exit(0);
    }

    return true;
}


int RandomPartition::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    int unsigned i, j, k, n, d, num, nlabel, clabel = 0;
    vector<unsigned> tmpcluster;

    this->D1    = new float[this->ndim];
    this->D2    = new float[this->ndim];
    this->C1    = new float[this->ndim];
    this->C2    = new float[this->ndim];
    this->tmpC1 = new float[this->ndim];
    this->tmpC2 = new float[this->ndim];


    this->clusters = new vector<unsigned>[this->clnumb];
    for(i = 0; i < this->clnumb; i++)
    {
        this->clusters[i].clear();
    }
    for(i = 0; i < this->count; i++)
    {
        this->clusters[0].push_back(i);
    }
    for(j = 0; j < clust_num; j++)
    {
        this->infoMap[j].E = 0.0f;
        this->infoMap[j].n = 0;
    }

    clabel = 0;
    nlabel = 1;
    d = 0;
    k = 0;
    num = 1;
    if(verbose)
        cout<<"Clustering ....................... on progress\n";
    while(nlabel < clnumb)
    {

        if(!this->_REFINE_)
            this->I2P(clabel, nlabel);
        else
            this->random(clabel, nlabel);
        /**/
        tmpcluster.clear();
        tmpcluster.reserve(clusters[clabel].size());
        for(auto it = this->clusters[clabel].begin(); it != this->clusters[clabel].end(); it++)
        {
            tmpcluster.push_back(*it);
        }
        clusters[clabel].clear();
        clusters[clabel].reserve(this->Ns[clabel]);
        clusters[nlabel].reserve(this->Ns[nlabel]);
        for(auto it = tmpcluster.begin(); it != tmpcluster.end(); it++)
        {
            clusters[labels[*it]].push_back(*it);
        }
        /**/
        k++;
        if(k == num)
        {
            d++;
            k = 0;
            num = num*2;
        }
        clabel = k;
        nlabel++;
    }

    for(i = 0; i < count; i++)
    {
        clabel = labels[i];
        for(j = 0; j < this->ndim; j++)
        {
            this->arrayD[clabel*this->ndim+j] += this->data[i*this->ndim+j];
        }
    }

    this->calAVGDist(this->arrayD, this->Ns, clust_num, infoMap);
    if(!this->_REFER_)
    {
        this->saveCenters(dstfn, 0);
    }



    for(i = 0; i < this->clnumb; i++)
    {
        this->clusters[i].clear();
    }
    delete [] this->clusters;
    this->clusters = NULL;
    delete []D1;
    delete []D2;
    delete []C1;
    delete []C2;
    delete []tmpC1;
    delete []tmpC2;
    D1 = NULL;
    D2 = NULL;
    C1 = NULL;
    C2 = NULL;
    tmpC1 = NULL;
    tmpC2 = NULL;

    return 0;
}

int RandomPartition::I2P(unsigned clabel, unsigned nwlbl)
{
    unsigned int _loc1 = 0, _loc2 = 0, loc = 0;
    unsigned int i = 0, j = 0, r = 0;
    unsigned int cpysize = sizeof(float)*this->ndim;
    float tmp1 = 0, tmp2 = 0, len = 0, E1, E2, n1, n2;
    float innDcts[2];

    vector<unsigned> crntClust = this->clusters[clabel];

    unsigned int numb = crntClust.size();
    int *rl = new int[numb];
    for(i = 0; i < numb; i++)
    {
        rl[i] = i;
    }

    if(numb == 2)
    {
        i  = crntClust[1];
        this->labels[i] = nwlbl;
        this->Ns[clabel] = 1;
        this->Ns[nwlbl] = 1;
        return true;
    }
    else if(numb < 2)
    {
        cout<<"unseperatable\n";
    }

    int *nwlabels = new int[numb];
    memset(nwlabels, 0, sizeof(int)*numb);

    int v = 0, tmpn1 = 0, tmpn2 = 0;
    unsigned int sed1 = 0, sed2 = 0, iter = 0;
    float tmpEg = 0, optEg = 0, dst1 = 0,  dst2 = 0;
    bool _UPDATE_ = true, F2S = false;
    n1 = n2 = 0;
    /**** Initial assigment ****************/


    memset(D1, 0, cpysize);
    memset(D2, 0, cpysize);
    innDcts[0] = innDcts[1] = 0;

    if(this->seed != _non_)
    {
        random_shuffle(rl, rl+numb);
        sed1 = crntClust[rl[0]];
        sed2 = crntClust[rl[1]];
        _loc1 = sed1*this->ndim;
        _loc2 = sed2*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            tmpC1[j] = this->data[_loc1+j];
            tmpC2[j] = this->data[_loc2+j];
        }
        for(i = 0; i < numb; i++)
        {
            v = crntClust[i];
            mvs[0]++;
            dst1 = PQMath::l2f(tmpC1, 0, this->data, v, this->ndim);
            dst2 = PQMath::l2f(tmpC2, 0, this->data, v, this->ndim);
            loc  = v*this->ndim;

            if(dst1 < dst2)
            {
                for(j = 0; j < this->ndim; j++)
                {
                    this->D1[j] += this->data[loc+j];
                }
                tmpn1++;
                nwlabels[i] = -1;
            }
            else
            {
                for(j = 0; j < this->ndim; j++)
                {
                    D2[j] += this->data[loc+j];
                }
                tmpn2++;
                nwlabels[i] = v;
            }
        }
    }
    else
    {
        random_shuffle(rl, rl+numb);
        loc = rand()%numb;
        if(loc == 0)
            loc = 1;
        for(i = 0; i < loc; i++)
        {
            v = crntClust[rl[i]];
            for(j = 0; j < this->ndim; j++)
            {
                this->D1[j] += this->data[v*ndim+j];
            }
            tmpn1++;
            nwlabels[rl[i]] = -1;
        }
        while(i < numb)
        {
            v = crntClust[rl[i]];
            for(j = 0; j < this->ndim; j++)
            {
                this->D2[j] += this->data[v*ndim+j];
            }
            tmpn2++;
            nwlabels[rl[i]] = v;
            i++;
        }
    }


    innDcts[0] = PQMath::fvec_norm(D1, this->ndim, 2);
    innDcts[1] = PQMath::fvec_norm(D2, this->ndim, 2);
    innDcts[0] = innDcts[0]*innDcts[0];
    innDcts[1] = innDcts[1]*innDcts[1];

    optEg = innDcts[0]/tmpn1 +  innDcts[1]/tmpn2;


    /******** Incremental optimization *****/
    do
    {

        // random_shuffle(rl, rl + numb);
        _UPDATE_ = false;
        for(iter = 0; iter < numb; iter++)
        {
            mvs[0]++;
            r   = rl[iter];
            v   = crntClust[r];
            loc = v*this->ndim;
            F2S = false;
            if(nwlabels[r] == -1)
            {
                if(tmpn1 > 1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                    F2S = true;
                }
                else
                {
                    continue;
                }
            }
            else
            {
                if(tmpn2 > 1)
                {

                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                    F2S = false;
                }
                else
                {
                    continue;
                }
            }
            //tts1 = clock();
            tmp1 = tmp2 = 0;
            for(j = 0; j < this->ndim; j++)
            {
                tmp1 += D1[j]*data[loc+j];
                tmp2 += D2[j]*data[loc+j];
            }

            len = this->lens[v];
            if(F2S)
            {
                tmp1 = innDcts[0] - 2*tmp1 + len;
                tmp2 = innDcts[1] + 2*tmp2 + len;

            }
            else
            {
                tmp1 = innDcts[0] + 2*tmp1 + len;
                tmp2 = innDcts[1] - 2*tmp2 + len;
            }

            tmpEg = tmp1/tmpn1 + tmp2/tmpn2;
            //  t1 += clock() - tts1;
            //  tts2 = clock();
            if(tmpEg > (optEg+0.00001))
            {
                if(F2S)
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        D1[j] -= this->data[loc+j];
                        D2[j] += this->data[loc+j];
                    }
                }
                else
                {
                    for(j = 0; j < this->ndim; j++)
                    {
                        D1[j] += this->data[loc+j];
                        D2[j] -= this->data[loc+j];
                    }
                }
                optEg = tmpEg;
                innDcts[0] = tmp1;
                innDcts[1] = tmp2;
                _UPDATE_ = true;
                mvs[1]++;
            }
            else
            {
                ///roll-back if being updated
                if(nwlabels[r] == -1)
                {
                    tmpn1--;
                    tmpn2++;
                    nwlabels[r] = v;
                }
                else
                {
                    tmpn1++;
                    tmpn2--;
                    nwlabels[r] = -1;
                }
            }
            //  t2 += clock() - tts2;
        }
        iter++;
        if(iter > this->NIter)
            break;
    }
    while(_UPDATE_);

    float vr;
    pair<float, unsigned> *vlmap = new pair<float, unsigned>[numb];
    for(i = 0; i < ndim; i++)
    {
        D2[i] -= D1[i];
        D2[i] /= numb;
    }

    for(i = 0; i < numb; i++)
    {
        r = crntClust[i];
        vr = 0;
        for(j = 0; j < ndim; j++)
        {
            vr += data[r*ndim+j]*D2[j];
        }
        vlmap[i].first = vr;
        vlmap[i].second = r;
    }

    sort(vlmap, vlmap+numb, mapLs);
    n1 = (numb+1)/2;
    n2 = numb-n1;
    this->Ns[clabel] = n1;
    this->Ns[nwlbl] = n2;
    for(i = 0; i < n1; i++)
    {
        v = vlmap[i].second;
        labels[v] = clabel;
    }

    while(i < numb)
    {
        v = vlmap[i].second;
        labels[v] = nwlbl;
        i++;
    }

    delete [] vlmap;
    delete [] nwlabels;
    delete [] rl;
    nwlabels = NULL;

    return true;
}

int RandomPartition::RP(unsigned clabel, unsigned nlabel)
{
    unsigned i, j, k,  numb, n1, n2, r, v, m1, m2;
    float mid1, mid2, v1, v2, val, vmin = -1;//RAND_MAX;
    vector<unsigned > crntClust = this->clusters[clabel];
    pair<float, unsigned> *vlmap;
    numb = crntClust.size();
    float vr;
    bool iszero = true;

    vlmap = new pair<float, unsigned>[numb];

    unsigned *rl = new unsigned[numb];
    for(i = 0; i < numb; i++)
    {
        rl[i] = i;
    }
    for(k = 0; k < this->NIter; k++)
    {
        ///random direction generate
        random_shuffle(rl, rl+numb);
        r = crntClust[rl[0]];
        v = crntClust[rl[1]];
        for(i = 0; i < ndim; i++)
        {
            D2[i] = data[r*ndim+i];
            D1[i] = data[v*ndim+i];
        }
        for(i = 0; i < ndim; i++)
        {
            D2[i] = D2[i] - D1[i];
            if(D2[i] != 0)
                iszero = false;
        }
        while(iszero)
        {
            for(i = 0; i < ndim; i++)
            {
                D2[i] = (rand()+0.0)/RAND_MAX;
                if(D2[i] != 0)
                    iszero = false;
            }
        }
        /**
            for(i = 0; i < ndim; i++)
            {
                D2[i] = (rand()-0.0)/RAND_MAX-0.5;
            }
        /**/
        for(i = 0; i < numb; i++)
        {
            r = crntClust[i];
            vr = this->distance->dot(data+r*ndim, D2, ndim);
            vlmap[i].first = vr;
            vlmap[i].second = r;
        }

        sort(vlmap, vlmap+numb, mapLs);

        n1 = (numb+1)/2;
        n2 = numb-n1;
        m1 = n1/2;
        m2 = (numb+n1)/2;

        val = vlmap[m2].first - vlmap[m1].first;

///updata points if find the min erro direction
        if(val > vmin)
        {
            //cout<<numb<<endl;
            vmin = val;
            for(i = 0; i < n1; i++)
            {
                v = vlmap[i].second;
                labels[v] = clabel;
            }
            while(i < numb)
            {
                v = vlmap[i].second;
                labels[v] = nlabel;
                i++;
            }
            this->Ns[clabel] = n1;
            this->Ns[nlabel] = n2;
        }
    }
    ///cout<<"n1 "<<n1<<" n2 "<<n2<<endl;

    delete [] vlmap;
    delete [] rl;
    return true;
}

int RandomPartition::KD(unsigned clabel, unsigned nlabel)
{
    unsigned i, j, k,  numb, n1, n2, r, v, m1, m2;
    float mid1, mid2, v1, v2, val, vmin = -1;//RAND_MAX;
    vector<unsigned > crntClust = this->clusters[clabel];
    pair<float, unsigned> *vlmap;
    numb = crntClust.size();
    float vr;
    bool iszero = true;

    vlmap = new pair<float, unsigned>[numb];

    ///random direction generate
    k = rand()%this->ndim;
    for(i = 0; i < numb; i++)
    {
        r = crntClust[i];
        vr = 0;
        vr += data[r*ndim+k]*D2[k];
        vlmap[i].first = vr;
        vlmap[i].second = r;
    }

    sort(vlmap, vlmap+numb, mapLs);

    n1 = (numb+1)/2;
    n2 = numb-n1;
    m1 = n1/2;
    m2 = (numb+n1)/2;

    val = vlmap[m2].first - vlmap[m1].first;

///updata points if find the min erro direction
    if(val > vmin)
    {
        //cout<<numb<<endl;
        vmin = val;
        for(i = 0; i < n1; i++)
        {
            v = vlmap[i].second;
            labels[v] = clabel;
        }
        while(i < numb)
        {
            v = vlmap[i].second;
            labels[v] = nlabel;
            i++;
        }
        this->Ns[clabel] = n1;
        this->Ns[nlabel] = n2;
    }
    ///cout<<"n1 "<<n1<<" n2 "<<n2<<endl;

    delete [] vlmap;
    return true;
}

int RandomPartition::var(unsigned label, float *v)
{

}

int RandomPartition::random(unsigned clabel, unsigned nlabel)
{
    unsigned i, j, k,  numb, n1, n2, r, v;
    float v1, v2;//RAND_MAX;
    vector<unsigned > crntClust = this->clusters[clabel];
    pair<float, unsigned> *vlmap;
    numb = crntClust.size();
    float vr;
    bool iszero = true;

    vlmap = new pair<float, unsigned>[numb];

    unsigned *rl = new unsigned[numb];
    for(i = 0; i < numb; i++)
    {
        rl[i] = i;
    }
    ///random direction generate
    random_shuffle(rl, rl+numb);
    r = crntClust[rl[0]];
    v = crntClust[rl[1]];
    for(i = 0; i < ndim; i++)
    {
        D2[i] = data[r*ndim+i];
        D1[i] = data[v*ndim+i];
    }
    for(i = 0; i < ndim; i++)
    {
        D2[i] = D2[i] - D1[i];
        if(D2[i] != 0)
            iszero = false;
    }
    while(iszero)
    {
        for(i = 0; i < ndim; i++)
        {
            D2[i] = (rand()+0.0)/RAND_MAX;
            if(D2[i] != 0)
                iszero = false;
        }
    }

    for(i = 0; i < numb; i++)
    {
        r = crntClust[i];

        vr = this->distance->dot(D2, this->data+r*ndim, ndim);
        vlmap[i].first = vr;
        vlmap[i].second = r;
    }

    sort(vlmap, vlmap+numb, mapLs);

    n1 = (numb+1)/2;
    n2 = numb-n1;

///updata points if find the min erro direction

    for(i = 0; i < n1; i++)
    {
        v = vlmap[i].second;
        labels[v] = clabel;
    }
    while(i < numb)
    {
        v = vlmap[i].second;
        labels[v] = nlabel;
        i++;
    }
    this->Ns[clabel] = n1;
    this->Ns[nlabel] = n2;
    ///cout<<"n1 "<<n1<<" n2 "<<n2<<endl;

    delete [] vlmap;
    delete [] rl;
    return true;
}

void RandomPartition::saveCenters(const char *dstfn, bool append)
{
    unsigned int clabel = 0, j, loc, rCNum = 0;
    ofstream *outStrm  = NULL;

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->infoMap[clabel].n > 0)
        {
            rCNum++;
        }
    }

    if(!this->_INIT_||rCNum == 0)
    {
        return ;
    }

    if(append)
    {
        outStrm = new ofstream(dstfn, ios::app);
    }
    else
    {
        outStrm = new ofstream(dstfn, ios::out);
        (*outStrm)<<rCNum<<" "<<this->ndim<<endl;;
    }
    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        loc  = clabel*this->ndim;

        if(this->infoMap[clabel].n <= 0)
            continue;

        for(j = 0; j < this->ndim; j++)
        {
            (*outStrm)<<this->arrayD[loc+j]/this->infoMap[clabel].n<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
    return ;
}

int RandomPartition::fetchCenters(float *centers)
{
    unsigned int clabel = 0, j = 0, loc = 0, idxi = 0, rCNum = 0;
    assert(centers);

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->infoMap[clabel].n > 0)
        {
            rCNum++;
        }
    }

    if(!this->_INIT_||rCNum == 0)
    {
        memset(centers, 0, this->clnumb*this->ndim*sizeof(float));
        return 0;
    }

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        loc   = clabel*this->ndim;

        if(this->infoMap[clabel].n <= 0)
            continue;

        for(j = 0; j < this->ndim; j++)
        {
            centers[idxi + j] = this->arrayD[loc+j]/this->infoMap[clabel].n;
        }
        idxi += this->ndim;
    }
    return rCNum;
}

RandomPartition::~RandomPartition()
{
    if(this->Ns != NULL)
    {
        free(this->Ns);
        this->Ns = NULL;
    }
    if(this->arrayD != NULL)
    {
        free(this->arrayD);
        this->arrayD = NULL;
    }
}

bool vtest()
{
    /**/
    const char *srcfn = "/home/chenghaod/data/sift/raw/sift_learn.txt";
    const char *dstfn = "/home/chenghaod/dataset/vlad/ox5kbowhe.out";
    unsigned i, j, k, id;
    unsigned clnumb;
    unsigned cllv = 10;
    unsigned it = 10;
    unsigned num = 100000;
    unsigned rt = 0;
    clnumb = pow(2, cllv);
    int *label = new int[num];

    unsigned tmpid[200];
    unsigned *rtop = new unsigned[200*num];
    ifstream is("nn.txt");
    for(i = 0; i < num; i++)
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

    RandomPartition *mykm;
    for(i = 0; i < it; i++)
    {
        mykm = new RandomPartition();
        mykm->buildcluster(srcfn, dstfn, "non", "large", "i2", clnumb, false);
        mykm->getLabel(label);
        delete mykm;

        for(j = 0; j < num; j++)
        {
            if(label[j] == label[rtop[j*200+99]])
                rt++;
        }
    }

    cout<<rt<<" rate "<<(rt+0.0)/num/it<<endl;

    delete [] label;
    delete [] rtop;
    /**/

}

void RandomPartition::test()
{
    // vtest();
    //  exit(0);
    const char *srcMatFn1 = "/home/chenghaod/datasets/bignn/mat/sift1m/sift_base.txt";
    const char *srcMatFn2 = "/home/chenghaod/datasets/bignn/mat/sift_learn.txt";
    const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.itm";
    const char *dstfn = "/home/chenghaod/datasets/vlad/ox5kbowhe.out";
    //const char *srcfn = "/home/chenghaod/datasets/bow/ox5k_nrsift_65k_hesaff_bow.txt";
    //const char *srcfn = "/home/chenghaod/src/finished/minhash/etc/data/holidays_nrsift_20k_hess_bow.mat";
    const char *srcfn = "/home/chenghaod/data/sift/raw/sift_base.fvecs";
    ///const char *srcfn = "/home/chenghaod/data/sift/raw/sift_learn.txt";

    // do{
    RandomPartition *mykm = new RandomPartition();
    mykm->buildcluster(srcfn, dstfn, "non", "large", "i2", 1000, false);
    delete mykm;
    //cout<<"t1 "<<t1<<" \n"<<"t2 "<<t2<<"\n"<<"t3 "<<t3<<endl;
    //    i = i*2;
    // }while(i <= 8192);
}

