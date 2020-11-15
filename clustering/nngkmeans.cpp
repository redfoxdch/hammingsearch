#include "randompartition.h"
#include "nngkmeans.h"
#include "util/ioagent.h"
#include "util/vstring.h"
#include "util/pqmath.h"
#include "util/timer.h"
#include "nn/nn.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <map>

using namespace std;

static unsigned mx = 50;
static unsigned rvn = 20;
static unsigned Check = 10;
static unsigned topRange = 0;
static unsigned p = 50;

const unsigned NNGKMeans::RangeFirst = mx+p;
const unsigned NNGKMeans::RangeSecond = topRange+p;

const unsigned int NNGKMeans::NIter   = 30;
/**/
const unsigned int NNGKMeans::TreeNum = 10;
static unsigned ENum = 0;
const float NNGKMeans::Err0           = 0;
const unsigned long NNGKMeans::Treesize    = 50;
const unsigned long NNGKMeans::TopK        = 50;
const unsigned long NNGKMeans::NNTop = 25;

const unsigned NNGKMeans::NBIter = 1;
/**/


//static unsigned TreeNum = 10, TopK = 25, NBIter = 1, Treesize = 50;
static unsigned t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0, t8 = 0, titer = 0, titer0 = 0;
static double rc, dt;
static ofstream os;

bool static inline fuLs(const pair<float, unsigned> &a, const pair<float, unsigned> &b)
{
    if(a.first == b.first)
        return a.second < b.second;
    return a.first < b.first;
}

NNGKMeans::NNGKMeans()
{
    _INIT_ = false;
    this->data    = NULL;
    this->Ds      = NULL;
    this->Ns = NULL;
    this->Es = NULL;
    this->arrayD  = NULL;
    this->infoMap = NULL;
    this->visitList = NULL;
    this->radius = NULL;
    this->topDist = NULL;
    this->topList = NULL;
    kmMtd   = _xtkmn_;
    strcpy(mthStr, "_xtk_");
    cout<<"Method ........................... Neiborhood  K-means\n";
}

bool NNGKMeans::init(const char *srcfn)
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

    cout<<"clust num\t"<<this->clnumb<<endl;
    return true;
}

bool NNGKMeans::init(float *mat, const int row, const int dim)
{
    assert(mat);
    if(this->data != NULL)
    {
        cout<<"warning : reuse data before\n";
        return -1;
    }
    this->data   = mat;
    this->count  = row;
    this->ndim   = dim;

    cout<<row<<"x"<<dim<<endl;

    this->_INIT_  = true;

    cout<<"clust num\t"<<this->clnumb<<endl;
    return true;
}

bool NNGKMeans::config(const char *_seed_, const char *lg_first, const char *crtrn, int verbose)
{


    return true;
}

int NNGKMeans::updateNNGraph()
{
    unsigned sz = 0;
    vector<unsigned> *clusters;
    clusters = new vector<unsigned>[this->clnumb];
    unsigned long i = 0, j = 0, k = 0, id, r, col, r1 = 0, r2 = 0;
    float dst = 0;

    this->label2cluster(clusters);
    cout<<"update top\n";

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
                for(col = 0; col < ndim; col++)
                {
                    dst += (data[r2+col] - data[r1+col])*(data[r2+col] - data[r1+col]);
                }
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
    // cout<<"s "<<(s+0.0)/this->count<<endl;
    this->adjustNN();

    delete [] cdst;
    cdst = NULL;
    delete [] clusters;
    clusters = NULL;
    cout<<"update top finished\n";
    return 0;
}

int NNGKMeans::adjustNN()
{
    unsigned long i, j;
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

int NNGKMeans::initNNGraph()
{
    unsigned long i, j, k, id, rid, r1 = 0, r2 = 0;
    float dst = 0, mxdst = 0;
    unsigned *rl = new unsigned[this->count];
    for(i = 0; i < this->count; i++)
    {
        rl[i] = i;
    }
    random_shuffle(rl, rl+this->count);

    for(i = 0; i < this->count; i++)
    {
        rid = rand()%this->count;
        mxdst = 0;
        for(j = 0; j < TopK; j++)
        {
            id = rl[(j+rid)%this->count];
            if(id == i)
            {
                id = rl[(rid+TopK+1)%this->count];
            }

            topList[i*TopK+j] = id;
            dst = 0;
            r1  = i*this->ndim;
            r2  = id*this->ndim;
            for(k = 0; k < this->ndim; k++)
            {
                dst += (data[r1 + k]-data[r2+k])*(data[r1+k]-data[r2+k]);
            }
            topDist[i*TopK+j] = dst;
            if(dst > mxdst)
                mxdst = dst;
        }
        radius[i] = mxdst;
    }

    delete [] rl;
    rl = NULL;
    return 0;
}

int NNGKMeans::buildNNGraph()
{
    unsigned clustNum, i, j, ct = true;
    /// init
    clustNum = this->clnumb;
    cout<<"cl "<<this->clnumb<<endl;
    this->clnumb = this->count/Treesize;
    this->clnumb = pow(2, floor(log(this->clnumb)/log(2)));
    cout<<this->clnumb<<endl;

    if(this->Ds != NULL)
    {
        free(this->Ds);
        this->Ds = NULL;
    }
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
    if(this->Es != NULL)
    {
        free(this->Es);
        this->Es = NULL;
    }

    for(i = 0; i < TreeNum; i++)
    {
        cout<<"initialize cluster begin .......\n";
        RandomPartition *rpTree = new RandomPartition;
        rpTree->buildcluster(this->data, this->count, this->ndim, "tmpfile.txt",  "rnd", "large", "i2", this->clnumb, 0);
        rpTree->getLabel(this->labels);
        delete rpTree;
        cout<<"Initialize cluster finished\n";

        Timer *itm = new Timer;
        itm->start();
        if(i > 0)
        {
            this->Es     = (double*)calloc(this->clnumb, sizeof(double));
            this->arrayD = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
            this->Ds     = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
            this->Ns     = (int*)calloc(this->clnumb, sizeof(int));

            this->optzI2(1, NNTop, NBIter);
            free(this->Ds);
            this->Ds = NULL;
            free(this->Ns);
            this->Ns = NULL;
            free(this->arrayD);
            this->arrayD = NULL;
            free(this->Es);
            this->Es = NULL;

        }
        cout<<"i time";
        itm->end(true);
        Timer *utm = new Timer;
        utm->start();
        updateNNGraph();
        cout<<"u time";
        utm->end(true);
    }

    this->clnumb = clustNum;


    return true;
}

int NNGKMeans::recall()
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
        NNGKMeans::appRecall();
        return 0;
    }

    if(!is.is_open())
    {
        cout<<"cant not open groundturth , use approximate ground truth \n";
        NNGKMeans::appRecall();
        return 0;
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

    /*ofstream os("sift_grd.txt");
    os<<this->count<<" "<<200<<endl;
    for(i = 0; i < this->count; i++)
    {
        os<<i<<" ";
        for(j = 0; j < 200; j++)
            os<<rtop[i*200+j]<<" ";
    }
    os.close();*/

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
    vector<unsigned> s1;
    vector<unsigned> s2;
    vector<unsigned> s3;
    s3.resize(TopK);
    s1.resize(TopK);
    s2.resize(TopK);
    for(i = 0; i < this->count; i++)
    {
        for(unsigned d = 0; d < 5; d++)
        {
            if(c[d] > TopK)
                break;
            for(j = 0; j < TopK; j++)
            {
                s1[j] = topList[i*TopK+j];
                s2[j] = rtop[i*200+j];
            }
            sort(s1.begin(), s1.begin()+c[d]);
            sort(s2.begin(), s2.begin()+c[d]);
            auto it = set_intersection(s1.begin(), s1.begin()+c[d], s2.begin(), s2.begin() + c[d], s3.begin());
            cr[d] += it - s3.begin();
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
    rc = top1;
    return 0;
}

double NNGKMeans::appRecall()
{
    const unsigned num = 100;
    double rc;
    unsigned row = this->count;
    unsigned col = this->ndim;
    unsigned *seq = new unsigned[row];

    for(unsigned i = 0; i < row; i++)
    {
        seq[i] = i;
    }

    random_shuffle(seq, seq+row);
    unsigned *sid = new unsigned[num];

    for(unsigned i = 0; i < num; i++)
    {
        sid[i] = seq[i];
    }
    delete [] seq;

    unsigned *rtop = new unsigned[num*TopK];

    unsigned i, j, k, l, id;
    double dst;
    multimap<double, unsigned> * mp;
    mp = new multimap<double, unsigned>[num];
    multimap<double, unsigned>::iterator it;
    ///get the nearst neighbor map;
    cout<<"get nearst neighbor start\n";
    for(i = 0; i < num; i++)
    {
        id = sid[i];
        mp[i].clear();
        for(j = 0; j < row; j++)
        {
            if(id == j)
                continue;
            dst = 0;
            for(k = 0; k < col; k++)
            {
                dst += (data[id*col+k]-data[j*col+k])*(data[id*col+k]-data[j*col+k]);
            }
            if(mp[i].size() < TopK)
            {
                mp[i].insert(pair<double, unsigned>(dst, j));
            }
            else
            {
                it = mp[i].end();
                it--;
                if(it->first > dst)
                {
                    mp[i].erase(it);
                    mp[i].insert(pair<double, unsigned>(dst, j));
                }
            }
        }
    }
    ///get nearst neighbor;
    for(i = 0; i < num; i++)
    {
        it = mp[i].begin();
        for(j = 0; j < TopK; j++)
        {
            rtop[i*TopK+j] = it->second;
            it++;
        }
    }
    for(i = 0; i < num; i++)
    {
        mp[i].clear();
    }
    delete [] mp;

    vector<unsigned> s1;
    vector<unsigned> s2;
    vector<unsigned> s3;
    s3.resize(TopK);
    s1.resize(TopK);
    s2.resize(TopK);
    unsigned inter = 0;
    unsigned cid[5] = {1,10,25,50,100};
    double rt[5] = {0};
    set<unsigned> tmpset;
    for(i = 0; i < num; i++)
    {
        id = sid[i];
        tmpset.clear();
        for(j = 0; j < TopK; j++)
        {
            s1[j] = topList[id*TopK+j];
            s2[j] = rtop[i*TopK+j];
            tmpset.insert(s1[j]);
        }
        if(tmpset.size() < TopK)
        {
            cout<<"erro same id \n";
            return -1;
        }

        for(unsigned d = 0; d < 5; d++)
        {
            if(cid[d] > TopK)
                break;
            sort(s1.begin(), s1.begin()+cid[d]);
            sort(s2.begin(), s2.begin()+cid[d]);
            auto it = set_intersection(s1.begin(), s1.begin()+cid[d], s2.begin(), s2.begin()+cid[d], s3.begin());
            rt[d] += it - s3.begin();
        }
    }

    for(unsigned d = 0; d < 5; d++)
    {
        rt[d] = (rt[d])/num/cid[d];
    }
    cout<<"recall: "<<endl;
    for(unsigned d = 0; d < 5; d++)
    {
        cout<<"top "<<cid[d]<<" : "<<rt[d]<<" ";
    }
    cout<<endl;

    delete []sid;
    delete []rtop;

    return rc;
}

int NNGKMeans::ITNNGraph()
{
    Timer *toptime = new Timer;
    toptime->start();

    this->topList = (unsigned*)calloc(this->count, TopK*sizeof(unsigned));
    /// init
    this->topDist = (float*)calloc(this->count, TopK*sizeof(float));
    this->radius = (float*)calloc(this->count, sizeof(float));

    /**/
    cout<<"init NNG\n";
    initNNGraph();
    cout<<"init end\n";
    Timer *ext0 = new Timer;
    ext0->start();
    cout<<"build NNG\n";
    buildNNGraph();
    cout<<"build time\t";
    ext0->end(true);
    toptime->end(true);

    free(topDist);
    free(radius);
    this->topDist = NULL;
    this->radius = NULL;

    return 0;
}

int NNGKMeans::clust(const unsigned int clust_num, const char *dstfn, const int verbose)
{
    unsigned long clustNum = this->clnumb, i, j, id;

    ITNNGraph();


    this->Es     = (double*)calloc(this->clnumb, sizeof(double));
    this->arrayD = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
    this->Ds     = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
    this->Ns     = (int*)calloc(this->clnumb, sizeof(int));
    for(i = 0; i < this->clnumb; i++)
    {
        memset(Ds+i*this->ndim, 0, sizeof(double)*this->ndim);
        memset(Ns+i, 0, sizeof(int));
    }

    cout<<"Initialize cluster begin .......\n";
    RandomPartition *xbk = new RandomPartition;
    xbk->buildcluster(this->data, this->count, this->ndim, "tmpfile.txt",  "rnd", "large", "i2", this->clnumb, 0);
    xbk->getLabel(this->labels);
    delete xbk;
    cout<<"initialize cluster finished all\n";
    Timer *mytm2 = new Timer();
    mytm2->start();
    this->optzI2(1, this->TopK, NIter);
    mytm2->end(true);
    this->calAVGDist(this->arrayD, this->Ns, clust_num, this->infoMap);
    //this->save_clust(dstfn);
    cout<<"sc "<<this->clnumb<<endl;

    return clustNum;
}

int NNGKMeans::NNI2(bool verbose, const unsigned iter0)
{
    if(iter0 < 1)
        return 0;
    double delta, delta0, optEg, tmpEg, tmpEg1 = 0, tmpEg2 = 0, len, allEg;
    unsigned long i = 0, j = 0, k = 0, iter = 0, id = 0, row = 0;
    int label, nlabel = 0;
    float *tmpData = NULL;
    bool UPDATE = false;
    unordered_set<unsigned> tmptop;
    allEg = 0;
    for(i = 0; i < this->count; i++)
    {
        allEg += lens[i];
        label = this->labels[i];
        Ns[label]++;
        row   = i*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            Ds[label*this->ndim+j] += data[row+j];
        }
    }
    for(i = 0; i < this->clnumb; i++)
    {
        Es[i] = 0;
        row   = i*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            Es[i] += this->Ds[row+j]*this->Ds[row+j];
        }
    }
    optEg = getI2(this->Ds, this->clnumb, this->ndim, Ns);
    iter = 0;
    if(verbose)
        cout<<"iter "<<iter<<"\toptEg "<<optEg<<"\tMSE "<<(allEg -optEg)/this->count<<endl;
    ///circle

    do
    {
        UPDATE = false;
        for(i = 0; i < this->count; i++)
        {
            ///only compare with active points
            id = i;
            label = labels[id];
            nlabel = label;
            /**/
            tmptop.clear();

            for(auto it = visitList[i].begin(); it != visitList[i].end(); it++)
            {
                if(label != labels[*it])
                    tmptop.insert(labels[*it]);
            }

            if(tmptop.empty())
                continue;
            /**/

            if(Ns[label] <= 1)
            {
                continue;
            }
            len = this->lens[id];
            tmpData = this->data + id*this->ndim;
            ///only compare with neiborhood cluster;
            delta0 = 0;
            tmpEg1 = 0;
            for(j = 0; j < ndim; j++)
            {
                tmpEg1 += Ds[label*this->ndim+j]*tmpData[j];
            }
            tmpEg1 = Es[label] - 2*tmpEg1 + len;
            for(auto it = tmptop.begin(); it != tmptop.end(); it++)
            {
                k = *it;
                if(k == label)
                    continue;
                tmpEg2 = 0;
                for(j = 0; j < ndim; j++)
                {
                    tmpEg2 += Ds[k*this->ndim+j]*tmpData[j];
                }
                tmpEg2 = Es[k] + 2*tmpEg2 + len;
                delta  = tmpEg2/(Ns[k]+1) - Es[k]/Ns[k] + tmpEg1/(Ns[label]-1) - Es[label]/Ns[label];
                if(delta > delta0)
                {
                    delta0 = delta;
                    nlabel = k;
                    tmpEg = tmpEg2;
                }
            }
            if(delta0 > 0)
            {
                UPDATE = true;
                ///update points belong to new cluster
                Ns[label]--;
                Ns[nlabel]++;
                for(j = 0; j < this->ndim; j++)
                {
                    Ds[label*this->ndim+j] -= tmpData[j];
                    Ds[nlabel*this->ndim+j] += tmpData[j];
                }
                Es[label]  = tmpEg1;
                Es[nlabel] = tmpEg;
                labels[id] = nlabel;
            }
            //  t3 += clock()-tts3;
        }

        iter++;
        optEg = getI2(this->Ds, this->clnumb, this->ndim, Ns);

        if(verbose)
            cout<<"iter "<<iter<<"\toptEg "<<optEg<<"\tMSE "<<(allEg -optEg)/this->count<<endl;


    }
    while(iter < iter0 && UPDATE);

    for(i = 0; i < this->clnumb; i++)
    {
        memcpy(this->arrayD+i*this->ndim, Ds + i*this->ndim, sizeof(double)*this->ndim);
    }

    return clnumb;
}

int NNGKMeans::optzI2(bool verbose, const unsigned NNTop, const unsigned iter0)
{
    if(iter0 < 1)
        return 0;
    double delta, delta0, optEg, tmpEg, tmpEg1 = 0, tmpEg2 = 0, len, allEg;
    unsigned long i = 0, j = 0, k = 0, iter = 0, id = 0, row = 0;
    int label, nlabel = 0;
    float *tmpData = NULL;
    bool UPDATE = false;
    unordered_set<unsigned> tmptop;
    allEg = 0;
    for(i = 0; i < this->count; i++)
    {
        allEg += lens[i];
        label = this->labels[i];
        Ns[label]++;
        row   = i*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            Ds[label*this->ndim+j] += data[row+j];
        }
    }
    for(i = 0; i < this->clnumb; i++)
    {
        Es[i] = 0;
        row   = i*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            Es[i] += this->Ds[row+j]*this->Ds[row+j];
        }
    }
    optEg = getI2(this->Ds, this->clnumb, this->ndim, Ns);
    iter = 0;
    if(verbose)
        cout<<"iter "<<iter<<"\toptEg "<<optEg<<"\tMSE "<<(allEg -optEg)/this->count<<endl;

    ///circle
    /**
    for(i = 0; i < this->clnumb; i++)
    {
        tmptop.insert(i);
    }
    /**/

    t1 = clock();
    t2 = clock() - t1;

    // os<<iter<<" "<<(allEg-optEg)/this->count<<" "<<t2/CLOCKS_PER_SEC<<endl;

    unsigned *rl = new unsigned[this->count];
    for(i = 0; i < this->count; i++)
    {
        rl[i] = i;
    }
    random_shuffle(rl, rl+this->count);

    do
    {

        UPDATE = false;
        for(i = 0; i < this->count; i++)
        {
            ///only compare with active points
            id = rl[i];
            label = labels[id];
            nlabel = label;
            /**/
            tmptop.clear();
            for(j = 0; j < NNTop; j++)
            {
                if(label != labels[topList[id*TopK+j]])
                    tmptop.insert(labels[topList[id*TopK+j]]);
            }
            if(tmptop.empty())
                continue;
            /**/

            if(Ns[label] <= 1)
            {
                continue;
            }
            len = this->lens[id];
            tmpData = this->data + id*this->ndim;
            ///only compare with neiborhood cluster;
            delta0 = 0;
            tmpEg1 = 0;
            for(j = 0; j < ndim; j++)
            {
                tmpEg1 += Ds[label*this->ndim+j]*tmpData[j];
            }
            tmpEg1 = Es[label] - 2*tmpEg1 + len;
            for(auto it = tmptop.begin(); it != tmptop.end(); it++)
            {
                k = *it;
                if(k == label)
                    continue;
                tmpEg2 = 0;
                for(j = 0; j < ndim; j++)
                {
                    tmpEg2 += Ds[k*this->ndim+j]*tmpData[j];
                }
                tmpEg2 = Es[k] + 2*tmpEg2 + len;
                delta  = tmpEg2/(Ns[k]+1) - Es[k]/Ns[k] + tmpEg1/(Ns[label]-1) - Es[label]/Ns[label];
                if(delta > delta0)
                {
                    delta0 = delta;
                    nlabel = k;
                    tmpEg = tmpEg2;
                }
            }
            if(delta0 > 0)
            {
                UPDATE = true;
                ///update points belong to new cluster
                Ns[label]--;
                Ns[nlabel]++;
                for(j = 0; j < this->ndim; j++)
                {
                    Ds[label*this->ndim+j] -= tmpData[j];
                    Ds[nlabel*this->ndim+j] += tmpData[j];
                }
                Es[label]  = tmpEg1;
                Es[nlabel] = tmpEg;
                labels[id] = nlabel;
            }
        }

        iter++;
        optEg = getI2(this->Ds, this->clnumb, this->ndim, Ns);

        if(verbose)
            cout<<"iter "<<iter<<"\toptEg "<<optEg<<"\tMSE "<<(allEg -optEg)/this->count<<endl;
        t2 = clock() - t1;


        //  os<<iter<<" "<<(allEg-optEg)/this->count<<" "<<t2/CLOCKS_PER_SEC<<endl;

    }
    while(iter < iter0 && UPDATE);


    for(i = 0; i < this->clnumb; i++)
    {
        memcpy(this->arrayD+i*this->ndim, Ds + i*this->ndim, sizeof(double)*this->ndim);
    }

    delete []rl;

    return clnumb;
}

int NNGKMeans::NNKM(bool verbose, const unsigned NNTop, const unsigned iter0)
{
    if(iter0 < 1)
        return 0;
    double delta, delta0, optEg, dst, tmpdst, len, allEg;
    unsigned long i = 0, j = 0, k = 0, iter = 0, id = 0, row = 0;
    int label, nlabel = 0;
    float *tmpData = NULL;
    bool UPDATE = false;
    unordered_set<unsigned> tmptop;
    allEg = 0;

    float *Cs = new float[this->clnumb*ndim];
    memset(Cs, 0, sizeof(float)*this->clnumb*ndim);

    for(i = 0; i < this->count; i++)
    {
        allEg += lens[i];
        label = this->labels[i];
        Ns[label]++;
        row   = i*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            Ds[label*this->ndim+j] += data[row+j];
        }
    }
    for(i = 0; i < this->clnumb; i++)
    {
        Cs[i] = 0;
        row   = i*this->ndim;
        if(this->Ns[i] == 0)
            continue;
        for(j = 0; j < this->ndim; j++)
        {
            Cs[row+j] = this->Ds[row+j]/this->Ns[i];
        }
    }
    optEg = getI2(this->Ds, this->clnumb, this->ndim, Ns);
    iter = 0;
    if(verbose)
        cout<<"iter "<<iter<<"\toptEg "<<optEg<<"\tMSE "<<(allEg -optEg)/this->count<<endl;


    t1 = clock();
    t2 = clock()-t1;
    os<<iter<<" "<<(allEg -optEg)/this->count<<t2/CLOCKS_PER_SEC<<endl;

    do
    {
        UPDATE = false;
        for(i = 0; i < this->count; i++)
        {
            ///only compare with active points
            id = i;
            label = labels[id];
            nlabel = label;
            /**/
            tmptop.clear();
            for(j = 0; j < NNTop; j++)
            {
                if(label != labels[topList[id*TopK+j]])
                    tmptop.insert(labels[topList[id*TopK+j]]);
            }
            if(tmptop.empty())
                continue;
            /**/
            if(Ns[label] <= 1)
            {
                continue;
            }
            tmpData = data + i*ndim;
            dst = 0;
            dst = PQMath::l2f(tmpData, 0, Cs, nlabel, ndim);
            ///only compare with neiborhood cluster;
            for(auto it = tmptop.begin(); it != tmptop.end(); it++)
            {
                k = *it;
                if(k == label)
                    continue;

                tmpdst = PQMath::l2f(tmpData, 0, Cs, k, ndim);

                if(dst > tmpdst)
                {
                    dst = tmpdst;
                    nlabel = k;
                    UPDATE = true;
                }
            }
            labels[i] = nlabel;
        }

        memset(Ds, 0, sizeof(double)*ndim*this->clnumb);
        memset(Cs, 0, sizeof(float)*ndim*clnumb);
        memset(Ns, 0, sizeof(int)*clnumb);
        for(i = 0; i < this->count; i++)
        {
            label = this->labels[i];
            Ns[label]++;
            row   = i*this->ndim;
            for(j = 0; j < this->ndim; j++)
            {
                Ds[label*this->ndim+j] += data[row+j];
            }
        }
        for(i = 0; i < this->clnumb; i++)
        {
            Cs[i] = 0;
            row   = i*this->ndim;
            if(this->Ns[i] == 0)
                continue;
            for(j = 0; j < this->ndim; j++)
            {
                Cs[row+j] = this->Ds[row+j]/this->Ns[i];
            }
        }

        iter++;
        optEg = getI2(this->Ds, this->clnumb, this->ndim, Ns);

        if(verbose)
            cout<<"iter "<<iter<<"\toptEg "<<optEg<<"\tMSE "<<(allEg -optEg)/this->count<<endl;

        t2 = clock()-t1;
        os<<iter<<" "<<(allEg -optEg)/this->count<<t2/CLOCKS_PER_SEC<<endl;

    }
    while(iter < iter0 && UPDATE);


    for(i = 0; i < this->clnumb; i++)
    {
        memcpy(this->arrayD+i*this->ndim, Ds + i*this->ndim, sizeof(double)*this->ndim);
    }

    delete [] Cs;
    Cs = NULL;

    return clnumb;
}


int NNGKMeans::buildNN(const char* srcfn, const char* dstfn)
{
    this->init(srcfn);
    this->initMemry(this->count, 0);
    this->ITNNGraph();
    this->saveKNNGraph(dstfn, this->TopK);

    return 0;
}

int NNGKMeans::buildCluster(unsigned* NN, float* data, unsigned row, unsigned ndim, unsigned clustnumb)
{
    unsigned i, j;
    this->init(data, row, ndim);
    this->clnumb = clustnumb;
    this->initMemry(row, clustnumb);
    this->copyFromTopList(NN, 25);
    this->Es     = (double*)calloc(this->clnumb, sizeof(double));
    this->arrayD = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
    this->Ds     = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
    this->Ns     = (int*)calloc(this->clnumb, sizeof(int));

    for(i = 0; i < this->clnumb; i++)
    {
        memset(Ds+i*this->ndim, 0, sizeof(double)*this->ndim);
        memset(Ns+i, 0, sizeof(int));
    }

    cout<<"Initialize cluster begin .......\n";
    RandomPartition *xbk = new RandomPartition;
    xbk->buildcluster(this->data, this->count, this->ndim, "tmpfile.txt",  "rnd", "large", "i2", this->clnumb, 0);
    xbk->getLabel(this->labels);
    delete xbk;
    cout<<"initialize cluster finished all\n";
    Timer *mytm2 = new Timer();
    mytm2->start();
    this->optzI2(1, this->TopK, NIter);
    mytm2->end(true);
    this->calAVGDist(this->arrayD, this->Ns, clnumb, this->infoMap);

    free(Es);
    free(arrayD);
    free(Ds);
    free(Ns);
    Es = NULL;
    arrayD = NULL;
    Ds = NULL;
    Ns = NULL;

    return 0;
}

int NNGKMeans::buildCluster(const char *srcfn, const char *nnfn, const char *dstfn, unsigned clustnumb)
{
    unsigned i;
    this->clnumb = clustnumb;
    this->init(srcfn);
    this->initMemry(this->count, clustnumb);
    this->Es     = (double*)calloc(this->clnumb, sizeof(double));
    this->arrayD = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
    this->Ds     = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
    this->Ns     = (int*)calloc(this->clnumb, sizeof(int));

    unsigned row;
    this->topList = IOAgent::loadNN(nnfn, row, TopK);

    if(row != this->count)
    {
        cout<<"erro nn file\n";
        exit(0);
    }

    //this->recall();

    /**/

    t1 = clock();

    for(i = 0; i < this->clnumb; i++)
    {
        memset(Ds+i*this->ndim, 0, sizeof(double)*this->ndim);
        memset(Ns+i, 0, sizeof(int));
    }

    cout<<"Initialize cluster begin .......\n";
    RandomPartition *xbk = new RandomPartition;
    xbk->buildcluster(this->data, this->count, this->ndim, "tmpfile.txt",  "non", "large", "i2", this->clnumb, 0);
    xbk->getLabel(this->labels);
    delete xbk;
    cout<<"initialize cluster finished all\n";
    Timer *mytm2 = new Timer();
    mytm2->start();
    this->optzI2(1, TopK, NIter);
    ///this->NNKM(1, TopK, NIter);
    mytm2->end(true);
    this->calAVGDist(this->arrayD, this->Ns, clnumb, this->infoMap);
    //this->save_clust(dstfn);
    /**/

    free(Es);
    free(arrayD);
    free(Ds);
    free(Ns);
    Es = NULL;
    arrayD = NULL;
    Ds = NULL;
    Ns = NULL;

    return 0;
}

int NNGKMeans::copy2TopList(unsigned *top)
{
    for(unsigned i = 0; i < this->count; i++)
    {
        memcpy(top+i*TopK, this->topList+i*TopK, sizeof(unsigned)*TopK);
    }
    return 0;
}

int NNGKMeans::copyFromTopList(unsigned *top, unsigned topNum)
{
    this->topList = (unsigned*)calloc(this->count, TopK*sizeof(unsigned));
    for(unsigned i = 0; i < this->count; i++)
    {
        memcpy(this->topList+i*TopK, top+i*TopK, sizeof(unsigned)*TopK);
    }
    return 0;
}

void NNGKMeans::saveCenters(const char *dstfn, bool append)
{
    unsigned long clabel = 0, i, j, loc, rCNum = 0;
    bool isNULL = false;
    ofstream *outStrm  = NULL;

    if(this->arrayD == NULL)
    {
        cout<<"empty arrayD!\n";
        isNULL = true;
        this->arrayD = (double*)calloc(this->clnumb, this->ndim*sizeof(double));
    }
    cout<<"save Centers1\n";
    for(i = 0; i < count; i++)
    {
        clabel = labels[i];
        for(j = 0; j < ndim; j++)
        {
            arrayD[clabel*ndim+j] += data[i*ndim+j];
        }
    }
    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        if(this->infoMap[clabel].n > 0)
        {
            rCNum++;
        }
    }

    if(rCNum == 0)
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

    if(isNULL)
    {
        free(arrayD);
        this->arrayD = NULL;
    }

    cout<<"done\n";
    return ;
}

void NNGKMeans::saveKNNGraph(const char *dstfn, unsigned int k0)
{
    ofstream *outStrm  = NULL;
    unsigned int i = 0, j = 0;
    outStrm = new ofstream(dstfn, ios::out);
    unsigned n = this->count;
    if(k0 > this->TopK)
        k0 = this->TopK;
    (*outStrm)<<n<<"\t"<<k0<<endl;
    for(i = 0; i < n; i++)
    {
        (*outStrm)<<i;
        for(j = 0; j < k0; j++)
        {
            (*outStrm)<<" "<<topList[i*TopK+j];
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
    cout<<"done\n";
}

void NNGKMeans::saveKNNGraph(const char *dstfn, const unsigned int n, const unsigned int k0, vector<unsigned int> &ids)
{
    ofstream *outStrm  = NULL;
    unsigned int i = 0, j = 0, id;
    outStrm = new ofstream(dstfn, ios::out);
    (*outStrm)<<n<<"\t"<<k0<<endl;
    for(i = 0; i < n; i++)
    {
        (*outStrm)<<ids[i];
        for(j = 0; j < k0; j++)
        {
            id = topList[i*TopK+j];
            (*outStrm)<<" "<<ids[id];
        }
        (*outStrm)<<endl;
    }

    outStrm->close();
    cout<<"done\n";
}

int NNGKMeans::fetchCenters(float *centers)
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

NNGKMeans::~NNGKMeans()
{
    if(this->arrayD != NULL)
    {
        free(arrayD);
        this->arrayD = NULL;
    }
    if(this->Ds != NULL)
    {
        free(Ds);
        Ds = NULL;
    }
    if(this->Ns != NULL)
    {
        free(this->Ns);
        this->Ns = NULL;
    }
    if(this->Es != NULL)
    {
        free(this->Es);
        this->Es = NULL;
    }
    if(this->topDist != NULL)
    {
        free(topDist);
        topDist = NULL;
    }
    if(!this->outNN&&this->topList != NULL)
    {
        free(topList);
        topList = NULL;
    }
}

void NNGKMeans::test()
{
    const char *srcMatFn1 = "/home/chenghaod/datasets/bignn/mat/sift1m/sift_base.txt";
    const char *srcMatFn2 = "/home/chenghaod/datasets/bignn/mat/sift_learn.txt";
    const char *srcMatFn3 = "itm_vlad_hesaff_flickr10m.itm";
    const char *dstfn = "/home/chenghaod/sift_learn_test.out";
    //const char *srcfn = "/home/chenghaod/dataset/bow/ox5k_nrsift_65k_hesaff_bow.txt";
    //const char *srcfn = "/home/chenghaod/src/finished/minhash/etc/data/holidays_nrsift_20k_hess_bow.mat";
    const char *srcfn = "sift_base.fvecs";
    //const char *srcfn = "/home/chenghaod/data/sift/raw/sift_base.fvecs";
    ///const char *srcfn = "/home/chenghaod/data/sift/raw/sift_learn.txt";
    ///const char *srcfn = "/home/chenghaod/datasets/gist/gist_base.fvecs";
    ///const char *srcfn = "/home/chenghaod/datasets/gist/gist_learn.fvecs";

    //paraTest("/home/chenghaod/data/sift/raw/sift_learn.txt", "tmpfile.txt");
    //exit(0);

    NNGKMeans *mykm = new NNGKMeans();
    mykm->buildcluster(srcfn, dstfn, "non", "large", "i2", 10000, false);
    delete mykm;
    cout<<"time\n";
    cout<<"t2 "<<t2<<" t3 "<<t3<<" t4 "<<t4<<" t5 "<<t5<<" t6 "<<t6<<" t7 "<<t7<<" t8 "<<t8<<endl;

}
