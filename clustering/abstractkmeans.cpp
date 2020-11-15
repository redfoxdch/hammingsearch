#include "abstractkmeans.h"
#include "util/pqmath.h"
#include "util/timer.h"

#include <iostream>
#include <iomanip>
#include <cassert>
#include <fstream>
#include <cstring>
#include <cmath>

#include "util/pqmath.h"


using namespace std;

const unsigned int AbstractKMeans::paraBound = 1024;
const float AbstractKMeans::smallVal0 = 0.0000001f;

AbstractKMeans::AbstractKMeans()
{
    this->ndim    = 0;
    this->labels  = NULL;
    this->data    = NULL;
    this->lens    = NULL;
    this->dataType = 0;
    this->seed    = _kpp_;
    this->infoMap = NULL;
    this->sdata.col = NULL;
    this->sdata.index = NULL;
    this->sdata.data = NULL;
    this->_INIT_  = false;
    this->_REFINE_= false;
    this->dstort  = 0;
    //  this->nThrd   = HardWareSetup::getCoreNum()/4;
    this->allEg   = 0;
    //  assert(this->nThrd > 0);
    strcpy(srcmatfn, "");
}

bool AbstractKMeans::initMemry(const unsigned int dNum, const unsigned int clustNum)
{
    unsigned int i = 0;
    this->lens   = (double*)calloc(this->count, sizeof(double));
    this->labels = (int*)calloc(this->count, sizeof(int));
    this->infoMap = new CLSInfo[this->clnumb];
    cout<<"init \n";

    if(kmMtd != _prc_ && kmMtd != _pqc_)
        if(kmMtd == _rbkmn_ /**/ || kmMtd == _hek_/**/ || kmMtd == _hec_ /**|| kmMtd == _xtkmn_ || /** kmMtd == _xbkmn_ /**|| kmMtd == _tkmn_**/)
        {
            if(!dataType)
            {
                this->normVects(this->data, this->ndim, this->count, this->lens);
            }
            else
                this->normVects(this->sdata, this->ndim, this->count, this->lens);
            for(i = 0; i < this->count; i++)
            {
                if(lens[i] == 0)
                {
                    labels[i] = -1;
                }
                else
                {
                    labels[i] = 0;
                }
            }
            cout<<"Normalize input vectors .......... ";
            cout<<"yes\n";

        }
        else
        {
            if(!dataType)
                allEg = this->normVects(this->data, this->lens, this->ndim, this->count);
            else
                this->normVects(this->sdata, this->lens, this->ndim, this->count);
            cout<<"Normalize input vectors .......... ";
            cout<<"no\n";
        }

    return true;
}

bool AbstractKMeans::label2cluster(vector<unsigned> *clusters)
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

    return true;
}

bool AbstractKMeans::getLabel(int *labels)
{
    memcpy(labels, this->labels, sizeof(int)*this->count);
    return 0;
}

unsigned int AbstractKMeans::buildcluster(const char *srcfn, const char *dstfn, const char *_seed_, const char *lg_first,
        const char *crtrn, const int num, bool _refine_)
{
    if(kmMtd == _kppkmn_)
        cout<<"Threads for calcu. distance ...... "<<this->nThrd<<endl;
    this->config(_seed_, lg_first, crtrn, 1);
    this->clnumb   = num;
    this->_INIT_   = this->init(srcfn);

    this->_REFINE_ = _refine_;
    char msgStr[64];
    Timer *mytm = new Timer();
    mytm->start();
    if(!this->_INIT_ || this->count == 0)
    {
        cout<<"erro : no row\n";
        return -1;
    }
    else if(this->count <= this->clnumb)
    {
        cout<<"too many cluster num:"<<" row = "<<this->count<<" clnumb = "<<this->clnumb<<endl;
        return -1;
    }
    else
    {
        this->initMemry(this->count, clnumb);
        cout<<"init \n";
        this->clnumb  = this->clust(clnumb, dstfn, 1);

    }
    sprintf(msgStr, "%s\t%s\t%s\t%s\t%d\t%lf\t", mthStr, _seed_, lg_first, crtrn, num, this->dstort);
    mytm->end(true);
    //mytm->end(msgStr, "./recordnw.txt");
    delete mytm;
    return this->clnumb;
}

unsigned int AbstractKMeans::buildcluster(float *mat, const int row, const int dim, const char *dstfn,
        const char *_seed_, const char *lg_first, const char *crtrn, const int num, bool _refine_)
{
    assert(mat);
    assert(crtrn);
    assert(lg_first);
    this->_REFER_ = true;
    this->config(_seed_, lg_first, crtrn, 0);
    this->clnumb = num;
    this->_INIT_ = this->init(mat, row, dim);
    this->_REFINE_ = _refine_;

    Timer *mytm = new Timer();
    mytm->start();
    if(!this->_INIT_ || this->count == 0)
    {
        cout<<"erro : no row\n";
        return -1;
    }
    else if(this->count <= this->clnumb)
    {
        cout<<"too many cluster num:"<<" row = "<<this->count<<" clnumb = "<<this->clnumb<<endl;
        return -1;
    }
    else
    {
        this->initMemry(this->count, clnumb);
        this->clnumb  = this->clust(clnumb, dstfn, 0);
    }

    mytm->end(true);

    delete mytm;

    return this->clnumb;
}

float *AbstractKMeans::getCluster(const unsigned int clabel0, unsigned int &row, unsigned int &dim)
{
    vector<int>::iterator it;
    vector<int> vects;
    unsigned int i = 0;
    int loc1 = 0, loc2 = 0, v = 0;
    unsigned int cpysize = sizeof(float)*this->ndim;

    for(i = 0; i < this->count; i++)
    {
        if(clabel0 == this->labels[i])
        {
            vects.push_back(i);
        }
    }
    row = vects.size();
    float *cl_data = new float[this->infoMap[clabel0].n*this->ndim];

    dim = this->ndim;
    assert(row == this->infoMap[clabel0].n);
    for(it = vects.begin(), i = 0; it != vects.end(); it++)
    {
        v = *it;
        loc1 = v*this->ndim;
        loc2 = i*this->ndim;
        memcpy(cl_data+loc2, this->data+loc1, cpysize);
        i++;
    }
    vects.clear();

    return cl_data;
}

double AbstractKMeans::getI2(const double *Ds, const unsigned int k, const unsigned int dim, const int *tmpns)
{
    double sumDst = 0, tmpE = 0;
    unsigned int i = 0, j = 0, loc = 0;

    for(j = 0; j < k; j++)
    {
        if(tmpns[j] <= 0)
            continue;

        loc = dim*j;
        tmpE = 0;

        for(i = 0; i < dim; i++)
        {
            tmpE += Ds[loc+i]*Ds[loc+i];
        }
        sumDst += tmpE/tmpns[j];
    }
    return sumDst;
}

double AbstractKMeans::calAVGDist(const double *Ds, const int *Ns, const unsigned int clustNum, CLSInfo *infos)
{
    double *dcts = new double[clustNum], tmpE = 0, sumEg = 0, tmpEg = 0;
    double sumDcts = 0;
    double *C    = new double[this->ndim];
    memset(dcts, 0, sizeof(double)*clustNum);
    unsigned int i, j = 0, loc = 0;
    int c = 0;
    for(i = 0; i < this->count; i++)
    {
        c   = labels[i];
        if(kmMtd == _rbkmn_)
        {
            dcts[c] += lens[i]*lens[i];
            sumDcts += lens[i]*lens[i];
        }
        else
        {
            dcts[c] += lens[i];
            sumDcts += lens[i];
        }
    }
    for(i = 0; i < clustNum; i++)
    {
        infos[i].n = Ns[i];

        if(Ns[i] < 1)
        {
            infos[i].E = 0;
            continue;
        }
        loc = i*this->ndim;
        tmpEg = 0;
        for(j = 0; j < this->ndim; j++)
        {
            C[j] = Ds[loc+j]/Ns[i];
            tmpEg += Ds[loc+j]*Ds[loc+j];
        }
        sumEg += tmpEg/Ns[i];
        tmpE   = PQMath::dvec_norm(C, this->ndim, 2);
        tmpE   = dcts[i] - infos[i].n*tmpE*tmpE;
        if(infos[i].n > 1)
        {
            infos[i].E = 2*(tmpE/(infos[i].n - 1));
        }
        else
        {
            infos[i].E = 0;
        }
    }

    this->dstort = (sumDcts - sumEg)/this->count;
    cout<<"E: "<<sumDcts<<"\t"<<sumEg<<"\t"<<"Distortion: "<<this->dstort<<endl;
    delete [] dcts;
    delete [] C;
    dcts = C = NULL;
    return this->dstort;
}

double AbstractKMeans::calAVGDist(const double *Ds, const int *Ns, const unsigned int clustNum)
{
    double *dcts = new double[clustNum], tmpE = 0, sumEg = 0, tmpEg = 0;
    double sumDcts = 0;
    double *C    = new double[this->ndim];
    memset(dcts, 0, sizeof(double)*clustNum);
    unsigned int i, j = 0, loc = 0;
    int c = 0;
    for(i = 0; i < this->count; i++)
    {
        c   = labels[i];
        if(kmMtd == _rbkmn_)
        {
            dcts[c] += lens[i]*lens[i];
            sumDcts += lens[i]*lens[i];
        }
        else
        {
            dcts[c] += lens[i];
            sumDcts += lens[i];
        }
    }
    for(i = 0; i < clustNum; i++)
    {
        if(Ns[i] < 1)
        {
            continue;
        }
        loc = i*this->ndim;
        tmpEg = 0;
        for(j = 0; j < this->ndim; j++)
        {
            C[j] = Ds[loc+j]/Ns[i];
            tmpEg += Ds[loc+j]*Ds[loc+j];
        }
        sumEg += tmpEg/Ns[i];
        tmpE   = PQMath::dvec_norm(C, this->ndim, 2);
        tmpE   = dcts[i] - Ns[i]*tmpE*tmpE;
    }

    this->dstort = (sumDcts - sumEg)/this->count;
    cout<<"E: "<<sumDcts<<"\t"<<sumEg<<"\t"<<"Distortion: "<<this->dstort<<endl;
    delete [] dcts;
    delete [] C;
    dcts = C = NULL;
    return this->dstort;
}


void AbstractKMeans::printClusters(const char *dstdir)
{
    unsigned int clabel = 0, i,  j, loc, counter = this->count;
    CLSInfo &crntinfo    = infoMap[0];
    FILE *fp = NULL;
    char matfn[1024];
    assert(dstdir);

    for(j = 0; j < this->clnumb; j++)
    {
        crntinfo = infoMap[j];
        sprintf(matfn, "%s%d.mat", dstdir, j);
        fp = fopen(matfn, "w");
        fprintf(fp,"%d %d\n", crntinfo.n, this->ndim);
        fclose(fp);
    }

    for(i = 0; i < this->count; i++)
    {
        clabel = labels[i];
        loc = i*this->ndim;

        sprintf(matfn, "%s%d.mat", dstdir, clabel);
        fp = fopen(matfn, "a");
        for(j = 0; j < this->ndim; j++)
        {
            fprintf(fp,"%f ", data[loc+j]);
        }
        fprintf(fp,"\n");
        fclose(fp);
        counter--;
        cout<<"\r\r"<<setw(8)<<counter;
    }
    cout<<endl;
}

void AbstractKMeans::printCluster(const char *dstfn)
{
    FILE *fp = fopen(dstfn, "w");
    if(fp == NULL)
    {
        cout<<"File '"<<dstfn<<"' cannot open for write!\n";
        return;
    }
    unsigned int i;

    for(i = 0; i < this->count; i++)
    {
        fprintf(fp, "%d\n", this->labels[i]);
    }

    fclose(fp);
    return ;
}

int AbstractKMeans::idx2clabel(const int i)
{
    if(this->labels == NULL)
        return 0;

    if(i < 0 ||i >= this->count || this->clnumb == 0)
    {
        return 0;
    }
    else
    {
        return this->labels[i];
    }
}

double AbstractKMeans::normVects(float *vects, double *lens,
                                 const unsigned int d0, const unsigned int n0)
{
    assert(vects);
    assert(lens);
    unsigned int i = 0, j = 0, loc = 0;
    double len = 0, sumEg1 = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = 0; j < d0; j++)
        {
            len += vects[loc+j]*vects[loc+j];
        }
        lens[i] = len;
        sumEg1  += lens[i];
        loc += d0;
    }
    return sumEg1;
}

void AbstractKMeans::normVects(const SparseMatrix &svects, double *lens,
                               const unsigned int d0, const unsigned int n0)
{
    assert(svects.data);
    assert(svects.index);
    assert(svects.col);
    assert(lens);
    unsigned int i = 0, j = 0;
    double len = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = svects.col[i]; j < svects.col[i+1]; j++)
        {
            len += svects.data[j]*svects.data[j];
        }
        lens[i] = len;
    }
}


void AbstractKMeans::normVects(SparseMatrix &svects, const unsigned int d0,
                               const unsigned int n0, double *lens)
{
    assert(svects.data);
    assert(svects.index);
    assert(svects.col);
    assert(lens);
    unsigned int i = 0, j = 0, loc = 0;
    float len = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = svects.col[i]; j < svects.col[i+1]; j++)
        {
            len += svects.data[j]*svects.data[j];
        }
        lens[i] = len;
        len     = sqrt(len);
        lens[i] = 1;
        if(len > AbstractKMeans::smallVal0)
            for(j = svects.col[i]; j < svects.col[i+1]; j++)
            {
                svects.data[j] = svects.data[j]/len;
            }
        else
        {
            lens[i] = 0;
        }

        loc += d0;
    }
    return ;
}

void AbstractKMeans::normVects(float *vects, const unsigned int d0,
                               const unsigned int n0, double *lens)
{
    assert(vects);
    assert(lens);
    cout<<"snorm\n";
    unsigned int i = 0, j = 0, loc = 0;
    float len = 0;
    for(i = 0; i < n0; i++)
    {
        len = 0;
        for(j = 0; j < d0; j++)
        {
            len += vects[loc+j]*vects[loc+j];
        }
        lens[i] = len;
        len     = sqrt(len);
        lens[i] = 1;
        if(len > AbstractKMeans::smallVal0)
        {
            for(j = 0; j < d0; j++)
            {
                vects[loc+j] = vects[loc+j]/len;
            }
        }
        else
        {
            lens[i] = 0;
        }

        loc += d0;
    }

    return ;
}

void AbstractKMeans::saveCenters(const char *dstfn, double *arrayD, bool append, bool _norm_)
{
    unsigned int clabel = 0, j, loc;
    ofstream *outStrm  = NULL;
    if(append)
    {
        outStrm = new ofstream(dstfn, ios::app);
    }
    else
    {
        outStrm = new ofstream(dstfn, ios::out);
        (*outStrm)<<this->clnumb<<" "<<this->ndim<<endl;;
    }

    for(clabel = 0; clabel < this->clnumb; clabel++)
    {
        loc  = clabel*this->ndim;
        for(j = 0; j < this->ndim; j++)
        {
            arrayD[loc+j] = arrayD[loc+j]/this->infoMap[clabel].n;
            ///cout<<this->arrayD[loc+j]<<" ";
        }

        /**/
        if(_norm_)
        {
            PQMath::l2_norm(arrayD+loc, ndim);
        }
        /**/

        for(j = 0; j < this->ndim; j++)
        {
            (*outStrm)<<arrayD[loc+j]<<" ";
        }
        (*outStrm)<<endl;
    }

    outStrm->close();

    return ;
}

bool AbstractKMeans::save_clust(const char *dstfn)
{
    if(dstfn == NULL || !strcmp(dstfn,""))
        return false;

    ofstream *outStrm = new ofstream(dstfn, ios::out);
    if(!outStrm->is_open())
    {
        cout<<"\nDestine file '"<<dstfn<<"' cannot open for write!\n";
        return false;
    }

    for(unsigned int i = 0; i < this->count; i++)
    {
        (*outStrm)<<this->labels[i]<<endl;
    }
    outStrm->close();
    return true;
}


void AbstractKMeans::resetlabels(unsigned int &n)
{
    unsigned int i;
    vector<unsigned > realid;
    map<int, unsigned> lbmap;
    unsigned int *lb = new unsigned[this->clnumb];
    memset(lb, 0, sizeof(unsigned)*(this->clnumb));
    for(i = 0; i < this->count; i++)
    {
        lb[labels[i]] = 1;
    }
    for(i = 0; i < this->clnumb; i++)
    {
        if(lb[i] != 0)
            realid.push_back(i);
    }
    for(i = 0; i < realid.size(); i++)
    {
        lbmap.insert(pair<unsigned, unsigned>(realid[i], i));
    }
    for(i = 0; i < this->count; i++)
    {
        labels[i] = lbmap[labels[i]];
    }
    n = lbmap.size();

    delete [] lb;
    lb = NULL;
}

AbstractKMeans::~AbstractKMeans()
{
    if(this->_REFER_)
    {
        this->data = NULL;
        this->sdata.data = NULL;
        this->sdata.index = NULL;
        this->sdata.col = NULL;
    }
    else
    {
        if(this->data != NULL)
        {
            free(this->data);
            this->data = NULL;
        }
        if(this->sdata.index != NULL)
        {
            delete [] this->sdata.index;
            delete [] this->sdata.col;
            this->sdata.index = NULL;
            this->sdata.col = NULL;

        }
        if(this->sdata.data != NULL)
        {
            delete [] this->sdata.data;
            this->sdata.data = NULL;
        }
    }
    if(this->infoMap != NULL)
    {
        delete [] this->infoMap;
        this->infoMap = NULL;
    }
    if(this->lens != NULL)
    {
        free(this->lens);
        this->lens = NULL;
    }
    if(this->labels != NULL)
    {
        free(this->labels);
        this->labels = NULL;
    }
}
