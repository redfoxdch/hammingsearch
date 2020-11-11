#include "missionagent.h"
#include "loadSift.h"
#include "clustering/nngkmeans.h"
#include "hammingEmbedding/he.h"
#include "nn/bruteforcenn.h"
#include "cnpy/cnpy.h"
#include "nn/geohenn.h"
#include "info.h"
#include <fstream>
#include <iostream>
#include <tuple>

using namespace std;


bool MissionAgent::genCenters(const char *name, const char * centers, int clustNum)
{

    float *mat = NULL;
    long row = 0;
    int dim = 0;
    loadSift *mload = new loadSift();
    mat = mload->loadFromNameList(name,row,dim,10000000);
    ///mat = mload->loadFromFile(sift_npy,row,dim);
    ///cout<<mat[0]<<endl;
    AbstractKMeans * myKmeans = new NNGKMeans();
    myKmeans->buildcluster(mat, row, dim, "tmp.txt", "non", "large", "i2", clustNum, false);
    myKmeans->saveCenters(centers, 0);

    delete mload;
    return 0;
}

bool MissionAgent::findGEOTop(const char *siftfile,const char *geofile, const char *shapefile, const char * projectionfile, const char * medianfile,
                              const char * centersfile, const char * wordsfile, const char * hesfile, const char *geosfile, const char *szfile, int TOPK)
{
    float *sift;
    float *projection;
    double *median;
    float *centers;
    long row;
    int col;
    float *geo = NULL;
    vector<unsigned>shape;
    vector<unsigned>shapes;
    float *centerMat = NULL;
    long crow = 0;
    DISTANCE::Distance<float> * d;
    //d = new DISTANCE::FastL2Distance<float>;
    d = new DISTANCE::L2DistanceAVXr4<float>;


    int tmp1,tmp2;
    loadSift *mload = new loadSift();
    sift = mload->loadFromDir(siftfile,crow,col);
    geo = mload->loadGeoDir(geofile,crow,tmp1);
    shape = mload->loadShape(shapefile,crow);


    ///load centers end

    ///load sift feature

    cout<<"load all end\n";
    cout<<shape.size()<<endl;

    cout<<"shape file\n";
    GEOHENN *mgnn = new GEOHENN(projectionfile, medianfile, centersfile, "filelist.csv",szfile, wordsfile, hesfile, geosfile);
    ifstream is("filelist.csv");
    vector<string > filemap;
    string tmp;
    while(is>>tmp)
    {
        is>>tmp;
        is>>tmp;
        filemap.push_back(tmp);
    }
    is.close();
    string *mtop;
    int id = 0;
    ofstream os("result.txt");
    for(int i= 0 ; i< 128*2; i++)
    {
        //   sift[i] = 1;
    }
    for(int i= 0; i < 2*2; i++)
    {
        //   geo[i] = 1.2;
    }
    //shape[0] = 2;
    for(int j = 0; j < 1; j++)
    {
        cout<<"start match \n";
        mtop = matchANDupdate(j,j,"sp(2).jpg",sift+id*128,geo+id*2,shape[j],mgnn,TOPK);
        //id += shape[j];
        cout<<"id "<<id<<endl;
        cout<<"find id "<<j<<".........."<<endl;
        for(int i = 0; i < TOPK; i++)
        {
            cout<<mtop[i]<<endl;
            os<<i<<" "<<mtop[i]<<endl;
        }
        cout<<mtop[TOPK]<<endl;
        os<<endl;
        delete []mtop;
    }

    delete mgnn;

    os.close();
    return 0;
}

string * MissionAgent::matchANDupdate(unsigned uid, unsigned entry_id, string url,
                                      float *sift,float *geo, unsigned row,GEOHENN *mgnn, int &TOPK)
{
    pair<long, float> *mtop;
    cout<<"start find topk ...\n";
    float value;
    mtop = mgnn->findTop(sift,geo,row,TOPK, value);
    cout<<"find topk finished\n";

    unsigned *word = mgnn->getWord();
    unsigned *he = mgnn->getHE();
    //cout<<"url "<<url<<endl;
    mgnn->storeInfo->add(uid,entry_id,url.c_str(),row,word,he,geo);
    tuple<unsigned,unsigned,string, float> *topInfo = new tuple<unsigned,unsigned,string, float>[TOPK];
    string *s = new string[TOPK+1] ;
    std::ostringstream ss;

    for(int i = 0; i < TOPK; i++)
    {
        get<0>(topInfo[i]) = mgnn->storeInfo->get_uid(mtop[i].first);
        get<1>(topInfo[i]) = mgnn->storeInfo->get_entry_id(mtop[i].first);
        get<2>(topInfo[i]) = mgnn->storeInfo->get_url(mtop[i].first);
        cout<<" top "<<i<<" "<<mtop[i].first<<" "<<mgnn->storeInfo->get_url(mtop[i].first)<<"  "<<mtop[i].second<<endl;
        get<3>(topInfo[i]) = mtop[i].second;
        ss.clear();
        ss.str("");
        ss<<mgnn->storeInfo->get_uid(mtop[i].first);
        s[i] += ss.str();
        ss.clear();
        ss.str("");
        ss<<mgnn->storeInfo->get_entry_id(mtop[i].first);
        s[i] += " "+ ss.str();
        ss.clear();ss.str("");
        ss<<mgnn->storeInfo->get_url(mtop[i].first);
        s[i] += " "+ss.str();
        ss.clear();ss.str("");
        ss<<mtop[i].second;
        s[i] += " "+ss.str();
        //s[i]+= to_string(mgnn->storeInfo->get_uid(mtop[i].first));
        //s[i]+=" "+ to_string(mgnn->storeInfo->get_entry_id(mtop[i].first));
        //s[i]+= " "+mgnn->storeInfo->get_url(mtop[i].first);
        //s[i]+=" ";
        //s[i]+= to_string(mtop[i].second);
    }
    ss.clear();
    ss.str("");
    ss<<value;
    s[TOPK] = ss.str();
    delete []mtop;



    return s;
}


bool MissionAgent::HEQuantization(const char *sift_npy,const char *projectionfile,const char *medianfile, const char *wordfile,const char *HEWords)
{
    float *projection;
    double *tmppro;
    double *median;
    unsigned *words;
    unsigned *he;
    long row;
    int col;
    DISTANCE::Distance<float> * d;
    //d = new DISTANCE::FastL2Distance<float>;
    d = new DISTANCE::L2DistanceAVXr4<float>;
    float *sift = NULL;
    long crow = 0;
    int tmp1,tmp2;

    ///load sift feature
    loadSift *mload = new loadSift();
    projection = mload->loadProjection(projectionfile,tmp1,tmp2);
    cout<<"load projection\n";

    words = mload->loadWords(wordfile,crow);
    cout<<"w\n";

    HE *gm = new HE(projection,d);
    cout<<"gen median .....\n";
    gm->genMedian(sift_npy,wordfile,medianfile);
    cout<<"median gened \n";
    delete gm;

    median = mload->loadMedian(medianfile,tmp1,tmp2);
    cout<<"m\n";

    HE *myHE = new HE(projection, median,d);
    myHE->enCoding(sift_npy,wordfile,HEWords);



    delete mload;
    delete myHE;
    return 0;
}

bool MissionAgent::nn(const char *sift_npy, const char *centers, const char *ALLwords, const char *ALLHEWords)
{
    return 0;
}

bool MissionAgent::BruteForceSearch(const char *name, const char * centers, const char *words)
{
    DISTANCE::Distance<float> * d;
    //d = new DISTANCE::FastL2Distance<float>;
    d = new DISTANCE::L2DistanceAVXr4<float>;
    float *centerMat = NULL;
    long row = 0;
    int col = 0, crow = 0;
    ifstream is(centers);
    float ct;
    while(is>>ct)
    {
        crow = ct;
        is>>ct;
        col = ct;
        centerMat = new float[crow*col];
        cout<<"crow "<<crow<<" col "<<col<<endl;
        for (int i = 0; i < crow; i++)
        {
            for (int j = 0; j < col; j++)
            {
                is>>ct;
                centerMat[i*col+j] = ct;
            }
        }
    }
    is.close();
    loadSift *mload = new loadSift();
    BruteForceNN *bfnn;
    bfnn = new BruteForceNN(centerMat,crow,d);
    bfnn->BOWQuantization(name,words);
    delete bfnn;

    //cnpy::npy_save(words,wordsMat, {row},"w");

    delete []centerMat;
    centerMat = NULL;
    delete mload;
}




