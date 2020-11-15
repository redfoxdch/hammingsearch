#include "geohenn.h"
#include "bruteforcenn.h"
#include "he.h"
#include "loadSift.h"
#include <algorithm>
#include <cstring>
#include <thread>
#include <set>
#include <map>

using namespace std;

float *projection;
float *median;
float *centers;

/***descend order***/
bool pairCompare(pair<long,float> &a, pair<long, float> &b)
{
    return a.second>b.second;
}

GEOHENN::GEOHENN(string projectionfile, string medianfile, string centersfile, string pic_info, string szfile,
                 string wordsfile, string hesfile, string geosfile)
{
    cout<<"init geohenn......\n";

    this->distance = new DISTANCE::L2DistanceAVXr4<float>;
    unsigned NTHREADS = 0;

    int tmp1, tmp2;
    loadSift *mload = new loadSift();
    cout<<"proj\n";
    this->projection = mload->loadProjection(projectionfile.c_str(),tmp1,tmp2);
    cout<<"ct\n";
    this->centers = mload->loadCenters(centersfile.c_str(),tmp1,tmp2);
    cout<<"md\n";
    this->median = mload->loadMedian(medianfile.c_str(),tmp1,tmp2);
    cout<<"en\n";
    delete mload;
    cout<<"projection file , centersfile, median file : "<<projectionfile<<" "<<centersfile<<" "<<medianfile<<" \n";

    this->storeInfo = new INFO(pic_info.c_str(),szfile.c_str(), wordsfile.c_str(), hesfile.c_str(), geosfile.c_str());
    if(NTHREADS == 0)
        this->NTHREADS = thread::hardware_concurrency()-1;
    if(this->NTHREADS < 1)
        this->NTHREADS = 1;
    this->word = NULL;
    this->he = NULL;
    cout<<"threads num "<<this->NTHREADS<<endl;
    cout<<"file keep in memory "<<this->storeInfo->words_npy<<" w "<<this->storeInfo->he_npy<<" h "<<this->storeInfo->geos_npy<<"  g "<<this->storeInfo->pic_info<<" sp "<<this->storeInfo->shape_file<<endl;
    cout<<"init end\n";
}

GEOHENN::GEOHENN(float *projection, double *median, float *centers, const DISTANCE::Distance<float> *distance,unsigned NTHREADS)
{
    this->projection = projection;
    this->median = median;
    this->centers = centers;
    this->distance = distance;

    if(NTHREADS == 0)
        this->NTHREADS = thread::hardware_concurrency()-1;
    if(this->NTHREADS < 1)
        this->NTHREADS = 1;
    this->word = NULL;
    this->he = NULL;
    this->storeInfo = NULL;
}
GEOHENN::GEOHENN(float *projection, double *median, float *centers, const DISTANCE::Distance<float> *distance,const char *szfile,
                 const char * wordsfile, const char * hesfile, const char * geosfile,unsigned NTHREADS)
{
    this->projection = projection;
    this->median = median;
    this->centers = centers;
    this->distance = distance;

    this->storeInfo = new INFO("filelist.csv",szfile, wordsfile, hesfile, geosfile);
    this->NTHREADS = NTHREADS;
    if(NTHREADS == 0)
        this->NTHREADS = thread::hardware_concurrency()-1;
    if(this->NTHREADS < 1)
        this->NTHREADS = 1;
    this->word = NULL;
    this->he = NULL;
}

GEOHENN::GEOHENN()
{

    this->word = NULL;
    this->he = NULL;

}
void GEOHENN::miniTOPK( long st,const  unsigned *word, const  unsigned *he,const  float *geo,const   unsigned *wmap,  unsigned row,
                        const  unsigned *words, const  unsigned *hes, const  float *geos, const  unsigned * sz,  long num, float *score)
{
    long id = 0;
    for(long i = 0; i < num; i++)
    {
        //cout<<"id "<<st<<" "<<i<<" "<<(st+i)<<" rid "<<*(sz+i)<<endl;
        score[i] = PairScore(word,he,geo,wmap,row,words+id,hes+id*2,geos+id*2,*(sz+i));
        id += *(sz+i);
    }
}

pair<long, float> * GEOHENN::findTop(const float *sift,const float *geo,unsigned row, const unsigned *words,
                                     const unsigned *hes, const float *geos, const unsigned *sz, long num, int &TOPK, float &value)
{
    //cout<<"feature "<<sift[0]<<" "<<geo[0]<<" "<<row<<" "<<words[0]<<" "<<hes[0]<<" "<<geos[0]<<" "<<sz[0]<<endl;
    if(this->word != NULL)
    {
        delete []word;
        this->word = NULL;
    }
    if(this->he != NULL)
    {
        delete []he;
        this->he = NULL;
    }
    pair<long, float> *top;

    BruteForceNN myNN(this->centers,WORD_LEN,this->distance,this->NTHREADS);
    this->word = myNN.BOWQuantization(sift,row,COL,this->centers,WORD_LEN);
    HE myHE(this->projection,this->median, this->distance,this->NTHREADS);
    this->he = myHE.enCoding(sift, word, row);
    top = GEOHENN::findTop(word, he, geo, row, words, hes, geos, sz, num, TOPK, value,this->NTHREADS);

    return top;
}

pair<long, float> *GEOHENN::findTop(const float *sift, const float *geo, unsigned row, INFO *storeInfo, int &TOPK, float &value)
{

    BruteForceNN myNN(this->centers,WORD_LEN,this->distance,this->NTHREADS);
    cout<<"word0\n";
    this->word = myNN.BOWQuantization(sift,row,COL,this->centers,WORD_LEN);
    cout<<"word1\n";
    HE myHE(this->projection,this->median, this->distance,this->NTHREADS);
    this->he = myHE.enCoding(sift, word, row);
//    for(int i = 270; i < row; i++)
//    {
//    cout<<"word "<<word[i]<<endl;
//        cout<<"he "<<he[i*2]<<" "<<he[i*2+1]<<endl;
//        cout<<"geo "<<geo[i*2]<<" "<<geo[i*2+1]<<endl;
//    }
    cout<<"here\n";
    unsigned *words = storeInfo->get_words();
    unsigned *hes = storeInfo->get_hes();
    float *geos = storeInfo->get_geos();
    unsigned * sz = storeInfo->get_shapes();
    long num = storeInfo->shapes.size();
    cout<<"num "<<num<<endl;
    pair<long, float> *topScores;
    topScores = GEOHENN::findTop(word, he, geo, row, words, hes, geos, sz, num, TOPK,value,this->NTHREADS);

    return topScores;
}

pair<long, float> *GEOHENN::findTop(const  unsigned *word, const  unsigned *he, const float *geo, unsigned row, INFO *storeInfo, int &TOPK, float &value)
{
    unsigned *words = storeInfo->get_words();
    unsigned *hes = storeInfo->get_hes();
    float *geos = storeInfo->get_geos();
    unsigned * sz = storeInfo->get_shapes();
    long num = storeInfo->shapes.size();
    pair<long, float> *topScores;
    topScores = GEOHENN::findTop(word, he, geo, row, words, hes, geos, sz, num, TOPK,value,thread::hardware_concurrency());

    return topScores;
}

pair<long, float> * GEOHENN::findTop(const  unsigned *word, const  unsigned *he, const float *geo, unsigned row, const  unsigned *words,
                                     const  unsigned *hes, const  float *geos,  const unsigned * sz,  long num,  int &TOPK, float &value,unsigned NTHREADS)
{
    unsigned *wmap = new unsigned[WORD_LEN];
    memset(wmap,0,sizeof(unsigned)*WORD_LEN);
    for(int i = 0; i < row; i++)
    {
        wmap[word[i]] += 1;
    }

    long id = 0;

    long *vid = new long[num];
    float *score = new float[num];
    for(long i = 0; i < num; i++)
    {
        vid[i] = id;
        id += sz[i];
    }
    if(NTHREADS == 0)
        NTHREADS = thread::hardware_concurrency();
    //NTHREADS = 1;
    long len = num/NTHREADS;
    vector<std::thread> threads;
    for(id = 0; id < num; )
    {
        if (num - id < len)
            break;
        threads.push_back(thread(miniTOPK,id,word,he,geo,wmap,row, words+vid[id],hes+vid[id]*2,geos+vid[id]*2,sz+id, len, score+id));
        id += len;
    }
    for (auto& th : threads) th.join();

    miniTOPK(id,word, he, geo, wmap, row, words+vid[id], hes+vid[id]*2, geos+vid[id]*2, sz+id, num-id, score+id);

    pair<long, float> *scores = new pair<long, float>[num];
    for(long i = 0; i < num; i++)
    {
        scores[i].first = i;
        scores[i].second = score[i];
    }

    sort(scores, scores + num, pairCompare);
    value = scores[num/2].second;
    for(int i = 0; i < num; i++)
    {
        if(scores[i].second < value * 10)
            break;
        if(i >= TOPK)
            TOPK = i+1;

    }
    //cout<<"topk "<<TOPK<<" "<<value<<endl;

    pair<long, float> *topScores = new pair<long, float>[TOPK];

    for(int i = 0; i < TOPK; i++)
    {
        topScores[i].first = scores[i].first;
        topScores[i].second = scores[i].second;
    }

    delete []wmap;
    delete []vid;
    delete []score;
    delete []scores;
    return topScores;
}


float GEOHENN::PairScore(const unsigned *word, const unsigned *he,const float *geo, const unsigned *wmap, const int row, const unsigned *words, const unsigned *hes, const float *geos, const int srow)
{
//cout<<"2\n";
    set<unsigned> hit;
    unsigned w;
    for(int i = 0; i < srow; i++)
    {
        w = words[i];
        if (wmap[w] < 1 || wmap[w] > SPARSE_THR)
            continue;
        //cout<<"w1 "<<" "<<srow<<" "<<w<<endl;
        hit.insert(w);
    }
    map<unsigned, vector<unsigned> >hes1,hes2;
    // cout<<"start he\n";
    for(int i = 0 ; i < row; i++)
    {
        w = word[i];
        //cout<<"i"<<i<<" "<<he[i*2]<<" "<<he[i*2+1]<<endl;
        if(hit.find(w) != hit.end())
            hes1[w].push_back(i);
    }
    //cout<<"nhe\n\n";
    // cout<<"4\n";
    for(int i = 0 ; i < srow; i++)
    {
        w = words[i];
        // cout<<"w "<<w<<endl;
        //cout<<"i "<<i<<" "<<srow<<" "<<hes[i*2]<<" "<<hes[i*2+1]<<endl;
        if(hit.find(w) != hit.end())
            hes2[w].push_back(i);
    }
    float score = 0, tmpScore = 0;
    int *as = new int[ANGLE_LEN];
    int *ss = new int[SCALE_LEN];
    memset(as,0,sizeof(int)*ANGLE_LEN);
    memset(ss,0,sizeof(int)*SCALE_LEN);
    for(auto &id : hit)
    {
        //cout<<"id "<<id<<endl;
        if (hes2[id].size() > SPARSE_THR)
        {
            continue;
        }
        for(auto &h1 : hes1[id])
        {
            for(auto &h2 : hes2[id])
            {
                tmpScore = HE::heScore(he+h1*2,hes+h2*2);
                as[Geo::angleScore(geo[h1*2],geos[h2*2])] += 1;
                ss[Geo::scaleScore(geo[h1*2+1],geos[h2*2+1])] += 1;
            }
        }
    }
    int am = 0;
    int sm = 0;
    int ct = 0;
    for(int i = 0; i< ANGLE_LEN; i++)
    {
        if (as[i] > ct)
        {
            am  = i;
            ct = as[i];
        }
    }
    ct = 0;
    for(int i = 0; i< SCALE_LEN; i++)
    {
        if (ss[i] > ct)
        {
            sm  = i;
            ct = ss[i];
        }
    }
    for(auto &id : hit)
    {
        if (hes2[id].size() > SPARSE_THR)
        {
            continue;
        }
        for(auto &h1 : hes1[id])
        {
            for(auto &h2 : hes2[id])
            {
                if(Geo::angleScore(geo[h1*2],geos[h2*2]) != am)
                    break;
                if(Geo::scaleScore(geo[h1*2+1],geos[h2*2+1])!= sm)
                    break;
                tmpScore = HE::heScore(he+h1*2,hes+h2*2);
                score += tmpScore;
            }
        }
    }
    // cout<<"bfs "<<score<<endl;
    score = score/row/srow;
    //cout<<row<<" "<<srow<<endl;
    // cout<<"sc "<<score<<endl;
    delete []as;
    delete []ss;
    //exit(0);
    return score;
}

GEOHENN::~GEOHENN()
{

    if(this->word != NULL)
    {
        delete []word;
        this->word = NULL;
    }
    if(this->he != NULL)
    {
        delete []he;
        this->he = NULL;
    }
    if(this->storeInfo != NULL)
        delete this->storeInfo;
}





