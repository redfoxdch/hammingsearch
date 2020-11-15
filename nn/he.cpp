#include "he.h"
#include "cnpy/cnpy.h"
#include "loadSift.h"
#include <bitset>
#include <thread>
#include <algorithm>
#include <cstring>

using namespace std;

float HE64_WGH[65] = {64.000000, 57.977632, 52.976939, 48.583169, 44.626691,
                      41.014689, 37.688568, 34.607987, 31.743503, 29.072698,
                      26.577953, 24.245070, 22.062373, 20.020098, 18.109963,
                      16.324856, 14.658598, 13.105768, 11.661558, 10.321661,
                      9.082177, 7.939533, 6.890407, 5.931667, 5.060308, 4.273387,
                      3.567956, 2.941002, 2.389365, 1.909665, 1.498218, 1.150956,
                      0.863353, 0.630372, 0.446452, 0.305550, 0.201274, 0.127087,
                      0.076601, 0.043900, 0.023831, 0.012212, 0.005889, 0.002664,
                      0.001128, 0.000445, 0.000164, 0.000056, 0.000018, 0.000005,
                      0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000
                     };

void HE::test()
{


}

double *HE::transpose(const double *median, int row, int col)
{
    double *mmedian = new double[col*row];
    for(int i = 0 ; i < row; i++)
    {
        for(int j = 0; j < col; j++)
        {
            mmedian[j*row+i] = median[i*col+j];
        }
    }


    return mmedian;
}
void HE::miniCoding(const float *sift, const unsigned *words, long row, unsigned *heCodes, const DISTANCE::Distance<float> * distance, const float *projection,const double *median)
{
    bitset<32> bs;
    /***projection***/
    float *psift = new float[LEN];
    long loc = 0,loc1= 0;
    for(long i = 0; i < row; i++)
    {
        loc = i*COL;
        for(int j = 0; j < LEN; j++)
        {
            psift[j] = distance->dot(sift+loc,projection+j*COL, COL);
        }
        for(int j = 0; j < 32; j++)
        {
            if(psift[j] > median[words[i]*LEN+j])
            {
                bs[31-j] = 1;
            }
            else
            {
                bs[31-j] = 0;
            }
        }
        heCodes[i*2] = bs.to_ulong();
        for(int j = 32; j < LEN; j++)
        {
            if(psift[j] > median[words[i]*LEN+j])
                bs[63-j] = 1;
            else
                bs[63-j] = 0;
        }
        heCodes[i*2+1] = bs.to_ulong();
        if(i%10000== 9999)
            cout<<"encoding in i = "<<i<<endl;
    }

    delete []psift;
}


void HE::enCoding(const char *filename, const char *word_file, const char *dst) const
{
    vector<string> files;
    string file;
    files = loadSift::loadFileName(filename);
    long nrow;
    unsigned *pwords;
    unsigned *words = loadSift::loadWords(word_file,nrow);
    long row = 0;
    int col = 0;
    unsigned *pheCodes;
    unsigned *heCodes = new unsigned[nrow*2];
    float *sift;
    pwords = words;
    pheCodes = heCodes;

    long id  = 0;
    long len = row/NTHREADS;
    vector<std::thread> threads;
    for(unsigned i = 0; i < files.size(); i++)
    {
        file = files[i];
        sift = loadSift::loadFromFile(file.c_str(),row,col);

        len = row/NTHREADS;
        threads.clear();
        id = 0;
        if (len > 0)
            for(id = 0; id < row; )
            {
                if (row - id < len)
                    break;
                threads.push_back(thread(miniCoding,sift+id*COL, pwords+id,len,pheCodes+id*2,this->distance,this->projection,this->median));
                id += len;
            }
        for (auto& th : threads) th.join();
        miniCoding(sift+id*COL, pwords+id,row-id, pheCodes+id*2, this->distance,this->projection,this->median);
        delete []sift;
        pwords = pwords + row;
        pheCodes = pheCodes + row*2;
    }
    delete []words;

    cnpy::npy_save(dst,heCodes,{nrow,2},"w");


    delete []heCodes;

}




void HE::miniAdd(const float *sift, const unsigned *words, long row, vector<double> *median, const DISTANCE::Distance<float> * distance, const float *projection)
{
    /***projection***/
    float *psift = new float[LEN];
    long loc = 0,loc1= 0;
    for(long i = 0; i < row; i++)
    {
        loc = i*COL;

        for(int j = 0; j < LEN; j++)
        {
            psift[j] = distance->dot(sift+loc,projection+j*COL, COL);
        }
        for(int j = 0; j < LEN; j++)
        {
            median[words[i]*LEN+j].push_back(psift[j]);
        }
        if(i%10000== 9999)
            cout<<"sum in i = "<<i<<endl;
    }

    delete []psift;
}


void HE::genMedian(const char *filename, const char *word_file, const char *dst) const
{
    vector<string> files;
    string file;
    files = loadSift::loadFileName(filename);
    long nrow;
    unsigned *pwords;
    unsigned *words = loadSift::loadWords(word_file,nrow);
    long row = 0;
    int n  = 0;
    int col = 0;
    double *median = new double[LEN*WORD_LEN];

    float *sift;
    pwords = words;
    long id  = 0;
    long len = row/NTHREADS;
    vector<double> *pmedian = new vector<double>[LEN*WORD_LEN*NTHREADS];
    int j = 0;
    vector<std::thread> threads;
    for(unsigned i = 0; i < files.size(); i++)
    {
        file = files[i];
        cout<<"file "<<file<<endl;
        sift = loadSift::loadFromFile(file.c_str(),row,col);

        len = row/NTHREADS;
        threads.clear();
        id = 0;
        j = 0;
        if (len > 0)
            for(id = 0; id < row; )
            {
                if (row - id < len)
                    break;
                threads.push_back(thread(miniAdd,sift+id*COL, pwords+id,len,pmedian+j*WORD_LEN*LEN,this->distance,this->projection));
                id += len;
                j += 1;
            }
        for (auto& th : threads) th.join();
        miniAdd(sift+id*COL, pwords+id,row-id, pmedian, this->distance,this->projection);
        delete []sift;
        pwords = pwords + row;
    }
    cout<<"combine\n";

    vector<double> tmp;
    for(j = 0; j < WORD_LEN; j++)
    {
        for(int dim = 0; dim < LEN; dim++)
        {
            tmp.clear();
            for(int s = 0; s < NTHREADS; s++)
            {
                tmp.insert(tmp.end(),pmedian[s*LEN*WORD_LEN+j*LEN+dim].begin(), pmedian[s*LEN*WORD_LEN+j*LEN+dim].end());
            }
            sort(tmp.begin(),tmp.end());
            median[j*LEN+dim] = tmp[tmp.size()/2];
        }
    }

    delete []words;
    delete []pmedian;


    cnpy::npy_save(dst,median, {WORD_LEN,LEN},"w");


    delete []median;

}





unsigned *HE::enCoding(const float *sift, const unsigned *words, long row) const
{
    unsigned *heCodes = new unsigned[row*2];

    long len = row/NTHREADS;
    vector<std::thread> threads;
    long id = 0;
    if (len > 0)
        for(id = 0; id < row; )
        {
            if (row - id < len)
                break;
            threads.push_back(thread(miniCoding,sift+id*COL, words+id,len,heCodes+id*2,this->distance,this->projection,this->median));
            id += len;
        }
    for (auto& th : threads) th.join();
    miniCoding(sift+id*COL, words+id,row-id, heCodes+id*2, this->distance,this->projection,this->median);

    return heCodes;
}

float HE::heScore(const unsigned *a, const unsigned *b)
{
    float score;
    int num;
    bitset<32> bit1(a[0]^b[0]);
    bitset<32> bit2(a[1]^b[1]);

    num = bit1.count()+bit2.count();
    score = HE64_WGH[num]/64.0;

    return score;
}

float HE::heScore(const pair<unsigned,unsigned> a, const pair<unsigned,unsigned> b)
{
    float score;
    int num;

    bitset<32> bit1(a.first^b.first);
    bitset<32> bit2(a.second^b.second);

    num = bit1.count()+bit2.count();
    score = HE64_WGH[num]/64.0;

    return score;
}
