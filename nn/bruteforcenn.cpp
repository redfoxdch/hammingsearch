#include "bruteforcenn.h"
#include "util/distance.hpp"
#include "cnpy/cnpy.h"
#include "loadSift.h"
#include <fstream>
#include <omp.h>
#include <thread>

using namespace std;
BruteForceNN::BruteForceNN()
{
}

BruteForceNN::BruteForceNN(const float *centers,unsigned crow,const DISTANCE::Distance<float>* d,unsigned NTHREADS)
{
    this->centers = centers;
    this->distance = d;
    this->NTHREADS = NTHREADS;
    this->crow = crow;
    if(NTHREADS == 0)
    {
        this->NTHREADS = thread::hardware_concurrency();
    }
    if(this->NTHREADS < 1)
    {
        this->NTHREADS = 1;
    }
}

BruteForceNN::~BruteForceNN()
{

}

pair<int, double> * BruteForceNN::buildNNG(float *data, long row, int col, int topK)
{

    pair<int, double> * KNNG;

    return KNNG;
}

void pal(const float *sift, long row, int col, const float *centers, int crow, unsigned *words,const DISTANCE::Distance<float>* distance)
{
    double dst = 0;
    double nst = RAND_MAX;
    long loc = 0;
    unsigned id = 0;
    //cout<<"row "<<row<<endl;
    for(unsigned i = 0; i < row; i++)
    {
        loc = col*i;
        nst = RAND_MAX;
        id = 0;
        for(unsigned j = 0; j < crow; j++)
        {
        //if(j > 19998)
        //cout<<"j "<<j<<endl;
            dst = distance->compare(sift+loc,centers+col*j, col);
            if (dst < nst)
            {
                nst = dst;
                id = j;
            }
        }
        //cout<<"words "<<i<<" "<<id<<endl;
        words[i] = id;
        //cout<<i<<" i"<<endl;
        if(i % 10000 == 9999)
            cout<<"finished compare i = "<<i<<endl;
    }
}

void BruteForceNN::BOWQuantization(const char *filename, const char *word_file) const
{
    vector<string> files;
    string file;
    files = loadSift::loadFileName(filename);
    long nrow;
    long row = 0;
    int col = 0;
    cnpy::NpyArray * arr;
    nrow = 0;
    for(int i = 0; i < files.size(); i++)
    {
        file = files[i];
        cout<<"load file ....... "<<file<<endl;
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        nrow += arr->shape[0];
        col = arr->shape[1];
        delete arr;
    }
    unsigned *pwords;
    unsigned *words = new unsigned[nrow];


    float *sift;
    pwords = words;

    long id  = 0;
    long len = row/NTHREADS;
    cout<<"start bow quantization \n";
    cout<<"compare threads num "<<NTHREADS<<endl;
    //cout<<row<<" "<<len<<endl;
    vector<std::thread> threads;
    for(unsigned i = 0; i < files.size(); i++)
    {
        file = files[i];
        sift = loadSift::loadFromFile(file.c_str(),row,col);

        threads.clear();
        len = row/NTHREADS;
        id = 0;
        //cout<<"row len nrow "<<row<<" "<<len<<" "<<nrow<<endl;
        if (len > 0)
        {
            for(id = 0; id < row; )
            {
                if (row - id < len)
                    break;
                threads.push_back(thread(pal,sift+id*col, len, col, this->centers, crow, pwords+id,this->distance));
                id += len;
            }
        }
        for (auto& th : threads) th.join();
        pal(sift+id*col, row-id, col, centers, crow, pwords+id,this->distance);
        pwords = pwords + row;
        delete []sift;
    }


    cnpy::npy_save(word_file,words, {nrow},"w");

    delete []words;
}

unsigned * BruteForceNN::BOWQuantization(const float *sift, long row, int col, const float *centers, int crow)const
{
    double dst = 0;
    double nst = RAND_MAX;
    long loc = 0;
    long id = 0;
    unsigned *words = new unsigned[row];
    cout<<"start bow quantization\n";

    int len = row/NTHREADS;
    cout<<"compare threads num "<<this->NTHREADS<<endl;
    vector<std::thread> threads;
    if (len > 0)
    {
        for(id = 0; id < row; )
        {
            if (row - id < len)
                break;
            threads.push_back(thread(pal,sift+id*col, len, col, centers, crow, words+id,this->distance));
            id += len;
        }
    }
    for (auto& th : threads) th.join();
    pal(sift+id*col, row-id, col, centers, crow, words+id,this->distance);
    cout<<"bow quantization end\n";
    return words;
}
void BruteForceNN::test()
{
    DISTANCE::Distance<float> * d;
    d = new DISTANCE::FastL2Distance<float>;
//    d = new DISTANCE::L2DistanceAVXr4<float>;
    float *mat = NULL;
    long row = 0;
    int col = 0, crow = 0;
    float * centers;
    ifstream is("data/centers.txt");
    float ct;
    while(is>>ct)
    {
        crow = ct;
        is>>ct;
        col = ct;
        centers = new float[crow*col];
        cout<<"crow "<<crow<<" col "<<col<<endl;
        for (int i = 0; i < crow; i++)
        {
            for (int j = 0; j < col; j++)
            {
                is>>ct;
                centers[i*col+j] = ct;
            }
        }
    }
    is.close();
    unsigned *words;
    const char *name = "data/siftfiles.txt";
    loadSift *mload = new loadSift();
    mat = mload->loadFromNameList(name,row,col);
    BruteForceNN *mnn = new BruteForceNN(centers,crow,d);
    words = mnn->BOWQuantization(mat, row, col, centers, crow);
    ofstream os("data/words.txt");
    for(int i = 0; i < row; i++)
    {
        os<<words[i]<<"\n";
    }
    os.close();
    delete []centers;
    delete []mat;
    delete mload;
}
