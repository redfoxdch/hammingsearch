#include "loadSift.h"
#include "cnpy/cnpy.h"
#include <iostream>
#include <cstring>
#include <sys/types.h>
#include <dirent.h>
#include <cstdio>
#include <fstream>



using namespace std;

int loadSift::scanFiles(string inputDirectory, vector<string> &fileList)
{

    DIR *p_dir;
    const char* str = inputDirectory.c_str();

    p_dir = opendir(str);
    if( p_dir == NULL)
    {
        cout<< "can't open :" << inputDirectory << endl;
    }

    struct dirent *p_dirent;

    while ( p_dirent = readdir(p_dir))
    {
        string tmpFileName = p_dirent->d_name;
        if( tmpFileName == "." || tmpFileName == "..")
        {
            continue;
        }
        else
        {
            cout<<tmpFileName<<endl;
            fileList.push_back(tmpFileName);
        }
    }
    closedir(p_dir);
    return fileList.size();
}

vector<string> loadSift::loadFileName(const char *fileName)
{
    ifstream is(fileName);
    string name;
    vector<string> filelist;
    while(is>>name)
    {
        filelist.push_back(name);
    }
    is.close();
    return filelist;
}

loadSift::loadSift()
{
    //std::cout<<"load sift\n";
}
loadSift::~loadSift()
{
    //std::cout<<"load sift end\n";
}

float* loadSift::loadFromDir(const char * dir, long &row, int &col)
{
    float* sift;
    vector<string> files;
    string file;
    scanFiles(dir, files);
    cnpy::NpyArray * arr;
    row = 0;
    for(int i = 0; i < files.size(); i++)
    {
        file = dir;
        file += files[i];
       // cout<<"load file ....... "<<file<<endl;
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        row += arr->shape[0];
        col = arr->shape[1];
        delete arr;
    }
    //cout<<"row "<<row<<endl;
    sift = new float[row*col];

    long pos = 0;
    int len = 0;
    long id = 0;
    float *tmp;
    for(int i = 0; i < files.size(); i++)
    {
        file = dir;
        file += files[i];
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        len = arr->shape[0];
        tmp = arr->data<float>();
        id = pos*col;
        for(unsigned j = 0; j < len; j++)
        {
            for(int dim = 0; dim < col; dim++)
                sift[id+j*col+dim] = tmp[j*col+dim];
        }
        arr = NULL;
        pos += len;
    }


   // std::cout<<"loaded\n";
    cout<<sift[12]<<endl;

    return sift;
}

float* loadSift::loadGeoDir(const char * dir, long &row, int &col)
{
    vector<string> files;
    string file;
    scanFiles(dir, files);
    cnpy::NpyArray * arr;
    row = 0;
    for(int i = 0; i < files.size(); i++)
    {
        file = dir;
        file += files[i];
       // cout<<"load file ....... "<<file<<endl;
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        row += arr->shape[0];
        col = arr->shape[1];
        delete arr;
    }
   // cout<<"row "<<row<<" "<<col<<endl;
    float *geo = new float[row*2];

    long pos = 0;
    int len = 0;
    long id = 0;
    double *sift;
    for(int i = 0; i < files.size(); i++)
    {
        file = dir;
        file += files[i];
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        len = arr->shape[0];
        sift = arr->data<double>();
        id = pos*2;
        for(unsigned j = 0; j < len; j++)
        {
            for(int dim = 0; dim < 2; dim++)
                geo[id+j*2+dim] = sift[j*col+dim+2];
        }
        delete arr;
        //delete []sift;
        arr = NULL;
        pos += len;
    }


    //std::cout<<"loaded\n";
    //delete []sift;
    return geo;
}

float* loadSift::loadGeoList(const char * filename, long &row, int &col)
{
    vector<string> files;
    string file;
    files = loadFileName(filename);
    cnpy::NpyArray * arr;
    row = 0;
    for(int i = 0; i < files.size(); i++)
    {
        ///file = "sift/";
        file = files[i];
       // cout<<"load geo file ....... "<<file<<endl;
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        row += arr->shape[0];
        col = arr->shape[1];
        delete arr;
    }
    //cout<<"row "<<row<<" "<<col<<endl;
    float *geo = new float[row*2];

    long pos = 0;
    int len = 0;
    long id = 0;
    double *sift;
    for(int i = 0; i < files.size(); i++)
    {
        file = files[i];
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        len = arr->shape[0];
        sift = arr->data<double>();
        id = pos*2;
        for(unsigned j = 0; j < len; j++)
        {
            for(int dim = 0; dim < 2; dim++)
                geo[id+j*2+dim] = sift[j*col+dim+2];
        }
        delete arr;
        //delete []sift;
        arr = NULL;
        pos += len;
    }


    //std::cout<<"loaded\n";
    //delete []sift;
    return geo;
}



float* loadSift::loadFromNameList(const char * filename, long &row, int &col,int MAX_LEN)
{
    float* sift;
    vector<string> files;
    string file;
    files = loadFileName(filename);
    cnpy::NpyArray * arr;
    row = 0;
    for(int i = 0; i < files.size(); i++)
    {
        ///file = "sift/";
        file = files[i];
       // cout<<"load sift file ....... "<<file<<endl;
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        row += arr->shape[0];
        col = arr->shape[1];
        delete arr;
        if (row > MAX_LEN)
            break;

    }
    //cout<<"size of pointer "<<sizeof(sift)<<endl;
    sift = new float[row*col];
    //cout<<"row "<<row<<endl;
    long pos = 0;
    long len = 0;
    float *tmp;
    long id = 0;
    for(int i = 0; i < files.size(); i++)
    {
        ///file = "sift/";
        file = files[i];
        arr = new cnpy::NpyArray;
        *arr = cnpy::npy_load(file);
        len = arr->shape[0];
        tmp = arr->data<float>();
        id = pos*col;
        for(unsigned j = 0; j < len; j++)
        {
            for(int dim = 0; dim < col; dim++)
                sift[id+j*col+dim] = tmp[j*col+dim];
        }
        delete arr;
        arr = NULL;
        pos += len;
        if (pos > MAX_LEN)
            break;
    }


    //std::cout<<"loaded\n";
    //cout<<sift[12]<<endl;
    return sift;
}


float* loadSift::loadFromFile(const char * file, long &row, int &col)
{
    float* sift;
    cnpy::NpyArray * arr = new cnpy::NpyArray;
    *arr = cnpy::npy_load(file);
    sift = arr->data<float>();
    row = arr->shape[0];
    col = arr->shape[1];
    //std::cout<<"loaded\n";
    //cout<<arr->word_size<<endl;
    return sift;
}

float* loadSift::loadGeo(const char *file, long &row, int &col)
{
    double* sift;
    cnpy::NpyArray * arr = new cnpy::NpyArray;
    *arr = cnpy::npy_load(file);
    sift = arr->data<double>();
    row = arr->shape[0];
    col = arr->shape[1];
    //cout<<"sp "<<sizeof(float)<<" "<<arr->word_size<<endl;
    float *geo = new float[row*2];
    //cout<<"col "<<col<<endl;
    for(long i = 0; i < row; i++)
    {
        for(int j= 0; j < 2; j++)
        {
            geo[i*2+j] = sift[i*col+j];
        }
    }
    //cout<<arr->word_size<<endl;
    delete arr;
    return geo;
}

void loadSift::saveGeo(const char *file,float *geo, unsigned &row)
{
    double* sift= new double[row*2];
    //cout<<"save geo\n";
    for(long i = 0; i < 2*row; i++)
    {
        sift[i] = geo[i];
    }
    cnpy::npy_save(file,sift,{row,2},"a");

    delete []sift;
}

float* loadSift::loadProjection(const char * file, int &row, int &col)
{
    float* projection;
    cnpy::NpyArray * arr = new cnpy::NpyArray;
    *arr = cnpy::npy_load(file);
    double *tmppro;
    tmppro = arr->data<double>();
    row = arr->shape[0];
    col  = arr->shape[1];
    projection = new float[arr->num_vals];
    for(int i = 0; i < row; i++)
    {
        for(int j= 0; j < col; j++)
        {
            projection[i*col+j] = tmppro[i*col+j];
        }
    }
    delete arr;
    return projection;
}
double* loadSift::loadMedian(const char * file, int &row, int &col)
{
    double* sift;
    cnpy::NpyArray * arr = new cnpy::NpyArray;
    *arr = cnpy::npy_load(file);
    sift = arr->data<double>();
    row = arr->shape[0];
    col = arr->shape[1];
    //std::cout<<"loaded\n";
    return sift;
}
float* loadSift::loadCenters(const char * file, int &row, int &col)
{
    float *centerMat = NULL;
    ifstream is(file);
    float ct;
    while(is>>ct)
    {
        row = ct;
        is>>ct;
        col = ct;
        centerMat = new float[row*col];
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                is>>ct;
                centerMat[i*col+j] = ct;
            }
        }
    }
    is.close();
    return centerMat;
}
unsigned* loadSift::loadWords(const char * file, long &row)
{
    unsigned* words;
    cnpy::NpyArray * arr = new cnpy::NpyArray;
    *arr = cnpy::npy_load(file);
    words = arr->data<unsigned>();
    row = arr->shape[0];
    //cout<<"words len "<<row<<endl;
    return words;
}

vector<unsigned> loadSift::loadShape(const char * file, long &row)
{
    vector<unsigned> sz ;
    ifstream is(file);
    unsigned ct;
    while(is>>ct)
    {
        sz.push_back(ct);
    }
    row = sz.size();
    is.close();

    return sz;
}

void loadSift::test()
{

    const char * file = "sift_last.npy";
    int row = 0, col = 0;
    loadSift *myload = new loadSift();
    cout<<"ld\n";
    myload->loadFileName("file.txt");
    cout<<"end\n";
    delete myload;


}
