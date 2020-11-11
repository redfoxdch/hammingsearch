#include <iostream>
#include "loadSift.h"
#include "cnpy/cnpy.h"
#include "missionagent.h"
using namespace std;

int main()
{
//
    string dir = "";
    const char *sift_npy = "sift/sift_last.npy";
    const char *sift_name = "siftfiles.txt";
    const char *geo_name = "kps.npy";
    const char *centersFile = "centers.txt";
    const char *wordsFile = "words.npy";
    const char *projectionf = "projection.npy";
    const char *medianf = "median.npy";
    const char *hefile = "he.npy";
    const char *geo_npy = "kps/kps_last.npy";
    const char *sz = "sift_shape.csv";

    string ms = dir + sift_name;
    string mz = dir + sz;
    string mg = dir + geo_name;
    string mc = dir + centersFile;
    string mw = dir + wordsFile;
    string sift_file = dir + sift_npy;
    string projectionfile = dir + projectionf;
    string medianfile = dir + medianf;
    string HEWords = dir + hefile;
//    double *sift ;
//    string file = "tmp.npy";
//    double b[2] = {0.5,0.3};
//    cnpy::NpyArray * arr = new cnpy::NpyArray;
//    //cnpy::npy_save(file,b, {2},"w");
//
//    double a[2] = {0.4,0.6};
//    cnpy::npy_save(file,a, {2},"a");
//    *arr = cnpy::npy_load(file);
//    sift = arr->data<double>();
//    for(int i = 0; i < arr->shape[0];i++)
//    {
//        cout<<sift[i]<<endl;
//    }
    //loadSift::test();
    MissionAgent::clusteringTest("data/sift_base.fvecs", "tmp.centers", 10000);
//    MissionAgent::genCenters(ms.c_str(),mc.c_str(), 20000);


    //MissionAgent::BruteForceSearch(ms.c_str(), mc.c_str(),mw.c_str());
//    MissionAgent::HEQuantization(ms.c_str(),projectionfile.c_str(),medianfile.c_str(),mw.c_str(),HEWords.c_str());
    //MissionAgent::findGEOTop("target/sift/","target/kps/", "target/sift_shape.csv", projectionfile.c_str(),medianfile.c_str(),mc.c_str(), mw.c_str(),HEWords.c_str() ,mg.c_str(), "data/sift_shape.csv", 11);
    cout << "Hello world!" << endl;
    return 0;
}
