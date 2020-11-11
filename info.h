#ifndef INFO_H
#define INFO_H

#include <vector>
#include <string>



using namespace std;


struct INFO
{
    string pic_info;
    string shape_file;
    string words_npy;
    string he_npy;
    string geos_npy;


    vector<unsigned > uids;
    vector<unsigned > entry_ids;
    vector<string> urls;
    vector<unsigned > shapes;

    vector<unsigned > words;
    vector<unsigned > heCodes;
    vector<float > geos;

    //INFO();
    INFO(const char *pic_info,const char *shape_file, const char *words_npy, const char *he_npy, const char *geos_npy);

    unsigned get_uid(long i);
    unsigned get_entry_id(long i);
    string get_url(long i);
    unsigned *get_words();
    unsigned *get_hes();
    float *get_geos();
    unsigned *get_shapes();

    bool save(unsigned uid, unsigned entry_id, const char *url, unsigned shape, unsigned *word, unsigned *heCode,float *geo);
    void add(unsigned uid, unsigned entry_id, const char *url, unsigned shape, unsigned *word, unsigned *heCode,float *geo);
    void add(unsigned uid, unsigned entry_id, const char *url,  vector<unsigned >word, vector<unsigned >heCode,vector<float >geo);

    virtual ~INFO();

};








#endif // INFO_H
