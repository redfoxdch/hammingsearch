#include "info.h"
#include "cnpy/cnpy.h"
#include "loadSift.h"
#include <fstream>
#include <cassert>


INFO::INFO(const char *pic_info, const char *shape_file, const char *words_npy, const char *he_npy, const char *geos_npy)
{



    this->pic_info = pic_info;
    this->shape_file = shape_file;
    this->words_npy = words_npy;
    this->he_npy = he_npy;
    this->geos_npy = geos_npy;



    long  row;
    long row1;
    int col;
    float *fgeo;
    unsigned *uwords;
    unsigned *uhe;

//// load picture info


    this->shapes = loadSift::loadShape(shape_file,row1);
    fgeo = loadSift::loadGeo(geos_npy,row,col);
    cout<<"row "<<row<<" "<<row1<<endl;
    //assert(row == row1);
    this->geos.insert(this->geos.end(),fgeo, fgeo+row*2);
    delete []fgeo;
    uwords = loadSift::loadWords(words_npy,row1);
    cout<<"row row1"<<row<<" "<<row1<<endl;
    assert(row == row1);
    this->words.insert(this->words.end(),uwords,uwords+row);
    delete []uwords;
    uhe = loadSift::loadWords(he_npy,row1);
    cout<<"row row1"<<row<<" "<<row1<<endl;
    assert(row == row1);
    this->heCodes.insert(this->heCodes.end(),uhe,uhe+row*2);
    delete []uhe;


    this->uids.reserve(row);
    this->entry_ids.reserve(row);
    this->urls.reserve(row);
    ifstream is(pic_info);
    int tmp1,tmp2;
    string url;
    const char *tmp;
    cout<<"pic info "<<pic_info<<endl;
    while(is>>tmp1)
    {
        is>>tmp2;
        is>>url;
        tmp = url.c_str();
        this->uids.push_back(tmp1);
        this->entry_ids.push_back(tmp2);
        this->urls.push_back(tmp);
        //cout<<tmp1<<" "<<tmp2<<" "<<url<<endl;
    }

    cout<<"row "<<this->shapes.size()<<" "<<uids.size()<<endl;
    assert(this->shapes.size() == uids.size());
    is.close();
}

INFO::~INFO()
{
}


bool INFO::save(unsigned uid, unsigned entry_id, const char *url, unsigned shape, unsigned *word, unsigned *heCode,float *geo)
{

    cnpy::npy_save(this->words_npy.c_str(),word, {shape},"a");
    cnpy::npy_save(this->he_npy.c_str(),heCode, {shape,2},"a");

    loadSift::saveGeo(this->geos_npy.c_str(),geo, shape);

    ofstream os;
    os.open(this->pic_info.c_str(),ios::out|ios::app);
    os<<uid<<" "<<entry_id<<" "<<url<<endl;
    os.close();
    ofstream os2;
    os2.open(this->shape_file.c_str(),ios::out|ios::app);
    os2<<shape<<endl;
    os2.close();
    cout<<"saved in memory succesfull"<<endl;
    return 0;
}

unsigned INFO::get_uid(long i)
{
    return uids[i];
}
unsigned INFO::get_entry_id(long i)
{
    return entry_ids[i];
}
string INFO::get_url(long i)
{
    return urls[i];
}
unsigned *INFO::get_words()
{
    return reinterpret_cast<unsigned*>(&(this->words[0]));
}
unsigned *INFO::get_shapes()
{
    return reinterpret_cast<unsigned*>(&(this->shapes[0]));
}
unsigned *INFO::get_hes()
{
    return reinterpret_cast<unsigned*>(&(this->heCodes[0]));
}
float *INFO::get_geos()
{
    return reinterpret_cast<float*>(&(this->geos[0]));
}

void INFO::add(unsigned uid, unsigned entry_id, const char *url, unsigned shape, unsigned *word, unsigned *heCode,float *geo)
{
    cout<<"save w"<<endl;
    //cout<<this->words_npy<<" w "<<this->he_npy<<" h "<<this->geos_npy<<"  g "<<this->pic_info<<" sp "<<this->shape_file<<endl;
    this->uids.push_back(uid);
    this->entry_ids.push_back(entry_id);
    this->urls.push_back(url);
    this->shapes.push_back(shape);
    this->geos.insert(this->geos.end(),geo, geo+shape*2);
    this->words.insert(this->words.end(),word,word+shape);
    this->heCodes.insert(this->heCodes.end(),heCode,heCode+shape*2);
    this->save(uid,entry_id,url,shape,word,heCode,geo);
}
void INFO::add(unsigned uid, unsigned entry_id, const char *url,  vector<unsigned >word, vector<unsigned >heCode,vector<float >geo)
{
    unsigned shape = word.size();
    this->uids.push_back(uid);
    this->entry_ids.push_back(entry_id);
    this->urls.push_back(url);
    this->shapes.push_back(shape);
    this->words.insert(this->words.end(),word.begin(),word.end());
    this->heCodes.insert(this->heCodes.end(),heCode.begin(),heCode.end());
    this->geos.insert(this->geos.end(),geo.begin(),geo.end());
}
