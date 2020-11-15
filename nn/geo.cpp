#include "geo.h"
#include <cmath>

using namespace std;

Geo::Geo()
{

}

int Geo::angleScore(float angle1, float angle2)
{
    int dif;
    dif = int((angle1-angle2)/(360.0/ANGLE_LEN))%ANGLE_LEN;
    dif = (dif+ANGLE_LEN)%ANGLE_LEN;
    return dif;
}
int Geo::scaleScore(float scale1, float scale2)
{
    int dif;
    dif = int(log(scale1/scale2)/log(2))%SCALE_LEN;
    dif = (dif+ANGLE_LEN)%SCALE_LEN;
    return dif;
}

int Geo::geoScore(float angle1, float angle2, float scale1, float scale2)
{
    int dif1,dif2,dif;
    dif1 = int((angle1-angle2)/(360.0/ANGLE_LEN))%ANGLE_LEN;
    dif2 = int(log(scale1/scale2)/log(2))%SCALE_LEN;
    dif = dif1 + dif2;
    return dif;
}
void Geo::test()
{

}
Geo::~Geo()
{

}


