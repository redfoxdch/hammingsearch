#ifndef GEO_H
#define GEO_H

#define ANGLE_LEN 8
#define SCALE_LEN 8


class Geo
{


public :

    Geo();

    static int angleScore(float angle1, float angle2);
    static int scaleScore(float scale1, float scale2);

    static int geoScore(float angle1, float angle2, float scale1, float scale2);

    static void test();
    virtual ~Geo();
};


#endif // GEO_H
