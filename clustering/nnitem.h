#ifndef NNITEM_H
#define NNITEM_H

#include <algorithm>
#include <vector>
#include <list>
#include <map>

using namespace std;

struct ValId
{
public :
    unsigned int idx;
    float val;
    static int LLcomparer(const ValId &a, const ValId &b)
    {
        return (a.val < b.val);
    }
    static int LGcomparer(const ValId &a, const ValId &b)
    {
        return (a.val > b.val);
    }
};

struct MiniNN
{
public:
    unsigned int idx;
    float val;

    /**ascending order**/
    static int LLcomparer(const MiniNN *a, const MiniNN *b)
    {
        return (a->val < b->val);
    }



    /**descending order**/
    static int LGcomparer(const MiniNN *a, const MiniNN *b)
    {
        return (a->val > b->val);
    }


};

class NNItem
{
public:
    NNItem(const unsigned int index, const double dist)
    {
        this->index = index;
        this->val   = dist;
        this->size  = 0;
        this->dvd   = true;
    }
    NNItem(const unsigned int index, const double dist, const int sz)
    {
        this->index = index;
        this->val   = dist;
        this->size  = sz;
        this->dvd   = true;
    }
    unsigned int index;
    double val;
    int   size;
    bool  dvd;

    /**ascending order**/
    static int LLcomparer(const NNItem *a,const NNItem *b)
    {
        return (a->val < b->val);
    }

    /**descending order**/
    static int LGcomparer(const NNItem *a, const NNItem *b)
    {
        return (a->val > b->val);
    }

    /**ascending order**/
    static int LLIDXcomparer(const NNItem *a, const NNItem *b)
    {

        return (a->index < b->index);
    }

    /**ascending order, consider the priority first**/
    static int LLVALcomparer(const NNItem *a, const NNItem *b)
    {
        if((a->dvd == false) && (b->dvd == false))
        {
            return false;
        }
        else if((a->dvd == false) && b->dvd)
        {
            return false;
        }
        else if(a->dvd && (b->dvd == false))
        {
            return true;
        }
        else
        {
            if(a->val < b->val)
            {
                return true;
            }
            else if(a->val == b->val)
            {
                if(a->size > b->size)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }

    /**descending order, consider the priority first**/
    static int LGVALcomparer(const NNItem *a, const NNItem *b)
    {
        if((a->dvd == false) && (b->dvd == false))
        {
            return true;
        }
        else if((a->dvd == false) && b->dvd)
        {
            return false;
        }
        else if(a->dvd && (b->dvd == false))
        {
            return true;
        }
        else
        {
            if(a->val > b->val)
            {
                return true;
            }
            else if(a->val == b->val)
            {
                if(a->size > b->size)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }

    /**descending order, consider the priority first**/
    static int LGSZcomparer(const NNItem *a, const NNItem *b)
    {
        if((a->dvd == false) && (b->dvd == false))
        {
            return true;
        }
        else if((a->dvd == false) && b->dvd)
        {
            return false;
        }
        else if(a->dvd && (b->dvd == false))
        {
            return true;
        }
        else
        {
            if(a->size > b->size)
            {
                return true;
            }
            else if(a->size == b->size)
            {
                if(a->val > b->val)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }
};


/******************************************/

#endif
