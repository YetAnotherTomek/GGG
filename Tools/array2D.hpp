// ~~~~~~~~~~~~~~~~~~~~~~
// a r r a y 2 D . h p p 
// ~~~~~~~~~~~~~~~~~~~~~~

#ifndef _ARRAY2D_HPP
#define _ARRAY2D_HPP

#include<vector>

using namespace std;

template <class type>
class array2D {
        public:
                array2D(void);
                inline void init(const unsigned int r,const unsigned int c,const type v);
                inline void free(void);
                inline unsigned int rows(void);
                inline unsigned int cols(void);
                inline unsigned int size(void);
                inline void add(const unsigned int i,const unsigned int j,const type v);
                inline void to(const unsigned int i,const unsigned int j,const type v);
                inline type at(const unsigned int i,const unsigned int j);
        private:
                unsigned int no_rows;
                unsigned int no_cols;
                vector<type> values;
};

template <class type>
inline array2D<type>::array2D(void) {}

template <class type>
inline void array2D<type>::init(const unsigned int r,const unsigned int c,const type v)
{
        no_rows = r;
        no_cols = c;
        values.clear();
        for(unsigned int i = 0; i < no_rows*no_cols; i++) values.push_back(v);
}

template <class type>
inline void array2D<type>::free(void)
{//test success
        values.~vector<type>;
}

template <class type>
inline unsigned int array2D<type>::rows(void)
{
        return no_rows;
}

template <class type>
inline unsigned int array2D<type>::cols(void)
{
        return no_cols;
}

template <class type>
inline unsigned int array2D<type>::size(void)
{
        return no_rows*no_cols;
}

template <class type>
inline void array2D<type>::add(const unsigned int i, const unsigned int j, const type v)
{
        values.at(i*no_cols + j) += v;
}

template <class type>
inline void array2D<type>::to(const unsigned int i, const unsigned int j, const type v)
{
        values.at(i*no_cols+j) = v;
}

template <class type>
inline type array2D<type>::at(const unsigned int i, const unsigned int j)
{
        return values.at(i*no_cols+j);
}


#endif // _ARRAY2D_HPP
