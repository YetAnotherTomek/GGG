// ~~~~~~~~~~~~~~~~~~~~~~
// h i s t o 2 D . h p p 
// ~~~~~~~~~~~~~~~~~~~~~~

#ifndef _HISTO2D_HPP
#define _HISTO2D_HPP

#include"array2D.hpp"

using namespace std;

class histo2D {
        public:
                histo2D(void) {};
                inline int no_cols(void);
                inline int no_rows(void);
                inline unsigned int events_at(const int r,const int c);
                inline double weight_at(const int r,const int c);
                inline double probab_at(const int r,const int c);
                inline double densty_at(const int r,const int c);
                void init(const int, const int, const int, const int, const double binsize=1.0);
                void add(const int r,const int c,const double w = 1);
                void normalise(void);
    public:
        double sum_probab(void);
        private:
                unsigned int sum_events(void);
                double sum_weight(void);
                int rows;
                int cols;
                int min_rows,max_rows;
                int min_cols,max_cols;
                double area;
                array2D<unsigned int> events;
                array2D<double> weight;
                array2D<double> probab;
};

inline int histo2D::no_cols(void)
{
        return cols;
}

inline int histo2D::no_rows(void)
{
        return rows;
}

inline unsigned int histo2D::events_at(const int r,const int c)
{
        return events.at(r-min_rows,c-min_cols);
}

inline double histo2D::weight_at(const int r,const int c)
{
        return weight.at(r-min_rows,c-min_cols);
}

inline double histo2D::probab_at(const int r,const int c)
{
        return probab.at(r-min_rows,c-min_cols);
}

inline double histo2D::densty_at(const int r,const int c)
{
        return area*probab.at(r-min_rows,c-min_cols);
}

#endif // _HISTO2D_HPP
