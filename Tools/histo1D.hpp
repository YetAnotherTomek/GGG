// ~~~~~~~~~~~~~~~~~~~~~~
// h i s t o 1 D . h p p 
// ~~~~~~~~~~~~~~~~~~~~~~

#ifndef _HISTO1D_HPP
#define _HISTO1D_HPP

#include<cmath>
#include<vector>

using namespace std;

class histo1D {
        public:
                histo1D(void){};
                void init(const int l,const int r,const double binsize = 1.0);
                void add(const int bin,const double w);
                void normalise(void);
                double average(void);
                double variance(void);
                double std_dev(void);
        private:
                int bins;
                int left_bin;
                int rght_bin;
                double wdth;
                vector<unsigned int> events;
                vector<double> weight;
                vector<double> probab;
                unsigned int sum_events(void);
                double sum_weight(void);
    public:
                double sum_probab(void);
        public:
                inline int no_bins(void);
                inline unsigned int events_at(const int bin);
                inline double weight_at(const int bin);
                inline double probab_at(const int bin);
                inline double densty_at(const int bin);
                inline void events_to(const int bin, const int e);
                inline void weight_to(const int bin, const double w);
};

inline int histo1D::no_bins(void)
{
        return bins;
}

inline unsigned int histo1D::events_at(const int bin)
{
        return events.at(bin-left_bin);
}

inline double histo1D::weight_at(const int bin)
{
        return weight.at(bin-left_bin);
}

inline double histo1D::probab_at(const int bin)
{
        return probab.at(bin-left_bin);
}

inline double histo1D::densty_at(const int bin)
{
        return wdth*probab.at(bin-left_bin);
}

inline void histo1D::events_to(const int bin, const int e)
{
        events.at(bin-left_bin) = e;
}

inline void histo1D::weight_to(const int bin, const double w)
{
        weight.at(bin-left_bin) = w;
}

#endif // _HISTO1D_HPP
