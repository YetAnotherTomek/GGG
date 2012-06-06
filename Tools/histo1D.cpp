// ~~~~~~~~~~~~~~~~~~~~~~
// h i s t o 1 D . c p p 
// ~~~~~~~~~~~~~~~~~~~~~~

#include"histo1D.hpp"

using namespace std;

void histo1D::init(const int l,const int r,const double d)
{
        left_bin = l;
        rght_bin = r;
        bins = rght_bin - left_bin + 1;
        wdth = d;
    
    events.clear();
    weight.clear();
    probab.clear();
        for(int i = 0; i <= bins; i++)
        {
                events.push_back(0);
                weight.push_back(0.0);
                probab.push_back(0.0);
        }
}

void histo1D::add(const int bin,const double w)
{
        if (bin >= left_bin && bin <= rght_bin)
        {
                events.at(bin-left_bin) += 1;
                weight.at(bin-left_bin) += w;
        }
}

void histo1D::normalise(void)
{
        double norm = sum_weight();
        if (norm > 0.0)
        {
                norm = 1.0/(wdth*norm);
        }
        else
        {
                norm = 0.0;
        }
        for (int bin = 0; bin < bins; bin++)
        {
                probab.at(bin) = norm*weight.at(bin);
        }
}

unsigned int histo1D::sum_events(void)
{
        unsigned int sum = 0;
        for (int bin = 0; bin < bins; bin++)
        {
                sum += events.at(bin);
        }
        return sum;
}

double histo1D::sum_weight(void)
{
        double sum = 0.0;
        for (int bin = 0; bin < bins; bin++)
        {
                sum += weight.at(bin);
        }
        return sum;
}

double histo1D::sum_probab(void)
{
        double sum = 0.0;
        for (int bin = 0; bin < bins; bin++)
        {
                sum += probab.at(bin);
        }
        return sum;
}

double histo1D::average(void)
{
        double av = 0.0;
        if (wdth == 1.0)
        {
                for(int bin = 0;bin < bins;bin++)
                {
                        av += (double(left_bin+bin))*probab.at(bin);
                }
        }
        else
        {
          for(int bin = 0;bin < bins;bin++)
          {
                        av += ((double(left_bin+bin)+0.5)*wdth)*probab.at(bin)*wdth;
          }
        }
        return av;
}

double histo1D::variance(void)
{
        double avg = average();
        double var = 0.0;
        if (wdth == 1.0)
        {
                for (int bin = 0; bin < bins; bin++)
                {
                        var += pow(double(left_bin+bin),2)*probab.at(bin);
                }
        }
        else
        {
          for (int bin = 0; bin < bins; bin++)
          {
                        var += pow(((double(left_bin+bin)+0.5)*wdth),2)*probab.at(bin)*wdth;
          }
        }
        return sqrt(var - avg*avg);
}

double histo1D::std_dev(void)
{
        return sqrt( this->variance() );
}
