// ~~~~~~~~~~~~~~~~~~~~~~
// h i s t o 2 D . c p p 
// ~~~~~~~~~~~~~~~~~~~~~~

#include"histo2D.hpp"

using namespace std;

void histo2D::init(const int min_r,const int max_r,const int min_c,const int max_c,const double d)
{
        min_cols = min_c;
        max_rows = max_r;
        min_rows = min_r;
        max_cols = max_c;
        rows = max_r - min_r + 1;
        cols = max_c - min_c + 1;
        area = d*d;
        events.init(rows,cols,0);
        weight.init(rows,cols,0.0);
        probab.init(rows,cols,0.0);
}

void histo2D::add(const int r,const int c,const double w)
{
        if (r >= min_rows && r <= max_rows && c >= min_cols && c <= max_cols)
        {
                events.add(r-min_rows,c-min_cols,1);
                weight.add(r-min_rows,c-min_cols,w);
        }
}

void histo2D::normalise(void)
{
        double norm = sum_weight();
        if (norm > 0.0)
                norm = 1.0/(area*norm);
        else
                norm = 0.0;
        for (int r = 0; r < rows; r++)
                for (int c = 0; c < cols; c++)
                        probab.to(r,c,norm*weight.at(r,c));
}

unsigned int histo2D::sum_events(void)
{
        unsigned int sum = 0;
        for (int r = 0; r < rows; r++)
                for (int c = 0; c < cols; c++)
                        sum += events.at(r,c);
        return sum;
}

double histo2D::sum_weight(void)
{
        double sum = 0.0;
        for (int r = 0; r < rows; r++)
                for (int c = 0; c < cols; c++)
                        sum += weight.at(r,c);
        return sum;
}

double histo2D::sum_probab(void)
{
        double sum = 0.0;
        for (int r = 0; r < rows; r++)
                for (int c = 0; c < cols; c++)
                        sum += probab.at(r,c);
        return sum;
}
