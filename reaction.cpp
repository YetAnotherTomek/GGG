// ~~~~~~~~~~~~~~~~~~~~~~~
// r e a c t i o n . c p p
// ~~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~
//     T.R. Sokolowski    
//     2011
// ~~~~~~~~~~~~~~~~~~~~~~~

#include "reaction.hpp"


reaction_class::reaction_class(void)
{
  
}

void reaction_class::init(int n, reac_type t, double r, int par1, int par2, int prod1, int prod2)
{
 
    number = n;
    type = t;
    rate = r;
    
    partner[0] = par1;
    partner[1] = par2;
    product[0] = prod1;
    product[1] = prod2;  
}

reaction_class& reaction_class::operator=(const reaction_class& r)
{
  
    number = r.number;
    type = r.type;
    rate = r.rate;
    
    partner[0] = r.partner[0];
    partner[1] = r.partner[1];
    product[0] = r.product[0];
    product[1] = r.product[1];
    
    return *this;
}
  
double reaction_class::prop(double* conf)
{
 
    double p;
    
    switch(type)
    {
      case PRODUCTION:
        
          p = rate * conf[partner[0]];
          break;
          
      case DECAY:
        
          p = rate * conf[partner[0]];
          break;
          
      case BINDING:
        
          p = rate * conf[partner[0]] * conf[partner[1]];
          break;
          
      case UNBINDING:
        
          p = rate * conf[partner[0]];
          break;                    
          
      case DIMERIZATION:
        
          // here by def. partner[0]==partner[1]
          p = rate * conf[partner[0]] * (conf[partner[0]] - 1.0);
          break;
    };
    
    return p;
}

void reaction_class::fire(double* conf)
{     
    switch(type)
    {
      case PRODUCTION:
        
          conf[product[0]] = conf[product[0]] + product[1];
          // product[1] holds the prod. burst size
          break;
          
      case DECAY:
        
          conf[partner[0]] = conf[partner[0]] - 1;
          break;

      case BINDING:
        
          conf[partner[0]] = conf[partner[0]] - 1;
          conf[partner[1]] = conf[partner[1]] - 1;
          conf[product[0]] = conf[product[0]] + 1;
          break;
          
      case UNBINDING:
        
          conf[partner[0]] = conf[partner[0]] - 1;
          conf[product[0]] = conf[product[0]] + 1;
          conf[product[1]] = conf[product[1]] + 1;
          break;
                    
      case DIMERIZATION:
        
          // here by def. partner[0]==partner[1]
          // dedimerization is a special case of unbinding
          conf[partner[0]] = conf[partner[0]] - 2;
          conf[product[0]] = conf[product[0]] + 1;
          break;       
    };    
}

void reaction_class::info(void)
{
    cout    << "Reaction " << number << " :" << endl
            << " type = " << type << endl
            << " rate = " << rate << endl
            << " partner[0] = " << partner[0] << endl
            << " partner[1] = " << partner[1] << endl
            << " product[0] = " << product[0] << endl
            << " product[1] = " << product[1] << endl
            << endl;
}
