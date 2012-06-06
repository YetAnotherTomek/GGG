#ifndef _OBSERVABLE_HPP
#define _OBSERVABLE_HPP

#include<stdio.h>
#include<string.h>

class observable
{
  
  public:
    
        int number;
        char name[30];
        int no_components;
        int component[100];
        double value;
        double threshold;                
        
        observable(void);
        observable& operator=(const observable&);
        void info(void);
};

#endif // _OBSERVABLE_HPP

