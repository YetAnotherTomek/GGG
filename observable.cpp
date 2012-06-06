// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
// o b s e r v a b l e . c p p
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~

//  ~~~~~~~~~~~~~~~~~~~~~~~~~~
//        T.R. Sokolowski    
//        2011
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~


#include "observable.hpp"

observable::observable(void)
{

}

observable& observable::operator=(const observable& o)
{
    number = o.number;
    strcpy(name, o.name);
    no_components = o.no_components;
    for(int i=0; i<no_components; i++)  component[i] = o.component[i];
    value = o.value;
    threshold = o.threshold;
}

void observable::info(void)
{
    printf("Observable %u = %s \n", number, name);
    printf("no_components = %u \n", no_components);
    for(int i=0; i<no_components; i++)
        printf("\t component %u = %u \n", i, component[i]);
    printf("threshold = %g \n", threshold);
    printf("value = %g \n\n", value);    
}
