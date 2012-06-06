// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// b i n a r y _ t r e e . h p p
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef _BINARY_TREE_HPP
#define _BINARY_TREE_HPP

#include<iostream>
#include<cstdlib>

#include "tree_node.hpp"

using namespace std;

class binary_tree
{
        public:
                binary_tree(void) {};
                int init(const unsigned int size,const double * values);
                void tree_sort(void);
                void tree_updt(const unsigned int position,const double old_value,const double new_value);
                void tree_rsrt(const unsigned int position);
                bool is_this_tree(void);
        private:
                void sift_down(const unsigned int l,const unsigned int r);
        public:
                inline void set_index(const unsigned int i,const int idx);
                inline void set_value(const unsigned int i,const double val);
                inline void set_size(const unsigned int i);
                inline int get_index(const unsigned int i);
                inline double get_value(const unsigned int i);
                inline unsigned int get_size(void);
        private:
                unsigned int tree_size;
                idxNode<double> * tree;
};

inline void binary_tree::set_index(const unsigned int i,const int idx)
{
        tree[i].index = idx;
}

inline void binary_tree::set_value(const int unsigned i,const double v)
{
        tree[i].value = v;
}

inline void binary_tree::set_size(const unsigned int i)
{
        tree_size = i;
}

inline int binary_tree::get_index(const unsigned int i)
{
        return tree[i].index;
}

inline double binary_tree::get_value(const unsigned int i)
{
        return tree[i].value;
}

inline unsigned int binary_tree::get_size(void)
{
        return tree_size;
}

#endif // _BINARY_TREE_HPP
