// ~~~~~~~~~~~~~~~~~~~~~~~~~
// t r e e _ n o d e . h p p
// ~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef _TREE_NODE_HPP
#define _TREE_NODE_HPP

using namespace std;

template <class type>
class idxNode
{
        friend class binary_heap;
        public:
                idxNode(void);
                void init(const int i, const type v);
                void set_index(const int i);
                void set_value(const type v);
                int get_index(void);
                type get_value(void);
        public:
                int index;
                type value;
        public:
                void operator=(idxNode & rt);
                bool operator<(idxNode rt);
                bool operator<=(idxNode rt);
                bool operator==(idxNode rt);
                bool operator>=(idxNode rt);
                bool operator>(idxNode rt);
};

template <class type>
inline idxNode<type>::idxNode(void)
{
        index = 0;
        value = 0.0;
}

template <class type>
inline void idxNode<type>::init(const int i, const type v)
{
        index = i;
        value = v;
}

template <class type>
inline void idxNode<type>::set_index(const int i)
{
        index = i;
}

template <class type>
inline void idxNode<type>::set_value(const type v)
{
        value = v;
}

template <class type>
inline int idxNode<type>::get_index(void)
{
        return index;
}

template <class type>
inline type idxNode<type>::get_value(void)
{
        return value;
}

template <class type>
inline void idxNode<type>::operator=(idxNode & rt)
{
        if (this == &rt) return;
        index = rt.index;
        value = rt.value;
        return;
}

template <class type>
inline bool idxNode<type>::operator<(idxNode rt)
{
        if (this == &rt) return false;
        return (value < rt.value);
}

template <class type>
inline bool idxNode<type>::operator<=(idxNode rt)
{
        if (this == &rt) return true;
        return (value <= rt.value);
}

template <class type>
inline bool idxNode<type>::operator==(idxNode rt)
{
        if (this == &rt) return false;
        return (value == rt.value);
}

template <class type>
inline bool idxNode<type>::operator>=(idxNode rt)
{
        if (this == &rt) return true;
        return (value >= rt.value);
}

template <class type>
inline bool idxNode<type>::operator>(idxNode rt)
{
        if (this == &rt) return false;
        return (value > rt.value);
}

#endif // _TREE_NODE_HPP
