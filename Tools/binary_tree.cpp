// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// b i n a r y _ t r e e . c p p
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include"binary_tree.hpp"

using namespace std;

int binary_tree::init(const unsigned int size, const double * values)
{
        tree_size = size;
        tree = new idxNode<double>[tree_size];
        for(unsigned int i = 0; i < tree_size; i++) {
                tree[i].init(i,values[i]);
        }
        tree_sort();

        if(!is_this_tree())
        {
                cerr << "Binary tree is not properly sorted in function binary_tree::init()!" << endl;
                exit(EXIT_FAILURE);
        }
        return(EXIT_SUCCESS);
}

// insert free(void)
void binary_tree::tree_sort(void)
{
        for(int i = tree_size/2-1; i >= 0; i--) {
                sift_down(i,tree_size-1);
        }
}

void binary_tree::sift_down(const unsigned int lft, const unsigned int rgt)
{
        idxNode<double> tmp_node = tree[lft];
        unsigned int jold = lft;
        unsigned int jnew = 2*lft+1;
        while(jnew <= rgt)
        {
                if((jnew < rgt) && (tree[jnew].value > tree[jnew+1].value)) jnew++;
                if(tmp_node.value <= tree[jnew].value) break;
                tree[jold] = tree[jnew];
                jold = jnew;
                jnew = 2*jnew+1;
        }
        tree[jold] = tmp_node;
}

void binary_tree::tree_updt(const unsigned int position,const double old_value,const double new_value)
{
        if(new_value < old_value)
        {
                cerr << "Attempted to update to smaller time in function binary_tree::tree_updt!" << endl;
                exit(EXIT_FAILURE);
        }
        tree[position].set_value(new_value);
}

void binary_tree::tree_rsrt(const unsigned int position)
{
        int i = (int)position;
        while(i > 0)
        {
                sift_down(i,tree_size-1),
                i = (i-1)/2;
        }
        sift_down(0,tree_size-1);
        if(!is_this_tree())
        {
                cerr << "Binary tree is not properly sorted in function binary_tree::tree_updt()!" << endl;
                exit(EXIT_FAILURE);
        }
}

bool binary_tree::is_this_tree(void)
{
        bool this_is_tree = true;
        for(unsigned int i = 1; i < tree_size; i++)
        {
                if(tree[(i-1)/2] > tree[i])
                {
                        this_is_tree = false;
                        break;
                }
        }
        return this_is_tree;
}
