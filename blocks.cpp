/*
the problem is such: determine the number of ways you can create a panel of 48" long, 10" tall using only blocks length 3.0" and 4.5", blocks are 1" tall. endpoints between blocks cannot line up (except for in the end). AKA there cannot be "cracks" that are more than 2" tall. Therefore this is not just a simple permutation problem. Additionally (for simplicity sake) assume that blocks can only be placed laterally.

Do what you want with this code I don't care. HOWEVER, do not claim that I stole (ANY) of this code because the entire project was created, designed, coded by me, I have log files...

The problem basically breaks down to two different subset sum problems (NP-hard), first we need to solve the number of configurations of length 48.

I am aware that there is a pseudo-polynomial time dynamic programming algorithm that solves this in time 0(nb), but that seems like overkill for just a sum of 48 inches so it might not neccessarily be all that much faster. In addition the dynamic programming method does not solve the problem of finding permutations while the brutish (DFS based) algorithm
exhaustively lists them

Additionally there is an optimization that can be done at generate_adj_matrix. It involves stepping down a decision tree but I'm not going to describe because I'm too lazy. Actually I'll elaborate a little bit more -- basically it boils down to the fact if you have a descendant node from a reference node (+1 or more degrees away), every sub-descendant will also be INCOMPATIBLE if the descendent is INCOMPATIBLE with the reference node. The implementation would have been a pain in the ass to complete, but a place to start would be to look into the reasoning behind why I created sides A and B of the bipartite graph. (HINT: everything in A is incompatible with everything in A... same thing with B)

For the most part the code should be fairly self documenting.

-a.liang
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

int length = 48; //length of the panel is 48 inches

using namespace std;

/* a row is described as such: we will have a bit for each 1.5 inch that is taken up by a block, a 0 denotes where each block ends: 
 * 3" block: 01
 * 4.5" block: 011
 * for example: 0110101101101 would be a panel consisting of 4.5", 3", 4.5", 4.5", 3" blocks
 */

//clearly, the number of bits that will be taken up is 48/1.5 -> 32 bits -> each row can be represented in a 4 bytes

class Row {

        //making everything public, lol
        public:
                unsigned int blocks;
                double inches; //length of the block in inches, probably not needed
        
                Row (){
                        
                        blocks = 0;
                        length = 0;
                }

                bool add_block(double length){ 
                
                        //length needs to be a double of either 3.0 or 4.5
                        //we need to make space for this block, we need to shift left by the number of bits

                        int num_bits = (int)(length/1.5);
                        int bits;
        
                        //unable to add a block!
                        if((length+inches) > 48)
                                return false; 

                        if(length == 3.0)
                                bits = 1;
                        else if(length == 4.5)
                                bits = 3;
                        else{
                                printf("invalid length = %f\n", length);
                                exit(1);
                                //then we have an invalid block length, shouldn't come up at all
                        }

                        int shifted = bits << (int) (inches/1.5);
                        blocks = blocks | shifted;

                        inches += length;

                        return true;
                }

                //for those of you who can't do decimal to binary quickly :)
                void print_row(){
                        int cnt, mask = 1 << 31;
                        int temp = blocks;
        
                        for(cnt=1;cnt<=32;++cnt)
                        {
                                putchar(((temp & mask) == 0) ? '0' : '1');
                                temp <<= 1;
                                if(cnt % 8 == 0 && cnt !=32)
                                        putchar(' ');
                                if(cnt == 32)
                                        putchar('\n');
                        }
                }

                //i'm not going to explain this because I'm too lazy. read the description in the beginning
                unsigned int check_holes(Row *b){
                        //checks if there is a hole
                        int bits = 0;
                        unsigned int ret = (~blocks)&(~b->blocks);

                        //the smaller row (if any) determines the shift value
                        if(b->inches < inches){
                                bits = (int)((48-b->inches)/(1.5));
                        }
                        else{
                                bits = (int)((48-inches)/1.5);
                        }
                        ret = ret << (bits +1); //got to shift out the meaningless 1s!
                        ret = ret >> (bits +1);
                        return ret;
                
                }
};

//prints the int in binary form
void int_to_binary(int x)
{
        int cnt, mask = 1 << 31;
        
        for(cnt=1;cnt<=32;++cnt)
        {
                putchar(((x & mask) == 0) ? '0' : '1');
                x <<= 1;
                if(cnt % 8 == 0 && cnt !=32)
                        putchar(' ');
                if(cnt == 32)
                        putchar('\n');
        }
}

//basic linked list structure

template <class T> class LinkedList{

        public:

                struct node{
                        T element;
                        node *next;
                };

                node *head;
                node *tail;
                int length;
                //each node in the linked list is a class of row types.
          
                LinkedList(){
                        head = NULL;
                        tail = NULL;
                        length = 0;
                }

                ~LinkedList(){
                        node *temp = head;
                        while(temp != NULL){
                                node *next = temp->next;
                                delete temp;
                                temp = next;
                        }
                        head = NULL;
                        tail = NULL;
                }

                void insert(T el){
                        if(head == NULL){
                                head = new node;
                                head->element = el;
                                head->next = NULL;
                                tail = head;
                        }
                        else{
                                node *newest;
                                newest = new node;
                                newest->element = el;
                                newest->next = NULL;
                                tail->next = newest;
                                tail = newest;
                        }
                        length++;
                }
             
};

//its unfortunate that C++ doesnt support local types in a containers, declaring this here is really the only solution.(I use it in MagicGraph)
struct tuple{
        int index;
        int value;
};

class MagicGraph{

        struct gnode{
                Row row;
                LinkedList <gnode*> adj;
        };

        public:
                //A, B are two sides of a bipartite graph. each node in A is a Row and a list of pointers to nodes (adjacency lists)

                LinkedList <gnode> A;
                LinkedList <gnode> B;

                //adjacency matrix representation is used to calculate the number of paths of n-length
                //since the matrix is huge we should go for a sparse matrix representation

                //basically we need an array for the rows and an array for the columns, each respective array cell stores a linked list of tuples containing the index and the value of the non-zero elements of the adjacency matrix
                        
                //tuple stores the index in the row/column of a non-zero element

                struct s_matrix{
                        LinkedList<tuple> * s_row; //ok they're pointers but they're also (correctly  allocated) arrays, I swear!
                        LinkedList<tuple> * s_col;
                } mat;

                int N;

                void set_N(int length){
                        N = length;
                }

                int counter;

                MagicGraph(){
                        counter = 0;
                }

                /* force a "bipartiteness" quality on our generated data so on side A we have the nodes that are trivially uncompatible with each other and same thing on side B - most obvious choice is to put all rows that                      * begin with a 3.0 in side A, and rows that begin with a 4.5 in side B, It is possible to enforce this quality throughout each subtree but I don't think I'm going to do it because it might be overkill...
                 */
                 
                void generate_decision_tree(){
                        Row base_a;
                        base_a.add_block(3.0);
                        generate_side_a(base_a, 0);
                        Row base_b;
                        base_b.add_block(4.5);
                        generate_side_a(base_b, 1);
        
                }
                

                //"blindly" builds side A - blindly A.head->element.adj.head->element->row.print_row();in the sense that it doesnt check for compatibility between blocks simply because it doesn't need to

                void generate_side_a(Row row, int control){
Draft comment:Edit
admittedly this is poorly named. it should just be named generate_side. the control int could simply just be a control char, '0' for side a '1' for side a
                        //printf("blocks: %d\n", row.blocks);
                        //row.print_row();
                        if(row.inches == 48){
                                //row.print_row();

                                gnode wrapper;
                                wrapper.row = row;

                                if(control == 0)
                                        A.insert(wrapper);
                                else
                                        B.insert(wrapper);
                        }
                        Row left = row;
                        if(left.add_block(3.0) == false)
                                return;
                        generate_side_a(left, control);
                        Row right = row;
                        if(right.add_block(4.5) == false)
                                return;
                        generate_side_a(right, control);
                }

                //makes the edges that describe a bipartite mapping of A to B

                //note: this is approximately a O((n/2)^2) algorithm, so O(n^2) really...

                void generate_adj_matrix(){

                        //it was getting pretty annoying to type it out so i'll just use macros

                        #define A_ADJ temp_a->graph_node.adj 
                        #define B_ADJ temp_b->graph_node.adj

                        #define A_ROW temp_a->graph_node.row
                        #define B_ROW temp_b->graph_node.row
                        
                        /*needed to take stuff off of a linked list, admittedly I could have done this within the LinkedList class but it would have complicated things especially in the DFS style recursion*/

                        struct node{
                                gnode graph_node;
                                node* next;
                        };

                        int N = A.length + B.length;
                        set_N(N);

                        mat.s_row = new LinkedList <tuple>[N];
                        mat.s_col = new LinkedList <tuple>[N];

                        //the coordinates of the matrix (they aren't true coordinates since I exploited the fact that the matrix is symmetrical)
                        int i = 0;
                        int j = 0;
                        int offset_mat_i;

                        node* temp_a = (node*)A.head;
                        while(temp_a != NULL){

                                node* temp_b = (node*) B.head;
                                while(temp_b != NULL){

                                        if(A_ROW.check_holes(&(B_ROW)) == 0){
                                                //if zero we are good!

                                                //populating a (sparse!) adjacency matrix
                                                //adj[i][A.length + j] = 1;                                             
                                                //adj[A.length + j][i] = 1;
                                                tuple t_row, t_col;

                                                t_row.index = A.length + j;
                                                t_row.value = 1;
                                                mat.s_row[i].insert(t_row);
                                                t_row.index = i;
                                                mat.s_row[A.length + j].insert(t_row);
                                                t_col.index = i;
                                                t_col.value = 1;
                                                mat.s_col[A.length +j].insert(t_col);
                                                t_col.index = A.length + j;
                                                mat.s_col[i].insert(t_col);
                                        
                                        }
                                        else{
                                                //do nothing
                                                //adj[i][A.length + j] = 0;
                                                //adj[A.length + j][i] = 0;
                                                
                                        }
                                        temp_b = temp_b->next;
                                        j++;
                                   
                                }

                                temp_a = temp_a->next;
                                i++;
                                j = 0;
                        }
        
                        printf("N: %d\n", N);

                }
                        
                //generates panels that are a specified height
                int generate_panels(int height){
                        struct d_node{
                                tuple s_tuple;
                                d_node* next;
                        };

                        d_node *temp;

                        //using dynamic programming
                        int N = A.length + B.length;

                        int64_t *tableau = (int64_t *)malloc(N*height*sizeof(int64_t));
                        //values of the table at position i, j of the table means the number of ways to stack a panel of height (i+1), with row type j at the very top

                        printf("sizeof int64_t: %d\n", sizeof(int64_t));

                        for(int j = 0; j < N; j++){
                                tableau[j] = 1; 
                                //trivial case, simply 1 way to stack a row of type j and height 1

                                
                        }

                        //now it gets interesting...
                        for(int i = 1; i < 10; i++){
                                for(int j = 0; j < N; j++){
                                        //consult the (sparse) adjacency matrix, look up the row for type j
                                        temp = (d_node *) mat.s_row[j].head;
                                        int64_t sum = 0;
                                        while(temp){
                                                int t_index = temp->s_tuple.index;
                                                
                                                sum = tableau[(i-1)*N+t_index] + sum;
                                                temp = temp->next;
                                        
                                                
                                        }                                       
                                        tableau[i*N+j] = sum;
                                        printf("sum: %lld\n", sum); 
                                }
                        }
                        printf("done!\n");

                        //sum the top row, and we're done!
                        int64_t total_sum = 0;
                        for(int j = 0; j < N; j++){
                                //printf("printing: %lld\n", tableau[(9)*N+j]);
                                total_sum = tableau[(height-1)*N+j] + total_sum;
                        }
                        printf("total paneling: %lld\n", total_sum);

                }
                
                //raises a sparse matrix to the nth power, not used
                void sparse_power(s_matrix a, s_matrix *c, int n){
                        //node for dequeuing off the linked list
                        struct d_node{
                                tuple s_tuple;
                                d_node* next;
                        };

                        a.s_row;
                        a.s_col;                        

                }
                
                //used in debug (which isn't used)
                void sparse_multiply(s_matrix a, s_matrix b, s_matrix *c, int length){
                        struct d_node{
                                tuple s_tuple;
                                d_node* next;
                        };

                        d_node* t_row_a;
                        d_node* t_col_b;
                        for(int j = 0; j < length; j++){
                                for(int i = 0; i < length; i++){
                                        t_col_b = (d_node *) b.s_col[j].head;
                                        t_row_a = (d_node *) a.s_row[i].head;
                                        int sum = 0;
                                        //dot product the lists
                                        //while neither is not NULL
                                        while(t_col_b && t_row_a){
                                                if(t_col_b->s_tuple.index == t_row_a->s_tuple.index){
                                                        sum += (t_col_b->s_tuple.value * t_row_a->s_tuple.value);
                                                        t_col_b = t_col_b->next;
                                                        t_row_a = t_row_a->next;
                                                }
                                                else if(t_col_b->s_tuple.index > t_row_a->s_tuple.index)
                                                        t_row_a = t_row_a->next;
                                                else
                                                        t_col_b = t_col_b->next;
                                        }
                                        //REMEMBER TO ADD THE SUM INTO s_matrix *c
                                        tuple c_element;
                                        c_element.value = sum;
                                        if(sum != 0){
                                                printf("sum @ [%d, %d]: %d\n", i, j, sum);
                                        }
                                        c_element.index = i;
                                        c->s_col[j].insert(c_element);
                                        c_element.index = j;
                                        c->s_row[i].insert(c_element);
                                }
                        }
                }
                
                //not used
                void debug(){
                        s_matrix test;

                        test.s_row = new LinkedList <tuple>[N];
                        test.s_col = new LinkedList <tuple>[N];

                        s_matrix test2;

                        test2.s_row = new LinkedList <tuple>[N];
                        test2.s_col = new LinkedList <tuple>[N];

                        sparse_multiply(mat, mat, &test, N);
                        sparse_multiply(mat, test, &test2,N);

                }
};


int main(){

        MagicGraph tree;
        tree.generate_decision_tree();
        tree.generate_adj_matrix();
        tree.generate_panels(10);
        
}