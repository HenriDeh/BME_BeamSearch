#include "fastme.h"
#include <iostream>
#include <omp.h>
#include <sched.h>
#include "results.h"


extern "C" {

    int* test (double* d, int* init_adj, int n_taxa, int m) {

        double ** D = new double*[n_taxa];
        int i,j;

        for(i = 0; i< n_taxa; i++) D[i] = &d[i*n_taxa];
 
        int ** A = new int*[m];
        for(i = 0; i< m; i++) A[i] = &init_adj[i*m];
                
        int* solution_mat = new int[m*m];
        double obj_val;
        int nni_count;
        int spr_count;
        run(D, A, solution_mat, n_taxa, m, obj_val, nni_count, spr_count);

        if (D!= nullptr) delete[] D;
        if (A!= nullptr) delete[] A;

        return solution_mat;

        }




    results* test_parallel (double* d, int* init_adj, int n_taxa, int m, int population_size, int num_procs) {

        
        double ** D = new double*[n_taxa];
        int i,j,t, mat_size;
        mat_size = m*m;

        for(i = 0; i< n_taxa; i++) D[i] = &d[i*n_taxa];

        int *** A = new int**[population_size];
        int* solution_mat = new int[population_size*mat_size];

        

        double* obj_vals = new double[population_size];
        int * nni_counts = new int[population_size];
        int * spr_counts = new int[population_size];



        omp_set_num_threads(num_procs);
        #pragma omp parallel for schedule(dynamic) shared(A, d)
        for(t = 0; t < population_size; t++) {
            A[t]= new int*[m];
            for(i = 0; i< m; i++) A[t][i] = &init_adj[t*mat_size + i*m];
                
            run(D, A[t], &solution_mat[t*mat_size], n_taxa, m, obj_vals[t], nni_counts[t], spr_counts[t]);
        }

        results* res = new results[1];
        res -> nni_counts = 0;
        res -> spr_counts = 0;
        res -> solution_adjs = solution_mat;
        res -> objs = obj_vals;

        for(t = 0; t < population_size; t++){
            delete[] A[t];
            res -> nni_counts += nni_counts[t];
            res -> spr_counts += spr_counts[t];
            // std::cout<<obj_vals[t]<<" ";
        }
        // std::cout<<std::endl;
        delete[] D;
        delete[] A; 
        delete[] nni_counts;
        delete[] spr_counts;

        return res;
        }

    void free_result(results* res){
        delete[] res -> solution_adjs;
        delete[] res -> objs;
        res -> solution_adjs = nullptr;
        res -> objs = nullptr;
        res  = nullptr;
    }


    results* test_obj(){

        int * adj = new int[10];
        double* obj_vals = new double[10];

        for(int i=0; i< 10; i ++ ) {
            adj[i] = 2*i;
            obj_vals[i] = 0.2*i;
        }
        results* obj = new results[1];
        obj ->objs = obj_vals; obj->nni_counts=4; obj->spr_counts=3; obj->solution_adjs=adj;
        std::cout<<"fjfjfjfjf"<< obj->spr_counts<<std::endl;

        return obj;
    }

}



