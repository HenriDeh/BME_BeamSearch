#include "fastme.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fcntl.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <libgen.h>
#include <assert.h>
#include <iostream>
#include <fstream>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>


using string = std::string;

int** fill_int_matrix(int m) {
    int** matrix = new int*[m];
    std::fstream file("init_mat");
    std::string line;
    // std::cout<<" ***\n";
    int row=0;
    while (getline(file, line) and row < m){
        std::stringstream ss( line );                     
        std::string data;
        int col = 0;
        matrix[row] = new int[m];
        while (getline( ss, data, ' ' ) )           
        {
            matrix[row][col]=stoi(data);
            // std::cout<<matrix[row][col]<<" ";
            col++;
            
        }
        // std::cout<<"\n";
        row++;
    }
    return matrix;
}


double** fill_matrix(int n_taxa) {
    double** matrix = new double*[n_taxa];
    std::fstream file("mat");
    std::string line;
    // std::cout<<" ***\n";
    int row=0;
    while (getline(file, line) and row < n_taxa){
        std::stringstream ss( line );                     
        std::string data;
        int col = 0;
        matrix[row] = new double[n_taxa];
        while (getline( ss, data, ' ' ) )           
        {
            matrix[row][col]=stod(data);
            col++;
            // std::cout<<stod(data)<<" ";
        }
        // std::cout<<"\n";
        row++;
    }
    return matrix;
}

int get_n_taxa() {

    std::fstream file("n_taxa");
    std::string line;
    int n_taxa;
    while (getline(file, line)){
        std::stringstream ss( line );                     
        std::string data;

        while (getline( ss, data, ' ' ) )   n_taxa = stoi(data);
        }

    std::cout<<"nuuuum taxa "<<n_taxa<<std::endl;
    return n_taxa;
}

void save_mat(int* adj_mat, int mat_size){
    std::fstream myfile;

    myfile.open("../result_adj_mat.txt",std::fstream::out);

    for (int i=0; i< mat_size;i++) //This variable is for each row below the x 
    {        

        for (int j=0; j<mat_size;j++)
        {                      
            myfile << adj_mat[i*mat_size + j] << "\t";
        }
        myfile<<std::endl;
    }
    myfile.close();
}

int main () {
    // std::cout.precision(17);
    // std::cout<<"from main"<<std::endl;
    int n_taxa = get_n_taxa();
    int m = (n_taxa * 2) -2;
    double** D = fill_matrix(n_taxa);
    int** init_adj = fill_int_matrix(m);

    clock_t start, end;
    start = clock();
    int* solution_mat = new int[m*m];
    double obj_val;
    int nni_count;
    int spr_count;
    run(D, init_adj, solution_mat, n_taxa, m, obj_val, nni_count, spr_count);
    end = clock();

    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    std::cout << "Time taken by program is : " << time_taken << std::endl;
    save_mat(solution_mat, m);

    for(int i=0; i<n_taxa; i++){
        delete[] D[i];
    }

    for(int i=0; i<m; i++){
        delete[] init_adj[i];
    }

    delete[] D; delete[] init_adj; delete[] solution_mat;
    std::cout<<"done"<<std::endl;
    
}
