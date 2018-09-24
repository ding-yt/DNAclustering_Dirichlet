//
//  main.cpp
//  simulation_dirichlet
//
//  Created by Yuantong Ding on 4/21/13.
//  Copyright (c) 2013 Yuantong Ding. All rights reserved.
//

#include <vector>
#include <algorithm>
#include <iostream>
extern "C" {
    #include <gsl/gsl_sf_bessel.h>
    #include <gsl/gsl_rng.h>
    #include <gsl/gsl_randist.h>
    #include <gsl/gsl_sf_pow_int.h>
    #include <gsl/gsl_sf_log.h>
    #include <gsl/gsl_statistics_double.h>
}

// ###############################################
// The fuction below simulates a base from a multinomial distribution prob[]
// ###############################################
int simulate_base(const double prob[]){
    const gsl_rng_type *T;
    gsl_rng *r;
    size_t t = 4;
    int base;
    const unsigned int n = 1;
    unsigned int theta[4];
    
    T = gsl_rng_ranlxs0;
    gsl_rng_default_seed = (rand());
    r = gsl_rng_alloc(T);
    gsl_ran_multinomial(r, t, n, prob, theta);
    for (int i = 0; i < t; i++) {
        if (theta[i]==1) {
            base = i;
        }
    }
    return base;
}

int simulate_cluster(int size, const double prob[]){
    const gsl_rng_type *T;
    gsl_rng *r;
    size_t t = size;
    int cluster;
    const unsigned int n = 1;
    unsigned int theta[t];
    
    T = gsl_rng_ranlxs0;
    gsl_rng_default_seed = (rand());
    r = gsl_rng_alloc(T);
    gsl_ran_multinomial(r, t, n, prob, theta);
    for (int i = 0; i < t; i++) {
        if (theta[i]==1) {
            cluster = i;
        }
    }
    return cluster;
}


int count_match(int size,std::vector<int>& seq_1,std::vector<int>& seq_2){
    int count=0;
    for (int i=0; i<size; i++) {
        if (seq_1[i]==seq_2[i]) {
            count++;
        }
    }
    return count;
}







int main (){
    // ###############################################
    //set up parameters for simulation (the length of H0, mutation rate gamma, number of original haplotypes K, error rate, number of reads per original haplotypes).
    // ###############################################
    int L_h0 = 200;
    int K = 3;
    int h[3][100];
    int n = 100;
    double gamma = 0.6;
    double theta = 0.95;
    int r[3][n][100];
    int N_reads = K*n;
    int reads[N_reads][100];
    
    
    // ###############################################
    //simulate H0
    // ###############################################
    std::vector<int>h0(L_h0);
    double base_distribution[] = {0.25,0.25,0.25,0.25};
    double* p = base_distribution;
    std::cout << ">master_haplotype\n";
    for (int i=0; i < L_h0; i++){
        h0[i] = simulate_base(p);
        std::cout <<h0[i];
    }
    std::cout <<"\n";
//    int equal1 = std::count(h0.begin(), h0.end(), 1);
//    std::cout <<"number of 1 is "<<equal1<<"\n";
//    std::cout << std::count(h0.begin(),h0.end(),1)<<"\n";


    
    // ###############################################
    // simulate original haplotypes h1-hk
    // ###############################################
    double error_rate = (1-gamma)/3;
    for (int j=0; j < L_h0; j++){
        double prob[] = {error_rate,error_rate,error_rate,error_rate};
        prob[h0[j]] = gamma;
        double* pp = prob;
        for (int i=0; i < K; i++){
            h[i][j] = simulate_base(pp);
        }
    }
    
    for (int i=0; i < K; i++) {
        std::cout <<"<h"<<i<<"\n";
        for (int j=0; j < L_h0; j++) {
            std::cout << h[i][j];
        }
        std::cout <<"\n";
    }
    
    // ###############################################
    // simulate reads
    // ###############################################
    double error_rate2 = (1-theta)/3;
    for (int i = 0; i<K; i++) {
        for (int j=0; j<L_h0; j++) {
             double prob2[] = {error_rate2,error_rate2,error_rate2,error_rate2};
            prob2[h[i][j]] = theta;
            //std::cout << "\n"<<prob2[0]<<" "<<prob2[1]<<" "<<prob2[2]<<" "<<prob2[3]<<"\n";
            double* rp = prob2;
            for (int m=0; m<n; m++) {
                r[i][m][j] = simulate_base(rp);
            }
        }
    }
    
    int index = 0;
    for (int i=0; i<K; i++) {
        for (int j=0; j<n; j++) {
            std::cout << ">r" << index <<"\n";
            
            for (int m=0; m<L_h0; m++) {
                reads[index][m]=r[i][j][m];
                std::cout << r[i][j][m];
            }
            index++;
            std::cout <<"\n";
        }
    }
    
    //compair sequence
    std::vector<int> temp(L_h0);
        for (int i=0; i<K; i++) {
            for (int j=0; j<L_h0; j++) {
                temp[j] = h[i][j];
            }
            int match = count_match(L_h0, h0, temp);
            std::cout << "master haplotype to h" <<i <<":  "<< match<<"\n";
        }
    std::vector<int> hap1(L_h0);
    for (int j=0; j<L_h0; j++) {
        hap1[j] = h[1][j];
    }

    std::vector<int> temp_seq(L_h0);
    std::vector<long> P_r_master(L_h0);
    for (int i=0; i<N_reads; i++) {
        for (int j=0; j<L_h0; j++) {
            temp_seq[j] = reads[i][j];
        }
// check reads similarity to original haplotypes
//        for (int j=0; j<L_h0; j++) {
//            std::cout<<h0[j];
//        }
//        std::cout <<"\n";
//        for (int j=0; j<L_h0; j++) {
//            std::cout<<temp_seq[j];
//        }
//        std::cout <<"\n";
        
// check reads similarity to master haplotypes
        P_r_master[i]=count_match(L_h0, temp_seq, h0);
//        std::cout <<P_r_master[i]<<"\n";
    }
    
    // ###############################################
    // ###############################################
    // Dirichlet prepare
    // ###############################################
    // ###############################################
    
    
    // prior for master haplotype
    std::cout << "what's prior master?\n";
    std::vector<int> prior_master(L_h0);
    for (int i=0; i<L_h0; i++) {
        prior_master[i]=h0[i];
        std::cout <<prior_master[i];
    }
     std::cout << "what's prior master?\n";
    
    // starting value for original haplotypes (random sequence generated from (0.25,0.25,0.25,02.5))
    int prior_K = 5;
    std::vector<std::vector<int> > prior_org_hap(10000,std::vector<int>(L_h0));
    for (int i=0; i<prior_K; i++) {
        std::cout << "\n hap " << i << "\n";
        for (int j=0; j<L_h0; j++) {
            prior_org_hap[i][j] = simulate_base(p);
            std::cout <<prior_org_hap[i][j];
        }
    }
    
    for (int i=0; i<L_h0; i++) {
        prior_org_hap[0][i]=h[0][i];
    }
    std::cout << "\nbegin Gibbs\n";
    // prior for theta, gamma,alpha
    double prior_theta = 0.8;
    double prior_gamma = 0.7;
    double alpha = 1.5;
    
    //starting value for cluster assign
    std::cout << "\n This is starting assign for cluster of each reads. \n";
    std::vector<int> cluster_Ci(N_reads);
    double prob_cluster[] = {0.5,0.5,0.5,0.5,0.5};
    double* pc = prob_cluster;
    for (int i=0; i<N_reads; i++) {
        cluster_Ci[i]=simulate_cluster(prior_K, pc);
        std::cout <<cluster_Ci[i]<<"\t";
    }
    
    
    
    
    // ###############################################
    // Gibbs sampler
    // ###############################################
    int generation = 500;
    int k_star = prior_K;
    double theta_star = prior_theta;
    double gamma_star = prior_gamma;
    double prob_in_cluster[k_star+1];
    double* pointer_prob = prob_in_cluster;
//    std::vector<double> trans_prob_in_cluster(k_star+1);
    double trans_prob_in_cluster[k_star+1];
    int N_in_cluster[100000];
    size_t index_min;
    double min;
    
    //caculate P(read|prior master haplotype) for each read
    std::vector<int> P_r_count(N_reads);
    double a1 = theta_star*gamma_star+(1-gamma_star)*(1-theta_star)/3;
    double a2 = (theta_star+gamma_star+3*(1-gamma_star*theta_star)-2)/9;
    for (int i=0; i<N_reads; i++) {
//        std::cout << "This is read " <<i<<"\n";
        for (int j=0; j<L_h0; j++) {
            temp_seq[j] = reads[i][j];
        }
            P_r_count[i]=count_match(L_h0, temp_seq, h0);
            P_r_master[i] = P_r_count[i]*gsl_sf_log(a1)+ (L_h0-P_r_count[i])*gsl_sf_log(a2);
            std::cout <<P_r_count[i]<<"\t"<<P_r_master[i]<<"\n";
    }
    
    
    
    for (int g=0; g<generation; g++) {
        std::cout <<"####################### This is generation "<<g<<"\n";
        
        // ###############################################
        //update original haplotypes
        // ###############################################
        int cluster_match_table[k_star][L_h0][4];
        double log_prob_base[4];
        double prob_base[4];
        double* pointer_base;
        double min_log_prob;
        int min_index;
        for (int i=0; i<k_star; i++) {
            if (std::count(cluster_Ci.begin(),cluster_Ci.end(),i)!=0) {
                // find the index of reads which belong to this cluster i
                for (int j=0; j<N_reads; j++) {
                    if (cluster_Ci[j]==i) {
                        for (int m=0; m<L_h0; m++) {
                            switch (reads[j][m]) {
                                case 0:
                                    cluster_match_table[i][j][0]++;
                                    break;
                                case 1:
                                    cluster_match_table[i][j][1]++;
                                    break;
                                case 2:
                                    cluster_match_table[i][j][2]++;
                                    break;
                                case 3:
                                    cluster_match_table[i][j][3]++;
                                default:
                                    break;
                            }
                        }
                    }
                }
                
                for (int j=0; j<L_h0; j++) {
                    for (int m=0; m<4; m++) {
                        log_prob_base[m] = cluster_match_table[i][j][m]*gsl_sf_log(theta_star)+(L_h0-cluster_match_table[i][j][m])*gsl_sf_log((1-theta_star)/3);
                    }
                    index_min = gsl_stats_min_index(log_prob_base, 1, 4);
                    min_log_prob = gsl_stats_min(log_prob_base, 1, 4);
                    for (int m=0; m<4; m++) {
                        prob_base[m]=log_prob_base[m]-min_log_prob;
                    }
                    for (int m=0; m<4; m++) {
                        if (prob_base[m]==0){
                            prob_base[m]=1;
                        }
                    }
                    
                    pointer_base = prob_base;
                    prior_org_hap[i][j]=simulate_cluster(4, pointer_base);
                }
                
                
            }
        }
        std::cout << "Updated haplotypes number: " << k_star<<"\n";
//        for (int i=0; i<k_star; i++) {
//            if (std::count(cluster_Ci.begin(),cluster_Ci.end(),i)!=0) {
//                std::cout << "hap"<<i<<"\t";
//                for (int j=0;j<L_h0; j++) {
//                    std::cout<<prior_org_hap[i][j];
//                }
//                std::cout << "\n";
//            }
//            
//        }
//        std::cout << "\n";

        // ###############################################
        // update cluster_Ci
        // ###############################################
        for (int i=0; i<N_reads; i++) {
//            std::cout <<"\n\t This is read "<<i<<"\n";
            
            for (int j=0; j<k_star; j++) {
                if (std::count(cluster_Ci.begin(),cluster_Ci.end(),j)==0) {
                    for (int m=j+1; m<k_star; m++) {
                        for (int l=0; l<L_h0; l++) {
                            prior_org_hap[m-1][l]=prior_org_hap[m][l];
                        }
                    }
                    k_star--;
                    for (int m=0; m<N_reads; m++) {
                        if (cluster_Ci[m]>j){
                            cluster_Ci[m]--;
                        }
                    }
                }
                
            }

//            std::cout <<"\n Modified cluster assignment\n ";
//            for (int j=0; j<N_reads; j++) {
//                std::cout <<cluster_Ci[j]<<"\t";
//            }
//            std::cout <<"\n";
            //count how many reads in each cluster
            for (int j=0; j<k_star; j++) {
                N_in_cluster[j]=std::count(cluster_Ci.begin(),cluster_Ci.end(),j);
            }
            N_in_cluster[cluster_Ci[i]]--;
            for (int j=0; j<k_star; j++) {
//                std::cout << N_in_cluster[j]<<"\t";
            }
            
            //caculate probability of read i come from each cluster
//            std::cout <<"\n";
            
            std::vector<int> temp_seq1(L_h0);
            std::vector<int> temp_seq2(L_h0);
            std::vector<int> match_in_cluster(k_star);
            for (int m=0; m<L_h0; m++) {
                temp_seq1[m]=reads[i][m];
            }
            for (n=0; n<k_star; n++) {
                for (int l=0; l<L_h0; l++) {
                    temp_seq2[l]=prior_org_hap[n][l];
                }
                match_in_cluster[n]=count_match(L_h0, temp_seq1,temp_seq2);
//                std::cout << match_in_cluster[n]<<"\t";
            }
                
//            std::cout <<"\n";
            
            //assign new cluster for read i
            for (int j=0; j<k_star; j++) {
                if (N_in_cluster[j]==0) {
                    prob_in_cluster[j]=0;
                }else{
                prob_in_cluster[j]=gsl_sf_log(N_in_cluster[j]/(N_reads-1+alpha))+match_in_cluster[j]*gsl_sf_log(theta_star)+(L_h0-match_in_cluster[j])*gsl_sf_log((1-theta_star)/3);
                }
                
            }
            prob_in_cluster[k_star] = gsl_sf_log(alpha/(N_reads-1+alpha))+P_r_master[i];
            for (int j=0; j<k_star+1; j++) {
 //               std::cout << prob_in_cluster[j]<<"\t";
            }
        
 //            std::cout <<"\n";
            index_min = gsl_stats_min_index(prob_in_cluster, 1, k_star+1);
            min = gsl_stats_min(prob_in_cluster, 1, k_star+1);
            for (int j=0; j<k_star+1; j++) {
                if (prob_in_cluster[j]==0) {
                    trans_prob_in_cluster[j]=-1;
                }else{
                trans_prob_in_cluster[j]=prob_in_cluster[j]-min;
                }
            }
            for (int j=0; j<k_star+1; j++) {

                if (trans_prob_in_cluster[j]==0) {
                    trans_prob_in_cluster[j]=1;
                }
                if (trans_prob_in_cluster[j]==-1) {
                    trans_prob_in_cluster[j]=0;
                }
//                std::cout << trans_prob_in_cluster[j]<<"\t";
            }
//            std::cout <<"\n";
            cluster_Ci[i]=simulate_cluster(k_star+1, trans_prob_in_cluster);
//            std::cout <<"new class is "<<cluster_Ci[i]<<"\n";
            
            //get new cluster number k_star
            if (cluster_Ci[i]==k_star) {
                k_star++;
                for (int j=0; j<L_h0; j++) {
                    prior_org_hap[k_star-1][j] = reads[i][j];
                }
            }
//            std::cout <<"\n original cluster assignment\n ";
//            for (int j=0; j<N_reads; j++) {
//                std::cout <<cluster_Ci[j]<<"\t";
//            }
//            std::cout <<"\n";
            
        }
        
        std::cout << "Updated cluster assignment is \n";
        for (int i=0; i<N_reads; i++) {
            std::cout <<cluster_Ci[i]<<"\t";
        }
        std::cout <<"\n";
        
                
        // ###############################################
        //update theta, gamma
        // ###############################################
        int total_match=0;
        int actural_k=0;
        int hap_match=0;
        std::vector<int> temp_hap(L_h0);
        for (int i=0; i<N_reads; i++) {
            std::vector<int> temp_seq1(L_h0);
            std::vector<int> temp_seq2(L_h0);
            for (int j=0; j<L_h0; j++) {
                temp_seq1[j]=reads[i][j];
                temp_seq2[j]=prior_org_hap[cluster_Ci[i]][j];
            }
//            std::cout << count_match(L_h0, temp_seq1, temp_seq2)<<"\n";
            total_match = total_match + count_match(L_h0, temp_seq1, temp_seq2);
        }
//        std::cout << "total match is"<<total_match<<"\n";
//        std::cout << "\n N_reads*L_h0 is "<< N_reads <<"\t"<<L_h0 <<"\n";
        theta_star = double(total_match)/ double(N_reads*L_h0);
        for (int i=0; i<k_star; i++) {
            
            if (std::count(cluster_Ci.begin(),cluster_Ci.end(),i)!=0) {
                actural_k++;
                for (int j=0; j<L_h0; j++) {
                    temp_hap[j]=reads[i][j];
                }
                hap_match = hap_match + count_match(L_h0, temp_hap, h0);
            }
        }
        gamma_star = double(hap_match)/double(actural_k*L_h0);
//        std::cout << "cluster number is " << actural_k<<"\n";
//        std::cout << "hap match is "<<hap_match<<"\n";
        actural_k = 0;
        std::cout <<"Updated theta is "<<theta_star<<"\n";
        std::cout <<"Updated gamma is "<<gamma_star<<"\n";

        
        
    }
    



 }

