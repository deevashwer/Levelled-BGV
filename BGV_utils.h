#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define PARI_OLD_NAMES

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <pari/pari.h>
//#include <pari.h>
#include <time.h>
#include <vector>
#include <string.h>
#include <sys/time.h>

struct timeval tv;

void print(GEN x){
    std::cout << GENtostr(x) << std::endl;
}

double Uniform(void) {
    return ((double) rand() + 1.0) / ((double) RAND_MAX + 2.0);
}

double Normal(void) {
    return sqrt(-log(Uniform())*2.0) * sin(2.0 * M_PI * Uniform());
}

double Gauss(double mu, double sigma) {
    double z = sqrt(-2.0 * log(Uniform())) * sin(2.0 * M_PI * Uniform());
    return mu + sigma*z;
}

GEN Sample(int n, double sigma) {
    GEN ret = cgetg(n + 1, t_COL);
    double z;
    int i;
    
    for (i = 1; i <= n; i++) {
        z = Gauss(0, sigma);
        z = abs(round(z)); /*absolute value of Gaussian distribution */
        ret[i] = (long) stoi((long) z);
    }
    
    return ret;
}

GEN generate_random(int bit_length){
    gettimeofday(&tv, NULL);
    setrand(stoi(tv.tv_usec + tv.tv_sec*1000000));
    GEN r = randomi(gshift(gen_1, bit_length));
    return r;
}

GEN decimal_to_binary(GEN x, int bit_size){
    GEN binary_m = cgetg(bit_size + 1, t_COL);
    if(gcmp(x, gen_0) == 0){
        for(int i = 0; i < bit_size; i++)
            gel(binary_m, i + 1) = gen_0;
        return binary_m;
    }
    std::vector<bool> bits;
    GEN temp, m = x;
    while(mpcmp(m, gen_1) > 0){
        if(mpodd(m) == 1)
            bits.push_back(true);
        else
            bits.push_back(false);
        m = gdivexact(m, gen_2);
    }
    bits.push_back(true);
    int i = 1;
    std::vector<bool>::iterator it;
    for(it = bits.begin(); it != bits.end(); it++)
        if(*it == true)
            gel(binary_m, i++) = gen_1;
        else
            gel(binary_m, i++) = gen_0;
    bits.clear();
    for(i; i < bit_size + 1; i++)
        gel(binary_m, i) = gen_0;
    return binary_m;
}

GEN bit_decomposition_utility(GEN x, int bit_size, int size){
    GEN binary_x = cgetg(bit_size + 1, t_VEC);
    GEN temp = cgetg(size + 1, t_VEC);
    GEN vec = gtovecrev(lift(lift(x)));
    int k = lg(vec) - 1;
    for(int i = 0; i < k; i++)
        gel(temp, i + 1) = gel(vec, i + 1);
    for(int i = k; i < size; i++)
        gel(temp, i + 1) = gen_0;
    GEN tmp = cgetg(size + 1, t_MAT);
    for(int i = 0; i < size; i++)
        gel(tmp, i + 1) = decimal_to_binary(gel(temp, i + 1), bit_size);
    tmp = shallowtrans(tmp);
    for(int i = 0; i < bit_size; i++)
        gel(binary_x, i + 1) = gtopolyrev(gel(tmp, i + 1), -1);
    return binary_x;
}

GEN bit_decomposition(GEN x, int log_q, int n, int d){
    GEN u = cgetg(log_q*n + 1, t_COL);
    GEN tmp;
    for(int i = 0; i < n; i++){
        tmp = bit_decomposition_utility(gel(x, i + 1), log_q, d);
        for(int j = 0; j < log_q; j++)
            gel(u, j*(n) + i + 1) = gel(tmp, j + 1);
    }
    return u;
}

GEN PowersOf2(GEN x, int log_q, int n){
    GEN u = cgetg(log_q*n + 1, t_COL);
    GEN tmp = x;
    for(int i = 0; i < log_q; i++){
        for(int j = 0; j < n; j++)
            gel(u, i*n + j + 1) = gel(tmp, j + 1);
        tmp = gmul(gen_2, tmp);
    }
    return u;
}








