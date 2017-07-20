#include "BGV_utils.h"

#define sigma 2400//8
#define mu 60
#define var 26

pari_sp ltop, lbot;

using namespace std;

class cryptosystem{
public:
    int n, d, L, N;
    GEN q, f, s, t, A, T, tensor_product_s;
    
    cryptosystem(bool b, int l){
        if(b == 0)
            d = 1, n = 20;
        else
            n = 1, d = 20;
        L = l;
        f = cgetg(d + 2, t_VEC);
        gel(f, 1) = gen_1;
        gel(f, d + 1) = gen_1;
        for(int i = 1; i < d; i++)
            gel(f, i + 1) = gen_0;
        f = gtopolyrev(f, -1);
        N = 2*n*mu;
        q = cgetg(L + 1, t_VEC);
        s = cgetg(L + 1, t_VEC);
        t = cgetg(L + 1, t_VEC);
        A = cgetg(L + 1, t_VEC);
        T = cgetg(L, t_VEC);
        tensor_product_s = cgetg(L + 1, t_VEC);
        for(int i = L; i > 0; i--){
            gel(q, i) = gshift(gen_1, mu + var*i);
            gel(q, i) = nextprime(gel(q, i));
            generate_secret_key(i);
            generate_public_key(i);
            generate_tensor_products(i);
            if(i != L)
                switch_key_gen(i, (n + 1)*(n + 1), mu + (i + 1)*var);
        }
    }
    
    GEN sample_error_polynomial(long variance, int j){
        ltop = avma;
        GEN tmp = Sample(d, variance);
        tmp = gmodulo(tmp, gel(q, j));
        //tmp = gmodulo(tmp, gel(q, L));
        tmp = gtopolyrev(tmp, -1);
        tmp = gmodulo(tmp, f);
        //tmp = lift(lift(tmp));
        lbot = avma;
        tmp = gerepilecopy(ltop, tmp);
        return tmp;
    }
    
    GEN sample_polynomial(int j){
        ltop = avma;
        GEN tmp = cgetg(d + 1, t_COL);
        for(int i = 0; i < d; i++)
            gel(tmp, i + 1) = generate_random(mu + var*j);
        tmp = gmodulo(tmp, gel(q, j));
        //tmp = gmodulo(tmp, gel(q, L));
        tmp = gtopolyrev(tmp, -1);
        tmp = gmodulo(tmp, f);
        lbot = avma;
        tmp = gerepilecopy(ltop, tmp);
        return tmp;
    }
    
    void generate_secret_key(int j){
        gel(s, j) = cgetg(n + 2, t_COL);
        gel(gel(s, j), 1) = gen_1;
        for(int i = 0; i < n; i++)
            gel(gel(s, j), i + 2) = sample_error_polynomial(sigma, j);
        gel(t, j) = cgetg(n + 1, t_COL);
        for(int i = 0; i < n; i++)
            gel(gel(t, j), i + 1) = gel(gel(s, j), i + 2);
        //print(gel(s, j));
    }
    
    void generate_public_key(int j){
        GEN B = cgetg(n + 1, t_MAT);
        for(int i = 0; i < n; i++){
            gel(B, i + 1) = cgetg(N + 1, t_COL);
            for(int k = 0; k < N; k++){
                //cout << i << " " << j << endl;
                gel(gel(B, i + 1), k + 1) = sample_polynomial(j);
            }
        }
        GEN e = cgetg(N + 1, t_COL);
        for(int i = 0; i < N; i++)
            gel(e, i + 1) = sample_error_polynomial(sigma, j);
        e = gmul(gen_2, e);
        GEN b = gmul(B, gel(t, j));
        b = gadd(b, e);
        gel(A, j) = cgetg(n + 2, t_MAT);
        for(int i = 0; i < n + 1; i++)
            gel(gel(A, j), i + 1) = cgetg(N + 1, t_COL);
        for(int i = 0; i < N; i++)
            gel(gel(gel(A, j), 1), i + 1) = gel(b, i + 1);
        for(int i = 1; i < n + 1; i++)
            gel(gel(A, j), i + 1) = gmul(gen_m1, gel(B, i));
    }
    
    void generate_tensor_products(int j){
        gel(tensor_product_s, j) = cgetg((n + 1)*(n + 1) + 1, t_COL);
        for(int i = 0; i < n + 1; i++)
            for(int k = 0; k < n + 1; k++)
                gel(gel(tensor_product_s, j), i*(n + 1) + k + 1) = gmul(gel(gel(s, j), i + 1), gel(gel(s, j), k + 1));
    }
    
    void switch_key_gen(int j, int n_1, int log_q){
        gel(T, j) = cgetg(n + 2, t_MAT);
        int N_custom = n_1*log_q;
        GEN B = cgetg(n + 1, t_MAT);
        for(int i = 0; i < n; i++){
            gel(B, i + 1) = cgetg(N_custom + 1, t_COL);
            for(int k = 0; k < N_custom; k++)
                gel(gel(B, i + 1), k + 1) = sample_polynomial(j + 1);
        }
        GEN e = cgetg(N_custom + 1, t_COL);
        for(int i = 0; i < N_custom; i++)
            gel(e, i + 1) = sample_error_polynomial(sigma, j + 1);
        e = gmul(gen_2, e);
        GEN b = gmul(B, lift(lift(gel(t, j))));
        b = gadd(b, e);
        for(int i = 0; i < n + 1; i++)
            gel(gel(T, j), i + 1) = cgetg(N_custom + 1, t_COL);
        for(int i = 0; i < N_custom; i++)
            gel(gel(gel(T, j), 1), i + 1) = gel(b, i + 1);
        for(int i = 1; i < n + 1; i++)
            gel(gel(T, j), i + 1) = gmul(gen_m1, gel(B, i));
        //GEN test = gmul(gel(T, j), lift(lift(gel(s, j))));
        //print(lift(lift(test)));
        GEN tmp = PowersOf2(gel(tensor_product_s, j + 1), log_q, n_1);
        //GEN tmp = PowersOf2(gel(s, j + 1), log_q, n_1);
        //print(tmp);
        gel(gel(T, j), 1) = gadd(gel(gel(T, j), 1), tmp);
        return;
    }
    
    GEN switch_key(GEN ct_1, int j, int n_1, int log_q){
        GEN ct_2 = gmul(shallowtrans(gel(T, j)), bit_decomposition(ct_1, log_q, n_1, d));
        return ct_2;
    }
    
    GEN switch_moduli(GEN ct_1, int j){
        GEN ct_2 = lift(lift(ct_1));
        GEN tmp_1, tmp_2, tmp_3, tmp_4;
        for(int i = 0; i < n + 1; i++){
            tmp_1 = gtovecrev(gel(ct_2, i + 1));
            for(int k = 0; k < d; k++){
                gel(tmp_1, k + 1) = mulii(gel(tmp_1, k + 1), gel(q, j));
                gel(tmp_1, k + 1) = diviiround(gel(tmp_1, k + 1), gel(q, j + 1));
            }
            gel(ct_2, i + 1) = gtopolyrev(tmp_1, -1);
        }
        for(int i = 0; i < n + 1; i++){
            tmp_1 = gmod(gtovecrev(gel(ct_2, i + 1)), gel(q, j));
            tmp_2 = gtovecrev(gel(lift(lift(ct_1)), i + 1));
            tmp_4 = gmod(tmp_2, gen_2);
            tmp_3 = gmod(tmp_1, gen_2);
            for(int k = 0; k < d; k++)
                if(mpcmp(gel(tmp_3, k + 1), gel(tmp_4, k + 1)) != 0)
                    gel(tmp_1, k + 1) = gadd(gel(tmp_1, k + 1), gen_1);
            tmp_1 = gmodulo(tmp_1, gel(q, j));
            tmp_1 = gtopolyrev(tmp_1, -1);
            gel(ct_2, i + 1) = gmodulo(tmp_1, f);
        }
        return ct_2;
    }
    
    GEN encrypt(GEN m, int level){
        ltop = avma;
        GEN ct = cgetg(n + 2, t_COL);
        gel(ct, 1) = m;
        for(int i = 1; i < n + 1; i++)
            gel(ct, i + 1) = gen_0;
        GEN r = cgetg(N + 1, t_COL);
        for(int i = 0; i < N; i++)
            gel(r, i + 1) = sample_error_polynomial(1, level);
        GEN tmp = gmul(shallowtrans(gel(A, level)), r);
        ct = gadd(ct, tmp);
        ct = gerepilecopy(ltop, ct);
        return ct;
    }
    
    GEN decrypt(GEN ct, int level, bool show_error){
        GEN m = cgetg(n + 2, t_COL);
        m = gmul(shallowtrans(ct), lift(lift(gel(s, level))));
        m = lift(lift(m));
        m = gtovecrev(m);
        for(int i = 0; i < d; i++){
            //if(level < L && gcmp(gel(m, i + 1), gdiv(gel(q, level + 1), gen_2)) > 0)
            //    gel(m, i + 1) = gsub(gel(m, i + 1), gel(q, level + 1));
            if(gcmp(gel(m, i + 1), gdiv(gel(q, level), gen_2)) > 0)
                gel(m, i + 1) = gsub(gel(m, i + 1), gel(q, level));
            if(!show_error)
                gel(m, i + 1) = gmod(gel(m, i + 1), gen_2);
        }
        return m;
    }
    
    GEN addition(GEN ct_1, GEN ct_2){
        GEN ct = gadd(ct_1, ct_2);
        return ct;
    }
    
    GEN multiply(GEN ct_1, GEN ct_2){
        GEN ct = cgetg((n + 1)*(n + 1) + 1, t_COL);
        for(int i = 0; i < n + 1; i++)
            for(int j = 0; j < n + 1; j++)
                gel(ct, i*(n + 1) + j + 1) = gmul(gel(ct_1, i + 1), gel(ct_2, j + 1));
        return ct;
    }
};

class ciphertext{
public:
    GEN value;
    int level;
    cryptosystem* pkc;
    
    ciphertext(){};
    
    ciphertext(cryptosystem* PKC){
        pkc = PKC;
    }
    
    ciphertext(cryptosystem* PKC, GEN m){
        level = PKC->L;
        value = PKC->encrypt(m, level);
        pkc = PKC;
    }
    
    ciphertext(cryptosystem* PKC, GEN m, int lev){
        level = lev;
        value = PKC->encrypt(m, level);
        pkc = PKC;
    }
    
    GEN decrypt(){
        return pkc->decrypt(value, level, false);
    }
    
    ciphertext operator+(const ciphertext& ct_2){
        ciphertext result(pkc);
        if(this->level == ct_2.level){
            result.value = pkc->addition(this->value, ct_2.value);
            result.level = this->level;
        }
        else
            cout << "Error: The ciphertexts belong to different levels" << endl;
        return result;
    }
    
    ciphertext operator*(const ciphertext& ct_2){
        ciphertext result(pkc);
        if(this->level == ct_2.level){
            GEN temp = pkc->multiply(this->value, ct_2.value);
            temp = pkc->switch_key(temp, this->level - 1, (pkc->n + 1)*(pkc->n + 1), mu + (this->level)*var);
            result.value = pkc->switch_moduli(temp, this->level - 1);
            result.level = this->level - 1;
        }
        else
            cout << "Error: The ciphertexts belong to different levels" << endl;
        return result;
    }
    
    void initialize(cryptosystem* PKC, GEN m){
        level = PKC->L;
        value = PKC->encrypt(m, level);
        pkc = PKC;
    }
    
    void initialize(cryptosystem* PKC, GEN m, int lev){
        level = lev;
        value = PKC->encrypt(m, level);
        pkc = PKC;
    }
    
    void custom_setup(GEN val, int lev, cryptosystem* PKC){
        value = val;
        level = lev;
        pkc = PKC;
    }
    
    void print(){
        cout << "Error (Level:" << level << "): " <<  GENtostr(pkc->decrypt(value, level, true)) << endl;
        cout << "Plaintext: " <<  GENtostr(pkc->decrypt(value, level, false)) << endl;
    }
};

int main(int argc, const char * argv[]) {
    pari_init(20000000000, 2);
    srand(time(NULL));
    cryptosystem pkc(true, 10);
    ciphertext ct[9];
    ciphertext result(&pkc, gen_1);
    for(int i = 0; i < 9; i++){
        ct[i].initialize(&pkc, gen_1, 10 - i);
        result = result * ct[i];
        result.print();
    }
    pari_close();
    return 0;
}
