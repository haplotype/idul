/////////////////////////////////////////////////////////////////////////////////
//The MIT License (MIT)
//
//Copyright (c) 2023 Yongtao Guan
//  Bug report: ytguan@gmail.com 
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

#include <random> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <algorithm> 
#include <vector> 
#include <stdio.h> 
#include <stdlib.h> 
#include <cstring> 
#include <string.h> 
#include <unistd.h>
//#include <sys/stat.h>
//#include <sys/types.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"              
#include <map> 
#include <cmath>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <zlib.h>
#include "thpool.h"
//#include <omp.h>

#define EIGEN_USE_OPENMP
//#define EIGEN_USE_LAPACKE

#include <Eigen/Dense>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
//#include <unsupported/Eigen/SelfAdjointEigenSolver>

#include <boost/iostreams/filtering_stream.hpp> 
#include <boost/iostreams/device/file_descriptor.hpp> 
#include <boost/iostreams/filter/gzip.hpp> 

//#include <boost/random/uniform_real_distribution.hpp>
//#include <boost/random/mersenne_twister.hpp>

#ifdef EIGEN_USE_LAPACKE
#include <lapacke.h>
#endif


using namespace std; 

#define BUFFER_SIZE  4096 
#define VERSION "0.81"
//6 Apr 2023  0.5
//29 July 2023  0.51
//15 Dec 2023 0.71  //bf and conditioning analysis. option changes -a -b -q 
//14 Mar 2024 0.81  //revised IDUL; 

typedef struct Stats {
    double pv; 
    int niter; 
    double eta; 
    double beta; 
    double sigma; 
    double l0; 
    double l1; 
    double logbf; 
} Stats; 

typedef struct KData {
    int ctrl; 
    int ni;      //number of individuals
    int ns;      //number of snps. 
    int nth;     //number of threads
    Eigen::MatrixXd evec; 
    Eigen::VectorXd eval; 
    Eigen::MatrixXd ww; 
    Eigen::MatrixXd wtw; 
    Eigen::MatrixXd mph; 

    int nph; //number of phenotypes; 
    vector<float> v8; 
    vector<float> vu; //untransmitted genotypes; 
    int flag_kk; 
    int m_lrt; 
    int m_wald; 
    int flag_untransmitted;

    double * veta0; 
    double * vloglike0; //these are sets of values are for likelihood ratio test. 
    
    Stats * pstats; 
} KData; 

KData dadu; 

typedef void (*FuncPtr)(void *); 

void print_progress_num(int last, int p)
{
    if(!last) {
	printf("##processed variants = %d \r", p); 
	fflush(stdout); 	
    } else {
	printf("##processed variants = %d \n", p); 
	fflush(stdout); 	
    }
}

void print_progress_bar(int last, long p, long total)
{
    char str[] = "##processed pairs:"; 
    int progress = (int) (100.0 * p / total); 
//    int barsize = (int) (progress / 2.0); 
//    char bar[100];
//    memset(bar, '\0', 100); 
    if(!last) {
//	    for (int i = 0; i < barsize; i++)
//		    bar[i] = '>'; 
	    printf("%s %ld or %d%%\r", str, p, progress); 
	    fflush(stdout); 	
    } else {
//	    for (int i = 0; i < barsize; i++)
//		    bar[i] = '>'; 
//	    printf("%s [%-50s] 100%%\n", str, bar); 
	    printf("%s %ld or %d%%\n", str, p, progress); 
    }
}


void read_from_gzipped_file(string fn, vector<double> &val, int &rows) {
    gzFile file = gzopen(fn.c_str(), "rb");
    if (file == NULL) {
        std::cerr << "Failed to open file: " << fn << std::endl;
        exit(1);
    }

    rows = 0;
    char buffer[BUFFER_SIZE];
    std::string decompressed_data;

    // Read and decompress the entire file into a string
    int bytes_read;
    while ((bytes_read = gzread(file, buffer, BUFFER_SIZE - 1)) > 0) {
        buffer[bytes_read] = '\0';
        decompressed_data.append(buffer);
    }
    gzclose(file);

    // Process the decompressed data line by line
    std::stringstream ss(decompressed_data);
    std::string line;
    while (std::getline(ss, line)) {
        // Ignore lines starting with #
        if (!line.empty() && line[0] == '#') {
            continue;
        }

	std::stringstream line_stream(line);
	double value;
	while (line_stream >> value) {
	    val.push_back(value); 
	}
        rows++;
    }
}

void read_from_gzipped_file(string fn, vector<double> &val) {
    gzFile file = gzopen(fn.c_str(), "rb");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file %s\n", fn.c_str());
        exit(1);
    }

    char buffer[BUFFER_SIZE + 1];
//    int value_count = 0;
    int buffer_offset = 0;
    int bytes_read;

//    int ii = 0; 
//    int jj = 0; 
    while ((bytes_read = gzread(file, buffer + buffer_offset, BUFFER_SIZE - buffer_offset)) > 0) {
        buffer[bytes_read + buffer_offset] = '\0';
        char *current = buffer;
        char *next = NULL;

	float temp; 
        while ((next = strpbrk(current, " \t\n")) != NULL) {
            *next = '\0';
	    if (strlen(current) > 0) {
		sscanf(current, "%f", &temp);
		val.push_back(temp); 
	    }
	    current = next + 1;
        }
        buffer_offset = strlen(current);
        memmove(buffer, current, buffer_offset);
    }
    gzclose(file);
}

using namespace boost::math;


void herit(int flag, string fn) 
{
    int ni = dadu.ni;    
    int cols = dadu.ww.cols(); 

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd d1(ni); 
    Eigen::VectorXd d2(ni);
    Eigen::MatrixXd X(ni, dadu.ww.cols()); 
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    
    Eigen::MatrixXd xtx;
    Eigen::VectorXd beta1;    

    Eigen::MatrixXd X2(ni, 2); 
    Eigen::MatrixXd xtx2;

    Eigen::MatrixXd X3(ni, dadu.ww.cols()); 

    gzFile fp1 = gzopen(fn.c_str(), "w");
    if (fp1 == NULL) {
	    fprintf(stderr, "can't open file %s\n", fn.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs pvals etc. 

    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	Eigen::VectorXd curph = dadu.mph.col(p); 
//	double pval; 
	int steps = 0; 
	int max_iter = 100;    
	double eta0=1, eta1=0.1; 
    //    double beta0, beta1; 
	double diff = 1; 
	while (fabs(diff)>1e-6 && steps<max_iter) 
	{
	    steps++;  
	    eta0 = eta1; 
	    //y=x beta + e; 
	    v1.setOnes(); 
//	    double h = eta0/(1+eta0); 
	    Eigen::VectorXd dd = eta0 * dadu.eval + v1; 
//	    dd = dd / (eta0 + 1); 
	    d2 = dd.cwiseInverse(); 
	    d1 = d2.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    for (int i = 0; i < ww2.cols(); i++)
		ww2.col(i).cwiseProduct(d1); 
	    X << ww2; 
	    xtx = X.transpose() * X; 
	    Eigen::VectorXd y1 = curph.cwiseProduct(d1); 
	    Eigen::VectorXd xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    //here e1 needs to scale back by multiplying 1/d1; 
	    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  

	    
	    if(flag == 2) { //mle 
		//r2 = (1, D) beta2 + err;                  
	//	Eigen::VectorXd y2 = r2.cwiseProduct(d2);  
		v1.setOnes(); 
		X2.col(0) = v1.cwiseProduct(d2); 
		X2.col(1) = dadu.eval.cwiseProduct(d2); 
		Eigen::MatrixXd xtx2 = X2.transpose() * X2; 
		Eigen::VectorXd xty2 = X2.transpose() * r2; 
		Eigen::VectorXd beta2 = xtx2.ldlt().solve(xty2);
	//	Eigen::VectorXd beta2 = xtx2.inverse() * xty2;
		double t1 = beta2(1) / max(1e-20, beta2(0));  
		double t2 = beta2(1) / (r2.sum() / ni);  
		eta1 = (t1 + t2)/2.0; 
		if(eta1 >= 1e5) eta1 = 1e5; 
		else if (eta1 <= 1e-5) eta1 = 1e-5;
	    }
	    else { //reml
		double rh1r = r2.sum(); 
		double rh2r = r2.cwiseProduct(d2).sum(); 
		double trh1 = d2.sum(); 
		double trh2 = d2.cwiseProduct(d2).sum();
	
		double mh = trh1/ni; 
		double vh = trh2/ni - mh*mh; 
	
		X3 << dadu.ww.cwiseProduct(d2); 
		xtx2 = X3.transpose() * X3; 
		Eigen::MatrixXd mat = xtx.inverse() * xtx2; 
		mat.trace(); 
		double fp1 = (trh1 - mat.trace() - (ni-cols)*rh2r/rh1r); 
	
		eta1 = eta0 + fp1 * eta0 /ni / (vh); 
		if(eta1 >= 1e5) eta1 = 1e5; 
		else if (eta1 <= 1e-5) eta1 = 1e-5;
	    }
	    diff = eta1/(1+eta1) - eta0/(1+eta0); 
	}
	gzprintf(fp1, "%.6e \n", eta1/(1+eta1)); 
    }
    gzclose(fp1); 
}

void loglike0(void) 
{
    int ni = dadu.ni;    
//    int cols = dadu.ww.cols(); 

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd d1(ni); 
    Eigen::VectorXd d2(ni);
    Eigen::MatrixXd X(ni, dadu.ww.cols()); 
    Eigen::MatrixXd X2(ni, 2); 
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    
    Eigen::MatrixXd xtx;
    Eigen::VectorXd beta1;    
    Eigen::VectorXd r1(ni);    
    Eigen::VectorXd r0(ni);    

    Eigen::MatrixXd X3(ni, dadu.ww.cols()); 
    Eigen::MatrixXd xtx2;

    dadu.veta0 = new double[dadu.mph.cols()]; 
    dadu.vloglike0 = new double[dadu.mph.cols()]; 
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	Eigen::VectorXd curph = dadu.mph.col(p); 
//	double pval; 
	int steps = 0; 
	int max_iter = 100;    
	double eta0=1, eta1=0.1; 
    //    double beta0, beta1; 
	double diff = 1; 
	while (fabs(diff)>1e-6 && steps<max_iter) 
	{
	    steps++;  
	    eta0 = eta1; 
	    //y=x beta + e; 
	    v1.setOnes(); 
	    Eigen::VectorXd dd = eta0 * dadu.eval + v1; 
	    d2 = dd.cwiseInverse(); 
	    d1 = d2.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << ww2.array().colwise() * d1.array(); 
	    xtx = X.transpose() * X; 
	    Eigen::VectorXd y1 = curph.cwiseProduct(d1); 
	    Eigen::VectorXd xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    //here e1 needs to scale back by multiplying 1/d1; 
	    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  

	    
	    //r2 = (1, D) beta2 + err;                  
	    v1.setOnes(); 
	    X2.col(0) = v1.cwiseProduct(d2); 
	    X2.col(1) = dadu.eval.cwiseProduct(d2); 
	    Eigen::MatrixXd xtx2 = X2.transpose() * X2; 
	    Eigen::VectorXd xty2 = X2.transpose() * r2; 
	    Eigen::VectorXd beta2 = xtx2.ldlt().solve(xty2);
	    double t1 = beta2(1) / max(1e-20, beta2(0));  
	    double t2 = beta2(1) / (r2.sum() / ni);  
	    eta1 = (t1 + t2)/2.0; 
	    if(eta1 >= 1e5) eta1 = 1e5; 
	    else if (eta1 <= 1e-5) eta1 = 1e-5;
	    diff = eta1/(1+eta1) - eta0/(1+eta0); 
	}
        dadu.veta0[p] = eta1; 
	dadu.vloglike0[p] = 0; 
	for (int i = 0; i < ni; i++)
	    dadu.vloglike0[p] += log(d1(i)); 
	dadu.vloglike0[p] -=  log(r2.sum()) * ni / 2.0;   
    }
}

void logpost0(void) 
{
    int ni = dadu.ni;    
//    int cols = dadu.ww.cols(); 

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd d1(ni); 
    Eigen::VectorXd d2(ni);
    Eigen::MatrixXd X(ni, dadu.ww.cols()); 
    Eigen::MatrixXd X2(ni, 2); 
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    
    Eigen::MatrixXd xtx;
    Eigen::VectorXd beta1;    
    Eigen::VectorXd r1(ni);    
    Eigen::VectorXd r0(ni);    

    Eigen::MatrixXd X3(ni, dadu.ww.cols()); 
    Eigen::MatrixXd xtx2;

    Eigen::VectorXd prior1(dadu.ww.cols());    
    prior1.setOnes(); 
    prior1 *= 1e-10; 

    dadu.veta0 = new double[dadu.mph.cols()]; 
    dadu.vloglike0 = new double[dadu.mph.cols()]; 
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	Eigen::VectorXd curph = dadu.mph.col(p); 
//	double pval; 
	int steps = 0; 
	int max_iter = 100;    
	double eta0=1, eta1=0.1; 
    //    double beta0, beta1; 
	double diff = 1; 
	while (fabs(diff)>1e-6 && steps<max_iter) 
	{
	    steps++;  
	    eta0 = eta1; 
	    //y=x beta + e; 
	    v1.setOnes(); 
	    Eigen::VectorXd dd = eta0 * dadu.eval + v1; 
	    d2 = dd.cwiseInverse(); 
	    d1 = d2.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << ww2.array().colwise() * d1.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior1.matrix().asDiagonal(); 
	    Eigen::VectorXd y1 = curph.cwiseProduct(d1); 
	    Eigen::VectorXd xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    //here e1 needs to scale back by multiplying 1/d1; 
	    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  

	    
	    //r2 = (1, D) beta2 + err;                  
	    v1.setOnes(); 
	    X2.col(0) = v1.cwiseProduct(d2); 
	    X2.col(1) = dadu.eval.cwiseProduct(d2); 
	    Eigen::MatrixXd xtx2 = X2.transpose() * X2; 
	    Eigen::VectorXd xty2 = X2.transpose() * r2; 
	    Eigen::VectorXd beta2 = xtx2.ldlt().solve(xty2);
	    double t1 = beta2(1) / max(1e-20, beta2(0));  
	    double t2 = beta2(1) / (r2.sum() / ni);  
	    eta1 = (t1 + t2)/2.0; 
	    if(eta1 >= 1e5) eta1 = 1e5; 
	    else if (eta1 <= 1e-5) eta1 = 1e-5;
	    diff = eta1/(1+eta1) - eta0/(1+eta0); 
	}
        dadu.veta0[p] = eta1; 
	dadu.vloglike0[p] = 0; 
	dadu.vloglike0[p] -= 0.5 * log(xtx.determinant()); 
	dadu.vloglike0[p] -=  log(r2.sum()) * ni / 2.0;   
    }
}

void jacquard_pv2(void *par)    //parent of origin joint test;
{
    int ni = dadu.ni;    
    long tiktok = (long) par; 
    long beg = tiktok * ni * 3; 
    Eigen::VectorXd tx1(ni); 
    Eigen::VectorXd tx2(ni); 
    for (int i = 0; i < ni; i++)
    {
	tx1(i) = dadu.v8[beg + i + ni]; 
	tx2(i) = dadu.v8[beg + i + ni + ni]; 
    }
    tx1.array() -= tx1.mean();
    tx2.array() -= tx2.mean();
    Eigen::VectorXd x1 = dadu.evec.transpose() * tx1; 
    Eigen::VectorXd x2 = dadu.evec.transpose() * tx2; 
    //snp; 
    int cols = dadu.ww.cols()+2;

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd hh(ni);
    Eigen::VectorXd h2inv(ni);   //old d2; 
    Eigen::VectorXd h1inv(ni);   //old d1; 
    //h2inv=1/hh; h1inv=sqrt(h2inv); 
    Eigen::MatrixXd X(ni, dadu.ww.cols()+2); 
    Eigen::MatrixXd xtx(cols,cols);
    Eigen::VectorXd xty(ni); 
    Eigen::VectorXd y1(ni); 
    Eigen::VectorXd beta1;    
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    

    Eigen::VectorXd t1(cols);    
    t1.setOnes(); 
    t1 *= 1e-10; 

    Eigen::MatrixXd X2(ni, 2); 
    Eigen::MatrixXd xtx2(2,2);
    Eigen::VectorXd xty2(ni); 
    Eigen::VectorXd beta2(ni);

//    Eigen::VectorXd r1(ni);    
//    Eigen::VectorXd r0(ni);    

    long shift = tiktok*dadu.nph; 
    int max_iter = 100;    
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	Eigen::VectorXd curph = dadu.mph.col(p); 
	double pval=1; 
	int steps = 0; 
    	double eta0=1, eta1 = 0.1; 
//	if(dadu.m_lrt == 1) eta1 = dadu.veta0[p]; 
//	else eta1 = 0.1; 
	double diff = 1; 
//	double like_old = -1e100; 
//	int hybrid = 1; 
	while (fabs(diff)>1e-6 && steps<max_iter) 
	{
	    steps++;  
	    eta0 = eta1; 
	    //y=x beta + e; 
	    v1.setOnes(); 
	    hh = eta0 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x1.cwiseProduct(h1inv), x2.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += t1.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    //here e1 needs to scale back by multiplying 1/d1; 
	    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  

//	    double like_new = 0;
//	    for(int i = 0; i < ni; i++)
//		like_new -= log(hh(i)) / 2.0; 
//	    like_new -= log(r2.sum()) * ni / 2.0; 
//	    if(like_new < like_old) {
//		cout << "likelhood decrease detected, switch to idle IDUL." << endl; 
//		hybrid = 1;   
//	    }
//	    like_old = like_new; 
	    
	    if(dadu.m_lrt == 1) {
		//r2 = (1, D) beta2 + err;                  
		v1.setOnes(); 
		X2.col(0) = v1.cwiseProduct(h2inv); 
		X2.col(1) = dadu.eval.cwiseProduct(h2inv); 
		xtx2 = X2.transpose() * X2; 
		xty2 = X2.transpose() * r2; 
		beta2 = xtx2.ldlt().solve(xty2);
	//	beta2 = xtx2.inverse() * xty2;
//		if(hybrid == 0) 
//        	    eta1 = beta2(1) / max(1e-20, beta2(0)); 
//		else 
		{
		    double t1 = beta2(1) / max(1e-20, beta2(0));  
		    double t2 = beta2(1) / (r2.sum() / ni);  
		    eta1 = (t1 + t2)/2.0; 
		}
	        if(eta1 >= 1e5) eta1 = 1e5; 
	        else if (eta1 <= 1e-5) eta1 = 1e-5;
	    } 
//	    else {
//		double rh1r = r2.sum(); 
//		r0 = r2.cwiseProduct(h2inv); 
//		double rh2r = r0.sum(); 
//		double trh1 = h2inv.sum(); 
//		r1 = h2inv.cwiseProduct(h2inv); 
//		double trh2 = r1.sum(); 
//	
//		double mh = trh1/ni; 
//		double vh = trh2/ni - mh*mh; 
//	
//		double fp1 = (trh1 - ni*rh2r/rh1r); 
//		double delta = fp1 * eta0 / ni / vh;  
//		eta1 = eta0 + delta; 
//		if(eta1 >= 1e5) eta1 = 1e5; 
//		else if (eta1 <= 1e-5) eta1 = 1e-5;
//	    }
	    else {   //if(dadu.m_wald == 1) 
		double rh1r = r2.sum(); 
		double rh2r = r2.cwiseProduct(h2inv).sum(); 
		double trh1 = h2inv.sum(); 
		double trh2 = h2inv.cwiseProduct(h2inv).sum(); 
	
		double mh = trh1/ni; 
		double vh = trh2/ni - mh*mh; 

		Eigen::MatrixXd ww2 = dadu.ww; 
                Eigen::MatrixXd X(ni, cols); 
	    	X << x1.cwiseProduct(h2inv), x2.cwiseProduct(h2inv), ww2.array().colwise() * h2inv.array(); 
		xtx2 = X.transpose() * X; 
		Eigen::MatrixXd mat = xtx.inverse() * xtx2; 
		double fp1 = (trh1 - mat.trace() - (ni-cols)*rh2r/rh1r); 
	
		eta1 = eta0 + fp1 * eta0 /ni / (vh); 
		if(eta1 >= 1e5) eta1 = 1e5; 
		else if (eta1 <= 1e-5) eta1 = 1e-5;
	    }

	    diff = eta1/(1+eta1) - eta0/(1+eta0); 
	}


//	double bvar = 1; 
//	Eigen::MatrixXd cov = xtx.inverse();
//	double tauinv = r2.sum() / (ni-cols);
//	bvar = tauinv * cov(0, 0);
	double loglike1 = 0; 
//	if(dadu.m_lrt == 1) 
	{
	    for (int i = 0; i < ni; i++)
		loglike1 += log(h1inv(i)); 
	    loglike1 -= log(r2.sum()) * ni / 2.0;   

	    double lrts = -2.0 * (dadu.vloglike0[p] - loglike1);   
	    if(lrts < 0 ) {
//		cout << "lrt test statistics = " << lrts << endl; 
		lrts = 0;
	    }
	    pval = cdf(complement(chi_squared(2), lrts));
	}
//	else {  //if(dadu.m_wald == 1) 
//	    double tstats = beta1(0) * beta1(0) / bvar;
//    //	    pval = 2 * cdf(complement(students_t(dof), std::abs(tstats)));
//	    //t-test; 
////	    pval = cdf(complement(chi_squared(1), tstats*tstats));
//	    //wald-test, chi_squared (1); 
//	    boost::math::fisher_f_distribution<double> dist(1, ni-cols); 
//	    pval = boost::math::cdf(complement(dist, tstats)); 
//	    //wald-test, f(1,n-c); 
//	}

//        long shift = tiktok*dadu.nph; 
	dadu.pstats[shift+p].pv = pval; 
	dadu.pstats[shift+p].eta = eta1/(1+eta1); 
	dadu.pstats[shift+p].niter = steps; 
//	dadu.pstats[shift+p].beta = beta1(0); 
//	dadu.pstats[shift+p].sigma = sqrt(bvar);  
	if (dadu.m_lrt == 1) {
	    dadu.pstats[shift+p].l1 = loglike1; 
	    dadu.pstats[shift+p].l0 = dadu.vloglike0[p]; 
	}
    }
}

void jacquard_pv1(void *par)    //parent of origin conditional test;
{
    int ni = dadu.ni;    
    long tiktok = (long) par; 
    long beg = tiktok * ni * 3; 
    Eigen::VectorXd tx1(ni); 
    Eigen::VectorXd tx2(ni); 
    for (int i = 0; i < ni; i++)
    {
	tx1(i) = dadu.v8[beg + i + ni]; 
	tx2(i) = dadu.v8[beg + i + ni + ni]; 
    }
    tx1.array() -= tx1.mean();
    tx2.array() -= tx2.mean();
    Eigen::VectorXd x1 = dadu.evec.transpose() * tx1; 
    Eigen::VectorXd x2 = dadu.evec.transpose() * tx2; 
    //snp; 
    if(dadu.ctrl == 0)  //compute bf( maternal | paternal)
    {
	double temp; 
	for(int i = 0; i < ni; i++)
	{
	    temp = x1(i); 
	    x1(i) = x2(i); 
	    x2(i) = temp; 
	}
    }

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd hh(ni);
    Eigen::VectorXd h2inv(ni);   //old d2; 
    Eigen::VectorXd h1inv(ni);   //old d1; 

    Eigen::VectorXd xty(ni); 
    Eigen::VectorXd y1(ni); 
    Eigen::VectorXd beta1;    
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    

    Eigen::MatrixXd X2(ni, 2); 
    Eigen::MatrixXd xtx2(2,2);
    Eigen::VectorXd xty2(ni); 
    Eigen::VectorXd beta2(ni);

      
    long shift = tiktok*dadu.nph; 
    int max_iter = 100;    
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	int steps = 0; 
	double eta0=1, eta1=0.1; 
//	if(dadu.m_lrt == 1) eta1 = dadu.veta0[p]; 
//	else eta1 = 0.1; 
	double loglike1 = 0; 
	{
	    int cols = dadu.ww.cols()+2;
	    //h2inv=1/hh; h1inv=sqrt(h2inv); 
	    Eigen::MatrixXd X(ni, cols); 
	    Eigen::MatrixXd xtx(cols,cols);

	    Eigen::VectorXd t1(cols);    
	    t1.setOnes(); 
	    t1 *= 1e-10; 

	    Eigen::VectorXd curph = dadu.mph.col(p); 
	    double diff = 1; 

	    while (fabs(diff)>1e-6 && steps<max_iter) 
	    {
		steps++;  
		eta0 = eta1; 
		//y=x beta + e; 
		v1.setOnes(); 
		hh = eta0 * dadu.eval + v1; 
		h2inv = hh.cwiseInverse(); 
		h1inv = h2inv.cwiseSqrt(); 
		Eigen::MatrixXd ww2 = dadu.ww; 
		X << x1.cwiseProduct(h1inv), x2.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
		xtx = X.transpose() * X; 
		xtx += t1.matrix().asDiagonal(); 
		y1 = curph.cwiseProduct(h1inv); 
		xty = X.transpose() * y1; 
		beta1 = xtx.ldlt().solve(xty);
	//	beta1 = xtx.inverse() * xty;
		e1 = (y1 - X * beta1); 
		r2 = e1.cwiseProduct(e1); 
		//here e1 needs to scale back by multiplying 1/d1; 
		//so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  
		
    //	    if(dadu.m_lrt == 1) 
		{
		    //r2 = (1, D) beta2 + err;                  
		    v1.setOnes(); 
		    X2.col(0) = v1.cwiseProduct(h2inv); 
		    X2.col(1) = dadu.eval.cwiseProduct(h2inv); 
		    xtx2 = X2.transpose() * X2; 
		    xty2 = X2.transpose() * r2; 
		    beta2 = xtx2.ldlt().solve(xty2);
		    {
			double t1 = beta2(1) / max(1e-20, beta2(0));  
			double t2 = beta2(1) / (r2.sum() / ni);  
			eta1 = (t1 + t2)/2.0; 
		    }
		    if(eta1 >= 1e5) eta1 = 1e5; 
		    else if (eta1 <= 1e-5) eta1 = 1e-5;
		} 
		diff = eta1/(1+eta1) - eta0/(1+eta0); 
	    }
	    loglike1 = 0; 
    //	    if(dadu.m_lrt == 1) 
	    {
		for (int i = 0; i < ni; i++)
		    loglike1 += log(h1inv(i)); 
		loglike1 -= log(r2.sum()) * ni / 2.0;   
	    }
	}

        //below is for the null; 
	double loglike0 = 0; 
	{
	    int cols = dadu.ww.cols()+1;
	    //h2inv=1/hh; h1inv=sqrt(h2inv); 
	    Eigen::MatrixXd X(ni, cols); 
	    Eigen::MatrixXd xtx(cols,cols);

	    Eigen::VectorXd t1(cols);    
	    t1.setOnes(); 
	    t1 *= 1e-10; 

	    Eigen::VectorXd curph = dadu.mph.col(p); 
	    int steps = 0; 
	    double eta0=1, eta1=0.1; 
	    double diff = 1; 

	    while (fabs(diff)>1e-6 && steps<max_iter) 
	    {
		steps++;  
		eta0 = eta1; 
		//y=x beta + e; 
		v1.setOnes(); 
		hh = eta0 * dadu.eval + v1; 
		h2inv = hh.cwiseInverse(); 
		h1inv = h2inv.cwiseSqrt(); 
		Eigen::MatrixXd ww2 = dadu.ww; 
		X << x2.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
		xtx = X.transpose() * X; 
		xtx += t1.matrix().asDiagonal(); 
		y1 = curph.cwiseProduct(h1inv); 
		xty = X.transpose() * y1; 
		beta1 = xtx.ldlt().solve(xty);
	//	beta1 = xtx.inverse() * xty;
		e1 = (y1 - X * beta1); 
		r2 = e1.cwiseProduct(e1); 
		//here e1 needs to scale back by multiplying 1/d1; 
		//so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  
		
    //	    if(dadu.m_lrt == 1) 
		{
		    //r2 = (1, D) beta2 + err;                  
		    v1.setOnes(); 
		    X2.col(0) = v1.cwiseProduct(h2inv); 
		    X2.col(1) = dadu.eval.cwiseProduct(h2inv); 
		    xtx2 = X2.transpose() * X2; 
		    xty2 = X2.transpose() * r2; 
		    beta2 = xtx2.ldlt().solve(xty2);
		    {
			double t1 = beta2(1) / max(1e-20, beta2(0));  
			double t2 = beta2(1) / (r2.sum() / ni);  
			eta1 = (t1 + t2)/2.0; 
		    }
		    if(eta1 >= 1e5) eta1 = 1e5; 
		    else if (eta1 <= 1e-5) eta1 = 1e-5;
		} 
		diff = eta1/(1+eta1) - eta0/(1+eta0); 
	    }
	    loglike0 = 0; 
	    {
		for (int i = 0; i < ni; i++)
		    loglike0 += log(h1inv(i)); 
		loglike0 -= log(r2.sum()) * ni / 2.0;   

	    }
	}
	double lrts = -2.0 * (loglike0 - loglike1);   
	if(lrts < 0 ) {
//		cout << "lrt test statistics = " << lrts << endl; 
	    lrts = 0;
	}
	double pval = cdf(complement(chi_squared(1), lrts));

//	long shift = tiktok*dadu.nph; 
	dadu.pstats[shift+p].pv = pval; 
	dadu.pstats[shift+p].eta = eta1/(1+eta1); 
	dadu.pstats[shift+p].niter = steps; 
//	dadu.pstats[shift+p].beta = beta1(0); 
//	dadu.pstats[shift+p].sigma = sqrt(bvar);  
	if (dadu.m_lrt == 1) {
	    dadu.pstats[shift+p].l1 = loglike1; 
	    dadu.pstats[shift+p].l0 = loglike0; 
	}
    }
}

void jacquard_pv(void *par) 
{
//    boost::random::mt19937 rng;
//    boost::random::uniform_real_distribution<> dist(0, 1.0);

    int ni = dadu.ni;    
    int flag = dadu.flag_untransmitted; 
    long tiktok = (long) par; 
    long beg = tiktok * ni * (1+flag); 
    Eigen::VectorXd tx1(ni); 
    for (int i = 0; i < ni; i++)
	tx1(i) = dadu.v8[beg + i]; 
    tx1.array() -= tx1.mean();
    Eigen::VectorXd x1 = dadu.evec.transpose() * tx1; 
    //snp; 
    int cols = dadu.ww.cols()+1;

//    if(tiktok == 1) {
//       for (int i = 0; i < ni; i++)
//	   cout << x1(i) << " "; 
//       cout << endl; 
//    }


    Eigen::VectorXd v1(ni);
    Eigen::VectorXd hh(ni);
    Eigen::VectorXd h2inv(ni);   //old d2; 
    Eigen::VectorXd h1inv(ni);   //old d1; 
    //h2inv=1/hh; h1inv=sqrt(h2inv); 
    Eigen::MatrixXd X(ni, dadu.ww.cols()+1); 
    Eigen::MatrixXd xtx(cols,cols);
    Eigen::VectorXd xty(ni); 
    Eigen::VectorXd y1(ni); 
    Eigen::VectorXd beta1;    
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    

    Eigen::VectorXd t1(dadu.ww.cols()+1);    
    t1.setOnes(); 
    t1 *= 1e-10; 

    Eigen::MatrixXd X2(ni, 2); 
    Eigen::MatrixXd xtx2(2,2);
    Eigen::VectorXd xty2(ni); 
    Eigen::VectorXd beta2(ni);

    Eigen::MatrixXd X3(ni, dadu.ww.cols()+1); 
//    Eigen::VectorXd r1(ni);    
//    Eigen::VectorXd r0(ni);    

    long shift = tiktok*dadu.nph; 
    int max_iter = 100;    
    //this uses eta1 under the null; 
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	Eigen::VectorXd curph = dadu.mph.col(p); 
	double pval=1; 
	int steps = 0; 
    	double eta0=1, eta1; 
	if(dadu.m_lrt == 1) eta1 = dadu.veta0[p]; 
	else eta1 = 0.1; 
	double diff = 1; 
	while (fabs(diff)>1e-6 && steps<max_iter) 
	{
	    steps++;  
	    eta0 = eta1; 
	    //y=x beta + e; 
	    v1.setOnes(); 
	    hh = eta0 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x1.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += t1.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    //here e1 needs to scale back by multiplying 1/d1; 
	    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  
	    double tauinv = r2.mean();  

	    if(dadu.m_lrt == 1) {
		//r2 = (1, D) beta2 + err;                  
		v1.setOnes(); 
		X2.col(0) = v1.cwiseProduct(h2inv); 
		X2.col(1) = dadu.eval.cwiseProduct(h2inv); 
		xtx2 = X2.transpose() * X2; 
		xty2 = X2.transpose() * r2; 
		beta2 = xtx2.ldlt().solve(xty2);
	//	beta2 = xtx2.inverse() * xty2;
		double mu = beta2(0); 
		double ga = beta2(1); 
//		eta1 = ga/tauinv+(1-mu/tauinv)*eta0; 
		    double t1 = ga/tauinv+(1-mu/tauinv)*eta0; 
//			double t1 = beta2(1) / max(1e-20, beta2(0));  
		    double t2 = ga / tauinv;  
		    eta1 = (t1 + t2)/2.0; 
	        if(eta1 >= 1e5) eta1 = 1e5; 
	        else if (eta1 <= 0) eta1 = 0;
	    } 
	    else {   //if(dadu.m_wald == 1) 
		double rh1r = r2.sum(); 
		double rh2r = r2.cwiseProduct(h2inv).sum(); 
		double trh1 = h2inv.sum(); 
		double trh2 = h2inv.cwiseProduct(h2inv).sum(); 
	
		double mh = trh1/ni; 
		double vh = trh2/ni - mh*mh; 

		Eigen::MatrixXd ww2 = dadu.ww; 
	    	X3 << x1.cwiseProduct(h2inv), ww2.array().colwise() * h2inv.array(); 
		xtx2 = X3.transpose() * X3; 
		Eigen::MatrixXd mat = xtx.inverse() * xtx2; 
		double fp1 = (trh1 - mat.trace() - (ni-cols)*rh2r/rh1r); 
	        double delta =  fp1 * eta0 /ni / vh;  
		double t1 = eta0 + delta; 
		double t2 = eta0 + delta * mh; 
		eta1 = (t1+t2)/2.0; 
//		eta1 = eta0 + fp1 * eta0 /ni / vh * (0.5+ 0.5*mh); 
		if(eta1 >= 1e5) eta1 = 1e5; 
		else if (eta1 <= 1e-5) eta1 = 1e-5;
	    }

	    diff = eta1/(1+eta1) - eta0/(1+eta0); 
//	    if(tiktok == 0) printf("%.8f %.8f \n", eta0, eta1); 
	}

	double bvar = 1; 
	Eigen::MatrixXd cov = xtx.inverse();
	double tauinv = r2.sum() / (ni-cols);
	bvar = tauinv * cov(0, 0);
	double loglike1 = 0; 
	if(dadu.m_lrt == 1) 
	{
	    for (int i = 0; i < ni; i++)
		loglike1 += log(h1inv(i)); 
	    loglike1 -= log(r2.sum()) * ni / 2.0;   

	    double lrts = -2.0 * (dadu.vloglike0[p] - loglike1);   
	    if(lrts < 0 ) {
//		cout << "lrt test statistics = " << lrts << endl; 
		lrts = 0;
	    }
	    pval = cdf(complement(chi_squared(1), lrts));
	}
	else {  //if(dadu.m_wald == 1)
	    double tstats = beta1(0) * beta1(0) / bvar;
    //	    pval = 2 * cdf(complement(students_t(dof), std::abs(tstats)));
	    //t-test; 
//	    pval = cdf(complement(chi_squared(1), tstats*tstats));
	    //wald-test, chi_squared (1); 
	    boost::math::fisher_f_distribution<double> dist(1, ni-cols); 
	    pval = boost::math::cdf(complement(dist, tstats)); 
	    //wald-test, f(1,n-c); 
	}

//	long shift = tiktok*dadu.nph; 
	dadu.pstats[shift+p].pv = pval; 
	dadu.pstats[shift+p].eta = eta1/(1+eta1); 
	dadu.pstats[shift+p].niter = steps; 
	dadu.pstats[shift+p].beta = beta1(0); 
	dadu.pstats[shift+p].sigma = sqrt(bvar);  
	if (dadu.m_lrt == 1) {
	    dadu.pstats[shift+p].l1 = loglike1; 
	    dadu.pstats[shift+p].l0 = dadu.vloglike0[p]; 
	}
    }
}

void jacquard_pv_plus(void *par)   //revised IDUL;
{
    int ni = dadu.ni;    
    int flag = dadu.flag_untransmitted; 
    long tiktok = (long) par; 
    long beg = tiktok * ni * (1+flag); 
    Eigen::VectorXd tx1(ni); 
    for (int i = 0; i < ni; i++)
	tx1(i) = dadu.v8[beg + i]; 
    tx1.array() -= tx1.mean();
    Eigen::VectorXd x1 = dadu.evec.transpose() * tx1; 
    //snp; 
    int cols = dadu.ww.cols()+1;

//    if(tiktok == 1) {
//       for (int i = 0; i < ni; i++)
//	   cout << x1(i) << " "; 
//       cout << endl; 
//    }

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd hh(ni);
    Eigen::VectorXd h2inv(ni);   //old d2; 
    Eigen::VectorXd h1inv(ni);   //old d1; 
    //h2inv=1/hh; h1inv=sqrt(h2inv); 
    Eigen::MatrixXd X(ni, dadu.ww.cols()+1); 
    Eigen::MatrixXd xtx(cols,cols);
    Eigen::VectorXd xty(ni); 
    Eigen::VectorXd y1(ni); 
    Eigen::VectorXd beta1;    
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    
    Eigen::VectorXd s2(ni);    

    Eigen::VectorXd t1(dadu.ww.cols()+1);    
    t1.setOnes(); 
    t1 *= 1e-10; 

    Eigen::MatrixXd X2(ni, 2); 
    Eigen::MatrixXd xtx2(2,2);
    Eigen::VectorXd xty2(ni); 
    Eigen::VectorXd beta2(ni);

    Eigen::MatrixXd X3(ni, dadu.ww.cols()+1); 
//    Eigen::VectorXd r1(ni);    
//    Eigen::VectorXd r0(ni);    

    long shift = tiktok*dadu.nph; 
//    int max_iter = 100;    
//    boost::random::mt19937 rng;
//    boost::random::uniform_real_distribution<> dist(0.1, 0.9);
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	Eigen::VectorXd curph = dadu.mph.col(p); 
	double pval=1; 
	int steps = 0; 
    	double eta0=1, eta1; 
	if(dadu.m_lrt == 1) eta1 = dadu.veta0[p]; 
	else eta1 = 0.1; 

	eta0 = eta1; 
	//y=x beta + e; 
	v1.setOnes(); 
	hh = eta0 * dadu.eval + v1; 
	h2inv = hh.cwiseInverse(); 
	h1inv = h2inv.cwiseSqrt(); 
	Eigen::MatrixXd ww2 = dadu.ww; 
	X << x1.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	xtx = X.transpose() * X; 
	xtx += t1.matrix().asDiagonal(); 
	y1 = curph.cwiseProduct(h1inv); 
	xty = X.transpose() * y1; 
	beta1 = xtx.ldlt().solve(xty);
//	beta1 = xtx.inverse() * xty;
	e1 = (y1 - X * beta1); 
	r2 = e1.cwiseProduct(e1); 
	//here e1 needs to scale back by multiplying 1/d1; 
	//so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  
	double tauinv = r2.mean(); 

	double tauinv1, like_eta0, like_eta1; 
	
	if(dadu.m_lrt == 1) {
	    like_eta0 = -log(tauinv) * ni;
	    for(int i = 0; i < ni; i++)
		like_eta0 -= log(hh(i)); 
	    //likelihood; 
	    while(1) {
	        //r2 = (1, D) beta2 + err;                  
		v1.setOnes(); 
		X2.col(0) = v1.cwiseProduct(h2inv); 
		X2.col(1) = dadu.eval.cwiseProduct(h2inv); 
		xtx2 = X2.transpose() * X2; 
		xty2 = X2.transpose() * r2; 
		beta2 = xtx2.ldlt().solve(xty2);

		double mu = beta2(0); 
		double ga = beta2(1); 
		eta1 = ga/tauinv+(1-mu/tauinv)*eta0/2.0; 
		//lazy update; 
		if(eta1 >= 1e5) eta1 = 1e5; 
		else if (eta1 <= 1e-5) eta1 = 1e-5;
		steps++;

		while(1) {
		    if(abs(eta1/(1+eta1) - eta0/(1+eta0)) < 1e-6) break;

		    hh = eta1 * dadu.eval + v1; 
		    h2inv = hh.cwiseInverse(); 
		    h1inv = h2inv.cwiseSqrt(); 
		    Eigen::MatrixXd ww2 = dadu.ww; 
		    X << x1.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
		    xtx = X.transpose() * X; 
		    xtx += t1.matrix().asDiagonal(); 
		    y1 = curph.cwiseProduct(h1inv); 
		    xty = X.transpose() * y1; 
		    beta1 = xtx.ldlt().solve(xty);
	    //	beta1 = xtx.inverse() * xty;
		    e1 = (y1 - X * beta1); 
		    r2 = e1.cwiseProduct(e1); 
		    //here e1 needs to scale back by multiplying 1/d1; 
		    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  
		    tauinv1 = r2.mean(); 
		    like_eta1 = -log(tauinv1) * ni;
		    for(int i = 0; i < ni; i++)
			like_eta1 -= log(hh(i)); 
		    if(like_eta1 >= like_eta0) break; 
		    eta1 = (eta1 + eta0)/2.0; 
		} //inner while loop; 
//	        if(tiktok == 0) cout << steps << "x\t" << eta0 << "\t" << eta1 << "\t" << abs(eta0-eta1) << endl; 
		if(abs(eta1/(1+eta1) - eta0/(1+eta0)) < 1e-6) break;
		eta0 = eta1; 
		tauinv = tauinv1; 
		like_eta0 = like_eta1;
	    } //outer while loop; 

	} 
	else {   //if(dadu.m_wald == 1) 
	    like_eta0 = -log(tauinv) * (ni-cols);
	    for(int i = 0; i < ni; i++)
		like_eta0 -= log(hh(i));        
	    like_eta0 -= log(xtx.determinant()); 
	    //reml-likelihood; 
	    while(1) {
		double rh1r = r2.sum(); 
		double rh2r = r2.cwiseProduct(h2inv).sum(); 
		double trh1 = h2inv.sum(); 
		double trh2 = h2inv.cwiseProduct(h2inv).sum(); 
	
		double mh = trh1/ni; 
		double vh = trh2/ni - mh*mh; 

		Eigen::MatrixXd ww2 = dadu.ww; 
		X3 << x1.cwiseProduct(h2inv), ww2.array().colwise() * h2inv.array(); 
		xtx2 = X3.transpose() * X3; 
		Eigen::MatrixXd mat = xtx.inverse() * xtx2; 
		double fp1 = (trh1 - mat.trace() - (ni-cols)*rh2r/rh1r); 
	
		eta1 = eta0 + fp1 * eta0 /ni / vh * (0.5+ 0.5*mh); 
//		eta1 = eta0 + fp1 * eta0 /ni / (vh); 
		if(eta1 >= 1e5) eta1 = 1e5; 
		else if (eta1 <= 1e-5) eta1 = 1e-5;
		steps++;

		while(1) {
		    if(abs(eta1/(1+eta1) - eta0/(1+eta0)) < 1e-6) break;

		    hh = eta1 * dadu.eval + v1; 
		    h2inv = hh.cwiseInverse(); 
		    h1inv = h2inv.cwiseSqrt(); 
		    Eigen::MatrixXd ww2 = dadu.ww; 
		    X << x1.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
		    xtx = X.transpose() * X; 
		    xtx += t1.matrix().asDiagonal(); 
		    y1 = curph.cwiseProduct(h1inv); 
		    xty = X.transpose() * y1; 
		    beta1 = xtx.ldlt().solve(xty);
	    //	beta1 = xtx.inverse() * xty;
		    e1 = (y1 - X * beta1); 
		    r2 = e1.cwiseProduct(e1); 
		    //here e1 needs to scale back by multiplying 1/d1; 
		    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  
		    tauinv1 = r2.mean(); 
		    like_eta1 = -log(tauinv1) * (ni-cols);
		    for(int i = 0; i < ni; i++)
			like_eta1 -= log(hh(i)); 
		    like_eta1 -= log(xtx.determinant()); 
		    if(like_eta1 >= like_eta0) break; 
		    eta1 = (eta1 + eta0)/2.0; 
		} //inner while loop; 
//	        if(tiktok == 0) cout << steps << "x\t" << eta0 << "\t" << eta1 << "\t" << abs(eta0-eta1) << endl; 
		if(abs(eta1/(1+eta1) - eta0/(1+eta0)) < 1e-6) break;
		eta0 = eta1; 
		tauinv = tauinv1; 
		like_eta0 = like_eta1;
	    }
	}

	Eigen::MatrixXd cov = xtx.inverse();
	double bvar = r2.sum() / (ni-cols) * cov(0, 0);
	double loglike1 = 0; 
	if(dadu.m_lrt == 1) 
	{
	    for (int i = 0; i < ni; i++)
		loglike1 += log(h1inv(i)); 
	    loglike1 -= log(r2.sum()) * ni / 2.0;   

	    double lrts = -2.0 * (dadu.vloglike0[p] - loglike1);   
	    if(lrts < 0 ) {
//		cout << "lrt test statistics = " << lrts << endl; 
		lrts = 0;
	    }
	    pval = cdf(complement(chi_squared(1), lrts));
	}
	else {  //if(dadu.m_wald == 1)
	    double tstats = beta1(0) * beta1(0) / bvar;
    //	    pval = 2 * cdf(complement(students_t(dof), std::abs(tstats)));
	    //t-test; 
//	    pval = cdf(complement(chi_squared(1), tstats*tstats));
	    //wald-test, chi_squared (1); 
	    boost::math::fisher_f_distribution<double> dist(1, ni-cols); 
	    pval = boost::math::cdf(complement(dist, tstats)); 
	    //wald-test, f(1,n-c); 
	}

//	long shift = tiktok*dadu.nph; 
	dadu.pstats[shift+p].pv = pval; 
	dadu.pstats[shift+p].eta = eta1/(1+eta1); 
	dadu.pstats[shift+p].niter = steps; 
	dadu.pstats[shift+p].beta = beta1(0); 
	dadu.pstats[shift+p].sigma = sqrt(bvar);  
	if (dadu.m_lrt == 1) {
	    dadu.pstats[shift+p].l1 = loglike1; 
	    dadu.pstats[shift+p].l0 = dadu.vloglike0[p]; 
	}
    }
}

void jacquard_bf(void *par)   //bayes factors
{
    int ni = dadu.ni;    
    int flag = dadu.flag_untransmitted; 
    long tiktok = (long) par; 
    long beg = tiktok * ni * (1+flag); 
    Eigen::VectorXd tx1(ni); 
    for (int i = 0; i < ni; i++)
	tx1(i) = dadu.v8[beg + i]; 
    tx1.array() -= tx1.mean();
    Eigen::VectorXd x1 = dadu.evec.transpose() * tx1; 
    //snp; 
    int cols = dadu.ww.cols()+1;

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd hh(ni);
    Eigen::VectorXd h2inv(ni);   //old d2; 
    Eigen::VectorXd h1inv(ni);   //old d1; 
    //h2inv=1/hh; h1inv=sqrt(h2inv); 
    Eigen::MatrixXd X(ni, dadu.ww.cols()+1); 
    Eigen::MatrixXd xtx(cols,cols);
    Eigen::VectorXd xty(ni); 
    Eigen::VectorXd y1(ni); 
    Eigen::VectorXd beta1;    
    Eigen::VectorXd beta0;    
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    

    Eigen::VectorXd prior1(dadu.ww.cols()+1);    
    prior1.setOnes(); 
    prior1 *= 1e-10; 
//    prior1(0) = 4.0; // 1/sigma^2 at simga=0.5; 

//    if(tiktok == 1) {
//	for (int i = 0; i < cols; i++)
//	    cout << prior1(i) << " "; 
//	cout << endl; 
//    }

    Eigen::MatrixXd X2(ni, 2); 
    Eigen::MatrixXd xtx2(2,2);
    Eigen::VectorXd xty2(ni); 
    Eigen::VectorXd beta2(ni);

    Eigen::MatrixXd X3(ni, dadu.ww.cols()+1); 
//    Eigen::VectorXd r1(ni);    
//    Eigen::VectorXd r0(ni);    

    long shift = tiktok*dadu.nph; 
    int max_iter = 100;    
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	Eigen::VectorXd curph = dadu.mph.col(p); 
	int steps = 0; 
    	double eta0=1, eta1; 
	eta1 = dadu.veta0[p]; 
	double diff = 1; 
//	double like_old = -1e100; 
//	int hybrid = 1; 
	while (fabs(diff)>1e-6 && steps<max_iter) 
	{
	    steps++;  
	    eta0 = eta1; 
	    //y=x beta + e; 
	    v1.setOnes(); 
	    hh = eta0 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x1.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior1.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    //here e1 needs to scale back by multiplying 1/d1; 
	    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  

//	    double like_new = 0;
//	    for(int i = 0; i < ni; i++)
//		like_new -= log(hh(i)) / 2.0; 
//	    like_new -= log(r2.sum()) * ni / 2.0; 
//	    if(like_new < like_old) {
//		cout << "likelhood decrease detected, switch to idle IDUL." << endl; 
//		hybrid = 1;   
//	    }
//	    like_old = like_new; 
	    
	    {
		//r2 = (1, D) beta2 + err;                  
		v1.setOnes(); 
		X2.col(0) = v1.cwiseProduct(h2inv); 
		X2.col(1) = dadu.eval.cwiseProduct(h2inv); 
		xtx2 = X2.transpose() * X2; 
		xty2 = X2.transpose() * r2; 
		beta2 = xtx2.ldlt().solve(xty2);
	//	beta2 = xtx2.inverse() * xty2;
//		if(hybrid == 0) 
//        	    eta1 = beta2(1) / max(1e-20, beta2(0)); 
//		else 
		{
		    double t1 = beta2(1) / max(1e-20, beta2(0));  
		    double t2 = beta2(1) / (r2.sum() / ni);  
		    eta1 = (t1 + t2)/2.0; 
		}
	        if(eta1 >= 1e5) eta1 = 1e5; 
	        else if (eta1 <= 1e-5) eta1 = 1e-5;
	    } 

	    diff = eta1/(1+eta1) - eta0/(1+eta0); 
	}
	double loglike1 = -log(0.5); //sigam=0.5, \beta \sim N(0,sigma^2);
	{
            prior1(0) = 4.0; // 1/sigma^2 at simga=0.5; 
	    v1.setOnes(); 
	    hh = eta1 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x1.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior1.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    loglike1 -= 0.5 * log(xtx.determinant()); 
	    loglike1 -= log(r2.sum()) * ni / 2.0;   
	}
	double bvar = 1; 
	Eigen::MatrixXd cov = xtx.inverse();
	double tauinv = r2.sum() / (ni-cols);
	bvar = tauinv * cov(0, 0);

	double loglikenull = 0; 
	{
	    int cols = dadu.ww.cols();
	    Eigen::MatrixXd X(ni, cols); 
	    Eigen::MatrixXd xtx(cols,cols);
	    Eigen::VectorXd prior0(cols);    
	    prior0.setOnes(); 
	    prior0 *= 1e-10; 
	    v1.setOnes(); 
	    hh = eta1 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior0.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta0 = xtx.ldlt().solve(xty);
	    e1 = (y1 - X * beta0); 
	    r2 = e1.cwiseProduct(e1); 
	    loglikenull -= 0.5 * log(xtx.determinant()); 
	    loglikenull -= log(r2.sum()) * ni /2.0; 
	}


	double logbf = loglike1 - loglikenull; 

//	long shift = tiktok*dadu.nph; 
	dadu.pstats[shift+p].logbf = logbf; 
	dadu.pstats[shift+p].eta = eta1/(1+eta1); 
	dadu.pstats[shift+p].niter = steps; 
	dadu.pstats[shift+p].beta = beta1(0); 
	dadu.pstats[shift+p].sigma = sqrt(bvar);  
//	if (dadu.m_lrt == 1) {
//	    dadu.pstats[shift+p].l1 = loglike1; 
//	    dadu.pstats[shift+p].l0 = dadu.vloglike0[p]; 
//	}
    }
}

void jacquard_bf2(void *par)   //joint bayes factors
{
    int ni = dadu.ni;    
    long tiktok = (long) par; 
    long beg = tiktok * ni * 3; 
    Eigen::VectorXd tx1(ni); 
    for (int i = 0; i < ni; i++)
	tx1(i) = dadu.v8[beg + ni+ i]; 
    tx1.array() -= tx1.mean();
    Eigen::VectorXd x1 = dadu.evec.transpose() * tx1; 
    Eigen::VectorXd tx2(ni); 
    for (int i = 0; i < ni; i++)
	tx2(i) = dadu.v8[beg + ni + ni+ i]; 
    tx2.array() -= tx2.mean();
    Eigen::VectorXd x2 = dadu.evec.transpose() * tx2; 
    //snp; 
    int cols = dadu.ww.cols()+2;

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd hh(ni);
    Eigen::VectorXd h2inv(ni);   //old d2; 
    Eigen::VectorXd h1inv(ni);   //old d1; 
    //h2inv=1/hh; h1inv=sqrt(h2inv); 
    Eigen::MatrixXd X(ni, dadu.ww.cols()+2); 
    Eigen::MatrixXd xtx(cols,cols);
    Eigen::VectorXd xty(ni); 
    Eigen::VectorXd y1(ni); 
    Eigen::VectorXd beta1;    
    Eigen::VectorXd beta0;    
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    

    Eigen::VectorXd prior1(dadu.ww.cols()+2);    
    prior1.setOnes(); 
    prior1 *= 1e-20; 
//    prior1(0) = 4.0; // 1/sigma^2 at simga=0.5; 

//    if(tiktok == 1) {
//	for (int i = 0; i < cols; i++)
//	    cout << prior1(i) << " "; 
//	cout << endl; 
//    }

    Eigen::MatrixXd X2(ni, 2); 
    Eigen::MatrixXd xtx2(2,2);
    Eigen::VectorXd xty2(ni); 
    Eigen::VectorXd beta2(ni);

    Eigen::MatrixXd X3(ni, dadu.ww.cols()+2); 
//    Eigen::VectorXd r1(ni);    
//    Eigen::VectorXd r0(ni);    

    long shift = tiktok*dadu.nph; 
    int max_iter = 100;    
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	Eigen::VectorXd curph = dadu.mph.col(p); 
	int steps = 0; 
    	double eta0=1, eta1=0.1; 
	double diff = 1; 
	while (fabs(diff)>1e-6 && steps<max_iter) 
	{
	    steps++;  
	    eta0 = eta1; 
	    //y=x beta + e; 
	    v1.setOnes(); 
	    hh = eta0 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x1.cwiseProduct(h1inv), x2.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior1.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    //here e1 needs to scale back by multiplying 1/d1; 
	    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  

	    {
		//r2 = (1, D) beta2 + err;                  
		v1.setOnes(); 
		X2.col(0) = v1.cwiseProduct(h2inv); 
		X2.col(1) = dadu.eval.cwiseProduct(h2inv); 
		xtx2 = X2.transpose() * X2; 
		xty2 = X2.transpose() * r2; 
		beta2 = xtx2.ldlt().solve(xty2);
	//	beta2 = xtx2.inverse() * xty2;
//		if(hybrid == 0) 
//        	    eta1 = beta2(1) / max(1e-20, beta2(0)); 
//		else 
		{
		    double t1 = beta2(1) / max(1e-20, beta2(0));  
		    double t2 = beta2(1) / (r2.sum() / ni);  
		    eta1 = (t1 + t2)/2.0; 
		}
	        if(eta1 >= 1e5) eta1 = 1e5; 
	        else if (eta1 <= 1e-5) eta1 = 1e-5;
	    } 

	    diff = eta1/(1+eta1) - eta0/(1+eta0); 
	}

	double loglike1 = -2*log(0.5); //sigam=0.5, \beta \sim N(0,sigma^2), times 2 because of joint test;
	{
            prior1(0) = 4.0; // 1/sigma^2 at simga=0.5; 
            prior1(1) = 4.0; // 1/sigma^2 at simga=0.5; 
	    v1.setOnes(); 
	    hh = eta1 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x1.cwiseProduct(h1inv),x2.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior1.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
//    	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    loglike1 -= 0.5 * log(xtx.determinant()); 
	    loglike1 -= log(r2.sum()) * ni / 2.0;   
	}

	Eigen::MatrixXd cov = xtx.inverse();
	double tauinv = r2.sum() / (ni-cols);
	double bvar = tauinv * (cov(0, 0) + cov(1,1) - 2.0 * cov(0,1));
        double db = beta1(0)-beta1(1); 
	double chisq1 = db * db / bvar; 
//	double pval = cdf(complement(chi_squared(1), chisq1));
        //this p-value is for H0: beta1(0) == beta1(1); H1: <> 

	double sig0 = sqrt(tauinv * cov(0,0)); 
	double sig1 = sqrt(tauinv * cov(1,1)); 

        double loglikenull = 0; 
	{
	    int cols = dadu.ww.cols();
	    Eigen::MatrixXd X(ni, cols); 
	    Eigen::MatrixXd xtx(cols,cols);
	    Eigen::VectorXd prior0(cols);    
	    prior0.setOnes(); 
	    prior0 *= 1e-20; 

	    v1.setOnes(); 
	    hh = eta1 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior0.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta0 = xtx.ldlt().solve(xty);
	    e1 = (y1 - X * beta0); 
	    r2 = e1.cwiseProduct(e1); 
	    loglikenull -= 0.5 * log(xtx.determinant()); 
	    loglikenull -= log(r2.sum()) * ni /2.0; 
	}

	double logbf = loglike1 - loglikenull; 

	dadu.pstats[shift+p].logbf = logbf; 
	dadu.pstats[shift+p].eta = eta1/(1+eta1); 
	dadu.pstats[shift+p].niter = steps; 
	dadu.pstats[shift+p].beta = beta1(0)/sig0; 
	dadu.pstats[shift+p].sigma = beta1(1)/sig1;  
	dadu.pstats[shift+p].pv = chisq1; //pval;   
//	if (dadu.m_lrt == 1) {
//	    dadu.pstats[shift+p].l1 = loglike1; 
//	    dadu.pstats[shift+p].l0 = dadu.vloglike0[p]; 
//	}
    }
}

void jacquard_bf1(void *par)   //conditioning bayes factors
{
    int ni = dadu.ni;    
    long tiktok = (long) par; 
    long beg = tiktok * ni * 3; 
    Eigen::VectorXd tx1(ni); 
    for (int i = 0; i < ni; i++)
	tx1(i) = dadu.v8[beg + ni+ i]; 
    tx1.array() -= tx1.mean();
    Eigen::VectorXd x1 = dadu.evec.transpose() * tx1; 
    Eigen::VectorXd tx2(ni); 
    for (int i = 0; i < ni; i++)
	tx2(i) = dadu.v8[beg + ni + ni+ i]; 
    tx2.array() -= tx2.mean();
    Eigen::VectorXd x2 = dadu.evec.transpose() * tx2; 

    //by default x1 is paternal and x2 is maternal. 
    if(dadu.ctrl == 0)  //compute bf( maternal | paternal)
    {
	double temp; 
	for(int i = 0; i < ni; i++)
	{
	    temp = x1(i); 
	    x1(i) = x2(i); 
	    x2(i) = temp; 
	}
    }

    Eigen::VectorXd v1(ni);
    Eigen::VectorXd hh(ni);
    Eigen::VectorXd h2inv(ni);   //old d2; 
    Eigen::VectorXd h1inv(ni);   //old d1; 
    //h2inv=1/hh; h1inv=sqrt(h2inv); 
    Eigen::VectorXd xty(ni); 
    Eigen::VectorXd y1(ni); 
    Eigen::VectorXd beta1;    
    Eigen::VectorXd beta0;    
    Eigen::VectorXd e1(ni);    
    Eigen::VectorXd r2(ni);    

//    prior1(0) = 4.0; // 1/sigma^2 at simga=0.5; 

//    if(tiktok == 1) {
//	for (int i = 0; i < cols; i++)
//	    cout << prior1(i) << " "; 
//	cout << endl; 
//    }

    Eigen::MatrixXd X2(ni, 2); 
    Eigen::MatrixXd xtx2(2,2);
    Eigen::VectorXd xty2(ni); 
    Eigen::VectorXd beta2(ni);

//    Eigen::MatrixXd X3(ni, dadu.ww.cols()+2); 
//    Eigen::VectorXd r1(ni);    
//    Eigen::VectorXd r0(ni);    

    long shift = tiktok*dadu.nph; 
    int max_iter = 100;    
    for (int p = 0; p < dadu.mph.cols(); p++) 
    {
	int cols = dadu.ww.cols()+2;
	Eigen::MatrixXd X(ni, dadu.ww.cols()+2); 
	Eigen::MatrixXd xtx(cols,cols);
	Eigen::VectorXd prior1(dadu.ww.cols()+2);    
	prior1.setOnes(); 
	prior1 *= 1e-20; 

	Eigen::VectorXd curph = dadu.mph.col(p); 
	int steps = 0; 
    	double eta0=1, eta1; 
	eta1 = dadu.veta0[p]; 
	double diff = 1; 
//	double like_old = -1e100; 
//	int hybrid = 1; 
	while (fabs(diff)>1e-6 && steps<max_iter) 
	{
	    steps++;  
	    eta0 = eta1; 
	    //y=x beta + e; 
	    v1.setOnes(); 
	    hh = eta0 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x1.cwiseProduct(h1inv),x2.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior1.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    //here e1 needs to scale back by multiplying 1/d1; 
	    //so r2 needs to multiply 1/d2 to scale back, which cancels with the multiplying d2 below.  

//	    double like_new = 0;
//	    for(int i = 0; i < ni; i++)
//		like_new -= log(hh(i)) / 2.0; 
//	    like_new -= log(r2.sum()) * ni / 2.0; 
//	    if(like_new < like_old) {
//		cout << "likelhood decrease detected, switch to idle IDUL." << endl; 
//		hybrid = 1;   
//	    }
//	    like_old = like_new; 
	    
	    {
		//r2 = (1, D) beta2 + err;                  
		v1.setOnes(); 
		X2.col(0) = v1.cwiseProduct(h2inv); 
		X2.col(1) = dadu.eval.cwiseProduct(h2inv); 
		xtx2 = X2.transpose() * X2; 
		xty2 = X2.transpose() * r2; 
		beta2 = xtx2.ldlt().solve(xty2);
	//	beta2 = xtx2.inverse() * xty2;
//		if(hybrid == 0) 
//        	    eta1 = beta2(1) / max(1e-20, beta2(0)); 
//		else 
		{
		    double t1 = beta2(1) / max(1e-20, beta2(0));  
		    double t2 = beta2(1) / (r2.sum() / ni);  
		    eta1 = (t1 + t2)/2.0; 
		}
	        if(eta1 >= 1e5) eta1 = 1e5; 
	        else if (eta1 <= 1e-5) eta1 = 1e-5;
	    } 

	    diff = eta1/(1+eta1) - eta0/(1+eta0); 
	}

	double loglike1 = -log(0.5); //sigam=0.5, \beta \sim N(0,sigma^2);
	{
            prior1(0) = 4.0; // 1/sigma^2 at simga=0.5; 
//            prior1(1) = 4.0; // 1/sigma^2 at simga=0.5; 
	    v1.setOnes(); 
	    hh = eta1 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x1.cwiseProduct(h1inv),x2.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior1.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta1 = xtx.ldlt().solve(xty);
    //	beta1 = xtx.inverse() * xty;
	    e1 = (y1 - X * beta1); 
	    r2 = e1.cwiseProduct(e1); 
	    loglike1 -= 0.5 * log(xtx.determinant()); 
	    loglike1 -= log(r2.sum()) * ni / 2.0;   
	}
	Eigen::MatrixXd cov = xtx.inverse();
	double tauinv = r2.sum() / (ni-cols);
	double bvar = tauinv * cov(0, 0);

	double loglikenull = 0; 
	{
	    int cols = dadu.ww.cols()+1;
	    Eigen::MatrixXd X(ni, cols); 
	    Eigen::MatrixXd xtx(cols,cols);
	    Eigen::VectorXd prior0(cols);    
	    prior0.setOnes(); 
	    prior0 *= 1e-20; 

	    v1.setOnes(); 
	    hh = eta1 * dadu.eval + v1; 
	    h2inv = hh.cwiseInverse(); 
	    h1inv = h2inv.cwiseSqrt(); 
	    Eigen::MatrixXd ww2 = dadu.ww; 
	    X << x2.cwiseProduct(h1inv), ww2.array().colwise() * h1inv.array(); 
	    xtx = X.transpose() * X; 
	    xtx += prior0.matrix().asDiagonal(); 
	    y1 = curph.cwiseProduct(h1inv); 
	    xty = X.transpose() * y1; 
	    beta0 = xtx.ldlt().solve(xty);
	    e1 = (y1 - X * beta0); 
	    r2 = e1.cwiseProduct(e1); 
	    loglikenull -= 0.5 * log(xtx.determinant()); 
	    loglikenull -= log(r2.sum()) * ni /2.0; 
	}

	double logbf = loglike1 - loglikenull; 

//	long shift = tiktok*dadu.nph; 
	dadu.pstats[shift+p].logbf = logbf; 
	dadu.pstats[shift+p].eta = eta1/(1+eta1); 
	dadu.pstats[shift+p].niter = steps; 
	dadu.pstats[shift+p].beta = beta1(0); 
	dadu.pstats[shift+p].sigma = sqrt(bvar);  
//	if (dadu.m_lrt == 1) {
//	    dadu.pstats[shift+p].l1 = loglike1; 
//	    dadu.pstats[shift+p].l0 = dadu.vloglike0[p]; 
//	}
    }
}

int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Version: %s\n", VERSION);
	fprintf(stderr, "Usage:   idul -i in.vcf.gz -p ph.txt -c cov.txt -k kinship.txt.gz [-abcevfghikoptwx]\n");
        fprintf(stderr, "Options: \n");
	fprintf(stderr, "         -0            compute maternal conditioning on paternal\n");
	fprintf(stderr, "         -1            compute paternal conditioning on maternal\n");
	fprintf(stderr, "         -a            output additional test statistics: beta, sd, etc.\n");
	fprintf(stderr, "         -b            compute bayes factor instead of p-values in the output\n");
        fprintf(stderr, "         -c str        covariate file (age, sex etc)\n");
//	fprintf(stderr, "         -d int        thin marker by a min neighor distance [1] \n");
	fprintf(stderr, "         -e str        file name for eigenvectors \n");
	fprintf(stderr, "         -v str        file name for eigenvalues \n");
	fprintf(stderr, "         -f flt        minor allele freq threshold [0.01]\n");    
	fprintf(stderr, "         -g str        read bimbam mgt format\n");    
	fprintf(stderr, "         -h            print usage\n");
	fprintf(stderr, "         -i str        (indexed) bgzipped vcf file\n");
	fprintf(stderr, "         -j int        joint (2) or conditional (1) test\n");
	fprintf(stderr, "         -k str        kinship matrix file\n");
	fprintf(stderr, "         -o str        output prefix [out]\n");
	fprintf(stderr, "         -p str        phenotype file\n");
//	fprintf(stderr, "         -q            compute logbf\n");
//	fprintf(stderr, "         -s            save bimbam mgt from vcf\n");
	fprintf(stderr, "         -t num        number of threads [16]\n");
	fprintf(stderr, "         -w            use wald test instead of lrt\n");
	fprintf(stderr, "         -x            use IDUL+ instead of IDUL\n");
//	fprintf(stderr, "         -y            output genotypes in bimbam mgt format\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Bug report: Yongtao Guan <ytguan@gmail.com>\n\n");
	return 1;
}

int main(int argc, char *argv[])
{
    time_t time0, time1, time2; 
    string fn_vcf; 
    string fn_ww; 
    string fn_ph; 
    string fn_kk; 
    string fn_evec; 
    string fn_eval; 
    string fn_gt012; 
    string pref("out");
    string af_tag("AF"); 
    string fn_untransmitted; 
    char c;
    int flag_vcf = 0; 
    int flag_kk = 0; 
    int flag_eigen = 0; 
    int flag_ww = 0; 
    int flag_gt012 = 0; 
    int flag_ph = 0; 
    int flag_help = 0; 
    int flag_pve = 0; 
    int snp_dist = 0; 
    int flag_save_mgt = 0; 
    int flag_print_additional = 0; 
    int flag_untransmitted = 0; 
    int flag_bf = 0; 
    int flag_joint_poe = 0;  //1:conditional analysis; 2: joint analysis. 

    dadu.nth = 16;
    double m_maf = 0.01; 
    double m_var = 1e-6; //2. * m_maf * (1-m_maf);  
    double m_eval_cutoff = 1e10; 
    double m_ctrl = 1; //paternal | maternal; 

    int flag_lrt = 1; 
    int flag_idul_plus = 0; 

    // Enable multi-threading
    Eigen::initParallel();

    while ((c = getopt(argc, argv, "01abc:d:e:f:g:hi:j:k:o:p:qr:st:u:v:wxy:z:")) >= 0) 
    {
	switch (c) {
	    case '0': 
		m_ctrl = 0; 
		break; 
	    case '1': 
		m_ctrl = 1; 
		break; 
	    case 'a': 
		flag_print_additional = 1; 
		break; 
		//af_tag disabled. 
	    case 'b': 
		flag_bf = 1; 
		flag_lrt = 0;  
		break;
	    case 'c': 
		flag_ww = 1; 
		fn_ww.assign(optarg); 
		break;
	    case 'd': 
		snp_dist = atoi(optarg); 
		break;
	    case 'e': 
		flag_eigen += 1; 
		fn_evec.assign(optarg); 
		break;
	    case 'v': 
		flag_eigen += 1; 
		fn_eval.assign(optarg); 
		break;
	    case 'f': 
		m_maf = atof(optarg); 
//		m_var = 2 * m_maf * (1-m_maf);
		break;
	    case 'g': 
		flag_gt012 = 1; 
		fn_gt012.assign(optarg); 
		break;
	    case 'h': 
		flag_help = 1; 
		break; 
	    case 'i': 
		flag_vcf = 1; 
		fn_vcf.assign(optarg); 
		break;
	    case 'j': 
		flag_joint_poe = atoi(optarg); 
		break;
	    case 'k': 
		flag_kk = 1; 
		fn_kk.assign(optarg); 
		break;
	    case 'o': 
		pref.assign(optarg); 
		break;
	    case 'p': 
		flag_ph = 1; 
		fn_ph.assign(optarg); 
		break;
	    case 'r': 
		flag_ph = 1; 
		fn_ph.assign(optarg); 
		break;
	    case 'q': 
		break;
	    case 's': 
		flag_save_mgt = 1; 
		break;
	    case 't': 
		dadu.nth = atoi(optarg); 
		break;
	    case 'u': 
		flag_untransmitted = 1; 
		fn_untransmitted.assign(optarg); 
		break;
	    case 'w': 
		flag_lrt = 0; 
		break;
	    case 'x': 
		flag_idul_plus = 1; 
		break;
	    case 'y':    //1 reml; 2 mle
		flag_pve = atoi(optarg); 
		break; 
	    case 'z': 
		m_eval_cutoff = atoi(optarg); 
		break;
	    default: 
		break; 
	}
    }

    dadu.ctrl = m_ctrl; 
    dadu.m_lrt = flag_lrt; 

    if (flag_help == 1 || flag_ph == 0) return usage();
    fprintf(stdout, "\n"); 
    fprintf(stdout, "IDUL %s by Yongtao Guan at Framingham Heart Study \n", VERSION); 
    fprintf(stdout, "  National Heart, Lung, and Blood Institute (C) 2023 \n"); 

    FuncPtr jacquard_ptr;  //function pointer;

    FILE * fplog = NULL;
    string buf; 
    buf.assign(pref); 
    buf.append(".log");  
    fplog = fopen(buf.c_str(), "w");
    if (fplog == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs. 
    for (int i = 0; i < argc; i++)
	fprintf(fplog, "%s ", argv[i]); 
    fprintf(fplog, "\n"); 

    dadu.flag_untransmitted = flag_untransmitted; 
    /////////////////////////////////////////////////////
    //read phenotypes; get ni; 
    Eigen::MatrixXd my1; 
    {
	int nrows = 0; 
	vector<double> vec_ph; 
	read_from_gzipped_file(fn_ph, vec_ph, nrows); 
	//a row begin with # will be ignored. 
	int ncols = vec_ph.size() / nrows; 
        my1 = Eigen::Map<Eigen::MatrixXd>(vec_ph.data(), ncols, nrows);
     	dadu.mph.conservativeResize(nrows, ncols); 
	dadu.mph = my1.transpose(); 
	Eigen::VectorXd mean = dadu.mph.colwise().mean(); 
	dadu.mph.rowwise() -= mean.transpose(); 
	Eigen::VectorXd std = (dadu.mph.array().square()).colwise().sum() / (dadu.mph.rows()-1.0); 
	std.array().sqrt(); 
	dadu.mph.array().rowwise() /= std.transpose().array(); 
        dadu.ni = nrows;
        dadu.nph = ncols;
	vector<double> ().swap(vec_ph); 
	cout << "number of samples = " << dadu.ni << endl;  
	cout << "number of phenotypes = " << dadu.nph << endl;  
    }

    int ni = dadu.ni; 
    int col = 1;   //always add a column of 1;
    Eigen::MatrixXd w2; 
    Eigen::VectorXd v1(ni); 
    v1.setOnes(ni); 
    if(flag_ww == 1) 
    {
    	vector<double> vec_ww; 
	read_from_gzipped_file(fn_ww, vec_ww); 
	int colw = vec_ww.size() / ni; 
	Eigen::MatrixXd w1 = Eigen::Map<Eigen::MatrixXd>(vec_ww.data(), colw, ni);
	Eigen::VectorXd mean = w1.rowwise().mean(); 
	w1.colwise() -= mean; 
	Eigen::VectorXd std = (w1.array().square()).rowwise().sum() / (w1.cols()-1.0); 
	std.array().sqrt(); 
	w1.array().colwise() /= std.array(); 
	col += colw; 
     	w2.conservativeResize(ni, col); 
	w2 << v1, w1.transpose() * std::sqrt(ni); 
    	vector<double> ().swap(vec_ww); 
    }
    else 
    {
     	w2.conservativeResize(ni, col); 
	w2.col(0) = v1; 
    }

    
    time0 = time(NULL); 

    /////////////////////////////////////////////////////
    if(flag_eigen == 2) 
    {
	vector<double> vec_evec; 
	read_from_gzipped_file(fn_evec, vec_evec); 
	int ni = sqrt(vec_evec.size()); 
//	dadu.ni = ni; 
	assert(ni == dadu.ni); 
        Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(vec_evec.data(), ni, ni);
	dadu.evec << A.transpose(); 
	vector<double> ().swap(vec_evec); 
	//need to test wether A.tranpose() is needed or not. 

	vector<double> vec_eval; 
	read_from_gzipped_file(fn_eval, vec_eval); 
	assert(ni == (int) vec_eval.size()); 
        Eigen::VectorXd B = Eigen::Map<Eigen::VectorXd>(vec_eval.data(), ni);
	dadu.eval << B; 
	vector<double> ().swap(vec_eval); 
    }
    else if(flag_kk == 1) 
    {   
	vector<double> vec_kk; 
	read_from_gzipped_file(fn_kk, vec_kk); 
	int ni = sqrt(vec_kk.size()); 
	assert(dadu.ni == ni); 

        Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(vec_kk.data(), ni, ni);
        std::cout << "kk.block(0, 0, 5, 5) =" << std::endl << A.block(0, 0, 5, 5) << std::endl;

      // Compute the eigen decomposition of A
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A);

        if (eigensolver.info() != Eigen::Success) {
          std::cerr << "Eigen decomposition failed!" << std::endl;
          return 1;
        }
	dadu.eval = eigensolver.eigenvalues();  
	dadu.evec = eigensolver.eigenvectors();  
//        std::cout << "dadu.evec.block(0, 0, 5, 5) =" << std::endl << dadu.evec.block(0, 0, 5, 5) << std::endl;
	vector<double> ().swap(vec_kk); 
    }
    else 
    {
	int ni = dadu.ni; 
	dadu.eval.setOnes(ni); 
	dadu.evec = Eigen::MatrixXd::Identity(ni, ni); 
    }
    for (int j = 0; j < dadu.ni; j++)
    {
	if(dadu.eval(j) < 0) 
	{
	    dadu.eval(j) = 1e-6; 
//	    dadu.evec.col(j).setZero(); 
	}
	if(dadu.eval(j) > m_eval_cutoff) 
	    dadu.eval(j) = m_eval_cutoff; 
    }
    //for negative eigenvalue, set corresponding eigenvector 0. 
//    cout << "eval: "; 
//    for (int j = 0; j < ni; j++)
//	cout << dadu.eval(j) << " "; 
//    cout << endl; 
    std::cout << "evec.block(0, 0, 5, 5) =" << std::endl << dadu.evec.block(0, 0, 5, 5) << std::endl;
    time1 = time(NULL); 
    cout << "eigendecomposition: " << difftime(time1, time0) << endl; 
	    
    /////////////////////////////////////////////////////
    dadu.ww = dadu.evec.transpose() * w2; 
//    dadu.wtw = dadu.ww.transpose() * dadu.ww; 
//    cout << "dadu.wtw = " << dadu.wtw << endl; 
//
//    for (int j = 0; j < dadu.nph; j++) {
//        Eigen::VectorXd y1 = dadu.mph.col(j); 
//	Eigen::VectorXd y2 = dadu.evec.transpose() * y1; 
//	Eigen::VectorXd wty = dadu.ww.transpose() * y2; 
//    //    Eigen::VectorXd beta = dadu.wtw.inverse() * wty; 
//	Eigen::VectorXd beta = dadu.wtw.ldlt().solve(wty); 
//	dadu.mph.col(j) = y2 - dadu.ww * beta;
//    }
    for (int j = 0; j < dadu.nph; j++) {
        Eigen::VectorXd y1 = dadu.mph.col(j); 
	dadu.mph.col(j) = dadu.evec.transpose() * y1;
    }
    /////////////////////////////////////////////////////
    if(flag_pve > 0) 
    {
	string buf1; 
	buf1.assign(pref); 
	buf1.append(".pve.gz");  
        herit(flag_pve, buf1); 
	exit(0); 
    }

    if(flag_bf == 1) 
	logpost0(); 
    else 
//    if(dadu.m_lrt == 1) 
    {
	loglike0(); 
    }
    
    if(flag_bf == 1) {
	if(flag_joint_poe == 2) 
	    jacquard_ptr = jacquard_bf2; 
	else if(flag_joint_poe == 1) 
	    jacquard_ptr = jacquard_bf1; 
	else 
	    jacquard_ptr = jacquard_bf; 
    }
    else {
	if(flag_joint_poe == 2) 
	    jacquard_ptr = jacquard_pv2; 
	else if(flag_joint_poe == 1) 
	    jacquard_ptr = jacquard_pv1; 
	else 
	{
	    if(flag_idul_plus == 0) 
	        jacquard_ptr = jacquard_pv; 
	    else
	        jacquard_ptr = jacquard_pv_plus; 
	}
    }

    /////////////////////////////////////////////////////
    string buf1; 
    buf1.assign(pref); 
    if(flag_bf == 1) 
    {
	if(flag_joint_poe == 2)
           buf1.append(".joint.logbf.gz");  
	else if(flag_joint_poe == 1) {
           if(dadu.ctrl == 0) 
              buf1.append(".0c1.logbf.gz");  
	   else 
              buf1.append(".1c0.logbf.gz");  
	}
	else 
           buf1.append(".logbf.gz");  
    }
    else {
	if(flag_joint_poe == 2)
           buf1.append(".joint.pval.gz");  
	else if(flag_joint_poe == 1) {
           if(dadu.ctrl == 0) 
              buf1.append(".0c1.pval.gz");  
	   else 
              buf1.append(".1c0.pval.gz");  
	}
	else 
           buf1.append(".pval.gz");  
    }
    gzFile fp1 = gzopen(buf1.c_str(), "w");
    if (fp1 == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf1.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs pvals etc. 

    if(flag_print_additional == 1) 
    {
	if(dadu.nph >= 1) {
	    if(flag_joint_poe == 2) {
		if(flag_bf == 1) {
		    for(int j = 1; j <= dadu.nph; j++) 
			gzprintf(fp1, "logbf_%d  beta1_%d  beta0_%d p_%d ", j, j, j, j); 
		    gzprintf(fp1, "\n"); 
		}
		else if(flag_lrt == 1) {
		    for(int j = 1; j <= dadu.nph; j++) 
			gzprintf(fp1, "p_lrt_%d  beta1_%d  beta0_%d ", j, j, j); 
		    gzprintf(fp1, "\n"); 
		}
		else {
		    for(int j = 1; j <= dadu.nph; j++) 
			gzprintf(fp1, "p_wald_%d  beta1_%d  beta0_%d ", j, j, j); 
		    gzprintf(fp1, "\n"); 
		}
	    }
	    else //if(flag_joint_poe == 1) 
	    {
		if(flag_bf == 1) {
		    for(int j = 1; j <= dadu.nph; j++) 
			gzprintf(fp1, "logbf_%d  beta_%d  sigma_%d  h_%d  iter_%d ", j, j, j, j, j); 
		    gzprintf(fp1, "\n"); 
		}
		else if(flag_lrt == 1) {
		    for(int j = 1; j <= dadu.nph; j++) 
			gzprintf(fp1, "p_lrt_%d  beta_%d  sigma_%d  l0_%d  l1_%d  h_%d  iter_%d ", j, j, j, j, j, j, j); 
		    gzprintf(fp1, "\n"); 
		}
		else {
		    for(int j = 1; j <= dadu.nph; j++) 
			gzprintf(fp1, "p_wald_%d  beta_%d  sigma_%d  h_%d  iter_%d ", j, j, j, j, j); 
		    gzprintf(fp1, "\n"); 
		}
	   }
	}
	else {
	    if(flag_joint_poe == 2) 
	    	gzprintf(fp1, "logbf  beta1  beta0  p  h  iter \n"); 
	    else if(flag_bf == 1) 
	    	gzprintf(fp1, "logbf  beta  sigma  h  iter \n"); 
	    else if(flag_lrt == 1) 
	    	gzprintf(fp1, "p_lrt  beta  sigma  l0  l1  h  iter \n"); 
	    else 
	    	gzprintf(fp1, "p_wald  beta  sigma  h  iter \n"); 
	}

    }
    else 
    {
	if(dadu.nph > 1) {
	    if(flag_bf == 1) {
		for(int j = 1; j <= dadu.nph; j++) 
		    gzprintf(fp1, "logbf_%d ", j); 
		gzprintf(fp1, "\n"); 
	    }
	    else if(flag_lrt == 1) {
		for(int j = 1; j <= dadu.nph; j++) 
		    gzprintf(fp1, "p_lrt_%d ", j); 
		gzprintf(fp1, "\n"); 
	    }
	    else {
		for(int j = 1; j <= dadu.nph; j++) 
		    gzprintf(fp1, "p_wald_%d ", j); 
		gzprintf(fp1, "\n"); 
	    }
	}
	else {
	    if(flag_bf == 1) 
	    	gzprintf(fp1, "logbf \n"); 
	    else if(flag_lrt == 1) 
	    	gzprintf(fp1, "p_lrt \n"); 
	    else 
	    	gzprintf(fp1, "p_wald \n"); 
	}
    }

    string buf2; 
    buf2.assign(pref); 
    buf2.append(".snpinfo.txt.gz");  
    gzFile fp2 = gzopen(buf2.c_str(), "w");
    if (fp2 == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf2.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs info. such as rsnumber, alleles, maf, etc.  
    gzprintf(fp2, "id  a  b  af \n"); 

    gzFile fp3 = NULL; 
    if(flag_save_mgt == 1) 
    {
	string buf3; 
	buf3.assign(pref); 
	buf3.append(".bimbam.mgt.gz");  
	fp3 = gzopen(buf3.c_str(), "w");
	if (fp3 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf3.c_str()); 
		exit(EXIT_FAILURE);
	}   // for SNPs pvals etc. 
    }

    dadu.flag_kk = flag_kk; 
    dadu.pstats = new Stats[dadu.nth*300*dadu.nph]; 
//    dadu.pv = new double[dadu.nth*100 * dadu.nph];
//    dadu.veta = new double[dadu.nth*100 * dadu.nph]; 
//    dadu.iter = new int[dadu.nth*100 * dadu.nph]; 

    fprintf(fplog, "##init thread pool n = %d\n", dadu.nth); 
    fprintf(stdout, "##init thread pool n = %d\n", dadu.nth); 
    threadpool thpool = thpool_init(dadu.nth);


    //////////////////////////////////////////////////////

    if(flag_gt012 == 1)  //this is to read bimbam mean genotype file. 
    {
	cout << fn_gt012 << endl; 
	std::ifstream file(fn_gt012, std::ios_base::in | std::ios_base::binary); 
	if(!file.is_open()) {
	    cout << "Failed to open file" << endl; 
	    return 1; 
	}
	                                      
	
	std::ifstream file2;
	if(flag_untransmitted == 1) {
	    cout << fn_untransmitted << endl;  
	    file2.open(fn_untransmitted, std::ios_base::in | std::ios_base::binary); 
	    if(!file2.is_open()) {
		cout << "Failed to open file" << endl; 
		return 1; 
	    }
	}

	try {
	    boost::iostreams::filtering_istream in; 
	    if(fn_gt012.substr(fn_gt012.find_last_of(".")) == ".gz") {
	    	in.push(boost::iostreams::gzip_decompressor()); 
	    }
	    in.push(file); 
	    std::string line; 
	    stringstream fss; 

	    boost::iostreams::filtering_istream in2; 
	    if(flag_untransmitted == 1) {
		if(fn_untransmitted.substr(fn_untransmitted.find_last_of(".")) == ".gz") {
		    in2.push(boost::iostreams::gzip_decompressor()); 
		}
	        in2.push(file2); 
	    }
	    std::string line2; 
	    stringstream fss2; 

	    dadu.ns = 0; 
	    int nlines = 0; 
	    while(std::getline(in, line)) {
		++dadu.ns;
		if(flag_untransmitted == 1) 
		    std::getline(in2, line2); 

		std::stringstream line_stream(line); 
		//Skip the first three columns (strings)
		//but write them into the output. 
		for (int i = 0; i < 3; ++i) {
		    std::string temp;
		    line_stream >> temp;
		    fss <<  temp << " "; 
		    
		}
		vector<float> vv; 

		// Read the remaining values and store them into the vector
		float value;
		int ii = 0; 
		double sum = 0; 
		double sum2 = 0; 
		while (line_stream >> value) {
//		    if(value > 2 || value < 0) 
//		    {
//			if(warnings == 1) {
//			    cout << "genotype out of the bounds " << value <<endl; 
//			    warnings == 0; 
//			}
//		    }
    //	        snpgt[ii++] = value; 
		    sum += value; 
		    sum2 += value * value; 
		    vv.push_back(value); 
		    ii++; 
		}
		sum /= ii; 
		sum2 /= ii; 
		sum2 -= sum*sum; 
//		if(sum/2.0 < m_maf || sum/2.0 > (1-m_maf) || sum2 < m_var)
//		{
//		    fss.clear(); 
//		    fss.str(""); 
//		    continue; 
//		}

		gzprintf(fp2, "%s %f \n", fss.str().c_str(), sum/2.0);   //allele frequency;
		fss.clear(); 
		fss.str(""); 
		
		for (unsigned i = 0; i < vv.size(); i++)
		    dadu.v8.push_back(vv.at(i));
		vector<float> ().swap(vv); 
		nlines++; 

		if(flag_untransmitted == 1) {
		    std::stringstream line_stream2(line2); 
		    while(line_stream2 >> value) {
			dadu.vu.push_back(value); 
		    }
		}


		if(nlines == dadu.nth*300) 
		{
		    assert((int) dadu.v8.size() == nlines*dadu.ni); 
		    long tiktok = 0; 
		    if(flag_joint_poe > 0) {
			for (int m = 0; m < nlines/3; m++)
			{
			    thpool_add_work(thpool, jacquard_ptr, (void*)tiktok);     //joint pv
			    tiktok++; 
			}
			thpool_wait(thpool); 
			print_progress_num(0, dadu.ns); 
			for (int j = 0; j < nlines/3; j++) 
			{                                                
			    int shift = j*dadu.nph; 
			    for (int p = 0; p < dadu.nph; p++) 
			    {
				Stats * ps = &dadu.pstats[p+shift]; 
				if(flag_print_additional == 1) {
				    if(dadu.nph >= 1 && flag_joint_poe == 2) {
					if(flag_bf == 1) 
					    gzprintf(fp1, "%.4e %f %f %f ", ps->logbf, ps->beta, ps->sigma, ps->pv); 
					else if(flag_lrt == 1) 
					    gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
					else 
					    gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				    }
				    else {
					if(flag_bf == 1) 
					    gzprintf(fp1, "%.4e %f %f %f %d ", ps->logbf, ps->beta, ps->sigma, ps->eta, ps->niter); 
					else if(flag_lrt == 1) 
					    gzprintf(fp1, "%.4e %f %f %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->l0, ps->l1, ps->eta, ps->niter); 
					else 
					    gzprintf(fp1, "%.4e %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->eta, ps->niter); 
				    }
				}
				else 
				{
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e  ", ps->logbf); 
				    else
					gzprintf(fp1, "%.4e  ", ps->pv); 
				}
			    }
			    gzprintf(fp1, "\n"); 
			}
			nlines = 0; 
			vector<float> ().swap(dadu.v8); 
		    }
		    else {
			for (int m = 0; m < nlines; m++)
			{
			    thpool_add_work(thpool, jacquard_ptr, (void*)tiktok);     
			    tiktok++; 
			}
			thpool_wait(thpool); 
			print_progress_num(0, dadu.ns); 
			for (int j = 0; j < nlines; j++) 
			{                                                
			    int shift = j*dadu.nph; 
			    for (int p = 0; p < dadu.nph; p++) 
			    {
				Stats * ps = &dadu.pstats[p+shift]; 
				if(flag_print_additional == 1) {
				    if(dadu.nph >= 1 && flag_joint_poe == 2) {
					if(flag_bf == 1) 
					    gzprintf(fp1, "%.4e %f %f %f ", ps->logbf, ps->beta, ps->sigma, ps->pv); 
					else if(flag_lrt == 1) 
					    gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
					else 
					    gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				    }
				    else {
					if(flag_bf == 1) 
					    gzprintf(fp1, "%.4e %f %f %f %d ", ps->logbf, ps->beta, ps->sigma, ps->eta, ps->niter); 
					else if(flag_lrt == 1) 
					    gzprintf(fp1, "%.4e %f %f %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->l0, ps-> l1, ps->eta, ps->niter); 
					else 
					    gzprintf(fp1, "%.4e %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->eta, ps->niter); 
				    }
				}
				else {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e  ", ps->logbf); 
				    else
					gzprintf(fp1, "%.4e  ", ps->pv); 
				}
			    }
			    gzprintf(fp1, "\n"); 
			}
			nlines = 0; 
			vector<float> ().swap(dadu.v8); 

		    }
		}
	    }
	    if(nlines > 0) {
		//process the remaing lines; 
		assert((int) dadu.v8.size() == nlines*dadu.ni); 
		long tiktok = 0; 
		if(flag_joint_poe > 0) {
		    for (int m = 0; m < nlines/3; m++)
		    {
			thpool_add_work(thpool, jacquard_ptr, (void*)tiktok);
			tiktok++; 
		    }
		    thpool_wait(thpool); 
		    print_progress_num(0, dadu.ns); 
		    for (int j = 0; j < nlines/3; j++) 
		    {
			int shift = j*dadu.nph; 
			for (int p = 0; p < dadu.nph; p++) 
			{
			    Stats * ps = &dadu.pstats[p+shift]; 

			    if(flag_print_additional == 1) {
				if(dadu.nph >= 1 && flag_joint_poe == 2) {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e %f %f %f ", ps->logbf, ps->beta, ps->sigma, ps->pv); 
				    else if(flag_lrt == 1) 
					gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				    else 
					gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				}
				else {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e %f %f %f %d ", ps->logbf, ps->beta, ps->sigma, ps->eta, ps->niter); 
				    else if(flag_lrt == 1) 
					gzprintf(fp1, "%.4e %f %f %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->l0, ps-> l1, ps->eta, ps->niter); 
				    else 
					gzprintf(fp1, "%.4e %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->eta, ps->niter); 
				}
			    }
			    else 
			    {
				if(flag_bf == 1) 
				    gzprintf(fp1, "%.4e  ", ps->logbf); 
				else
				    gzprintf(fp1, "%.4e  ", ps->pv); 
			    }
			}
			gzprintf(fp1, "\n"); 
		    }
		}
		else {
		    for (int m = 0; m < nlines; m++)
		    {
			thpool_add_work(thpool, jacquard_ptr, (void*)tiktok);
			tiktok++; 
		    }
		    thpool_wait(thpool); 
		    print_progress_num(0, dadu.ns); 
		    for (int j = 0; j < nlines; j++) 
		    {
			int shift = j*dadu.nph; 
			for (int p = 0; p < dadu.nph; p++) 
			{
			    Stats * ps = &dadu.pstats[p+shift]; 

			    if(flag_print_additional == 1) {
				if(dadu.nph >= 1 && flag_joint_poe == 2) {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e %f %f %f ", ps->logbf, ps->beta, ps->sigma, ps->pv); 
				    else if(flag_lrt == 1) 
					gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				    else 
					gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				}
				else {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e %f %f %f %d ", ps->logbf, ps->beta, ps->sigma, ps->eta, ps->niter); 
				    else if(flag_lrt == 1) 
					gzprintf(fp1, "%.4e %f %f %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->l0, ps-> l1, ps->eta, ps->niter); 
				    else 
					gzprintf(fp1, "%.4e %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->eta, ps->niter); 
				}
			    }
			    else {
				if(flag_bf == 1) 
				    gzprintf(fp1, "%.4e  ", ps->logbf); 
				else
				    gzprintf(fp1, "%.4e  ", ps->pv); 
			    }
			}
			gzprintf(fp1, "\n"); 
		    }
		}
	    }
	    cout << dadu.ni << endl; 
	    cout << dadu.ns << endl; 
	} catch(const boost::iostreams::gzip_error& e) {
	    std::cout << e.what() << std::endl; 
	}
    }

    ///////////////////////////////////
    int snps_by_chr[22]; 
    for (int i = 0; i < 22; i++)
	snps_by_chr[i] = 0; 
    int ns1 = 0; 
    if(flag_vcf == 1) 
    {
	htsFile *fpv = hts_open(fn_vcf.c_str(), "r");   
	bcf_hdr_t *hdr = bcf_hdr_read(fpv);  
	bcf1_t* line = bcf_init();   
	dadu.ni = bcf_hdr_nsamples(hdr); 
	fprintf(fplog, "##number of samples: %d \n", dadu.ni); 
	fprintf(stdout, "##number of samples: %d \n", dadu.ni); 
    //	if(hdr->samples != NULL) {
    //	    for (int i = 0; i < dadu.ni; i++)       
    //		cout << hdr->samples[i] << endl; 
    //	}
	dadu.ns = 0; 
//	int ns0 = 0;  
	int last_chr = -1; 
	int last_pos = 0; 


	float * snpgt = new float[dadu.ni]; 
	int nlines = 0; 
	while(bcf_read1(fpv, hdr, line) == 0) {   

	    if(ns1 % 1000 == 0) 
	        print_progress_num(0, ns1);  
	    ns1++; 

	    bcf_unpack(line,BCF_UN_STR); 
	    bcf_get_variant_types(line); 
	    if(line->d.var_type != VCF_SNP) continue; 
	    if(line->n_allele > 2) continue;       
	    if(snp_dist > 0) {
               if((last_chr != line->rid) || (line->pos - last_pos > snp_dist))
	       {
		   last_chr = line->rid; 
		   last_pos = line->pos; 
	       }
	       else 
		   continue; 
	    }  //this is to thin SNPs by bp distance.  simple, but not optimal. 

	    int32_t *gt_arr = NULL, ngt_arr = 0;
            int ngt;
	    ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
//	    cout << "ngt = " << ngt << endl; 
	    assert(ngt == dadu.ni*2); 
	    if ( ngt<=0 ) 
	    {
		cout << "no genotypes at " << bcf_seqname(hdr, line) << " " << line->pos << endl; 
		continue;
		//fiter by gt presence. 
	    }

	    for (int i=0; i<dadu.ni; i++)
	    {
		float gt = 3; 
		if (gt_arr[2*i+0] != bcf_gt_missing && gt_arr[2*i+1] != bcf_gt_missing)
		    gt = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]); 
		snpgt[i] = gt; 
	    }
	    int c3 = 0; 
	    float sum = 0; 
	    float sum2 = 0; 
	    for (int i = 0; i < dadu.ni; i++) 
	    {
		if(snpgt[i] > 2) c3++;
		else {
		    sum += snpgt[i]; 
		    sum2 += snpgt[i] * snpgt[i]; 
		}
	    }
	    sum /= (dadu.ni - c3); 
	    sum2 /= (dadu.ni - c3); 
	    sum2 -= sum*sum; 
       	    if(sum/2 < m_maf || sum/2 > 1-m_maf || sum2 < m_var)
       	       continue; 
	    
	    free(gt_arr);
	    dadu.ns++; 
            snps_by_chr[line->rid]++;  
	    //all satisfied, taking in the data; 

	    if(flag_save_mgt == 1) 
	    {
	        gzprintf(fp3, "%d-%d %s %s ", line->rid+1, line->pos+1, line->d.allele[0], line->d.allele[1]); 
		for (int i = 0; i < dadu.ni; i++) 
		{
		    if(snpgt[i] > 2) 
			gzprintf(fp3, " %.3f ", sum);  
		    else 
			gzprintf(fp3, " %.2f ", snpgt[i]);  
		}
		gzprintf(fp3, "\n");  
	    }
	              
	    gzprintf(fp2, "%d-%d %s %s %f\n", line->rid+1, line->pos+1, line->d.allele[0], line->d.allele[1], sum/2.0); 

	    for (int i = 0; i < dadu.ni; i++) 
	    {
		if(snpgt[i] > 2) 
		    dadu.v8.push_back(sum); 
		else 
		    dadu.v8.push_back(snpgt[i]); 
	    }
	    nlines++; 
	    //fill the missing gentoypes by the mean of the rest. 
	    if(nlines == dadu.nth*300) 
	    {
		assert((int) dadu.v8.size() == nlines*dadu.ni); 
		long tiktok = 0; 
		if(flag_joint_poe > 0) {
		    for (int m = 0; m < nlines/3; m++)
		    {
			thpool_add_work(thpool, jacquard_ptr, (void*)tiktok);
			tiktok++; 
		    }
		    thpool_wait(thpool); 
    //		print_progress_num(0, dadu.ns); 
		    for (int j = 0; j < nlines/3; j++) 
		    {
			int shift = j*dadu.nph; 
			for (int p = 0; p < dadu.nph; p++) 
			{
			    Stats * ps = &dadu.pstats[p+shift]; 
			    if(flag_print_additional == 1) {
				if(dadu.nph >= 1 && flag_joint_poe == 2) {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e %f %f %f ", ps->logbf, ps->beta, ps->sigma, ps->pv); 
				    else if(flag_lrt == 1) 
					gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				    else 
					gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				}
				else {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e %f %f %f %d ", ps->logbf, ps->beta, ps->sigma, ps->eta, ps->niter); 
				    else if(flag_lrt == 1) 
					gzprintf(fp1, "%.4e %f %f %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->l0, ps-> l1, ps->eta, ps->niter); 
				    else 
					gzprintf(fp1, "%.4e %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->eta, ps->niter); 
				}
			    }
			    else 
			    {
				if(flag_bf == 1) 
				    gzprintf(fp1, "%.4e  ", ps->logbf); 
				else
				    gzprintf(fp1, "%.4e  ", ps->pv); 
			    }
			}
			gzprintf(fp1, "\n"); 
		    }
		    nlines = 0; 
		    vector<float> ().swap(dadu.v8); 
		}
		else {
		    for (int m = 0; m < nlines; m++)
		    {
			thpool_add_work(thpool, jacquard_ptr, (void*)tiktok);
			tiktok++; 
		    }
		    thpool_wait(thpool); 
    //		print_progress_num(0, dadu.ns); 
		    for (int j = 0; j < nlines; j++) 
		    {
			int shift = j*dadu.nph; 
			for (int p = 0; p < dadu.nph; p++) 
			{
			    Stats * ps = &dadu.pstats[p+shift]; 
			    if(flag_print_additional == 1) {
				if(dadu.nph >= 1 && flag_joint_poe == 2) {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e %f %f %f", ps->logbf, ps->beta, ps->sigma, ps->pv); 
				    else if(flag_lrt == 1) 
					gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				    else 
					gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				}
				else {
				    if(flag_bf == 1) 
					gzprintf(fp1, "%.4e %f %f %f %d ", ps->logbf, ps->beta, ps->sigma, ps->eta, ps->niter); 
				    else if(flag_lrt == 1) 
					gzprintf(fp1, "%.4e %f %f %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->l0, ps-> l1, ps->eta, ps->niter); 
				    else 
					gzprintf(fp1, "%.4e %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->eta, ps->niter); 
				}
			    }
			    else {
				if(flag_bf == 1) 
				    gzprintf(fp1, "%.4e  ", ps->logbf); 
				else
				    gzprintf(fp1, "%.4e  ", ps->pv); 
			    }
			}
			gzprintf(fp1, "\n"); 
		    }
		    nlines = 0; 
		    vector<float> ().swap(dadu.v8); 
		}
	    }
        }
	if(nlines > 0) {
	    //process the remaing lines; 
	    assert((int)dadu.v8.size() == nlines*dadu.ni); 
	    long tiktok = 0; 
	    if(flag_joint_poe > 0) {
		for (int m = 0; m < nlines/3; m++)
		{
	            thpool_add_work(thpool, jacquard_ptr, (void*)tiktok);
		    tiktok++; 
		}
		thpool_wait(thpool); 
    //	    print_progress_num(0, dadu.ns); 
		for (int j = 0; j < nlines/3; j++) 
		{
		    int shift = j*dadu.nph; 
		    for (int p = 0; p < dadu.nph; p++) 
		    {
			Stats * ps = &dadu.pstats[p+shift]; 
			if(flag_print_additional == 1) {
	                    if(dadu.nph >= 1 && flag_joint_poe == 2) {
				if(flag_bf == 1) 
			            gzprintf(fp1, "%.4e %f %f %f", ps->logbf, ps->beta, ps->sigma, ps->pv); 
				else if(flag_lrt == 1) 
				    gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				else 
				    gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
			    }
			    else {
				if(flag_bf == 1) 
				    gzprintf(fp1, "%.4e %f %f %f %d ", ps->logbf, ps->beta, ps->sigma, ps->eta, ps->niter); 
				else if(flag_lrt == 1) 
				    gzprintf(fp1, "%.4e %f %f %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->l0, ps-> l1, ps->eta, ps->niter); 
				else 
				    gzprintf(fp1, "%.4e %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->eta, ps->niter); 
			    }
			}
			else 
			{
			    if(flag_bf == 1) 
				gzprintf(fp1, "%.4e  ", ps->logbf); 
			    else
				gzprintf(fp1, "%.4e  ", ps->pv); 
			}
		    }
		    gzprintf(fp1, "\n"); 
		}
	    }
	    else {
		for (int m = 0; m < nlines; m++)
		{
	            thpool_add_work(thpool, jacquard_ptr, (void*)tiktok);
		    tiktok++; 
		}
		thpool_wait(thpool); 
    //	    print_progress_num(0, dadu.ns); 
		for (int j = 0; j < nlines; j++) 
		{
		    int shift = j*dadu.nph; 
		    for (int p = 0; p < dadu.nph; p++) 
		    {
			Stats * ps = &dadu.pstats[p+shift]; 
			if(flag_print_additional) {
	                    if(dadu.nph >= 1 && flag_joint_poe == 2) {
				if(flag_bf == 1) 
			            gzprintf(fp1, "%.4e %f %f %f", ps->logbf, ps->beta, ps->sigma, ps->pv); 
				else if(flag_lrt == 1) 
				    gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
				else 
				    gzprintf(fp1, "%.4e %f %f ", ps->pv, ps->beta, ps->sigma); 
			    }
			    else {
				if(flag_bf == 1) 
				    gzprintf(fp1, "%.4e %f %f %f %d ", ps->logbf, ps->beta, ps->sigma, ps->eta, ps->niter); 
				else if(flag_lrt == 1) 
				    gzprintf(fp1, "%.4e %f %f %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->l0, ps-> l1, ps->eta, ps->niter); 
				else 
				    gzprintf(fp1, "%.4e %f %f %f %d ", ps->pv, ps->beta, ps->sigma, ps->eta, ps->niter); 
			    }
			}
			else {
			    if(flag_bf == 1) 
				gzprintf(fp1, "%.4e  ", ps->logbf); 
			    else
				gzprintf(fp1, "%.4e  ", ps->pv); 
			}
		    }
		    gzprintf(fp1, "\n"); 
		}
	    }
	}
	delete[] snpgt; 


	fprintf(fplog, "##number of biallelic SNPs: %d \n", dadu.ns); 
	fprintf(stdout, "##number of biallelic SNPs: %d \n", dadu.ns); 
//	fprintf(fplog, "##number of biallelic SNPs without tag: %s = %d \n", af_tag.c_str(), ns0); 
//	fprintf(stdout, "##number of biallelic SNPs without tag: %s = %d \n", af_tag.c_str(), ns0); 
	for (int i = 0; i < 22; i++)
	   fprintf(fplog, "%d %d \n", i+1, snps_by_chr[i]); 
	fflush(fplog); 
	hts_close(fpv); 
	if(dadu.ns == 0) 
	{
	   fprintf(fplog, "##no valid SNPs. abort. \n"); 
	   fprintf(stdout, "##no valid SNPs. abort. \n"); 
	   fflush(fplog); 
	   fclose(fplog); 
	   exit(0); 
	}
    }
    print_progress_num(1, ns1);  
    gzclose(fp1); 
    gzclose(fp2); 
    if(flag_save_mgt == 1) 
    	gzclose(fp3); 
///////////////////////////////////////////////////////
    thpool_wait(thpool);
    thpool_destroy(thpool);
    fprintf(fplog, "##drain thread pool. \n"); 
    fprintf(stdout, "##drain thread pool.\n"); 
///////////////////////////////////////////////////////

    time2 = time(NULL); 
    time_t dt1 = difftime (time1, time0); 
    time_t dt2 = difftime (time2, time1); 
    fprintf(fplog, "##read vcf and eigen decomposition < %ld seconds\n", dt1+1); 
    fprintf(fplog, "##compute pval < %ld seconds\n", dt2+1); 

    delete[] dadu.pstats; 
//    delete[] dadu.pv; 
//    delete[] dadu.veta; 
//    delete[] dadu.iter; 
    return 0;
}




