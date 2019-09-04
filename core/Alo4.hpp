#pragma once
#include "types.hpp"
#include "require.hpp"
#include "km_sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"
#include "counts_cache.hpp"
#include "termcolor.hpp"
#include "kmpc.hpp"
using namespace std;

double log2(int n) {
	return log(n) / log(2);
}

double logcnk(int n, int k) {
	double ans = 0;
	for (int i = n - k + 1; i <= n; i++) {
		ans += log(i);
	}
	for (int i = 1; i <= k; i++) {
		ans -= log(i);
	}
	return ans;
}

double pow2(double a) {
	return a*a;
}

double coverNum(const ugraph_t &graph,vector<ugraph_vertex_t>  A, KMSampler & sampler, bool flag,int n) {
	double sum=0;
	std::vector<double> fval_memo(n,0);
	for(size_t i=0; i<A.size(); i++) {
		std::vector<double> probabilities(n);
		sampler.connection_probabilities(graph, A[i], flag, probabilities);
		for(size_t j = 0; j < fval_memo.size(); j++) {
			if(probabilities[j] > fval_memo[j]) {
				fval_memo[j] = probabilities[j];
			}
		}
	}
	for(size_t i = 0; i < fval_memo.size(); i++) {
		sum+=fval_memo[i];
	}
	return sum/n;

}

vector<ugraph_vertex_t > D(vector<ClusterVertex> C) {
	vector<ugraph_vertex_t>  res;
	for(size_t i=0; i<C.size(); i++) {
		if(C[i].is_center()) {
			res.push_back(i);
		}
	}
	return res;
}

vector<ClusterVertex> Algorith4(const ugraph_t &graph, KMSampler & sampler, const int targetSize, const double epsilon,const double lambda, const double delta, bool celf, ConnectionCountsCache & cccache,bool random,bool topk) {
	cout<<"enter algorith4:"<<endl;
	const auto n = boost::num_vertices(graph);
	const double e = exp(1);
	double approx = 1 - 1.0 / e;
	if(topk) {
		approx=1;
	}
	cout<<"delta:"<<delta<<endl;
	const auto maxNumR = size_t(2*(2*e-1)*(e*epsilon+6*e-3)*n/(3*pow2(e)*epsilon*epsilon*targetSize)*(log(3)+logcnk(n,targetSize)+targetSize*log(n-targetSize)-log(delta))) +1;
	const auto numRbase = size_t(1000);
	LOG_INFO("maxNumR: "<<maxNumR);
	LOG_INFO("numRbase: "<<numRbase);
	const auto numIter = size_t(log2(maxNumR / numRbase)) + 1;
	const double a1 = log(numIter * 3.0 / delta);
	const double a2 = log(numIter * 3.0 / delta);
	LOG_INFO("a1 : "<<a1);
	double time1 = 0.0, time2 = 0.0;
	std::vector<ClusterVertex> clustering;
	for (auto idx = 0; idx < numIter; idx++) {
		const auto numR = numRbase << idx;
		cout<<"numR: "<<2*numR<<endl;
		auto start = std::chrono::steady_clock::now();
		sampler.sample_size(graph, numR);
		auto end = std::chrono::steady_clock::now();
		double elapsed = std::chrono::duration_cast< std::chrono::milliseconds >(end - start).count();
		LOG_INFO("sample used time : "<<elapsed);
		clustering = SearchKM(graph, sampler, cccache, targetSize, celf, random,topk, epsilon);
		vector<ugraph_vertex_t>  A;
		A.clear();
		A = D(clustering);
		double lambda2=0;
		for(ugraph_vertex_t v=0; v<n; v++) {
			lambda2+=clustering[v].probability();
		}
		lambda2=lambda2/n;
		sampler.sample_size(graph, 2*numR);
		const auto lambda1 = coverNum(graph,A, sampler, 0,n);
		cout<<"lambda2:"<<lambda2<<endl;
		cout<<"lambda1:"<<lambda1<<endl;
		auto theta = sampler.m_used_samples/2;
		const auto lowerSelect = pow2( sqrt(lambda1+2.0*a1/9.0/theta) - sqrt(a1 / (2.0*theta)) ) - a1 / (18.0*theta);
		auto upperOPT = pow2(sqrt(lambda2/approx + 8*a2 / (9.0*theta)) + sqrt(a2 / (2.0*theta)))-a2/(18.0*theta);
		if(upperOPT>1)
			upperOPT = 1.0;
		const auto approxOPIMC = lowerSelect / upperOPT;
		cout<<"lowerSelect:"<<lowerSelect<<endl;
		cout<<"upperOPT:"<<upperOPT<<endl;
		cout<<"approxOPIMC:"<<approxOPIMC<<endl;
		cout<<"approx - epsilon: "<<approx - epsilon<<endl;
		if (approxOPIMC >= approx - epsilon) {
			return clustering;
		}
	}
	return clustering;
}
