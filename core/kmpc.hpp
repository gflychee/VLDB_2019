#pragma once

#include "types.hpp"
#include "require.hpp"
#include "km_sampler.hpp"
#include "experiment_reporter.hpp"
#include "cluster_vertex.hpp"
#include "counts_cache.hpp"
#include "termcolor.hpp"
using namespace std;
struct node {
	ugraph_vertex_t v;
	double w;
};
bool comp(node x, node y) {
	return x.w>y.w;
}

class KM_Greedy {
	public:
		KM_Greedy(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache) {
			n = boost::num_vertices(graph);
			s = std::unordered_set<ugraph_vertex_t>();
			fval_memo = std::vector<double>(n,0.0);
			deltaf = std::vector<double>(n,0.0);
			center = std::vector<ugraph_vertex_t>(n);
			UB = std::vector<node>(n);
			randomNodes = std::vector<node>();
		}
		void greedy_lazy(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache , size_t k, double epsilon);
		void greedy(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache , size_t k);
		void greedy_celf(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache , size_t k);
		void greedy_lazycelf(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache , size_t k, double epsilon);
		void getRandomNodes(int k,int n);
		void addSeeds();
		void topk(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache,size_t k);
		void update_UB(const ugraph_t & graph, KMSampler & sampler,ugraph_vertex_t c,std::vector<node> &ub_upperbound);
		std::vector<double> get_prob() {
			return fval_memo;
		}

		std::vector<ugraph_vertex_t> get_center() {
			return center;
		}

		std::unordered_set<ugraph_vertex_t> get_s() {
			return s;
		}


	private:
		size_t n;
		std::unordered_set<ugraph_vertex_t> s;
		std::vector<double> deltaf;
		std::vector<double> fval_memo;
		std::vector<ugraph_vertex_t> center;
		std::vector<node> UB;
		std::vector<node> randomNodes;


		void fval_update(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache, ugraph_vertex_t xstar);
};

void
KM_Greedy::fval_update(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache, ugraph_vertex_t xstar) {
	std::vector<double> probabilities(n);
	sampler.connection_probabilities_cache(graph, xstar, cccache, probabilities);
	for(size_t i = 0; i < fval_memo.size(); i++) {
		if(probabilities[i] > fval_memo[i]) {
			fval_memo[i] = probabilities[i];
			center[i] = xstar;
		}
	}
}
void
KM_Greedy::update_UB(const ugraph_t & graph, KMSampler & sampler,ugraph_vertex_t c,std::vector<node> &ub_upperbound) {
	size_t c2n_size = sampler.c2n.size();
	const size_t n = boost::num_vertices(graph);
	std::vector<bool > flag(c2n_size,0);
	std::vector<int > cc(n,0);
	std::vector<double> ubb(n,0);
	int sum=0;
	double kmc=0;
	//std::vector<bool > flag315(c2n_size,0);
	for(size_t i=0; i<sampler.n2c[c].size(); i++) {
		flag[sampler.n2c[c][i]]=1;
	}

	for(size_t i = 0; i<n; i++) {
		for(size_t j=0; j<sampler.n2c[i].size(); j++) {
			if(flag[sampler.n2c[i][j]]==1) {
				cc[i]++;
			}
		}
	}

	for(size_t i = 0; i<n; i++) {
		sum+=cc[i];
	}
	sum-=cc[c];
	for(size_t i=0; i<n; i++) {
		ubb[i]=sum;
	}
	for(int i=0; i<sampler.c2n.size(); i++) {
		if(flag[i]==1)
			for(size_t j=0; j<sampler.c2n[i].size(); j++) {
				if(sampler.c2n[i][j]!=c) {
					ubb[sampler.c2n[i][j]]-=1;
				}

			}
		else {
			for(size_t j=0; j<sampler.c2n[i].size(); j++) {
				ubb[sampler.c2n[i][j]]+=sampler.c2n[i].size()-1;
			}
		}
	}
	for(size_t i=0; i<n; i++) {
		kmc+=fval_memo[i];
	}

	for(size_t i=0; i<n; i++) {
		ub_upperbound[i].w=ubb[ub_upperbound[i].v]/sampler.m_used_samples+1+1-kmc;
	}

}

void
KM_Greedy::topk(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache,size_t k) {
	LOG_INFO("topk***************************");
	const size_t n = boost::num_vertices(graph);
	std::vector<ugraph_vertex_t> u(n);
	std::vector<double> ub(n);
	std::ofstream out("time_topk.txt",std::ios::app);
	auto t0 = std::chrono::steady_clock::now();
	sampler.initialize_R(u, ub);
	auto t1 = std::chrono::steady_clock::now();
	double time1 = std::chrono::duration_cast< std::chrono::milliseconds >(t1 - t0).count();
	out<<"initialize_R:   "<<time1<<std::endl;
	//Initialize the structure
	for(size_t i = 0; i < n; i++) {
		UB[i].v = u[i];
		UB[i].w = ub[u[i]]/sampler.m_used_samples;
		//std::cout<<UB[i].w<<" ";
	}
	//std::cout<<std::endl;
	sort(UB.begin(),UB.end(),comp);
	//LOG_INFO("insert "<< UB[0].v <<" to s");
	for(int i=0; i<k; i++) {
		s.insert(UB[i].v);
		fval_update(graph, sampler, cccache, UB[i].v);
	}
	auto t2 = std::chrono::steady_clock::now();
	double time2 = std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t0).count();
	out<<"select topk used time :   "<<time2<<std::endl;
}


void
KM_Greedy::greedy_celf(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache,size_t k) {
	std::ofstream out("time_greedy_celf.txt",std::ios::app);
	std::ofstream outub("update_greedy_celf.txt",std::ios::app);
	auto t0 = std::chrono::steady_clock::now();
	const size_t n = boost::num_vertices(graph);
	std::vector<ugraph_vertex_t> u(n);
	std::vector<double> ub(n);
	sampler.initialize_R(u, ub);
	vector<node> test(n);
	//Initialize the structure
	for(size_t i = 0; i < n; i++) {
		UB[i].v = u[i];
		UB[i].w = ub[u[i]]/sampler.m_used_samples;
	}
	sort(UB.begin(),UB.end(),comp);
	LOG_INFO("insert "<< UB[0].v <<" to s");
	s.insert(UB[0].v);
	fval_update(graph, sampler, cccache, UB[0].v);
	UB[0].v=-1,UB[0].w=-1;
	sort(UB.begin(),UB.end(),comp);
	while(s.size() < k) {
		//int num=0;
		auto t1 = std::chrono::steady_clock::now();
		int sum=0;
		double time4=0;
		priority_queue<double> max_valve;
		for(ugraph_vertex_t x = 0; x < n; x++) {
			if(UB[x].v==-1) continue; //avoid to repeat calculate x
			//num++;
			double xx=UB[x].w;
			UB[x].w = 0;
			sum++;
			std::vector<double> probabilities(n);
			auto t5 = std::chrono::steady_clock::now();
			sampler.connection_probabilities_cache(graph, UB[x].v, cccache, probabilities);
			auto t6 = std::chrono::steady_clock::now();
			time4 += std::chrono::duration_cast< std::chrono::milliseconds >(t6 - t5).count();
			for(ugraph_vertex_t v = 0; v < n; v++) {
				if(probabilities[v] > fval_memo[v])
					UB[x].w += probabilities[v] - fval_memo[v];
			}
			if(sum%100==0) {
				out<<"compute  "<<sum<<" times probabilities used "<<time4<<std::endl;
				auto t2 = std::chrono::steady_clock::now();
				double time1 = std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t1).count();
				out<< sum <<" nodes : "<<time1<<std::endl;
			}
			max_valve.push(UB[x].w);
			if((x==n-1)||(max_valve.top() >= UB[x+1].w)) {
				// LOG_INFO("xx: "<<xx );
				std::cout<<"max_valve.top() >= UB[x+1].w "<<max_valve.top()<<" "<<UB[x+1].w<<std::endl;
				LOG_INFO("sum: "<<sum );
				cout<<"s.size(): "<<s.size()<<endl;
				//getchar();
				break;
			}
		}
		sort(UB.begin(),UB.end(),comp);
		if(UB[0].w >= 0) {
			LOG_INFO("insert "<< UB[0].v <<" to s");
			s.insert(UB[0].v);
		} else {
			std::cout<<"UB[0].w<0"<<std::endl;
		}
		fval_update(graph, sampler, cccache, UB[0].v);
		UB[0].v=-1,UB[0].w=-1;
	}
	auto t4 = std::chrono::steady_clock::now();
	double time3 = std::chrono::duration_cast< std::chrono::milliseconds >(t4 - t0).count();
	out<<"1 greedy_celf : "<<time3<<std::endl;
	out.close();
}

void
KM_Greedy::getRandomNodes(int k,int n) { //pick up k nodes from n nodes
	randomNodes.clear();
	int i;
	for(i=0; i<n; i++)
		if(rand()%(n -i) < k) {
			randomNodes.push_back(UB[i]);
			//std::cout<<i<<" ";
			k--;
		}
}

void
KM_Greedy::greedy_lazycelf(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache,size_t k, double epsilon) {
	std::cout<<"n/k*log(1/epsilon):"<<n/k*log(1/epsilon)<<std::endl;
	std::ofstream out("time_greedy_lazycelf.txt",std::ios::app);
	auto t0 = std::chrono::steady_clock::now();
	srand((unsigned)time(0));
	const size_t n = boost::num_vertices(graph);
	std::vector<ugraph_vertex_t> u(n);
	std::vector<double> ub(n);
	sampler.initialize_R(u, ub);
	//Initialize the structure
	for(size_t i = 0; i < n; i++) {
		UB[i].v = u[i];
		UB[i].w = ub[u[i]]/sampler.m_used_samples;
	}
	sort(UB.begin(),UB.end(),comp);
	LOG_INFO("insert "<< UB[0].v <<" to s");
	s.insert(UB[0].v);
	fval_update(graph, sampler, cccache, UB[0].v);
	UB[0].v=-1,UB[0].w=-1;
	while(s.size() < k) {
		auto t1 = std::chrono::steady_clock::now();
		int sum=0;
		double time4=0;
		getRandomNodes(n/k*log(1/epsilon),n);
		priority_queue<double> max_valve;
		node maxNode=randomNodes[0];
		sort(randomNodes.begin(),randomNodes.end(),comp);

		for(ugraph_vertex_t i = 0; i<randomNodes.size(); i++) {
			node x = randomNodes[i];
			if(x.v==-1) continue; //avoid to repeat calculate x
			x.w = 0;
			sum++;
			std::vector<double> probabilities(n);
			auto t5 = std::chrono::steady_clock::now();
			sampler.connection_probabilities_cache(graph, x.v, cccache, probabilities);
			auto t6 = std::chrono::steady_clock::now();
			time4 += std::chrono::duration_cast< std::chrono::milliseconds >(t6 - t5).count();
			for(ugraph_vertex_t v = 0; v < n; v++) {
				if(probabilities[v] > fval_memo[v])
					x.w += probabilities[v] - fval_memo[v];
			}
			for(int i=0; i<n; i++) {
				if(UB[i].v==x.v) {
					UB[i].w=x.w;
					break;
				}
			}
			if(max_valve.size()==0) {
				max_valve.push(x.w);
				maxNode=x;
			} else if(x.w > max_valve.top()) {
				maxNode=x;
			}
			max_valve.push(x.w);
			if(sum%100==0) {
				out<<"compute  "<<sum<<" times probabilities used "<<time4<<std::endl;
				auto t2 = std::chrono::steady_clock::now();
				double time1 = std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t1).count();
				out<< sum <<" nodes : "<<time1<<std::endl;
			}
			if((i==randomNodes.size()-1)||(max_valve.top() >= randomNodes[i+1].w)) {
				break;
			}
		}
		if(maxNode.w >= 0) {
			s.insert(maxNode.v);
		} else {
			std::cout<<"maxNode.w<0"<<std::endl;
		}
		fval_update(graph, sampler, cccache, maxNode.v);
		maxNode.v=-1,maxNode.w=-1;
	}
	auto t4 = std::chrono::steady_clock::now();
	double time3 = std::chrono::duration_cast< std::chrono::milliseconds >(t4 - t0).count();
	out<<"1 greedy_celf : "<<time3<<std::endl;
	out.close();
}

void
KM_Greedy::greedy(const ugraph_t & graph, KMSampler & sampler, ConnectionCountsCache & cccache,size_t k) {
	const size_t n = boost::num_vertices(graph);

	while(s.size() < k) {
		std::fill(deltaf.begin(), deltaf.end(), 0.0);
		ugraph_vertex_t xstar;
		double max_deltaf = 0.0;

		for(ugraph_vertex_t x = 0; x < n; x++) {
			if(s.count(x)) continue; //avoid to repeat calculate x
			deltaf[x]=0;
			std::vector<double> probabilities(n);
			sampler.connection_probabilities_cache(graph, x, cccache, probabilities);
			for(ugraph_vertex_t v = 0; v < n; v++) {
				if(probabilities[v] > fval_memo[v])
					deltaf[x] += probabilities[v] - fval_memo[v];
			}
			deltaf[x] /= n;
			if(max_deltaf < deltaf[x]) {
				max_deltaf = deltaf[x];
				xstar = x;
			}
		}
		LOG_INFO("insert "<< xstar <<" to s");
		s.insert(xstar);
		fval_update(graph, sampler, cccache, xstar);

	}
}

void
KM_Greedy::addSeeds() {
	ifstream readFile("seeds.txt");
	if(!readFile) {
		cout<<"No seeds.txt!";
	}
	ugraph_vertex_t a=0;
	while(readFile>>a) {
		s.insert(a);
	}
}

std::vector< ClusterVertex >
SearchKM(const ugraph_t & graph, KMSampler & sampler,
         ConnectionCountsCache & cccache, const size_t k, bool celf, bool random,bool topk, double epsilon) {

	const size_t n = boost::num_vertices(graph);
	KM_Greedy g(graph, sampler, cccache);

	if(topk) {
		g.topk(graph, sampler, cccache,k);
	} else if(random&&celf) {
		g.greedy_lazycelf(graph, sampler, cccache,k, epsilon);
	} else if(celf&&(!random)) {
		g.greedy_celf(graph, sampler, cccache,k);
	} else {
		g.greedy(graph, sampler, cccache,k);
	}

	std::vector< ClusterVertex > valid_clustering(n);
	std::unordered_set< ugraph_vertex_t > s = g.get_s();
	std::vector< ugraph_vertex_t > center = g.get_center();
	std::vector<double> prob = g.get_prob();

	for(ugraph_vertex_t v = 0; v < n; v++) {
		if(s.count(v)) {
			valid_clustering[v].make_center(v);
		} else {
			valid_clustering[v].cover(center[v],prob[v]);
		}
	}
	return valid_clustering;
}

