#include "km_sampler.hpp"
#include "logging.hpp"
#include "connected_components.hpp"

void KMSampler::min_probability(const ugraph_t & graph, probability_t prob) {
	m_min_probability = (prob < m_min_probability)? prob : m_min_probability;
}

void KMSampler::sample_size(const ugraph_t & graph, size_t total_samples ) {
	int n = boost::num_vertices(graph);
	Xorshift1024star rnd(seed);
	size_t new_samples = total_samples - m_used_samples;
	LOG_INFO("Using " << total_samples << " (taking " << new_samples << " new)");
	size_t start = m_used_samples;
	auto t1 = std::chrono::steady_clock::now();
	for (size_t i = start; i < start+new_samples; ++i) {
		if(i%100==0&&i!=0) {
			auto t2 = std::chrono::steady_clock::now();
			double elapsed = std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t1).count();
			LOG_INFO("num_samples "<<i<<" used time: "<<elapsed);
		}

		std::vector< int > components(n,-1);
		std::vector< size_t >  ranks(n,0);
		union_find(graph, rnd, ranks, components);
		std::vector< int > number(boost::num_vertices(graph),-1);
		std::vector< int > sum(boost::num_vertices(graph),0);
		std::vector< bool > a(boost::num_vertices(graph),0);


		for(size_t j=0; j<boost::num_vertices(graph); j++) {
			sum[components[j]]++;
		}
		for(size_t j=0; j<boost::num_vertices(graph); j++) {
			if(sum[j]>1) {
				a[j]=1;
			}
		}
		for(size_t j=0; j<boost::num_vertices(graph); j++) {
			if(a[components[j]]==1) {
				if(number[components[j]]==-1) {
					number[components[j]]=index;
					index++;
					std::vector < int > b;
					b.push_back(j);
					c2n.push_back(b);
				} else {
					c2n[number[components[j]]].push_back(j);
				}
				n2c[j].push_back(number[components[j]]);
			}
		}

	}
	m_used_samples = total_samples;
}

void KMSampler::connection_probabilities_cache(const ugraph_t & graph,
        const ugraph_vertex_t from,
        ConnectionCountsCache & cccache,
        std::vector< probability_t > & probabilities) {
	int n=boost::num_vertices(graph);
	for (auto & tstate : m_thread_states) {
		std::fill(tstate.connection_counts.begin(), tstate.connection_counts.end(), 0);
	}
	ConnectionCountsCacheElement & ccc_elem = cccache.get_or_new(from, n);
	// Accumulate, in parallel, the connection counts
	#pragma omp parallel for default(none) shared(graph,ccc_elem)
	for(int i = n2c[from].size()-1; i >=0 ; i--) {
		if(n2c[from][i]>=ccc_elem.num_samples) {
			for(int j = 0; j < c2n[n2c[from][i]].size(); j++) {
				auto tid = omp_get_thread_num();
				auto & connection_counts = m_thread_states[tid].connection_counts;
				connection_counts[c2n[n2c[from][i]][j]]++;
			}
		}
	}
	for (auto & tstate : m_thread_states) {
		for (size_t i=0; i< n; i++) {
			ccc_elem.counts[i] += tstate.connection_counts[i];
		}
	}

	ccc_elem.num_samples = index;
	for(size_t k = 0; k < n; k++) {
		probabilities[k]=ccc_elem.counts[k]/double(m_used_samples);

	}
	probabilities[from]=1;
}

void KMSampler::connection_probabilities(const ugraph_t & graph,
        const ugraph_vertex_t from,
        bool flag,
        std::vector< probability_t > & probabilities) {
	int n=boost::num_vertices(graph);
	std::vector< int > counts(n);
	if(flag) {
		for(int i = 0; i < n2c[from].size();  i++) {
			if(n2c[from][i]<index/2) {
				for(int j = 0; j < c2n[n2c[from][i]].size(); j++) {
					counts[c2n[n2c[from][i]][j]]++;
				}
			}

		}
	} else {
		for(int i = n2c[from].size()-1; i >=0 ; i--) {
			if(n2c[from][i]>index/2) {
				for(int j = 0; j < c2n[n2c[from][i]].size(); j++) {
					counts[c2n[n2c[from][i]][j]]++;
				}
			}

		}

	}

	for(size_t k = 0; k < n; k++) {
		probabilities[k]=counts[k]*2/double(m_used_samples);

	}
	probabilities[from]=1;

}

void KMSampler::initialize_R(std::vector<ugraph_vertex_t> & u, std::vector<double> & ub) {
	cc_count();

	for(size_t i = 0; i < cc_cnt.size(); i++) {
		ub[i] = 1.0 * cc_cnt[i];
	}
	std::vector<std::pair<int,int>> u_index(cc_cnt.size());
	for(size_t i = 0; i < cc_cnt.size(); i++) {
		u_index[i] = std::make_pair(cc_cnt[i],i);
	}

	sort(u_index.begin(), u_index.end(), [](std::pair<int,int> &a, std::pair<int,int> &b) -> bool {
		return a.first > b.first;
	});
	for(size_t i = 0; i < cc_cnt.size(); i++) {
		u[i] = u_index[i].second;
	}

}

void KMSampler::cc_count() {
	for(size_t i=0; i<n2c.size(); i++) {
		for(size_t j=0; j<n2c[i].size(); j++) {
			cc_cnt[i]+=c2n[n2c[i][j]].size();
		}
	}
}
