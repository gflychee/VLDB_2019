#pragma once
#include "prelude.hpp"
#include "types.hpp"
#include "rand.hpp"
#include "logging.hpp"
#include "counts_cache.hpp"
struct KMSamplerThreadState {

	KMSamplerThreadState(const ugraph_t &graph)
		: connection_counts(std::vector<size_t>(boost::num_vertices(graph))) {};

	std::vector< size_t > connection_counts;
};

namespace std {
	std::ostream & operator<<(std::ostream &os, KMSamplerThreadState & tstate);
}
class KMSampler {
	public:
		typedef std::vector< int > component_vector_t;
		// The minimum connection probability that is estimate reliably
		probability_t m_min_probability = 1.0;
		KMSampler(const ugraph_t & graph,
		          std::function<size_t()> prob_to_samples,
		          uint64_t seed,
		          size_t num_threads)
			: prob_to_samples(prob_to_samples),
			  m_thread_states(std::vector< KMSamplerThreadState >()),
			  seed(seed),
			  n2c(std::vector< std::vector< int > >(boost::num_vertices(graph))),
			  c2n(std::vector< std::vector< int > >()) {
			cc_cnt = std::vector<int>(boost::num_vertices(graph),0);
			for (size_t i=0; i<num_threads; ++i) {
				m_thread_states.emplace_back(graph);
			}

		}
		void min_probability(const ugraph_t & graph, probability_t prob);
		// void sample(const ugraph_t & graph);

		void sample_size(const ugraph_t & graph, size_t total_samples);
		void connection_probabilities(const ugraph_t & graph,
		                              const ugraph_vertex_t from,
		                              bool flag,
		                              std::vector< probability_t > & probabilities) ;
		void connection_probabilities_cache(const ugraph_t & graph,
		                                    const ugraph_vertex_t from,
		                                    ConnectionCountsCache & cccache,
		                                    std::vector<probability_t> & probabilities);

		void initialize_R(std::vector<ugraph_vertex_t> & u, std::vector<double> & ub);

		std::function<size_t()> prob_to_samples;

		std::vector< std::vector< int > > n2c;

		std::vector< std::vector< int > > c2n;

		uint64_t seed;

		int index=0;

		std::vector< KMSamplerThreadState > m_thread_states;

		size_t m_used_samples=0;

		std::vector< int > cc_cnt; // how many nodes in the same connected component, accumulate all samples

		void cc_count();

};
