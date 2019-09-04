#pragma once

#include "km_sampler.hpp"
#include "cluster_vertex.hpp"
#include "experiment_reporter.hpp"

probability_t min_probability(const std::vector< ClusterVertex > & vinfo);

probability_t sum_probability(const std::vector< ClusterVertex > & vinfo);

struct AVPR {
	double inner;
	double outer;

	AVPR(): inner(-1), outer(-1) {}
	AVPR(double inner, double outer): inner(inner), outer(outer) {
		REQUIRE(inner <= 1.0, "Inner AVPR should be less than 1.0");
		REQUIRE(outer <= 1.0, "Outer AVPR should be less than 1.0");
	}
};

/// Computes the Average Vertex Pairwise Reliability
AVPR average_vertex_pairwise_reliability(const ugraph_t & graph,
        const std::vector<ClusterVertex> & vinfo,
        KMSampler & sampler);

void add_scores(const ugraph_t & graph,
                const std::vector< ClusterVertex > & vinfo,
                KMSampler & sampler,
                ConnectionCountsCache & cccache,
                const bool with_acr,
                const bool with_avpr);
