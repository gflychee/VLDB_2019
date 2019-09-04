#include "prelude.hpp"
#include "rand.hpp"
#include "io.hpp"
#include "logging.hpp"
#include "git_info.hpp"
#include "Alo4.hpp"
#include "experiment_reporter.hpp"
#include "scores.hpp"
#include "counts_cache.hpp"
boost::program_options::variables_map
parse_args(int argc, char** argv) {
	using namespace boost;
	namespace po = boost::program_options;

	po::options_description desc("k-median clustering");
	desc.add_options()
	("help", "produce help message")
	("revision", "get the git revision of the code")
	("debug", "print debug output")
	("trace", "print trace output")
	("graph", po::value<std::string>(),
	 "input graph")
	("target,k", po::value<size_t>(),
	 "desired number of clusters")
	("epsilon", po::value<double>()->default_value(0.1),
	 "tolerated absolute error")
	("lambda", po::value<double>()->default_value(0.1),
	 "tolerated absolute error")
	("celf", po::value<bool>()->default_value(0),
	 "if false use greedy without celf")
	("random", po::value<bool>()->default_value(0),
	 "if false use greedy without lazy")
	("topk", po::value<bool>()->default_value(0),
	 "if true use topk")
	("delta", po::value<double>()->default_value(0.01),
	 "error probability")
	("seed", po::value<uint64_t>(),
	 "seed for random generator")
	("avpr", po::value<bool>()->default_value(0),
	 "if false use greedy without avpr")
	("acr", po::value<bool>()->default_value(0),
	 "if false use greedy without acr");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("revision")) {
		std::cout << g_GIT_SHA1 << std::endl;
		std::exit(0);
	}
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		std::exit(0);
	}
	if (vm.count("trace")) {
		logging::set_level(logging::Level::Trace);
	} else if (vm.count("debug")) {
		logging::set_level(logging::Level::Debug);
	}
	return vm;
}

void check_num_components(const ugraph_t & graph, const size_t target) {
	auto component_map = boost::make_vector_property_map<int>(boost::get(boost::vertex_index, graph));
	size_t num_components = boost::connected_components(graph, component_map);
	if (target < num_components) {
		LOG_ERROR("The target size ("
		          << target
		          << " + "
		          << ") is smaller than the number of connected components ("
		          << num_components << "): the algorithm can't terminate");
		throw std::logic_error("Target size too small");
	}
}

void add_clustering_info(const ugraph_t &graph,
                         const std::vector<ClusterVertex> &vinfo,
                         const std::string & table_name) {
	size_t n = vinfo.size();
	for (ugraph_vertex_t v = 0; v < n; v++) {
		ugraph_vertex_t center = vinfo[v].center();
		EXPERIMENT_APPEND(table_name, {{"id", v},
			{"center", center},
			{"label", graph[v].label},
			{"center label", graph[center].label},
			{"probability", vinfo[v].probability()}
		});
	}

}

double log_n_k(size_t n, size_t k) {
	double ans = 0;
	for (size_t j = 0; j < k; j++)
		ans += log(n - j) - log(k - j);
	return ans;
};

int main(int argc, char**argv) {
	auto args = parse_args(argc, argv);

	std::string graph_path(args["graph"].as<std::string>());

	uint64_t seed;
	if (args.count("seed")) {
		seed = args["seed"].as<uint64_t>();
	} else {
		std::random_device rd;
		seed = 7;
		LOG_INFO("Using random seed " << seed);
	}
	double
	epsilon = args["epsilon"].as<double>(),
	lambda = args["lambda"].as<double>(),
	delta = args["delta"].as<double>();
	bool
	celf = args["celf"].as<bool>();
	bool
	topk = args["topk"].as<bool>();
	bool
	random = args["random"].as<bool>();
	bool
	avpr = args["avpr"].as<bool>();
	bool
	acr = args["acr"].as<bool>();
	size_t
	k = args["target"].as<size_t>();
	auto omp_threads = omp_get_max_threads();
	LOG_INFO("Running with " << omp_threads << " threads");
	EXPERIMENT_TAG("algorithm", std::string("k-median"));
	EXPERIMENT_TAG("input", graph_path);
	EXPERIMENT_TAG("epsilon", epsilon);
	EXPERIMENT_TAG("lambda", lambda);
	EXPERIMENT_TAG("delta", delta);
	EXPERIMENT_TAG("seed", seed);
	EXPERIMENT_TAG("k", k);
	EXPERIMENT_TAG("num-threads", omp_threads);
	EXPERIMENT_TAG("celf", celf);
	EXPERIMENT_TAG("topk", topk);
	EXPERIMENT_TAG("random", random);
	EXPERIMENT_TAG("avpr", avpr);
	EXPERIMENT_TAG("acr", acr);
	ugraph_t graph;
	read_edge_list(graph, graph_path);
	LOG_INFO("Loaded graph with " << boost::num_vertices(graph) <<
	         " nodes and " << boost::num_edges(graph) << " edges");
	check_num_components(graph, k);
	delta = (double)1/boost::num_vertices(graph);
	const double e = exp(1);
	lambda = (2*(2*e-1)*(e*epsilon+6*e-3)*boost::num_vertices(graph)/(3*pow2(e)*epsilon*epsilon*k)*(log(3)+logcnk(boost::num_vertices(graph),k)+k*log(boost::num_vertices(graph)-k)-log(delta)) +1)/1000.0;
	LOG_INFO("lambda: "<< lambda);
	Splitmix64 seeder(seed);
	Xorshift1024star rnd(seeder.next());
	KMSampler sampler(graph, 0, seed, omp_threads);
	std::vector<ClusterVertex> clustering;
	auto start = std::chrono::steady_clock::now();
	const auto n = boost::num_vertices(graph);
	ConnectionCountsCache cccache(k);
	clustering = Algorith4(graph, sampler, k, epsilon, lambda, delta, celf,cccache,random,topk);
	auto end = std::chrono::steady_clock::now();
	double elapsed = std::chrono::duration_cast< std::chrono::milliseconds >(end - start).count();
	string filename;
	filename = "cluster_"+to_string(k)+".txt";
	std::ofstream outc(filename,std::ios::app);
	outc<<n<<endl;
	for (ugraph_vertex_t v = 0; v < n; v++) {
		ugraph_vertex_t center = clustering[v].center();
		outc<<v<<" "<<center<<" "<<clustering[v].probability()<<endl;
	}
	outc.close();
	EXPERIMENT_APPEND("performance", {{"time", elapsed},});
	std::ofstream out("resKMPC.txt",std::ios::app);
	out<< "\t" << epsilon <<
	   "\t " << delta <<
	   "\t" <<  lambda <<
	   "\t " << k <<
	   "\t" << elapsed/1000 <<","<<
	   "\t" << logging::vm_size() ;
	out.close();
	add_clustering_info(graph, clustering, "clustering");
	add_scores(graph, clustering, sampler,cccache, acr,avpr);
	EXPERIMENT_SAVE();
	LOG_INFO(elapsed << " ms elapsed.");

}
