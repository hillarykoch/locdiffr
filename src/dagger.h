#ifndef _DAGGER_H
#define _DAGGER_H

#include <lemon/list_graph.h>
using namespace lemon;

class Dagger {
    private:
        int _totalNodes = 0;
        double l;
        double m;
        int _L = 0;
        int _currLayer = 1;
    public:
        const lemon::ListDigraph& _gr;
        ListDigraph::NodeMap<double> _numOutArcs;
        ListDigraph::NodeMap<double> _numInArcs;
        ListDigraph::NodeMap<int> _layer;
        ListDigraph::NodeMap<int> _label;
        ListDigraph::NodeMap<double> _effectiveLeafNum;
        ListDigraph::NodeMap<double> _effectiveNodeNum;
        ListDigraph::NodeMap<bool> _isRejected;
        ListDigraph::NodeMap<bool> _isCandidateHypothesis;
        double _alpha;
        int _max_depth;
        arma::colvec _prob_theta_equals_zero;

        Dagger(const ListDigraph& gr, const ListDigraph::NodeMap<int>& layer, const ListDigraph::NodeMap<int>& label, arma::colvec prob_theta_equals_zero, const double alpha, const int max_depth) :_gr(gr), _numOutArcs(gr), _numInArcs(gr), _layer(gr), _label(gr), _effectiveLeafNum(gr), _effectiveNodeNum(gr), _isRejected(gr, false), _isCandidateHypothesis(gr, false) {
            // Create numInArcs and numOutArcs nodeMaps
            _alpha = alpha;
            _max_depth = max_depth;
            for (ListDigraph::NodeIt n(gr); n != INVALID; ++n) {
                _numOutArcs[n] = countOutArcs(gr, n);
                _numInArcs[n] = countInArcs(gr, n);
                _layer[n] = layer[n];
                _label[n] = label[n];
                _prob_theta_equals_zero = prob_theta_equals_zero;
                ++_totalNodes;

                if(_layer[n] == 1) {
                    _isCandidateHypothesis[n] = true;
                    _L++;
                }
            }

            if(_totalNodes != prob_theta_equals_zero.size()) {
                std::cout << "wrong number of probabilities passed" << std::endl;
            }
            
            // Calculate effective number of leaves and effective number of nodes nodeMaps
            for (ListDigraph::NodeIt n(gr); n != INVALID; ++n) {
                l = 0.0;
                m = 1.0;
                // At final layer, effective leaves and nodes are trivially 1
                if(layer[n] == max_depth) {
                    _effectiveLeafNum[n] = 1;
                    _effectiveNodeNum[n] = 1;
                } else { // Otherwise, pour the water
                    for (ListDigraph::OutArcIt a(gr, n); a != INVALID; ++a) {
                        l +=  (1/_numInArcs[gr.target(a)]) * _effectiveLeafNum[gr.target(a)];
                        m += (1/_numInArcs[gr.target(a)]) * _effectiveNodeNum[gr.target(a)];
                    }
                    _effectiveLeafNum[n] = l;
                    _effectiveNodeNum[n] = m;
                }
            }

            
        }

        int total_nodes() {
            return _totalNodes;
        }
        
        double out_arcs(const ListDigraph::Node& n) const {
            return _numOutArcs[n];
        }

        double in_arcs(const ListDigraph::Node& n) const {
            return _numInArcs[n];
        }

        double effective_leaf_num(const ListDigraph::Node& n) const {
            return _effectiveLeafNum[n];
        }

        double effective_node_num(const ListDigraph::Node& n) const {
            return _effectiveNodeNum[n];
        }

        int layer(const ListDigraph::Node& n) const {
            return _layer[n];
        }

        int label(const ListDigraph::Node& n) const {
            return _label[n];
        }

        int curr_layer() const {
            return _currLayer;
        }

        bool is_rejected(const ListDigraph::Node& n) const {
            return _isRejected[n];
        }

        bool is_candidate_hypothesis(const ListDigraph::Node& n) const {
            return _isCandidateHypothesis[n];
        }

        int num_currently_rejected() const {
            int nrej = 0;
            for(ListDigraph::NodeIt n(_gr); n != INVALID; ++n) {
                if(is_rejected(n)) {
                    nrej++;
                }
            }
            return nrej;
        }

        int num_currently_candidates() const {
            int ncan = 0;
            for(ListDigraph::NodeIt n(_gr); n != INVALID; ++n) {
                if(is_candidate_hypothesis(n)) {
                    ncan++;
                }
            }
            return ncan;
        }

        arma::colvec compute_alpha_for_fixed_r(const int r, const int ncan) const {
            int num_already_rejected = num_currently_rejected();
            arma::colvec candidate_alphas(ncan, arma::fill::none);
            int count = 0;
            for(ListDigraph::NodeIt n(_gr); n != INVALID; ++n) {
                if(is_candidate_hypothesis(n)) {
                    candidate_alphas(count) = (effective_leaf_num(n) / _L) * ((effective_node_num(n) + num_already_rejected + r - 1) / effective_node_num(n)) * _alpha;
                    count++;
                }
            }
            return candidate_alphas;
        }

        std::tuple<int, std::vector<double>> num_nodes_with_null_prob_below_alphar(const int r) const {
            std::tuple<int, std::vector<double>> out;
            int ncan = num_currently_candidates();
            arma::colvec curr_alpha = compute_alpha_for_fixed_r(r, ncan);
            std::vector<double> candidate_probabilities(ncan);

            int count = 0;
            for(ListDigraph::NodeIt n(_gr); n != INVALID; ++n) {
                if(is_candidate_hypothesis(n)) {
                    candidate_probabilities[count] = _prob_theta_equals_zero(label(n)-1);
                    count++;
                }
            }

            int num_below = 0;
            for(auto i = 0; i < ncan; i++) {
                if(curr_alpha(i) <= candidate_probabilities[i]) {
                    num_below++;
                }
            }
            out = std::make_tuple(num_below, candidate_probabilities);
            return out;
        }

        void add_candidate_hypothesis(ListDigraph::Node& n) {
            _isCandidateHypothesis[n] = true;
        }

        void remove_candidate_hypothesis(ListDigraph::Node& n) {
            _isCandidateHypothesis[n] = false;
        }

        void reject(ListDigraph::Node& n) {
            _isRejected[n] = true;
            remove_candidate_hypothesis(n);
        }

        void traverse() {
            int num_rejected_source_nodes;
            _currLayer++;
            // If we're at the end, stop the procedure
            if(curr_layer() > _max_depth) {
                // Do nothing
            } else {
                for (ListDigraph::NodeIt n(_gr); n != INVALID; ++n) {
                    // If all source nodes are rejected, make this a candidate hypothesis
                    if(layer(n) == curr_layer()) {
                        num_rejected_source_nodes = 0;
                        for(ListDigraph::InArcIt a(_gr, n); a != INVALID; ++a) {
                            if(is_rejected(_gr.source(a))) {
                                num_rejected_source_nodes++;
                            }
                        }
                        if(num_rejected_source_nodes == in_arcs(n)) {
                            add_candidate_hypothesis(n);
                        }
                    }
                }
            }
        }

        void test_hypothesis_at_current_layer() {
            int ncan = num_currently_candidates();
            if(ncan == 0) {
                // If no candidate hypotheses, just get out
                return;
            }

            double upper_bound;
            std::tuple<int, std::vector<double>> below_dat;
            for(auto r = 1; r <= ncan; r++) {
                below_dat = num_nodes_with_null_prob_below_alphar(r);

                if(std::get<0>(below_dat) < r) {
                    // If we can no longer make any more rejections, no more candidates
                    if(r == 1) {
                        for(ListDigraph::NodeIt n(_gr); n != INVALID; ++n) {
                            remove_candidate_hypothesis(n);
                        }
                    } else {
                        below_dat = num_nodes_with_null_prob_below_alphar(r-1);
                        std::sort(std::get<1>(below_dat).begin(), std::get<1>(below_dat).end());
                        upper_bound = std::get<1>(below_dat)[std::get<0>(below_dat)];

                        for(ListDigraph::NodeIt n(_gr); n != INVALID; ++n) {
                            if(is_candidate_hypothesis(n)) {
                                // if probability less than bound, reject
                                if(_prob_theta_equals_zero(label(n)-1) < upper_bound) {
                                    reject(n);
                                } else {
                                    // otherwise, don't reject and simply remove from
                                    //  list of candidate hypotheses
                                    remove_candidate_hypothesis(n);
                                }
                            }
                        }
                        // Move to the next layer (if available)
                        // For all nodes in new layer whose parents are not all rejected, do not make candidate hypotheses
                        traverse();
                    }
                    return;
                }
            }
            // if we made it this far, it's because we can reject all hypotheses, so do that
            for(ListDigraph::NodeIt n(_gr); n != INVALID; ++n) {
                if(is_candidate_hypothesis(n)) {
                    reject(n);
                }
            }
            traverse();
        }      
};

#endif