#pragma once

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <sarma/sarma.hpp>

#include <bitset>
#include <sstream>
#include <string>

namespace Zoltan2 {
    template<typename Adapter>
    class AlgSarma : public Algorithm<Adapter> {
        typedef typename Adapter::lno_t lno_t;       // local ids
        typedef typename Adapter::gno_t gno_t;       // global ids
        typedef typename Adapter::offset_t offset_t;
        typedef typename Adapter::scalar_t scalar_t; // scalars
        typedef typename Adapter::part_t part_t;     // part numbers
        typedef typename Adapter::user_t user_t;
        typedef typename Adapter::userCoord_t userCoord_t;

    private:
        const RCP<const Environment> env;
        const RCP<const Comm<int> > problemComm;
        const RCP<const typename Adapter::base_adapter_t> adapter;

        struct Parameters {
            std::string alg = "pal";
            sarma::Order order_type = sarma::Order::NAT;
            std::string order_str = "nat";
            offset_t p = 8, q = 0;
            scalar_t z = 0;
            int seed = 2147483647;
            double sparsify = 1.0;
            bool triangular = false, serialize = false, use_data = false;
        };

        std::map <std::string, std::pair<std::function < std::pair < std::vector < offset_t>, std::vector<offset_t>>(
        const sarma::Matrix <offset_t, scalar_t> &, const offset_t, const offset_t, const scalar_t, const int)>, bool>>
        algs;
        std::map <std::string, sarma::Order> orders;
        Parameters config;


    public:
        // Constructor
        AlgSarma(
                const RCP<const Environment> &env_,
                const RCP<const Comm<int> > &problemComm_,
                const RCP<const MatrixAdapter <user_t, userCoord_t> > &adapter_
        ) : env(env_), problemComm(problemComm_), adapter(adapter_),
            algs(sarma::get_algorithm_map<offset_t, scalar_t>()), orders(sarma::get_order_map()) {
            const Teuchos::ParameterList &pl = env->getParameters();
            const Teuchos::ParameterList &sparams = pl.sublist("sarma_params");

            // Set Params
            const Teuchos::ParameterEntry *pe = sparams.getEntryPtr("alg");
            if (pe)
                config.alg = pe->getValue(&config.alg);
            pe = sparams.getEntryPtr("order");
            if (pe) {
                std::string s;
                s = pe->getValue(&s);
                if (orders.find(s) == orders.end()) {
                    throw std::logic_error("Wrong order type.");
                }
                config.order_type = orders.at(s);
                config.order_str = s;
            }
            pe = sparams.getEntryPtr("p");
            if (pe)
                config.p = pe->getValue(&config.p);
            pe = sparams.getEntryPtr("q");
            if (pe)
                config.q = pe->getValue(&config.q);
            pe = sparams.getEntryPtr("z");
            if (pe)
                config.z = pe->getValue(&config.z);
            pe = sparams.getEntryPtr("sparsify");
            if (pe)
                config.sparsify = pe->getValue(&config.sparsify);
            pe = sparams.getEntryPtr("triangular");
            if (pe)
                config.triangular = pe->getValue(&config.triangular);
            pe = sparams.getEntryPtr("serialize");
            if (pe)
                config.serialize = pe->getValue(&config.serialize);
            pe = sparams.getEntryPtr("use_data");
            if (pe)
                config.use_data = pe->getValue(&config.use_data);
            pe = sparams.getEntryPtr("seed");
            if (pe)
                config.seed = pe->getValue(&config.seed);
        }

        // Partitioning method
        void partition(const RCP <PartitioningSolution<Adapter>> &solution) {
            env->debug(DETAILED_STATUS, std::string("Entering AlgSarma"));


            offset_t *offsets;
            gno_t *colIds;
            scalar_t *values;
            adapter->getCRSView(offsets, colIds, values);

            auto offset_vec = std::vector<offset_t>(offsets, offsets + adapter->getLocalNumRows());
            auto colids_vec = std::vector<gno_t>(colIds, colIds + adapter->getLocalNumEntries());
            auto values_vec = values ? std::vector<scalar_t>(values, values + adapter->getLocalNumEntries())
                                     : std::vector<scalar_t>();
            auto M = 1 + std::max_element(std::begin(colids_vec), std::end(colids_vec));

            auto mtx = std::make_shared < sarma::Matrix < offset_t, scalar_t> >(std::move(offset_vec),
                    std::move(colids_vec),
                    std::move(values_vec), M);

            const auto res = sarma::Run<offset_t, scalar_t>(algs.at(config.alg).first, std::cout, mtx,
                                                       config.order_type, config.p, config.q, config.z,
                                                       config.triangular,
                                                       config.serialize, config.sparsify, algs.at(config.alg).second,
                                                       config.use_data, config.seed);

            // Assign solutions to 1D part array
            ArrayRCP<part_t> parts = arcp(res.first.data(), 0, res.first.size());

            solution->setParts(parts);

            env->debug(DETAILED_STATUS, std::string("Exiting AlgSarma"));
        }
    };
}
