// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
 *   Created by mbenlioglu on Aug 31, 2020.
 */

#pragma once

#include <Teuchos_ParameterEntryValidator.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>

#ifndef HAVE_ZOLTAN2_SARMA

namespace Zoltan2 {
    template<typename Adapter>
    class AlgSarma : public Algorithm<Adapter> {
    public:
        using base_adapter_t = typename Adapter::base_adapter_t;

        AlgSarma(const RCP<const Environment> &/* env */,
                const RCP<const Comm<int> > &/* problemComm */,
                const RCP<const base_adapter_t> &/* adapter */) {
            throw std::runtime_error("SARMA requested but not enabled into Zoltan2\n"
                                     "Please set CMake flag TPL_ENABLE_SARMA:BOOL=ON");
        }
    };
}

#else

#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <sarma/sarma.hpp>

namespace Zoltan2 {
    template<typename Adapter>
    class AlgSarma : public Algorithm<Adapter> {
    public:
        using base_adapter_t = typename Adapter::base_adapter_t;
        using lno_t = typename Adapter::lno_t;
        using gno_t = typename Adapter::gno_t;
        using offset_t = typename Adapter::offset_t;
        using scalar_t = typename Adapter::scalar_t;
        using part_t = typename Adapter::part_t;
        using user_t = typename Adapter::user_t;
        using userCoord_t = typename Adapter::userCoord_t;

    private:
        const offset_t *offsets;
        const gno_t *colids;
        const scalar_t *vals;
        size_t offsize;
        size_t nnz;

    public:
        AlgSarma(const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
                 const RCP<const IdentifierAdapter <user_t> > &adapter)
                : env(env), comm(comm), adapter(adapter),
                  algs(sarma::get_algorithm_map<Ordinal, Value>()), orders(sarma::get_order_map()) {
            throw std::runtime_error("Cannot build SARMA from IdentifierAdapter");
        }

        AlgSarma(const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
                 const RCP<const VectorAdapter <user_t> > &adapter)
                : env(env), comm(comm), adapter(adapter),
                  algs(sarma::get_algorithm_map<Ordinal, Value>()), orders(sarma::get_order_map()) {
            throw std::runtime_error("Cannot build SARMA from VectorAdapter");

        }

        AlgSarma(const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
                 const RCP<const GraphAdapter <user_t, userCoord_t> > &adapter)
                : env(env), comm(comm), adapter(adapter),
                  algs(sarma::get_algorithm_map<Ordinal, Value>()), orders(sarma::get_order_map()) {
            throw std::runtime_error("Cannot build SARMA from GraphAdapter");
        }

        AlgSarma(const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
                 const RCP<const MeshAdapter <user_t> > &adapter)
                : env(env), comm(comm), adapter(adapter),
                  algs(sarma::get_algorithm_map<Ordinal, Value>()), orders(sarma::get_order_map()) {

            throw std::runtime_error("Cannot build SARMA from MeshAdapter");
        }

        AlgSarma(const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
                 const RCP<const MatrixAdapter <user_t, userCoord_t> > &adapter)
                : env(env), comm(comm), adapter(adapter),
                  offsize(1 + adapter->getLocalNumRows()), nnz(adapter->getLocalNumEntries()),
                  algs(sarma::get_algorithm_map<Ordinal, Value>()), orders(sarma::get_order_map()) {
            getParams();
            if (adapter->CRSViewAvailable()) {
                offsets = new offset_t[offsize]();
                colids = new gno_t[nnz]();
                vals = new scalar_t[nnz]();
                adapter->getCRSView(offsets, colids, vals);
            } else
                throw std::runtime_error("NO CRS Available");
        }

        static void getValidParameters(ParameterList &pl) {
            RCP<const Teuchos::ParameterEntryValidator> validator;
            pl.set("sarma_parameters", ParameterList(), "Sarma options", validator);
        }

        void partition(const RCP <PartitioningSolution<Adapter>> &solution);
    private:
        const RCP<const Environment> env;
        const RCP<const Comm<int> > comm;
        const RCP<const base_adapter_t> adapter;
        RCP<const CoordinateModel <base_adapter_t> > model;

        struct Parameters {
            std::string alg = "pal";
            sarma::Order order_type = sarma::Order::NAT;
            std::string order_str = "nat";
            Ordinal row_parts = 8, col_parts = 0;
            Value max_load = 0;
            int seed = 2147483647;
            double sparsify = 1.0;
            bool triangular = false, use_data = false;
        };

        std::map <std::string, std::pair<std::function < std::pair < std::vector < Ordinal>, std::vector<Ordinal>>(
        const sarma::Matrix <Ordinal, Value> &, const Ordinal, const Ordinal, const Value, const int)>, bool>>
        algs;
        std::map <std::string, sarma::Order> orders;
        Parameters config;

//        void buildModel(modelFlag_t &flags);
        void getParams();

    };

    template<typename Adapter>
    void AlgSarma<Adapter>::getParams() {
        const Teuchos::ParameterList &pl = env->getParameters();
        const Teuchos::ParameterList &sparams = pl.sublist("sarma_parameters");

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
        pe = sparams.getEntryPtr("row_parts"); // row_cut
        if (pe)
            config.row_parts = pe->getValue(&config.row_parts); 
        pe = sparams.getEntryPtr("col_parts"); //col_cut
        if (pe)
            config.col_parts = pe->getValue(&config.col_parts);  
        pe = sparams.getEntryPtr("max_load"); // max_load
        if (pe)
            config.max_load = pe->getValue(&config.max_load);
        pe = sparams.getEntryPtr("sparsify");
        if (pe)
            config.sparsify = pe->getValue(&config.sparsify);
        pe = sparams.getEntryPtr("triangular");
        if (pe)
            config.triangular = pe->getValue(&config.triangular);
        pe = sparams.getEntryPtr("use_data");
        if (pe)
            config.use_data = pe->getValue(&config.use_data);
        pe = sparams.getEntryPtr("seed");
        if (pe)
            config.seed = pe->getValue(&config.seed);
    }

    template<typename Adapter>
    void AlgSarma<Adapter>::partition(const RCP <PartitioningSolution<Adapter>> &solution) {

        std::ofstream null;
        null.setstate(std::ios_base::badbit);
        auto M = std::make_shared<sarma::Matrix<Ordinal, Value> >(std::move(std::vector<Ordinal>(offsets, offsets + offsize)),
                std::move(std::vector<Ordinal>(colids, colids + nnz)), std::move(std::vector<Value>(vals, vals + nnz)),
                1 + *std::max_element(colids, colids + nnz));
        auto parts = sarma::Run<Ordinal, Value>(algs.at(config.alg).first, null, M, config.order_type, config.row_parts,
                                                config.col_parts, config.max_load, config.triangular, false, config.sparsify,
                                                algs.at(config.alg).second, config.use_data, config.seed);

        unsigned result_size = parts.first.size() + parts.second.size();
        auto partl = ArrayRCP<int>(new int[result_size], 0, result_size, true);
        for (size_t i = 0; i < parts.first.size(); ++i) partl[i] = parts.first[i];
        for (size_t i = parts.first.size(); i < result_size; ++i)
            partl[i] = parts.second[i - parts.first.size()];

        solution->setParts(partl);
    }
}

#endif