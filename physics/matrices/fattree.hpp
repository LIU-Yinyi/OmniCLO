/**
 * @file fattree.hpp
 * @brief Fat-Tree matrices in silicon photonic
 * @author LIU-Yinyi
 * @version 1.0.2
 * @date 2020-04-08
 */

#ifndef OMNICLO_FATTREE_HPP
#define OMNICLO_FATTREE_HPP


// Custom Library Header
#include "matrices_base.hpp"

namespace matrices {

    /**
     * Fat-Tree Network Matrices
     */
    class FatTree : public MatricesBase {
    public:
        /**
         * Fat-Tree Constructor
         * @param in_ports: number of input ports
         * @param out_ports: number of output ports
         */
        FatTree(size_t in_ports, size_t out_ports): MatricesBase(in_ports, out_ports) {}

        ~FatTree() override = default;

        /**
         * Topology design
         * @note After initializng the constructor, please excute topology() for auto-generating crossbar networks
         */
        void topology() override {

        }

        /**
         * Control switch to conduct correctly from inputs to outputs
         * @param in_out_flag: tuple of input-index, output-index, on-off-flag
         */
        void switch_ports(const std::vector<std::tuple<size_t, size_t, bool>> &in_out_flag) override {
            if(!check_topology_ready()) return;
        }
    };
}


#endif //OMNICLO_FATTREE_HPP
