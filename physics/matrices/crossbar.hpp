/**
 * @file crossbar.hpp
 * @brief Crossbar matrices in silicon photonic
 * @author LIU-Yinyi
 * @version 1.0.2
 * @date 2020-04-08
 */

#ifndef OMNICLO_CROSSBAR_HPP
#define OMNICLO_CROSSBAR_HPP

// Custom Library Header
#include "matrices_base.hpp"

namespace matrices {

    /**
     * Crossbar Network Matrices
     */
    class Crossbar : public MatricesBase {
    public:
        /**
         * Crossbar Constructor
         * @param in_ports: number of input ports
         * @param out_ports: number of output ports
         */
        Crossbar(size_t in_ports, size_t out_ports): MatricesBase(in_ports, out_ports) {}

        ~Crossbar() override = default;

        /**
         * Topology design
         * @note After initializng the constructor, please excute topology() for auto-generating crossbar networks
         */
        void topology() override {
            cell_ptrs.clear();
            // device cell
            for(size_t index_in = 0; index_in < in_ports_size; index_in++) {
                for(size_t index_out = 0; index_out < out_ports_size; index_out++) {
                    add_cell<cells::MicroRingCross>({}, utils::naming_identifier("MicroRingCross"
                    + std::to_string(index_in) + "_" + std::to_string(index_out)));
                }
            }
            // terminator cell
            for(size_t index_in = 0; index_in < in_ports_size; index_in++) {
                add_cell<cells::InputPort>({}, utils::naming_identifier("InputPort" + std::to_string(index_in)));
            }
            for(size_t index_out = 0; index_out < out_ports_size; index_out++) {
                add_cell<cells::OutputPort>({}, utils::naming_identifier("OutputPort" + std::to_string(index_out)));
            }
            confirm_cell();

            // connect link
            size_t term_in_id = in_ports_size * out_ports_size;
            size_t term_out_id = term_in_id + in_ports_size;
            for(size_t index_in = 0; index_in < in_ports_size; index_in++) {
                add_link(term_in_id + index_in, 0, index_in * out_ports_size, 0);
                if(index_in < in_ports_size - 1) {
                    add_link(index_in * out_ports_size, 1, (index_in + 1) * out_ports_size, 1);
                } else {
                    add_link(index_in * out_ports_size, 1, term_out_id, 0);
                }
                for(size_t index_out = 1; index_out < out_ports_size; index_out++) {
                    size_t adjacent_hori = index_out + index_in * out_ports_size;
                    add_link(adjacent_hori - 1, 0, adjacent_hori, 0);
                    if(index_in < in_ports_size - 1) {
                        size_t adjacent_vert = index_out + (index_in + 1) * out_ports_size;
                        add_link(adjacent_hori, 1, adjacent_vert, 1);
                    } else {
                        add_link(adjacent_hori, 1, term_out_id + index_out, 0);
                    }
                 }
            }
            confirm_link();
        }

        /**
         * Control switch to conduct correctly from inputs to outputs
         * @param in_out_flag: tuple of input-index, output-index, on-off-flag
         */
        void switch_ports(const std::vector<std::tuple<size_t, size_t, bool>> &in_out_flag) override {
            if(!check_topology_ready()) return;
            for(const auto &itm : in_out_flag) {
                size_t _in = std::get<0>(itm), _out = std::get<1>(itm);
                bool _on_off = std::get<2>(itm);
                auto _id = _in * out_ports_size + _out;
                //TODO: add on control
                cell_ptrs[_id]->control({});
            }
        }
    };
}

#endif //OMNICLO_CROSSBAR_HPP
