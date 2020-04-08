/**
 * @file piloss.hpp
 * @brief PILOSS matrices in silicon photonic
 * @author LIU-Yinyi
 * @version 1.0.2
 * @date 2020-04-08
 */

#ifndef OMNICLO_PILOSS_HPP
#define OMNICLO_PILOSS_HPP

// Custom Library Header
#include "matrices_base.hpp"

namespace matrices {

    /**
     * Butterfly Network Matrices
     */
    class PILOSS : public MatricesBase {
    public:
        /**
         * PILOSS Constructor
         * @param stages: terminologically 2-ary N-fly topology, here N is the stage number, 2 is the radix
         * @note inputs and outputs are both \f$ 2n \f$ ports
         */
        explicit PILOSS(size_t stages): MatricesBase(2 * stages, 2 * stages) {}

        ~PILOSS() override = default;

        /**
         * Topology design
         * @note After initializng the constructor, please excute topology() for auto-generating crossbar networks
         */
        void topology() override {
            cell_ptrs.clear();

            // device cell
            for (size_t index_out = 0; index_out < out_ports_size / 2; index_out++) {
                for (size_t index_in = 0; index_in < in_ports_size / 2; index_in++) {
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
            size_t _stage = in_ports_size / 2;
            for (size_t index_out = 0; index_out < _stage - 1; index_out++) {
                for (size_t index_in = 0; index_in < _stage; index_in++) {
                    size_t _current = index_out * _stage + index_in;
                    size_t _next = (index_out + 1) * _stage + index_in;
                    size_t _next0, _next1, _p0, _p1;
                    if(index_in == 0) { _next0 = _next; _p0 = 0; } else { _next0 = _next - 1; _p0 = 1; }
                    if(index_in == (_stage - 1)) { _next1 = _next; _p1 = 1; } else { _next1 = _next + 1; _p1 = 0; }
                    add_link(_current, 0, _next0, _p0);
                    add_link(_current, 1, _next1, _p1);
                }
            }
            size_t _term_in_head = (in_ports_size * out_ports_size / 4);
            size_t _term_out_head = _term_in_head + in_ports_size;
            size_t _last_dev_head = _term_in_head - _stage;
            for (size_t index_in = 0; index_in < _stage; index_in++) {
                add_link(_term_in_head + index_in * 2, 0, index_in, 0);
                add_link(_term_in_head + index_in * 2 + 1, 0, index_in, 1);
                add_link(_last_dev_head + index_in, 0, _term_out_head + 2 * index_in, 0);
                add_link(_last_dev_head + index_in, 1, _term_out_head + 2 * index_in + 1, 0);
            }
            confirm_link();
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

#endif //OMNICLO_PILOSS_HPP
