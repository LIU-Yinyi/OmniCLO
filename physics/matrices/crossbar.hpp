/**
 * @file crossbar.hpp
 * @brief Crossbar matrices in silicon photonic
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-04-04
 */

#ifndef OMNICLO_CROSSBAR_HPP
#define OMNICLO_CROSSBAR_HPP

// Custom Library Header
#include "matrices_base.hpp"

namespace matrices {

    class Crossbar : public MatricesBase {
    public:
        Crossbar(size_t in_ports, size_t out_ports): MatricesBase(),
            in_ports_size(in_ports), out_ports_size(out_ports) {
            assert(in_ports > 0 && out_ports > 0);
        }

        virtual ~Crossbar() = default;

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

    protected:
        size_t in_ports_size;
        size_t out_ports_size;
    };
}

#endif //OMNICLO_CROSSBAR_HPP
