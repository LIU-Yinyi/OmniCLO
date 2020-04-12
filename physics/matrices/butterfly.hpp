/**
 * @file butterfly.hpp
 * @brief Butterfly matrices in silicon photonic
 * @author LIU-Yinyi
 * @version 1.0.2
 * @date 2020-04-08
 */

#ifndef OMNICLO_BUTTERFLY_HPP
#define OMNICLO_BUTTERFLY_HPP

// Custom Library Header
#include "matrices_base.hpp"

namespace matrices {

    /**
     * Butterfly Network Matrices
     */
    class Butterfly : public MatricesBase {
    public:
        /**
         * Butterfly Constructor
         * @param stages: stages of networks
         */
        explicit Butterfly(size_t stages): MatricesBase(pow(2, stages), pow(2, stages)), stages_value(stages) {}

        ~Butterfly() override = default;

        /**
         * Topology design
         * @note After initializng the constructor, please excute topology() for auto-generating butterfly networks
         */
        void topology() override {
            cell_ptrs.clear();
            size_t _rows = in_ports_size / 2;
            size_t _cols = stages_value;

            // device cell
            for (size_t index_col = 0; index_col < _cols; index_col++) {
                for (size_t index_row = 0; index_row < _rows; index_row++) {
                    add_cell<cells::MicroRingCross>({}, utils::naming_identifier("MicroRingCross"
                    + std::to_string(index_row) + "_" + std::to_string(index_col)));
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
            size_t _block = _rows, _thresh = _rows / 2;
            for (size_t index_col = 0; index_col < _cols - 1; index_col++) {
                for (size_t index_row = 0; index_row < _rows; index_row++) {
                    size_t _current = index_col * _rows + index_row;
                    size_t _next = (index_col + 1) * _rows + index_row;
                    size_t _next0, _next1, _p0, _p1;
                    if((index_row % _block) < _thresh) {
                        _next0 = _next; _p0 = 0;
                        _next1 = _next + _thresh; _p1 = 0;
                    } else {
                        _next0 = _next - _thresh; _p0 = 1;
                        _next1 = _next; _p1 = 1;
                    }
                    add_link(_current, 0, _next0, _p0);
                    add_link(_current, 1, _next1, _p1);
                }
                _thresh /= 2;
                _block /= 2;
            }
            size_t _term_in_head = _rows * _cols;
            size_t _term_out_head = _term_in_head + in_ports_size;
            size_t _last_dev_head = _term_in_head - _rows;
            for (size_t index_row = 0; index_row < _rows; index_row++) {
                add_link(_term_in_head + index_row * 2, 0, index_row, 0);
                add_link(_term_in_head + index_row * 2 + 1, 0, index_row, 1);
                add_link(_last_dev_head + index_row, 0, _term_out_head + 2 * index_row, 0);
                add_link(_last_dev_head + index_row, 1, _term_out_head + 2 * index_row + 1, 0);
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

    protected:
        size_t stages_value;
    };
}

#endif //OMNICLO_BUTTERFLY_HPP
