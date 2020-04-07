/**
 * @file matrices_base.hpp
 * @brief Base of physical matrices in silicon photonic
 * @author LIU-Yinyi
 * @version 1.0.1
 * @date 2020-04-06
 */

#ifndef OMNICLO_MATRICES_BASE_HPP
#define OMNICLO_MATRICES_BASE_HPP

// Scientific Library Header
#include <eigen3/Eigen/Sparse>

// Custom Library Header
#include "cells/cells.hpp"

namespace matrices {

    using Eigen::SparseMatrix;
    using Eigen::VectorXcd;


    class MatricesBase {
    public:
        MatricesBase(): link_graph(), cell_in_id(0), cell_out_id(0), inputs(), outputs() {}
        ~MatricesBase() = default;

        template <typename T>
        size_t add_cell(const std::map<std::string, std::any> &config_map, const std::string &name = std::string()) {
            static_assert(std::is_base_of<cells::CellBase, T>::value);
            size_t id = cell_ptrs.size();
            auto gen_name(name);
            if(gen_name.empty()) { gen_name = utils::naming_identifier(cells::cell_type_name<T>(), id); }
            auto gen_cell = std::make_unique<T>(gen_name, config_map);
            cell_ptrs.push_back(std::move(gen_cell));
            return id;
        }

        void add_link(size_t out_from, size_t out_port_index, size_t in_to, size_t in_port_index, double weight = 1.0) {
            assert(out_from < cell_ptrs.size() && in_to < cell_ptrs.size());
            //assert(out_port_index < cell_ptrs[out_from]->out_size());
            //assert(in_port_index < cell_ptrs[in_to]->in_size());
            link_graph.insert(cell_in_id[in_to] + in_port_index, cell_out_id[out_from] + out_port_index) = weight;
        }

        void confirm_cell() {
            size_t inport_total_num(0), outport_total_num(0);
            for(const auto &ptr: cell_ptrs) {
                cell_in_id.push_back(inport_total_num);
                cell_out_id.push_back(outport_total_num);
                inport_total_num += ptr->in_size();
                outport_total_num += ptr->out_size();
            }
            link_graph = SparseMatrix<double>(inport_total_num, outport_total_num);
            inputs = VectorXcd::Zero(inport_total_num);
            outputs = VectorXcd::Zero(outport_total_num);
        }

        void confirm_link() {
            std::cout << "---- Matrices Graph Link Info ----" << std::endl;
            for(size_t outIndex = 0; outIndex < link_graph.outerSize(); ++outIndex) {
                for(SparseMatrix<double>::InnerIterator iter(link_graph, outIndex); iter; ++iter) {
                    double link_length = iter.value();
                    size_t in_port = iter.row();
                    size_t out_port = iter.col();
                    auto [incell_index, inport_index] = cell_port_from_id(in_port, false);
                    auto [outcell_index, outport_index] = cell_port_from_id(out_port, true);
                    std::cout << "$ FROM [Cell: " << outcell_index << " ^ Port: " << outport_index << "] --> "
                        << " TO [Cell: " << incell_index << " ^ Port: " << inport_index << "]  WITH"
                        << " [Length(Weight) = " << link_length << "]." << std::endl;
                }
            }
        }

        void print(const std::vector<size_t> &cells_id = std::vector<size_t>()) const {
            if(cells_id.empty()) {
                for(const auto &itm: cell_ptrs) { itm->print(); }
            } else {
                for(const auto &id: cells_id) {
                    if(id < cell_ptrs.size()) { cell_ptrs[id]->print(); }
                }
            }
        }

        virtual void topology() {};

        size_t inputs_size() const { return inputs.size(); }
        size_t outputs_size() const { return outputs.size(); }
        size_t input_size(size_t cell_index) const { assert(cell_index < cell_ptrs.size()); return cell_ptrs[cell_index]->in_size(); }
        size_t output_size(size_t cell_index) const { assert(cell_index < cell_ptrs.size()); return cell_ptrs[cell_index]->out_size();}

        VectorXcd get_inputs() const { return inputs; }
        VectorXcd get_outputs() const { return outputs; }

        std::string get_cell_name(size_t cell_index) const { assert(cell_index < cell_ptrs.size()); return cell_ptrs[cell_index]->get_name(); }
        std::map<std::string, std::any> get_cell_vars(size_t cell_index) const { assert(cell_index < cell_ptrs.size()); return cell_ptrs[cell_index]->get_vars(); }

        void init_inputs(const VectorXcd &v) {
            assert(v.size() == inputs.size());
            inputs = v;
            load_inputs();
        }

        void update(double wavelength) {
            for(size_t idx = 0; idx < cell_ptrs.size(); idx++) {
                // update intra-cells system
                cell_ptrs[idx]->update(wavelength);
                // load outputs
                size_t _out_size = cell_ptrs[idx]->out_size();
                size_t _out_id = id_from_cell_port(idx, 0, true);
                outputs.segment(_out_id, _out_size) = cell_ptrs[idx]->get_E();
            }
            // update inter-cells system
            inputs = link_graph * outputs;
        }

        void update_iterate(double wavelength) {
            update(wavelength);
            load_inputs();
        }

    protected:
        std::tuple<size_t, size_t> cell_port_from_id(size_t port_id_of_graph, bool false_as_in_true_as_out = false) const {
            size_t idx(0);
            if(!false_as_in_true_as_out) {
                for (idx = 0; idx < cell_in_id.size() - 1; idx++) {
                    if(port_id_of_graph >= cell_in_id[idx] && port_id_of_graph < cell_in_id[idx + 1]) {
                        break;
                    }
                }
                return {idx, port_id_of_graph - cell_in_id[idx]};
            } else {
                for (idx = 0; idx < cell_out_id.size() - 1; idx++) {
                    if(port_id_of_graph >= cell_out_id[idx] && port_id_of_graph < cell_out_id[idx + 1]) {
                        break;
                    }
                }
                return {idx, port_id_of_graph - cell_out_id[idx]};
            }
        }

        size_t id_from_cell_port(size_t cell_index, size_t port_index, bool false_as_in_true_as_out = false) const {
            assert(cell_index < cell_ptrs.size());
            if(!false_as_in_true_as_out) {
                return cell_in_id[cell_index] + port_index;
            } else {
                return cell_out_id[cell_index] + port_index;
            }
        }

        void load_inputs() {
            for(size_t idx = 0; idx < cell_ptrs.size(); idx++) {
                // update inputs buffer
                size_t _in_size = cell_ptrs[idx]->in_size();
                size_t _in_id = id_from_cell_port(idx, 0, false);
                cell_ptrs[idx]->set_E(inputs.segment(_in_id, _in_size));
            }
        }


    protected:
        SparseMatrix<double> link_graph;    //!< Graph that indicates the links between cells
        VectorXcd inputs;
        VectorXcd outputs;
        std::vector<size_t> cell_in_id;    //!< ID List that indicates each head input ports of cells by index
        std::vector<size_t> cell_out_id;   //!< ID List that indicates each head output ports of cells by index
        std::vector<std::unique_ptr<cells::CellBase>> cell_ptrs;    //!< Cell Pointers for switches matrices
    };
}

#endif //OMNICLO_MATRICES_BASE_HPP
