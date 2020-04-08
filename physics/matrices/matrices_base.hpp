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

    /**
     * Matrices Network Base Class
     *
     * MatricesBase can be initialized directly, which is different from CellBase.\n
     * This feature is ready for further python interface that supported network topology customization.
     * @note The workflow: add_cell() -> topology() -> add_link() -> assign_inputs() -> update_iterate() -> fetch_outputs().
     */
    class MatricesBase {
    public:
        MatricesBase(size_t in_ports, size_t out_ports): link_graph(), cell_in_id(0), cell_out_id(0),
            inputs_global(), outputs_global(), in_ports_size(in_ports), out_ports_size(out_ports), is_topology_ready(false) {
            assert(in_ports > 0 && out_ports > 0);
        }
        ~MatricesBase() = default;

        /**
         * Add cell to network matrices
         * @tparam T: type should be derived from CellBase
         * @param config_map: configuration map of cell parameters in a form of std::map<std::string, std::any>
         * @param name: cell name, using utils::naming_identifier() to generate is also okay
         * @return current given index
         * @note not need to add in/out terminator, excuting topology() will auto generate
         * @sa cells, utils::naming_identifier()
         */
        template <typename T>
        size_t add_cell(const std::map<std::string, std::any> &config_map, const std::string &name = std::string()) {
            static_assert(std::is_base_of<cells::CellBase, T>::value);
            is_topology_ready = false;
            size_t id = cell_ptrs.size();
            auto gen_name(name);
            if(gen_name.empty()) { gen_name = utils::naming_identifier(cells::cell_type_name<T>(), id); }
            auto gen_cell = std::make_unique<T>(gen_name, config_map);
            cell_ptrs.push_back(std::move(gen_cell));
            return id;
        }

        /**
         * Add link to network matrices
         * @param out_from: output cell index of cell_ptrs (send)
         * @param out_port_index: output port index of given cell, port start as 0
         * @param in_to: input cell index of cell_ptrs (received)
         * @param in_port_index: input port index of given cell, port start as 0
         * @param weight: normally set as 1.0, it is a further feature for routing with consideration of line-length loss
         * @return whether manage or not
         */
        bool add_link(size_t out_from, size_t out_port_index, size_t in_to, size_t in_port_index, double weight = 1.0) {
            if(!check_topology_ready()) return false;
            std::function<void(void)> msg_func = [](){ std::cout << YELLOW << "[WARN] add_link wrong parameter." << RESET << std::endl; };
            if(!(out_from < cell_ptrs.size() && in_to < cell_ptrs.size())) { msg_func(); return false; }
            //auto &out_ref = *(cell_ptrs[out_from].get());
            //auto &in_ref = *(cell_ptrs[in_to].get());
            if(out_port_index >= cell_ptrs[out_from]->out_size()) { msg_func(); return false; }
            if(in_port_index >= cell_ptrs[in_to]->in_size()) { msg_func(); return false; }
            link_graph.insert(cell_in_id[in_to] + in_port_index, cell_out_id[out_from] + out_port_index) = weight;
            return true;
        }

        /**
         * After adding cell, excute confirm_cell().\n
         * This function will calculate the in/out size and allocate for the adjacent matrix for networks.
         */
        void confirm_cell() {
            size_t inport_total_num(0), outport_total_num(0);
            for(const auto &ptr: cell_ptrs) {
                cell_in_id.push_back(inport_total_num);
                cell_out_id.push_back(outport_total_num);
                inport_total_num += ptr->in_size();
                outport_total_num += ptr->out_size();
            }
            link_graph = SparseMatrix<double>(inport_total_num, outport_total_num);
            inputs_global = VectorXcd::Zero(inport_total_num);
            outputs_global = VectorXcd::Zero(outport_total_num);
            is_topology_ready = true;
        }

        /**
         * This function just visualizes the link of networks without any variables operation.
         */
        void confirm_link() {
            if(!check_topology_ready()) return;
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
            //std::cout << "---- Adjacent Matrix ----\n" << link_graph << std::endl;
        }

        /**
         * Print the information of cells.
         * @param cells_id: the id that you acquired from the return value of add_cell()
         */
        void print(const std::vector<size_t> &cells_id = std::vector<size_t>()) const {
            if(cells_id.empty()) {
                for(const auto &itm: cell_ptrs) { itm->print(); }
            } else {
                for(const auto &id: cells_id) {
                    if(id < cell_ptrs.size()) { cell_ptrs[id]->print(); }
                }
            }
        }

        /**
         * Once you finish adding cells, excute topology().\n
         * This function will automatically add in/out port for you, then you could add_link().
         * @note for override user, you should contain add_cell(), confirm_cell() and add_link() orderly.
         * @sa add_cell(), confirm_cell() and add_link()
         */
        virtual void topology() {
            // terminator cell
            for(size_t index_in = 0; index_in < in_ports_size; index_in++) {
                add_cell<cells::InputPort>({}, utils::naming_identifier("InputPort" + std::to_string(index_in)));
            }
            for(size_t index_out = 0; index_out < out_ports_size; index_out++) {
                add_cell<cells::OutputPort>({}, utils::naming_identifier("OutputPort" + std::to_string(index_out)));
            }
            confirm_cell();
        };

        /**
         * Get the size of matrices input ports.
         * @return
         */
        size_t inputs_size() const { return in_ports_size; }

        /**
         * Get the size of matrices output ports.
         * @return
         */
        size_t outputs_size() const { return out_ports_size; }

        /**
         * Get the size of global input ports, containing all the cells' input terminals, other than matrices input ports.
         * @return
         */
        size_t inputs_global_size() const { return inputs_global.size(); }

        /**
         * Get the size of global output ports, containing all the cells' output terminals, other than matrices output ports.
         * @return
         */
        size_t outputs_global_size() const { return outputs_global.size(); }

        /**
         * Get the size of each cell input ports
         * @param cell_index: the cell index that you acquired from return value of add_cell()
         * @return
         */
        size_t input_size(size_t cell_index) const { assert(cell_index < cell_ptrs.size()); return cell_ptrs[cell_index]->in_size(); }

        /**
         * Get the size of each cell output ports
         * @param cell_index: the cell index that you acquired from return value of add_cell()
         * @return
         */
        size_t output_size(size_t cell_index) const { assert(cell_index < cell_ptrs.size()); return cell_ptrs[cell_index]->out_size();}

        /**
         * Get the E-field vector of matrices global inputs, containing each cells' terminals
         * @return
         */
        VectorXcd get_inputs_global() const { return inputs_global; }

        /**
         * Get the E-field vector of matrices global outputs, containing each cells' terminals
         * @return
         */
        VectorXcd get_outputs_global() const { return outputs_global; }

        std::string get_cell_name(size_t cell_index) const { assert(cell_index < cell_ptrs.size()); return cell_ptrs[cell_index]->get_name(); }
        std::map<std::string, std::any> get_cell_vars(size_t cell_index) const { assert(cell_index < cell_ptrs.size()); return cell_ptrs[cell_index]->get_vars(); }

        /**
         * Set value for global inputs of each cells
         * @param v: vector of E-field components
         * @return whether succeeds or not, warns that size should be consistent with inputs_global_size.
         * @note DO NOT RECOMMEND USE THIS FUNCTION, this was used to debug
         */
        bool init_inputs_global(const VectorXcd &v) {
            if(v.size() != inputs_global.size()) { return false; }
            inputs_global = v;
            load_inputs_global();
            return true;
        }

        /**
         * Assign the initial value at inputs ports (not global)
         * @param v: vector to set
         * @return whether succeeds or not
         */
        bool assign_inputs(const VectorXcd &v) {
            if(v.size() != in_ports_size) { return false; }
            size_t _id = inputs_global.size() - in_ports_size - out_ports_size;
            for(size_t idx = 0; idx < in_ports_size; idx++) {
                inputs_global.segment(_id + idx, 1) = v.segment(idx, 1);
            }
            load_inputs_global();
            return true;
        }

        /**
         * Fetch the result from outputs ports (not global)
         * @return vector to fetch
         */
        VectorXcd fetch_outputs() {
            VectorXcd _outputs = VectorXcd::Zero(out_ports_size);
            size_t _id = cell_ptrs.size() - out_ports_size;
            for(size_t idx = 0; idx < out_ports_size; idx++) {
                _outputs.segment(idx, 1) = cell_ptrs[_id + idx]->get_E();
            }
            return _outputs;
        }

        /**
         * For arithmetical convenience, we use a local copy of inputs other than that saved in cells. \n
         * So when updated something in inputs_global or outputs_global, do not forget to call load_inputs global
         * in order to update parameters into cells.
         */
        void load_inputs_global() {
            if(!check_topology_ready()) return;
            for(size_t idx = 0; idx < cell_ptrs.size(); idx++) {
                // update inputs buffer
                size_t _in_size = cell_ptrs[idx]->in_size();
                size_t _in_id = id_from_cell_port(idx, 0, false);
                cell_ptrs[idx]->set_E(inputs_global.segment(_in_id, _in_size));
            }
        }

        /**
         * Update each intra-cells system and transfer to inter-cells without updating cell buffers. \n
         * For iteration, use update_iterate() instead.
         * @param wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         * @sa update_iterate(), load_inputs_global()
         */
        void update(double wavelength) {
            for(size_t idx = 0; idx < cell_ptrs.size(); idx++) {
                // update intra-cells system
                cell_ptrs[idx]->update(wavelength);
                // load outputs
                size_t _out_size = cell_ptrs[idx]->out_size();
                size_t _out_id = id_from_cell_port(idx, 0, true);
                outputs_global.segment(_out_id, _out_size) = cell_ptrs[idx]->get_E();
            }
            // update inter-cells system
            inputs_global = link_graph * outputs_global;
        }

        /**
         * Update function for iteration
         * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         */
        void update_iterate(double wavelength) {
            update(wavelength);
            load_inputs_global();
        }

    protected:
        /**
         * Get the cell and port index from id
         * @param port_id_of_graph: port id with respect to link_graph
         * @param false_as_in_true_as_out: false means inputs ports while true means outputs ports
         * @return std::tuple<size_t cell_index, size_t port_index>
         * @note above C++17, use auto [cell_index, port_index] = cell_port_from_id(...) for convenience
         * @sa id_from_cell_port
         */
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

        /**
         * Get the id from cell and port index
         * @param cell_index: cell index of cell_ptrs that once you acquired from return value of add_cell()
         * @param port_index: port index of each cell, starts from 0
         * @param false_as_in_true_as_out: false means inputs ports while true means outputs ports
         * @return id of link_graph
         * @sa cell_port_from_id
         */
        size_t id_from_cell_port(size_t cell_index, size_t port_index, bool false_as_in_true_as_out = false) const {
            assert(cell_index < cell_ptrs.size());
            if(!false_as_in_true_as_out) {
                return cell_in_id[cell_index] + port_index;
            } else {
                return cell_out_id[cell_index] + port_index;
            }
        }

        /**
         * A safety check to ensure you have done in correct workflow
         * @return
         */
        bool check_topology_ready() {
            if(is_topology_ready) { return true; }
            else {
                std::cout << YELLOW << "\n[WARN] Matrices/topology is not ready! Please check:"
                    << "\n\t" << "1. if you use MatricesBase for freestyle design, make sure excute "
                        << "add_cell() and " << CYAN << "topology()" << YELLOW << " instead of confirm_cell();"
                    << "\n\t" << "2. if you use Derived class, make sure you override the "
                        << CYAN << "topoloy()" << YELLOW << "." << RESET << std::endl;
                std::cout << MAGENTA << "[GUIDE] Workflow for topology is: "
                        << "\n\t 1. assign in/out port number at Constructor Function"
                        << "\n\t 2. add_cell() for device without consider in/out ports"
                        << "\n\t 3. topology() for MatricesBase and confirm_cell() for Derived Class overrided in topology()"
                        << RESET << std::endl;
                return false;
            }
        }


    protected:
        bool is_topology_ready;
        SparseMatrix<double> link_graph;    //!< Graph that indicates the links between cells
        size_t in_ports_size;       //!< Interface inputs size of matrices
        size_t out_ports_size;      //!< Interface outputs size of matrices
        VectorXcd inputs_global;    //!< Totally global inputs of each cells
        VectorXcd outputs_global;   //!< Totally global outputs of each cells
        std::vector<size_t> cell_in_id;    //!< ID List that indicates each head input ports of cells by index
        std::vector<size_t> cell_out_id;   //!< ID List that indicates each head output ports of cells by index
        std::vector<std::unique_ptr<cells::CellBase>> cell_ptrs;    //!< Cell Pointers for switches matrices
    };
}

#endif //OMNICLO_MATRICES_BASE_HPP
