/**
 * @file cells_base.hpp
 * @brief Base of physical cells in silicon photonic
 * @author LIU-Yinyi
 * @version 1.0.1
 * @date 2020-04-06
 */

#ifndef OMNICLO_CELLS_BASE_HPP
#define OMNICLO_CELLS_BASE_HPP

// C++ Built-in Header
#include <utility>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <any>
#include <utility>

// Custom Library Header
#include "display/console.hpp"
#include "effects/effects.hpp"

namespace utils {

    using Eigen::VectorXcd;

    /**
     * Calculate the Square-Root of Power from Electric field
     * @param E: Eletric field
     * @return square-root power of E-field
     */
    VectorXcd cal_P_sqrt(const Eigen::VectorXcd &E) { return abs(E.array()); }

    /**
     * Calculate the Value of Power from Electric field
     * @param E: Eletric field
     * @return power of E-field
     */
    VectorXcd cal_P(const Eigen::VectorXcd &E) { return pow(abs(E.array()), 2); }

    /**
     * Acquire the format string of VectorXcd
     * @param v: VectorXcd would like to print as string
     * @return string of format VectorXcd
     */
    std::string print_vectorXcd(const VectorXcd &v) {
        std::stringstream ss{};
        ss << "[";
        if(v.size() > 0) {
            for (size_t idx = 0; idx < v.size() - 1; idx++) {
                ss << v[idx].real();
                if (v[idx].imag() >= 0.0) { ss << " + " << v[idx].imag() << "i, "; }
                else { ss << " - " << -v[idx].imag() << "i, "; }
            }

            ss << v[v.size() - 1].real();
            if (v[v.size() - 1].imag() >= 0.0) { ss << " + " << v[v.size() - 1].imag() << "i"; }
            else { ss << " - " << -v[v.size() - 1].imag() << "i"; }
        }
        ss << "]";
        return ss.str();
    }

    /**
     * Acquire the format string of vector double
     * @param dv: std::vector<double> would like to print as string
     * @return string of format vector<double>
     */
    std::string print_vector_double(const std::vector<double> &dv) {
        std::stringstream ss{};
        ss << "[";
        const auto _separator = ", ";
        const auto* _sep = "";
        for (const auto &itm : dv) {
            ss << _sep << itm;
            _sep  = _separator;
        }
        ss << "]";
        return ss.str();
    }

    /**
     * Acquire the format string of complex double
     * @param cd: std::complex<double> would like to print as string
     * @return string of format complex<double>
     */
    std::string print_complex_double(const std::complex<double> &cd) {
        std::stringstream ss{};
        ss << cd.real();
        if(cd.imag() >= 0.0) { ss << " + " << cd.imag() << "i, "; }
        else { ss << " - " << -cd.imag() << "i"; }
        return ss.str();
    }

    std::string print_vector_complex_double(const std::vector<std::complex<double>> &v) {
        std::stringstream ss{};
        ss << "[";
        if(!v.empty()) {
            for (size_t idx = 0; idx < v.size() - 1; idx++) {
                ss << v[idx].real();
                if (v[idx].imag() >= 0.0) { ss << " + " << v[idx].imag() << "i, "; }
                else { ss << " - " << -v[idx].imag() << "i, "; }
            }

            ss << v[v.size() - 1].real();
            if (v[v.size() - 1].imag() >= 0.0) { ss << " + " << v[v.size() - 1].imag() << "i"; }
            else { ss << " - " << -v[v.size() - 1].imag() << "i"; }
        }
        ss << "]";
        return ss.str();
    }

    /**
     * Acquire the format string of value type in std::any
     * @param var: value type of std::any would like to print as string
     * @return string of format any value type or nil
     * @note supported type: size_t, int, double, std::complex<double>, std::vector<double>, std::vector<std::complex<double>>, Eigen::VectorXcd
     */
    std::string print_any_value(std::any var) {
        std::stringstream ss{};
        try { auto x = std::any_cast<size_t>(var); ss << x; return ss.str(); } catch(std::exception &ex) {}
        try { auto x = std::any_cast<int>(var); ss << x; return ss.str(); } catch(std::exception &ex) {}
        try { auto x = std::any_cast<double>(var); ss << x; return ss.str(); } catch(std::exception &ex) {}
        try { auto x = std::any_cast<std::complex<double>>(var); ss << print_complex_double(x); return ss.str(); } catch(std::exception &ex) {}
        try { auto x = std::any_cast<std::vector<std::complex<double>>>(var); ss << print_vector_complex_double(x); return ss.str(); } catch(std::exception &ex) {}
        try { auto x = std::any_cast<Eigen::VectorXcd>(var); ss << print_vectorXcd(x); return ss.str(); } catch(std::exception &ex) {}
        try { auto x = std::any_cast<std::vector<double>>(var); ss << print_vector_double(x); return ss.str(); } catch(std::exception &ex) {}
        std::cout << YELLOW << "\n[WARN] utils::print_any_value(std::any) was past an undefined value type." << RESET << std::endl;
        return "";
    }

    /**
     * Check if key exist in map
     * @param m: query map of {string: any}
     * @param key: query key in string
     * @return whether key exists in map as boolean
     */
    bool map_key_exist(const std::map<std::string, std::any> &m, const std::string &key) {
        return m.find(key) != m.end();
    }

    /**
     * Setup default value in map if key does not exist
     * @param m: map that tests and sets
     * @param key: query key in string
     * @param default_val: default value to set as std::any
     * @return the value with respect to key in map
     */
    std::any map_setup_default(std::map<std::string, std::any> &m, const std::string &key, const std::any &default_val) {
        if(!map_key_exist(m, key)) { m[key] = default_val; }
        return m[key];
    }

    /**
     * Naming Identifier utility
     * @param prefix: normally letter property name
     * @param index: normmally number ID
     * @return combined name in string
     */
    std::string naming_identifier(const std::string &prefix, int index = -1) {
        std::stringstream ss;
        ss << prefix;
        if(index >= 0) { ss << index; }
        return ss.str();
    };
}

namespace cells {

    using Eigen::VectorXcd;

    /**
     * Control function for switching on-off method
     * @sa CellBase::control()
     */
    typedef std::vector<std::any>(*CtrlFunc)(std::map<std::string, std::any>);

    /**
     * Base of cell that needs derived realization.
     * @note Derived class should realize update(), control() for development
     * @note after version 1.0.0+, print(), optimize() won't be overrided in derived class thanks to universal map
     */
    class CellBase {
    public:
        /**
         * Constructor with size of in/out port parameters
         * @param in_size: vector size of in-port
         * @param out_size: vector size of out-port
         * @param name: name of cell in std::string, use utils::naming_identifier()
         * @sa in_size(), out_size(), utils::naming_identifier()
         */
        explicit CellBase(size_t in_size, size_t out_size, std::string name): is_updated(true),
            ctrl_func(nullptr), device_vars{}, device_name(std::move(name)) {
            E_in = VectorXcd::Zero(in_size);
            E_out = VectorXcd::Zero(out_size);
        }
        virtual ~CellBase() = default;

        /**
         * Print the model properties of std::map<std::string, std::map> device_vars.
         * @param query_keys: vector string that queries from cell device variables
         */
        void print(const std::vector<std::string> &query_keys = std::vector<std::string>()) {
            std::cout << "/---------- Cell Parameters ----------\\" << std::endl;
            std::cout << "| \tDeviceName = " << device_name << std::endl;
            std::cout << "| \tIn-Ports = " << E_in.size() << ", Out-Ports = " << E_out.size() << std::endl;
            if(query_keys.empty()) {
                for(const auto &itm: device_vars) {
                    std::cout << "| \t" << itm.first << " = " << utils::print_any_value(itm.second) << std::endl;
                }
            } else {
                for(const auto &query_key: query_keys) {
                    if(utils::map_key_exist(device_vars, query_key)) {
                        std::cout << "| \t" << query_key << " = " << utils::print_any_value(device_vars[query_key]) << std::endl;
                    } else {
                        std::cout << YELLOW << "\n[WARN] cells::print() query key string <" << CYAN << query_key
                        << YELLOW << "> does not exists in map of cell device_vars." << std::endl;
                    }
                }
            }
            std::cout << "\\-------------------------------------/" << std::endl;
        }

        /**
         * Virtual function that needs realization in derived class.
         * *update* indicates a mapping relation between input and output.
         * In this function you should update output from input, even or use static memory for recursive method.
         *
         * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         * @note Sign *update* as mapping operator \f$ M_U \f$, *inputs* and *outputs* as \f$ V_I \f$ and \f$ V_O\f$, respectively. \n
         * Formula can be proposed: \f$ V_O = M_U V_I \f$
         * @sa update_with_flag()
         */
        virtual void update(double wavelength) = 0;

        /**
         * This is a wrapper for *update* function with alternating *is_update* flag when finished execution.
         * @sa update()
         */
        inline void update_with_flag(double wavelength) { update(wavelength); is_updated = true; }

        /**
         * Function pointer which indicates the control method of switch cell.
         * @param func: function pointer of on-off control
         * @note CtrlFunc: return control targets, normally for delta_beta [ \f$ \Delta \beta \f$ ]
         * @sa control()
         */
        void reg_ctrl_func(CtrlFunc func) { ctrl_func = func; }

        /**
         * Alternate the properties by changing control variables in order to achieve switch control.
         * @param ctrl_vars: control variables, for example, voltage, temperature, magnetic that affects optic
         * @sa reg_ctrl_func()
         */
        virtual void control(std::map<std::string, std::any> ctrl_vars) = 0;

        /**
         * Update the parameters of cell properties that availbale for optimization.
         * @param optim_map: former key for available parameter names, latter key for value
         * @note recommend use: size_t, int, double, std::complex<double>, std::vector<double>, Eigen::VectorXcd
         */
        void optimize(const std::map<std::string, std::any> &optim_map) {
            for(const auto &item : optim_map) {
                if(utils::map_key_exist(device_vars, item.first)) { device_vars[item.first] = item.second; }
                else { std::cout << YELLOW << "[WARN] {cells->" << device_name << "/optimize} undefined optimization parameter: "
                    << CYAN << item.first << RESET << std::endl; }
            }
        }

        /**
         * Set the values of \f$ E_{in} \f$
         * @param E: Electric field of inputs as a column vector
         */
        void set_E(const VectorXcd &E) { assert(E.size() == E_in.size()); E_in = E; is_updated = false; }

        /**
         * Get the values of \f$ E_{out} \f$
         * @return Electric field of outputs as a column vector
         */
        VectorXcd get_E() const { return E_out; }

        /**
         * Get the size of input vector, which indicates the number of input ports
         * @return number of input ports
         */
        size_t in_size() const { return E_in.size(); }

        /**
         * Get the size of output vector, which indicates the number of output ports
         * @return number of output ports
         */
        size_t out_size() const { return E_out.size(); }

        /**
         * Get the device name of cell saved as device_name
         * @return device name
         */
        std::string get_name() const { return device_name; }

        /**
         * Get the device variables of cell saved as device_vars
         * @return device variables/parameters
         */
        std::map<std::string, std::any> get_vars() const { return device_vars; }

        /**
         * For linear algebra, transfer matrix shape can be defined as \f$ N_{out} \times N_{in} \f$
         * @return std::tuple<size_t (RowNum), size_t (ColNum)>
         * @note for C++17 or above use: auto [rows, cols] = transfer_matrix_shape(); \n
         * for C++17 below use: std::tie(rows, cols) = transfer_matrix_shape();
         */
        std::tuple<size_t, size_t> transfer_matrix_shape() const { return std::make_tuple(out_size(), in_size()); }

        void set_E_in(const VectorXcd &E) { E_in = E; }
        void set_E_out(const VectorXcd &E) { E_out = E; }
        VectorXcd get_E_in() const { return E_in; }
        VectorXcd get_E_out() const { return E_out; }

    protected:
        bool is_updated;    //!< the flag indicates whether \f$ E_{out} \f$ updated or not after refreshing \f$ E_{in} \f$
        VectorXcd E_in;     //!< Electric field of inputs as column vector
        VectorXcd E_out;     //!< Electric field of outputs as column vector
        CtrlFunc ctrl_func; //!< Function pointer that controls the properties of cell
        std::string device_name;    //!< Name of cell device
        std::map<std::string, std::any> device_vars;    //!< Map of device cell's variables or parameters
    };

}

#endif //OMNICLO_CELLS_BASE_HPP
