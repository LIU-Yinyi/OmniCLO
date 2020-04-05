/**
 * @file cells_base.hpp
 * @brief Base of physical cells in silicon photonic
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-04-02
 */

#ifndef OMNICLO_CELLS_BASE_HPP
#define OMNICLO_CELLS_BASE_HPP

// C++ Built-in Header
#include <utility>
#include <vector>
#include <tuple>
#include <map>
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
        for (size_t idx = 0; idx < v.size() - 1; idx++) {
            ss << v[idx].real();
            if(v[idx].imag() >= 0.0) { ss << " + " << v[idx].imag() << "i, "; }
            else { ss << " - " << -v[idx].imag() << "i, "; }
        }
        ss << v[v.size() - 1].real();
        if(v[v.size() - 1].imag() >= 0.0) { ss << " + " << v[v.size() - 1].imag() << "i"; }
        else { ss << " - " << -v[v.size() - 1].imag() << "i"; }
        ss << "]";
        return ss.str();
    }
}

namespace cells {

    using Eigen::VectorXcd;

    typedef std::vector<double>(*CtrlFunc)(std::vector<double>);

    /**
     * Base of cell that needs derived realization.
     * @note Derived class should realize print(), update(), control(), optimize() for development.
     */
    class CellBase {
    public:
        /**
         * Constructor with size of in/out port parameters
         * @param in_size: vector size of in-port
         * @param out_size  vector size of out-port
         * @sa in_size(), out_size()
         */
        explicit CellBase(size_t in_size, size_t out_size): is_updated(true), ctrl_func(nullptr)
            { E_in = VectorXcd::Zero(in_size); E_out = VectorXcd::Zero(out_size); }
        virtual ~CellBase() = default;

        /**
         * Print the model properties and the protected or private parameters of the cell.
         */
        virtual void print() = 0;

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
         * @param ctrl_vars: control variables
         * @sa reg_ctrl_func()
         */
        virtual void control(std::vector<double> ctrl_vars) = 0;

        /**
         * Update the parameters of cell properties that availbale for optimization.
         * @param optim_map: former key for availableparameter names, latter key for value
         * @note a complex number or higher dimension should be divided into several double
         */
        virtual void optimize(const std::map<std::string, double> &optim_map) = 0;

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
    };

}

#endif //OMNICLO_CELLS_BASE_HPP
