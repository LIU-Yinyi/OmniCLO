/**
 * @file waveguide.hpp
 * @brief Relization of Waveguide Generalized Cell
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-04-03
 */

#ifndef OMNICLO_WAVEGUIDE_HPP
#define OMNICLO_WAVEGUIDE_HPP

// Custom Library Header
#include "cells_base.hpp"

namespace cells {

    using effects::constants::i;

    /**
     * Waveguide derived from CellBase.
     */
    class Waveguide : public CellBase {
    public:
        Waveguide(): CellBase(1, 1), refractive_index(1.0), delta_n(0.0), propagate_length(0.0) {}

        /**
         * Constructor with Waveguide's properties
         * @param n: refractive index of waveguide
         * @param length: propagation length of waveguide
         */
        Waveguide(std::complex<double> n, double length): CellBase(1, 1),
            refractive_index(n), delta_n(0.0), propagate_length(length) {}
        ~Waveguide() override = default;

        /**
         * Print the model properties and the protected or private parameters of the cell.
         */
        void print() override {

        }

        /**
         * Waveguide update strategy
         * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         */
        void update(double wavelength) override {
            using effects::constants::PI;
            std::complex<double> k = 2.0 * PI * (refractive_index + delta_n) / wavelength;
            E_out = exp(i * k * propagate_length) * E_in;
        }

        /**
         * Waveguide control strategy
         * @param ctrl_vars: control variables
         * @note CtrlFunc return value needs 1 element for delta_n [ \f$ \Delta n \f$ ]
         */
        void control(std::vector<double> ctrl_vars) override {
            if(!ctrl_func) {
                if(ctrl_vars.size() == 1) { delta_n = std::complex<double>(ctrl_vars[0], 0); }
                else if(ctrl_vars.size() == 2) { delta_n = std::complex<double>(ctrl_vars[0], ctrl_vars[1]); }
                else { delta_n = 0.0; }
            } else {
                auto ctrl_targets = ctrl_func(ctrl_vars);
                if(ctrl_targets.size() == 1) { delta_n = std::complex<double>(ctrl_targets[0], 0); }
                else if(ctrl_targets.size() == 2) { delta_n = std::complex<double>(ctrl_targets[0], ctrl_targets[1]); }
                else { delta_n = 0.0; }
            }
        }

        /**
         * Update the parameters of cell properties that availbale for optimization.
         * @param optim_map: former key for availableparameter names, latter key for value
         * @note a complex number or higher dimension should be divided into several double
         * @note available parameters:
         */
        void optimize(const std::map<std::string, double> &optim_map) override {
            for(const auto &item : optim_map) {
                if(item.first == "n_real") { refractive_index.real(item.second); }
                else if(item.first == "n_imag") { refractive_index.imag(item.second); }
                else if(item.first == "delta_n_real") { delta_n.real(item.second); }
                else if(item.first == "delta_n_imag") { delta_n.imag(item.second); }
                else if(item.first == "propagate_length") { propagate_length = item.second; }
                else { std::cout << YELLOW << "[WARN] {cells->Waveguide/optimize} undefined optimization parameter" << RESET << std::endl; }
            }
        }

    protected:
        std::complex<double> refractive_index;  //!< refractive index of waveguide [ \f$ n \f$ ]
        std::complex<double> delta_n;   //!< change of refractive index of waveguide [ \f$ \Delta n \f$ ]
        double propagate_length{};  //!< length of light propagation among media materials in unit: *micrometer* [ \f$ L \f$ ]
    };
}

#endif //OMNICLO_WAVEGUIDE_HPP
