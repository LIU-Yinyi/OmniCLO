/**
 * @file coupler.hpp
 * @brief Relization of Coupler Cell
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-04-03
 */

#ifndef OMNICLO_COUPLER_HPP
#define OMNICLO_COUPLER_HPP

// Custom Library Header
#include "cells_base.hpp"

namespace cells {

    using effects::constants::i;



    /**
     * Coupler derived from CellBase.
     *
     * When \f$ \kappa_1 = \kappa_2 \f$ means Symmetrical Coupler, \n
     * When \f$ \kappa_1 \neq \kappa_2 \f$ means Asymmetrical Coupler.
     */
    class Coupler : public CellBase {
    public:
        Coupler(): CellBase(2, 2), refractive_index(1.0), delta_n(0.0), propagate_length(10.0),
            kappa1(1.0), kappa2(1.0) {}

        /**
         * Constructor with Coupler's properties
         * @param kappa1: coupling coefficients 1 [ \f$ \kappa_1 \f$ ]
         * @param kappa2: coupling coefficients 2 [ \f$ \kappa_2 \f$ ]
         * @param length: length of light propagation among media materials in unit: *micrometer* [ \f$ L \f$ ]
         */
        Coupler(std::complex<double> n, double length, double kappa1, double kappa2): CellBase(2, 2),
            refractive_index(n), delta_n(0.0), propagate_length(length),
            kappa1(kappa1), kappa2(kappa2) {}
        ~Coupler() override = default;

        /**
         * Print the model properties and the protected or private parameters of the cell.
         */
        void print() override {

        }

        /**
         * Coupler update strategy
         * @param wavelength: wavelength of light at vacuum in unit: *micrometer* [ \f$ \lambda_0 \f$ ]
         */
        void update(double wavelength) override {
            using Eigen::Matrix2cd;
            using effects::constants::PI;
            Matrix2cd m;
            double delta_beta = 2.0 * PI * delta_n.real() / wavelength;
            double kappa = sqrt(kappa1 * kappa2);
            double K = sqrt(pow(delta_beta/2.0, 2.0) + pow(kappa, 2.0));
            double Kz = K * propagate_length;
            double cosKz = cos(Kz);
            double sinKz = sin(Kz);
            std::complex<double> exp_ipz = exp(i * delta_beta * propagate_length / 2.0);
            std::complex<double> exp_ipzn = exp(-i * delta_beta * propagate_length / 2.0);
            std::complex<double> d1 = exp_ipz * (cosKz - i * delta_beta * sinKz / (2 * K));
            std::complex<double> d2 = i * exp_ipzn * (kappa / K) * sinKz;
            m(0, 0) = d1;
            m(1, 1) = std::conj(d1);
            m(0, 1) = m(1, 0) = d2;
            E_out = m * E_in;
        }

        /**
         * Coupler control strategy
         * @param ctrl_vars: control variables
         * @note CtrlFunc return value needs 1 element for delta_beta [ \f$ \Delta \beta \f$ ]
         */
        void control(std::vector<double> ctrl_vars) override {
            if(!ctrl_func) {
                assert(ctrl_vars.size() == 1);
                delta_n = std::complex<double>(0.0, ctrl_vars[0]);
            } else {
                auto ctrl_targets = ctrl_func(ctrl_vars);
                assert(ctrl_targets.size() == 1);
                delta_n = std::complex<double>(0.0, ctrl_targets[0]);
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
                else if(item.first == "kappa1") { kappa1 = item.second; }
                else if(item.first == "kappa2") { kappa2 = item.second; }
                else { std::cout << YELLOW << "[WARN] {cells->Coupler/optimize} undefined optimization parameter" << RESET << std::endl; }
            }
        }

    protected:
        double kappa1, kappa2;  //!< coupling coefficients of two line waveguides
        std::complex<double> refractive_index;  //!< refractive index of waveguide [ \f$ n \f$ ]
        std::complex<double> delta_n;   //!< change of refractive index of waveguide [ \f$ \Delta n \f$ ]
        double propagate_length;    //!< length of light propagation among media materials in unit: *micrometer* [ \f$ L \f$ ]
    };
}

#endif //OMNICLO_COUPLER_HPP
