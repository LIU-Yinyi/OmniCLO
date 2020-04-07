/**
 * @file terminator.hpp
 * @brief Relization of Terminator Generalized Cell
 * @author LIU-Yinyi
 * @version 1.0.1
 * @date 2020-04-06
 */

#ifndef OMNICLO_TERMINATOR_HPP
#define OMNICLO_TERMINATOR_HPP

// Custom Library Header
#include "cells_base.hpp"

namespace cells {

    class InputPort: public CellBase {
    public:
        InputPort(const std::string &name, const std::map<std::string, std::any> &vars):
            CellBase(0, 1, name) { device_vars = vars; }

        ~InputPort() override = default;

        void update(double wavelength) override {}
        void control(std::map<std::string, std::any> ctrl_vars) override {}

        void set_E(const VectorXcd &E) override { assert(E.size() == E_out.size()); E_out = E; is_updated = false; }
        void set_E_in(const VectorXcd &E) override { E_out = E; }
        VectorXcd get_E_in() const override { return E_out; }
    };

    class OutputPort: public CellBase {
    public:
        OutputPort(const std::string &name, const std::map<std::string, std::any> &vars):
                CellBase(1, 0, name) { device_vars = vars; }

        ~OutputPort() override = default;

        void update(double wavelength) override {}
        void control(std::map<std::string, std::any> ctrl_vars) override {}

        VectorXcd get_E() const override { return E_in; }
        void set_E_out(const VectorXcd &E) override { E_in = E; }
        VectorXcd get_E_out() const override { return E_in; }
    };
}

#endif //OMNICLO_TERMINATOR_HPP
