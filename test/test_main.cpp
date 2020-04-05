/**
 * @file test_main.hpp
 * @brief Testify the functions
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-04-04
 */

// C++ Built-in Header
#include <iostream>
#include <cmath>
#include <utility>
#include <ctime>
#include <chrono>
#include <limits>
#include <numeric>
#include <any>

// Project Header
#include "matrices/matrices.hpp"
#include "problems/optimizer.hpp"

//-------------------------------------------------//
#define RUN_TIME_WRAPPER(test_case)\
    {\
        std::cout << "------- Func: <" #test_case "> -------" << std::endl;\
        std::cout << "[TEST] Start running..." << std::endl;\
        auto _t_start = std::chrono::high_resolution_clock::now();\
        test_case;\
        auto _t_end = std::chrono::high_resolution_clock::now();\
        auto _t_delta = _t_end - _t_start;\
        std::cout << "[TEST] Finished for " << std::chrono::duration<double, std::milli>(_t_delta).count() << " ms." << std::endl;\
    }

//-------------------------------------------------//
using namespace pagmo;

//-------------------------------------------------//
void random_testbench(int times = 1) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> urdis(-10.0, 10.0);

    std::cout << "Random Lists: [";
    for(int index = 0; index < times; index++) {
        std::cout << urdis(rd) << ", ";
    }
    std::cout << "]" << std::endl;
}


//-------------------------------------------------//
vector_double model_func(const vector_double &dv) {
    return { dv[0] * (dv[0] + dv[1] * dv[2]) + pow(sin(dv[3]), dv[4]) };
}

/**
 * @brief \f$ y = x_0 * (x_0 + x_1 * x_2) + x_3^{x_4} \f$
 * @note \f$ x_{i} \in  [-2, 2] \f$
 */
struct prob_v0 {

    vector_double fitness(const vector_double &dv) const {
        return {model_func(dv)};
    }

    std::pair<vector_double, vector_double> get_bounds() const {
        return {{-2, -2., -2, -2., -2.}, {2, 2., 2., 2., 2.}};
    }

    std::string get_name() const {
        return "problem v0";
    }
};

std::string printVector(const vector_double &dv) {
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

void pagmo_testbench() {

    // Define Problems
    std::cout << "--- Problem loading ---" << std::endl;
    problem prob{prob_v0{}};

    std::cout << "--- Compute value ---" << std::endl;
    std::cout << "at: (-2, -1, 1, -1, 2), value = " << prob.fitness({-2, -1, 1, -1, 2})[0] << std::endl;

    std::cout << "--- Get boundary ---" << std::endl;
    std::cout << "low bound: " << printVector(prob.get_lb()) << std::endl;
    std::cout << "up bound: " << printVector(prob.get_ub()) << std::endl;

    std::cout << "--- Problem Brief ---" << std::endl;
    std::cout << prob << std::endl;

    // Define Algorithm
    algorithm algo{pso(10000)};

    population pop{prob, 24};

    pop = algo.evolve(pop);

    std::cout << "--- Population ---" << std::endl;
    std::cout << pop << std::endl;
}

//-------------------------------------------------//

void cells_test() {
    cells::Waveguide wg(effects::cal_n(1, 0), 1);
    wg.update(1.55);

    cells::Coupler cp(1, 1, 2, 3);

    cells::CellBase *cell = &cp;
    cell->control(std::vector<double>(1));

    using effects::constants::i;
    Eigen::Vector2cd E(3.0 + i * 4.0, 4.0 - i * 3.0);
    std::cout << "[Power] " << utils::print_vectorXcd(utils::cal_P(E)) << std::endl;
}

void matrices_test() {
    matrices::MatricesBase mb;
    std::map<std::string, double> m;
    m["propagate_length"] = effects::constants::PI / 4.0;

    mb.add_cell<cells::Coupler>(m);
    mb.add_cell<cells::Coupler>(m);
    mb.add_cell<cells::Coupler>(m);
    mb.add_cell<cells::Coupler>(m);
    mb.confirm_cell();

    mb.add_link(0, 0, 1, 0);
    mb.add_link(0,1,3,0);
    mb.add_link(2,0,1,1);
    mb.add_link(2,1,3,1);
    mb.confirm_link();

    size_t in_num = mb.inputs_size();
    Eigen::VectorXcd vec_in = Eigen::VectorXcd::Zero(in_num);
    /*
    for(int i = 0; i < in_num; i++) {
        vec_in(i) = cos(i) + effects::constants::i * sin(i);
    }
     */
    vec_in(0) = vec_in(1) = 1.0;
    vec_in(5) = 1.0;

    std::cout << "Inputs: E = " << utils::print_vectorXcd(vec_in) << std::endl;
    mb.init_inputs(vec_in);
    mb.update_iterate(1.55);
    std::cout << "Inputs: E = " << utils::print_vectorXcd(mb.get_inputs()) << std::endl;
    mb.update_iterate(1.55);
    std::cout << "Outputs: E = " << utils::print_vectorXcd(mb.get_outputs()) << std::endl;
    std::cout << "Outputs: P = " << utils::print_vectorXcd(utils::cal_P(mb.get_outputs())) << std::endl;
}

//-------------------------------------------------//
template<typename T>
void print_any(std::any var) {
    std::cout << "[ANY] Type = " << var.type().name() << ", Value = " << std::any_cast<T>(var) << std::endl;
}

void print_value(std::any var) {
    try { auto x = std::any_cast<int>(var); std::cout << "Int = " << x << std::endl; } catch(std::exception &ex) {}
    try { auto x = std::any_cast<double>(var); std::cout << "Double = " << x << std::endl; } catch(std::exception &ex) {}
    try { auto x = std::any_cast<std::complex<double>>(var); std::cout << "Complex = " << x << std::endl; } catch(std::exception &ex) {}
    try { auto x = std::any_cast<std::vector<double>>(var); std::cout << "Vector = " << printVector(x) << std::endl; } catch(std::exception &ex) {}
}

void any_test() {
    std::map<std::string, std::any> univ_map;
    univ_map["int"] = int(1);
    univ_map["double"] = double(1.0);
    univ_map["complex"] = std::complex<double>(1.0, 2.0);
    univ_map["vector"] = std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0, 10.0};

    std::cout << std::any_cast<int>(univ_map["int"]) << std::endl;
    std::cout << std::any_cast<double>(univ_map["double"]) << std::endl;
    std::cout << std::any_cast<std::complex<double>>(univ_map["complex"]) << std::endl;

    print_any<int>(univ_map["int"]);
    print_any<double>(univ_map["double"]);
    print_any<std::complex<double>>(univ_map["complex"]);

    for(const auto &itm: univ_map) {
        std::cout << "Key = " << itm.first << " : ";
        print_value(itm.second);
    }

    auto &v = std::any_cast<std::vector<double>&>(univ_map["vector"]);
    std::cout << "Before: " << printVector(v) << std::endl;
    for(auto &itm: v) {
        itm *= 2.0;
    }
    std::cout << "After: " << printVector(v) << std::endl;
    for(const auto &itm: univ_map) {
        std::cout << "Key = " << itm.first << " : ";
        print_value(itm.second);
    }
}

//-------------------------------------------------//
int main() {

    /*
    RUN_TIME_WRAPPER(random_testbench(100));
    RUN_TIME_WRAPPER(pagmo_testbench());
    RUN_TIME_WRAPPER(cells_test());
    RUN_TIME_WRAPPER(matrices_test());
    */
    RUN_TIME_WRAPPER(any_test());

    return 0;
}