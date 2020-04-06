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
#include <regex>
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
    cells::Waveguide wg("waveguide", {});
    wg.update(1.55);

    cells::Coupler cp("coupler", {});

    cells::CellBase *cell = &cp;
    cell->control(std::map<std::string, std::any>{{"delta_n", std::complex<double>(1.0, 1.0)}});

    using effects::constants::i;
    Eigen::Vector2cd E(3.0 + i * 4.0, 4.0 - i * 3.0);
    std::cout << "[Power] " << utils::print_vectorXcd(utils::cal_P(E)) << std::endl;
}

void microring_test() {
    using effects::constants::i;
    std::map<std::string, std::any> config;
    config["ring_number"] = size_t(16);
    config["couple_length"] = std::vector<double>(17, 0.0);
    config["n"] = std::vector<std::complex<double>>(16, 1.4 + 6.04e-6 * i);
    cells::MicroRingSerial mrs("microring", config);
    mrs.print();
    mrs.set_E(Eigen::Vector2cd(1.0, 0.0));
    mrs.update(1.55);
    std::cout << "E_in = " << utils::print_vectorXcd(mrs.get_E_in()) << std::endl;
    std::cout << "E_out = " << utils::print_vectorXcd(mrs.get_E()) << std::endl;
    std::cout << "P_in = " << utils::print_vectorXcd(utils::cal_P(mrs.get_E_in())) << std::endl;
    std::cout << "P_out = " << utils::print_vectorXcd(utils::cal_P(mrs.get_E())) << std::endl;
}

void matrices_test() {
    matrices::MatricesBase mb;
    std::map<std::string, std::any> m;
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

    mb.print();
}

//-------------------------------------------------//
template<typename T>
void print_any(std::any var) {
    std::cout << "[ANY] Type = " << var.type().name() << ", Value = " << std::any_cast<T>(var) << std::endl;
}

void print_value(std::any var) {
    try { auto x = std::any_cast<size_t>(var); std::cout << "Size_t = " << x << std::endl; } catch(std::exception &ex) {}
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
struct gen_rand {
    double range;
public:
    explicit gen_rand(double r = 1.0) : range(r/(double)RAND_MAX) {}
    double operator()() {
        return rand() * range;
    }
};

std::vector<double> gen_rand_vec(size_t num, double range) {
    std::vector<double> vec(num);
    std::generate_n(vec.begin(), num, gen_rand(range));
    return vec;
}

std::string naming_utils(const std::string &prefix, int index) {
    std::stringstream ss;
    ss << prefix << index;
    return ss.str();
};

void regex_test() {
    std::regex re("([a-z]+)([0-9]+)");
    std::cmatch match_res;

    std::string test01{naming_utils("abc", 12)};
    if(std::regex_match(test01.c_str(), match_res, re)) {
        for(const auto &itm: match_res) {
            std::cout << itm.str() << ", ";
        }
        std::cout << std::endl;
    }

}

void any_dynamic_vector_test() {
    size_t N = 10;
    std::map<std::string, std::any> univ_map;
    for(size_t idx = 0; idx < N; idx++) {
        univ_map[naming_utils("kappa", idx)] = idx * 10;
    }
    for(size_t idx = N; idx < 2 * N; idx++) {
        univ_map[naming_utils("kappa", idx)] = gen_rand_vec(idx, 10);
    }
    for(const auto &itm: univ_map) {
        std::cout << "Key = " << itm.first << " : ";
        print_value(itm.second);
    }
}

void std_complex_vector_vs_eigen_vectorXcd() {
    std::cout << std::is_same<std::vector<std::complex<double>>, Eigen::VectorXcd>::value << std::endl;
}

//-------------------------------------------------//
int main() {

    RUN_TIME_WRAPPER(random_testbench(100));
    RUN_TIME_WRAPPER(pagmo_testbench());
    RUN_TIME_WRAPPER(cells_test());
    RUN_TIME_WRAPPER(matrices_test());

    RUN_TIME_WRAPPER(microring_test());
    RUN_TIME_WRAPPER(any_test());
    RUN_TIME_WRAPPER(any_dynamic_vector_test());
    RUN_TIME_WRAPPER(regex_test());
    RUN_TIME_WRAPPER(std_complex_vector_vs_eigen_vectorXcd());
    return 0;
}