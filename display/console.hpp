/**
 * @file console.hpp
 * @brief Console colorful display for command line
 * @author LIU-Yinyi
 * @version 0.1.0
 * @date 2020-04-04
 */

#ifndef OMNICLO_CONSOLE_HPP
#define OMNICLO_CONSOLE_HPP

#include <iostream>
#include <string>
#include <unistd.h>

#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define BLACK   "\033[30m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define WHITE   "\033[37m"

#define RED_WHITE       "\033[41;37m"
#define GREEN_WHITE     "\033[42;37m"
#define YELLOW_WHITE    "\033[43;37m"
#define BLUE_WHITE      "\033[44;37m"

#endif //OMNICLO_CONSOLE_HPP
