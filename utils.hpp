
#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <string>
#include <functional>
#include <utility>


template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& v)
{
    os << "array(";
    for (const auto& x : v)
    {
        os << x;
        if (&x != &v[v.size()-1])
            os << ", ";
    }
    os << ")";
    return os;
}

template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, std::array<T, N>&& v)
{
    os << "array(";
    for (const auto& x : v)
    {
        os << x;
        if (&x != &v[v.size()-1])
            os << ", ";
    }
    os << ")";
    return os;
}

template<typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::pair<A, B>& p)
{
    os << "pair(" << p.first << ", " << p.second << ")";
    return os;
}


template<typename Function, typename... Args>
double time_function(Function&& func, size_t times, Args&&... args) {
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < times; ++i) {
        std::invoke(func, std::forward<Args>(args)...);
    }
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

template<typename Function, typename... Args>
std::pair<double, double> timing_mean_std(Function&& func, size_t times, Args&&... args) {
    assert(times > 2);

    std::vector<double> times_vec;
    for (size_t i = 0; i < times; ++i) {
        times_vec.push_back(time_function(func, 1, std::forward<Args>(args)...));
    }

    double sum = 0;
    for (const auto& t : times_vec) {
        sum += t;
    }

    double mean = sum / static_cast<double>(times);

    double sq_sum = 0;
    for (const auto& t : times_vec) {
        sq_sum += (t - mean) * (t - mean);
    }

    double stdev = std::sqrt(sq_sum / static_cast<double>(times));

    return std::make_pair(mean, stdev);
}

# endif // UTILS_HPP
