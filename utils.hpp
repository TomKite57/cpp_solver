
#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <string>

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

template<typename T>
std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2)
{
    assert(v1.size() == v2.size());
    std::vector<T> v(v1.size());
    for (size_t i=0; i<v1.size(); ++i)
        v[i] = v1[i] + v2[i];
    return v;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2)
{
    assert(v1.size() == v2.size());
    std::vector<T> v(v1.size());
    for (size_t i=0; i<v1.size(); ++i)
        v[i] = v1[i] - v2[i];
    return v;
}

template<typename T>
std::vector<T> operator*(const std::vector<T>& v1, const std::vector<T>& v2)
{
    assert(v1.size() == v2.size());
    std::vector<T> v(v1.size());
    for (size_t i=0; i<v1.size(); ++i)
        v[i] = v1[i] * v2[i];
    return v;
}


template<typename F, typename... Args>
double time_function(F&& f, Args&&... args, size_t repeats=1)
{
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i=0; i<repeats; ++i)
        std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

# endif // UTILS_HPP
