
#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <string>

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::valarray<T>& v)
{
    os << "valarray(";
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
std::ostream& operator<<(std::ostream& os, std::valarray<T>&& v)
{
    os << "valarray(";
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

# endif // UTILS_HPP
