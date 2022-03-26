#ifndef UTILS_HPP
#define UTILS_HPP

#include <wup/wup.hpp>


template <typename T>
void
arraySquaredCollapse(T * data, size_t n)
{
    for (size_t i=0;i!=n;++i)
        data[i] = data[i] * data[i];
}

template <typename T>
T
arraySumCollapse(T * data, size_t n)
{
    size_t h = 0;

    while (n > 1)
    {
        if (n % 2)
        {
            data[0] += data[n-1];
            --n;
        }

        h = n >> 1;

        for (size_t k=0;k!=h;++k)
            data[k] += data[k+h];

        n = h;
    }

    return data[0];
}

template <typename T>
T
arraySumCollapse(std::vector<T> & data)
{
    return arraySumCollapse(data.data(), data.size());
}

template <typename T>
void
arraySquaredCollapse(std::vector<T> & data)
{
    return arraySquaredCollapse(data.data(), data.size());
}

#endif // UTILS_HPP
