/**
 * @file Distribution.h
 * @author JunPing Yuan
 * @brief Distribution1D Distribution2D
 * @version 0.1
 * @date 2022-09-29
 *
 * @copyright NJUMeta (c) 2022
 * www.njumeta.com
 *
 */

#pragma once

#include <memory>
#include <vector>
#include "CoreLayer/Math/Geometry.h"

#include <vector>

template <typename Predicate>
static int FindInterval(int size, const Predicate& pred) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        // Bisect range based on value of _pred_ at _middle_
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        }
        else
            len = half;
    }
    return std::clamp(first - 1, 0, size - 2);
}



struct Distribution1D {

    Distribution1D(const double* f, int n) : func(f, f + n), cdf(n + 1) {
        // Compute integral of step function at $x_i$
        cdf[0] = 0;
        for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;

        // Transform step function integral into CDF
        funcInt = cdf[n];
        if (funcInt == 0) {
            for (int i = 1; i < n + 1; ++i) cdf[i] = double(i) / double(n);
        }
        else {
            for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
        }
    }

    int SampleDiscrete(double sample, double* pdf) const {
        int offset = FindInterval((int)cdf.size(),
            [&](int index) { return cdf[index] <= sample; });
        if (pdf) *pdf = DiscretePDF(offset);
        return offset;
    }

    double SampleContinuous(double u, double* pdf, int* off = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindInterval((int)cdf.size(),
            [&](int index) { return cdf[index] <= u; });
        if (off) *off = offset;
        // Compute offset along CDF segment
        double du = u - cdf[offset];
        if ((cdf[offset + 1] - cdf[offset]) > 0) {
            du /= (cdf[offset + 1] - cdf[offset]);
        }

        // Compute PDF for sampled offset
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;
        if (pdf) *pdf = (funcInt > 0) ? cdf[offset + 1] - cdf[offset] : 0;

        // Return $x\in{}[0,1)$ corresponding to sample
        return (offset + du) / Count();
    }

    double DiscretePDF(int index) const {
        if (index >= func.size() || funcInt == 0) {
            return 0;
        }
        return func[index] / (funcInt * Count());
    }

    void warp(double& sample, int& index) {
        index = FindInterval((int)cdf.size(),
            [&](int index) { return cdf[index] <= sample; });
        sample = (sample - cdf[index]) / func[index];
    }

    int Count() const {
        return int(func.size());
    }
    std::vector<double> func, cdf;
    double funcInt;
};


