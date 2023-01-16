#include "EmpiricalModeDecomposition.hpp"

#include <cstdio>
#include <cmath>

EMD::EMD() {
    _emd = new EMDData();
}

EMD::~EMD() {
    if (_emd) {
        delete _emd;
    }
}

// public
int EMD::setup(int max_order, int max_iterations, int locality, float imf_thresh) {
    _emd->order = max_order;
    _emd->iterations = max_iterations;
    _emd->locality = locality;
    _emd->imf_thresh = imf_thresh;
    return 0;
}

bool EMD::check_peaks(const std::vector<float>& cur_imf) {
    for (int i = 0; i < _emd->max_size; ++i) {
        auto cur_max = cur_imf[_emd->max_points[i]];
        if (cur_max < 0) {
            return false;
        }
    }
    for (int i = 0; i < _emd->min_size; ++i) {
        auto cur_min = cur_imf[_emd->min_points[i]];
        if (cur_min > 0) {
            return false;
        }
    }
    return true;
}

std::vector<std::vector<float>> EMD::decompose(const std::vector<float>& signal) {
    std::vector<std::vector<float>> imfs;
    if (signal.empty()) {
        return imfs;
    }

    // temp init
    imfs.resize(1);
    imfs[0] = signal;
    _emd->residual = signal;
    _emd->min_points.resize(signal.size() / 2);
    _emd->max_points.resize(signal.size() / 2);
    _emd->min.resize(signal.size());
    _emd->max.resize(signal.size());

    int i = 0;
    while (i < _emd->order) {
        auto& cur_imf = imfs[i++];
        int j = 0;
        while (_emd->iterations == 0 || j++ < _emd->iterations) {
            make_extrema(cur_imf);
            printf("max_size: %d min_size: %d\n", _emd->max_size, _emd->min_size);

            if (_emd->min_size < 4 || _emd->max_size < 4) {
                break;  // can't fit splines
            }
            if (check_peaks(cur_imf)) {
                printf("peaks\n");
                break;
            }

            interpolate(cur_imf, _emd->min, _emd->min_points, _emd->min_size);
            interpolate(cur_imf, _emd->max, _emd->max_points, _emd->max_size);
            update_imf(cur_imf);
            printf("order: %d iter: %d\n", i - 1, j - 1);
        }
        make_residue(cur_imf);

        // calc residual sum
        auto sum = abs_sum(_emd->residual);
        printf("order: %d iter: %d\n", i - 1, j - 1);

        if (sum < _emd->imf_thresh) {
            printf("sums: %f\n", sum);
            break;  // stop calc
        }
        imfs.push_back(_emd->residual);
    }

    // clear temp vars
    emd_clear();

    return std::move(imfs);
}

// private
inline int EMD::mirror_index(int i, int size) {
    if (i < size) {
        if (i < 0) {
            return -i - 1;
        }
        return i;
    }
    return (size - 1) + (size - i);
}

inline float EMD::abs_sum(const std::vector<float>& vec) {
    float sum = 0.0f;
    for (auto& v : vec) {
        sum += std::abs(v);
    }
    return sum;
}

void EMD::make_extrema(const std::vector<float>& cur_imf) {
    int last_min = 0;
    int last_max = 0;
    _emd->min_size = 0;
    _emd->max_size = 0;
    for (int i = 1; i < cur_imf.size() - 1; i++) {
        if (cur_imf[i - 1] < cur_imf[i]) {
            if (cur_imf[i] > cur_imf[i + 1] && (i - last_max) > _emd->locality) {
                _emd->max_points[_emd->max_size++] = i;
                last_max = i;
            }
        } else {
            if (cur_imf[i] < cur_imf[i + 1] && (i - last_min) > _emd->locality) {
                _emd->min_points[_emd->min_size++] = i;
                last_min = i;
            }
        }
    }
}

void EMD::interpolate(
        const std::vector<float>& in,
        std::vector<float>& out,
        std::vector<int>& points,
        int points_size) {
    int size = in.size();
    int i, j, i0, i1, i2, i3, start, end;
    float a0, a1, a2, a3;
    float y0, y1, y2, y3, mu_scale, mu;
    for (i = -1; i < points_size; i++) {
        i0 = points[mirror_index(i - 1, points_size)];
        i1 = points[mirror_index(i, points_size)];
        i2 = points[mirror_index(i + 1, points_size)];
        i3 = points[mirror_index(i + 2, points_size)];

        y0 = in[i0];
        y1 = in[i1];
        y2 = in[i2];
        y3 = in[i3];

        a0 = y3 - y2 - y0 + y1;
        a1 = y0 - y1 - a0;
        a2 = y2 - y0;
        a3 = y1;

        // left boundary
        if (i == -1) {
            start = 0;
            i1 = -i1;
        } else {
            start = i1;
        }

        // right boundary
        if (i == points_size - 1) {
            end = size;
            i2 = size + size - i2;
        } else {
            end = i2;
        }

        mu_scale = 1.f / (i2 - i1);
        for (j = start; j < end; j++) {
            mu = (j - i1) * mu_scale;
            out[j] = ((a0 * mu + a1) * mu + a2) * mu + a3;
        }
    }
}

void EMD::update_imf(std::vector<float>& imf) {
    for (int i = 0; i < imf.size(); i++) {
        imf[i] -= (_emd->min[i] + _emd->max[i]) * .5f;
    }
}

void EMD::make_residue(const std::vector<float>& cur) {
    for (int i = 0; i < cur.size(); i++) {
        _emd->residual[i] -= cur[i];
    }
}

void EMD::emd_clear() {
    _emd->clear();
}
