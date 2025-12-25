#include <iostream>
#include <array>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <random>
#include <algorithm>
#include <numeric>
#include <bitset>
#include <chrono>
#include <string>
#include <sstream>

constexpr size_t n = 3;
constexpr size_t r = 19;
constexpr size_t g = 2;

constexpr size_t d3 = (r + 3) / 6;
constexpr size_t p3 = std::min(g, d3);
constexpr size_t s3 = g - p3;
constexpr size_t p1 = s3 ? 0 : r / 2 - 3 * p3;
constexpr size_t s1 = r - 3 * g - p1;

constexpr size_t factorial(size_t x) {
    size_t res = 1;
    for (size_t i = 2; i <= x; ++i) {
        res *= i;
    }
    return res;
}

constexpr size_t AUT = 2 * factorial(n);
constexpr size_t N = 3*n;

template<typename T, size_t k>
class cwr_it {
    const std::vector<T>& Src;
    std::array<decltype(Src.begin()), k> Its{};
public:
    cwr_it(const std::vector<T>& src) : Src(src) {
        Its.fill(Src.begin());
    }

    std::array<T, k> operator*() const noexcept {
        std::array<T, k> result;
        for (size_t i = 0; i < k; ++i)
            result[i] = *Its[i];
        return result;
    }

    cwr_it& operator++() noexcept {
        size_t i;
        for (i = k; i; --i)
            if (++Its[i-1] != Src.end())
                break;
        if (i)
            for (; i < k; ++i)
                Its[i] = Its[i-1];
        return *this;
    }

    bool operator==(const cwr_it& oth) const noexcept {
        return Its == oth.Its;
    }

    void exhaust() noexcept {
        Its.fill(Src.end());
    }
};

template<typename T>
class cwr_it<T, 0> {
    size_t K = 0;
public:
    cwr_it(const std::vector<T>& src) {
    }

    std::array<T, 0> operator*() const noexcept {
        return {};
    }

    cwr_it& operator++() noexcept {
        ++K;
        return *this;
    }

    bool operator==(const cwr_it& oth) const noexcept {
        return K == oth.K;
    }

    void exhaust() noexcept {
        K = 1;
    }
};

template<typename T, size_t k>
class cwr {
    const std::vector<T>& Src;
public:
    cwr(const std::vector<T>& src) : Src(src) {}
    
    cwr_it<T, k> begin() const noexcept {
        return cwr_it<T, k>(Src);
    }

    cwr_it<T, k> end() const noexcept {
        auto result = begin();
        result.exhaust();
        return result;
    }
};

template<typename T, size_t k1, size_t k2>
std::array<T, k1+k2> concatenate(const std::array<T, k1>& a1, const std::array<T, k2>& a2) {
    std::array<T, k1+k2> result;
    for (size_t i = 0; i < k1; ++i) {
        result[i] = a1[i];
    }
    for (size_t i = 0; i < k2; ++i) {
        result[k1+i] = a2[i];
    }
    return result; 
}

bool get(uint32_t m, uint32_t i, uint32_t j) {
    return m & (1u << (3*i+j));
}

void set(uint32_t& m, uint32_t i, uint32_t j, bool v) {
    if (v) {
        m |= (1u << (3*i+j));
    } else {
        m &= ~(1u << (3*i+j));
    }
}

uint32_t cycle(uint32_t m) {
    uint32_t res = 0;
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < 3; ++j) {
            set(res, i, j, get(m, i, (j+1)%3));
        }
    }
    return res;
}

std::array<uint32_t, 1u << N> normal_form;
std::array<uint32_t, (1u << N) + (1u << n)> unfold;
uint32_t RHS = 0;
std::array<std::array<uint32_t, (1u << N) + (1u << n)>, AUT> apply;
std::vector<uint32_t> ms;
std::vector<uint32_t> cs;

void init() {
    for (uint32_t code = 0; code < normal_form.size(); ++code) {
        uint32_t tmp = cycle(code);
        normal_form[code] = std::min(code, tmp);
        tmp = cycle(tmp);
        normal_form[code] = std::min(normal_form[code], tmp);
        if (normal_form[code] == code) {
            ms.push_back(code);
        }
    }
    for (uint32_t code = 0; code < normal_form.size(); ++code) {
        uint32_t res = 0;
        for (uint32_t i = 0; i < n; ++i) {
            for (uint32_t j = 0; j < n; ++j) {
                for (uint32_t k = 0; k < n; ++k) {
                    res ^= ((get(code, i, 0) & get(code, j, 1) & get(code, k, 2))
                          ^ (get(code, i, 1) & get(code, j, 2) & get(code, k, 0))
                          ^ (get(code, i, 2) & get(code, j, 0) & get(code, k, 1))) << (n*n*i + n*j + k);
                }
            }
        }
        unfold[code] = res;
    }
    for (uint32_t code = (1u << N); code < (1u << N) + (1u << n); ++code) {
        uint32_t res = 0;
        for (uint32_t i = 0; i < n; ++i) {
            for (uint32_t j = 0; j < n; ++j) {
                for (uint32_t k = 0; k < n; ++k) {
                    res ^= (1u & (code >> i) & (code >> j) & (code >> k)) << (n*n*i + n*j + k);
                }
            }
        }
        unfold[code] = res;
        cs.push_back(code);
    }
    for (uint32_t k = 0; k < n; ++k) {
        RHS ^= 1u << ((n*n + n + 1) * k);
    }
    for (uint32_t code = 0; code < normal_form.size(); ++code) {
        uint32_t k = 0;
        for (uint32_t swap = 0; swap <= 1; ++swap) {
            std::array<uint32_t, n> p;
            std::iota(p.begin(), p.end(), 0);
            do {
                uint32_t res = 0;
                for (uint32_t i = 0; i < n; ++i) {
                    for (uint32_t j = 0; j < 3; ++j) {
                        set(res, i, j, get(code, p[i], j == 2 ? j : (j^swap)));
                    }
                }
                apply[k][code] = normal_form[res];
                ++k;
            } while (std::ranges::next_permutation(p).found);
        }
    }
    for (uint32_t code = (1u << N); code < (1u << N) + (1u << n); ++code) {
        uint32_t k = 0;
        for (uint32_t swap = 0; swap <= 1; ++swap) {
            std::array<uint32_t, n> p;
            std::iota(p.begin(), p.end(), 0);
            do {
                uint32_t res = 1u << N;
                for (uint32_t i = 0; i < n; ++i) {
                    res |= ((code >> p[i]) & 1u) << i;
                }
                apply[k][code] = res;
                ++k;
            } while (std::ranges::next_permutation(p).found);
        }
    }

}

bool is_optimal(const std::array<uint32_t, r - 2*g> ds) {
    for (uint32_t f = 1; f < AUT; ++f) {
        std::array<uint32_t, r - 2*g> gds;
        for (size_t t = 0; t < ds.size(); ++t) {
            gds[t] = apply[f][ds[t]];
        }
        std::ranges::sort(gds);
        if (gds < ds) {
            return false;
        }
    }
    return true;
}

template <size_t k>
uint32_t tp(const std::array<uint32_t, k>& codes) {
    uint32_t s = 0;
    for (uint32_t code : codes)
        s ^= unfold[code];
    return s;
}

uint64_t glue(const std::array<uint32_t, r-2*g>& codes) {
    uint64_t result = 0;
    for (auto code : codes) {
        result <<= code & (1u << N) ? n : N;
        result |= code & ((1u << N) - 1u);
    }
    return result;
}

std::vector<uint64_t> generate_data() {
    init();
    std::unordered_map<uint32_t, std::vector<std::array<uint32_t, s3 + s1>>> halves;
    for (auto c3 : cwr<uint32_t, s3>(ms)) {
    	for (auto c1 : cwr<uint32_t, s1>(cs)) {
            auto codes = concatenate(c3, c1);
            halves[RHS ^ tp(codes)].push_back(codes);
        }
    }
    std::vector<uint64_t> result;
    for (auto c3 : cwr<uint32_t, p3>(ms)) {
        for (auto c1 : cwr<uint32_t, p1>(cs)) {
            auto h1 = concatenate(c3, c1);
            for (auto& h2 : halves[tp(h1)]) {
                if (h2.front() < h1.back())
                    continue;
                auto codes = concatenate(h1, h2);
                if (is_optimal(codes))
                    result.push_back(glue(codes));
            }
        }
    }
    return result;
}

int main(int argc, char** argv) {
    uint32_t signature = (n << 16u) | (r << 8u) | g;
    std::cout << "n = " << n << ", r = " << r << ", g = " << g << std::endl;
    std::cout << p1 << ' ' << p3 << ' ' << s1 << ' ' << s3 << std::endl;
    auto start = std::chrono::steady_clock::now();
    auto data = generate_data();
    auto end = std::chrono::steady_clock::now();
    uint64_t size = data.size();
    std::cout << "generated " << size << " in " <<
        std::chrono::duration<double>(end-start).count() << std::endl;
    if (argc < 2) {
        return 0;
    }
    std::mt19937 gen;
    std::ranges::shuffle(data, gen);
    if (argv[1] == std::string("-p")) {
        uint32_t chunksize = std::stoul(argv[2]);
        std::stringstream name;
        name << "p_" << chunksize << '_' << n << '_' << r << '_' << g;
        std::ofstream out(name.str(), std::ios::binary);
        out.write(reinterpret_cast<char*>(&signature), sizeof(signature));
        out.write(reinterpret_cast<char*>(&chunksize), sizeof(chunksize));
        for (uint32_t i = 0; i < chunksize; ++i) {
            out.write(reinterpret_cast<char*>(&data[i]), sizeof(data[i]));
        }
    } else if (argv[1] == std::string("-c")) {
        uint32_t nchunks = std::stoul(argv[2]);
        for (uint32_t k = 0; k < nchunks; ++k) {
            std::stringstream name;
            name << "c_" << k << "_of_" << nchunks << '_' << n << '_' << r << '_' << g;
            std::ofstream out(name.str(), std::ios::binary);
            uint32_t chunksize = size / nchunks + (k < (size % nchunks) ? 1 : 0);
            out.write(reinterpret_cast<char*>(&signature), sizeof(signature));
            out.write(reinterpret_cast<char*>(&chunksize), sizeof(chunksize));
            for (uint32_t i = k; i < size; i += nchunks) {
                out.write(reinterpret_cast<char*>(&data[i]), sizeof(data[i]));
            }
        }
    }
}

