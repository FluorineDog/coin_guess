#include <gmpxx.h>
#include <vector>
#include <cmath>
#include "doglib/common/common.h"
#include <random>

using namespace doglib::common;

constexpr int64_t Alpha_Up = 121;
constexpr int64_t Alpha_Down = 1889;
constexpr int precision = 10000;
constexpr int META_SIZEOF = 30;
const mpf_class Alpha = mpf_class(Alpha_Up, precision) / mpf_class(Alpha_Down, precision);
const mpf_class Beta = 1 - 3 * Alpha;

using SeqType = std::vector<bool>;
const auto logAlpha = -log2f128(Alpha.get_d());
const auto logBeta = -log2f128(Beta.get_d());

std::map<int, int> meta_map;
SeqType *global_b;
using std::vector;

std::vector<int> bits2states(std::vector<bool> bits, int transaction_size) {
    std::vector<int> states;
    auto indicator = mpf_class(0, precision);
    for(auto bit: bits) {
        indicator += bit;
        indicator /= 2;
    }
    // step3: calculate next states
    std::array<int, 4> statistics = {};
    _Float128 entropy = 0;
    states.clear();
    while (entropy < transaction_size + logAlpha) {
        int multiple = mpf_class(indicator / Alpha).get_ui();
        if (multiple < 3) {
            states.push_back(multiple);
            indicator -= multiple * Alpha;
            indicator /= Alpha;
            entropy += logAlpha;
            statistics[multiple]++;
        } else {
            states.push_back(3);
            indicator -= 3 * Alpha;
            indicator /= Beta;
            entropy += logBeta;
            statistics[3]++;
        }
    }
    return states;
}


std::pair<SeqType, SeqType> generate_a(int N, const SeqType &seq_ans) {
    SeqType seqA(N);
    SeqType seqB(N);
    auto fill_meta = [&seqA, &seqB](int offset, int meta) {
        meta_map[offset] = meta;
        for (auto i: Range(0, META_SIZEOF)) {
            auto index = offset + i;
            assert(meta < (1 << META_SIZEOF));
            bool bit = (meta >> i) & 1;
            seqA[index] = bit;
            seqB[index] = false;
        }
    };
    int transaction_size;
    std::vector<int> states(4096, 3);
    for (int step = 0;; step++) {
        // step 1: fill info, calculate next indicator
        transaction_size = states.size();
        vector<bool> bits;
        for (auto i: Range(states.size())) {
            auto index = N - 1 - i;
            auto state = states[i];
            auto bit_ans = seq_ans[index];
            auto bitA = !(state % 2) ^bit_ans;
            auto bitB = !(state / 2) ^bit_ans;
            seqA[index] = bitA;
            seqB[index] = bitB;
            bits.push_back(bitB);
        }
        N -= transaction_size;

        // step 2: fill meta
        N -= META_SIZEOF;
        fill_meta(N, transaction_size);

        states = bits2states(bits, transaction_size);

        if (N <= 2 * states.size() + 2 * META_SIZEOF) {
            break;
        }
    }

    // step 6 set meta and other bits
    fill_meta(0, N);

    for (auto i: Range(transaction_size)) {
        auto dst_index = i + META_SIZEOF;
        auto src_index = N + META_SIZEOF + i;
        seqA[dst_index] = seqB[src_index];
        seqB[dst_index] = false;
    }

    for (auto index: Range(META_SIZEOF + transaction_size, N)) {
        seqA[index] = seq_ans[index];
        seqB[index] = false;
    }
    return {seqA, seqB};
}

bool inc(std::vector<int> &states) {
    for (auto &x: states) {
        if (x < 3) {
            ++x;
            return false;
        } else {
            x = 0;
        }
    }
    return true;
}

SeqType verify_b(int N, const SeqType &seqA, const SeqType &seq_std) {
    SeqType seqB(N);
    int ack = 0;
    auto setB = [&](int offset, bool value) {
        assert(ack == offset);
        seqB[offset] = value;
        bool ref_bit = (*global_b)[offset];

        if(ref_bit != value) {
            auto str1 = std::to_string(ref_bit);
            auto str2 = std::to_string(value);
            throw std::runtime_error("fuck" + str1 + str2);
        }
        ++ack;
    };
    auto getA = [&](int offset) {
        assert(offset < ack);
        return seqA[offset];
    };
    auto getStd = [&](int offset) {
        assert(offset < ack);
        return seq_std[offset];
    };

    auto fetch_meta = [&](int offset) {
        int n = 0;
        assert(meta_map.count(offset));
        for (auto i: Range(0, META_SIZEOF)) {
            auto index = offset + i;
            setB(index, false);
            auto bit = getA(index);
            n |= bit << i;
        }
        auto expected = meta_map[offset];
        assert(meta_map[offset] == n);
        return n;
    };

    auto iter = fetch_meta(0);
    for (auto i: Range(META_SIZEOF, iter)) {
        setB(i, false);
    }
    std::vector<bool> candidateB;
    for (auto i: Range(META_SIZEOF, iter)) {
        // including paddings for the first time
        candidateB.push_back(getA(i));
    }

    while (iter < N) {
        auto transaction_size = fetch_meta(iter);
        iter += META_SIZEOF;
        std::vector<int> states;
        for (auto i: Range(0, transaction_size)) {
            auto dst_index = iter + i;
            bool bitB = candidateB[i];
            setB(dst_index, bitB);
            auto bitA = getA(dst_index);
            auto std_bit = getStd(dst_index);
            bitB ^= std_bit;
            bitA ^= std_bit;
            states.push_back(3 - bitA - bitB * 2);
        }
        iter += transaction_size;
        auto indicator = mpf_class(0, precision);
        if (false) {
            indicator += 1;
        } else {
            for (auto x: states) {
                if (x < 3) {
                    indicator *= Alpha;
                    indicator += x * Alpha;
                } else {
                    indicator *= Beta;
                    indicator += 3 * Alpha;
                }
            }
        }
        std::vector<bool> transaction;
        for (auto i: Range(0, transaction_size)) {
            indicator *= 2;
            auto bit = indicator.get_ui();
            indicator -= bit;
            transaction.push_back(bit);
        }

        candidateB.clear();
        for (bool bit: transaction) {
            candidateB.push_back(bit);
        }
    }
    assert(iter == N);
    assert(ack == N);
    return seqB;
}


int main() {
    SeqType expected;
    std::default_random_engine e(42);
    int N = 1 << 15;
    for (auto i: Range(0, N)) {
        bool bit = e() % 2;
        expected.push_back(bit);
    }
    auto[a, b] = generate_a(N, expected);
    global_b = &b;
    auto vb = verify_b(N, a, expected);

    std::vector<int> statistics(4, 0);
    if(b != vb) {
        throw std::runtime_error("fuck");
    }
    for(auto i: Range(N)) {
        auto a_yes = expected[i] == a[i];
        auto b_yes = expected[i] == b[i];
        auto state = a_yes + b_yes * 2;
        statistics[state]++;
    }
    for(auto i: Range(4)) {
        double count = statistics[i];
        std::cout << "state " << i << "->" << count << "->" << count / N << std::endl;
    }
}




