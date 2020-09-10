#include <gmpxx.h>
#include <vector>
#include <cmath>
#include "doglib/common/common.h"
#include <random>

using namespace doglib::common;

// TODO: modify these two parameters
constexpr int N = 1 << 22; // 4 * 1024 * 1024
constexpr int initial_transaction_size = 10000;

// transaction size should be much smaller than this
constexpr int precision = initial_transaction_size * 2;


constexpr int64_t Alpha_Up = 119;
constexpr int64_t Alpha_Down = 1889;
constexpr int META_SIZEOF = 16;
const mpf_class Alpha = mpf_class(Alpha_Up, precision) / mpf_class(Alpha_Down, precision) + 0.0001;
const mpf_class Beta = 1 - 3 * Alpha;
const mpf_class revAlpha = 1 / Alpha;
const mpf_class revBeta = 1 / Beta;

using SeqType = std::vector<bool>;
const auto logAlpha = -log2f128(Alpha.get_d());
const auto logBeta = -log2f128(Beta.get_d());

std::map<int, int> meta_map;

//std::map<int, std::vector<bool>> bitmap;
using std::vector;

std::vector<bool> states2bits(std::vector<char> states, int transaction_size);

std::vector<char> bits2states(const std::vector<bool> &bits) {
    std::vector<char> states;
    states.reserve(bits.size() * 1.1);
    auto indicator = mpf_class(0, precision);
    // bits is little endian
    for (auto bit: bits) {
        indicator += bit;
        indicator /= 2;
    }
    std::cout << "put" << bits.size() << "-" << indicator.get_d() << std::endl;
    // step3: calculate next states
    std::array<int, 4> statistics = {};
    _Float128 entropy = 0;
    states.clear();
    // states is big endian
    while (entropy < bits.size() + 1) {
        auto cache = mpf_class(indicator * revAlpha);
        int multiple = cache.get_ui();
        if (multiple < 3) {
            states.push_back(multiple);
            indicator = cache - multiple;
            entropy += logAlpha;
            statistics[multiple]++;
        } else {
            states.push_back(3);
            indicator -= 3 * Alpha;
            indicator *= revBeta;
            entropy += logBeta;
            statistics[3]++;
        }
    }
    std::reverse(states.begin(), states.end());
    assert(states2bits(states, bits.size()) == bits);
    return states;
}

std::vector<bool> states2bits(std::vector<char> states, int transaction_size) {
    auto indicator = mpf_class(1, precision);
    // states is little endian
    for (auto x: states) {
        if (x < 3) {
            indicator += x;
            indicator *= Alpha;
        } else {
            indicator *= Beta;
            indicator += 3 * Alpha;
        }
    }
    std::cout << "fetch" << states.size() << "-" << indicator.get_d() << std::endl;
    std::vector<bool> bits;
    bits.reserve(transaction_size);
    for (auto i: Range(0, transaction_size)) {
        indicator *= 2;
        auto bit = indicator.get_ui();
        indicator -= bit;
        bits.push_back(bit);
    }
    std::reverse(bits.begin(), bits.end());
    // bits is big endian
    return bits;
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
    std::vector<char> states(initial_transaction_size, 3);
    for (int step = 0;; step++) {
        // step 1: fill info, calculate next indicator
        transaction_size = states.size();
        vector<bool> bits;

        N -= transaction_size;
        for (auto i: Range(states.size())) {
            auto index = N + i;
            auto state = states[i];
            auto bit_ans = seq_ans[index];
            auto bitA = !(state % 2) ^bit_ans;
            auto bitB = !(state / 2) ^bit_ans;
            seqA[index] = bitA;
            seqB[index] = bitB;
            bits.push_back(bitB);
        }

        // step 2: fill meta
        N -= META_SIZEOF;
        fill_meta(N, transaction_size);

        assert(bits.size() == transaction_size);
        states = bits2states(bits);

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

SeqType verify_b(int N, const SeqType &seqA, const SeqType &seq_std) {
    SeqType seqB(N);
    int ack = 0;
    auto setB = [&](int offset, bool value) {
        assert(ack == offset);
        seqB[offset] = value;
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
    std::vector<bool> bits;
    bits.reserve(iter - META_SIZEOF);
    for (auto i: Range(META_SIZEOF, iter)) {
        // including paddings for the first time
        bits.push_back(getA(i));
    }

    std::vector<char> states;
    while (iter < N) {
        auto transaction_size = fetch_meta(iter);
        if (states.size()) {
            bits = states2bits(std::move(states), transaction_size);
        }
        iter += META_SIZEOF;
        states.clear();
        states.reserve(transaction_size);
        for (auto i: Range(0, transaction_size)) {
            auto dst_index = iter + i;
            bool bitB = bits[i];
            setB(dst_index, bitB);
            auto bitA = getA(dst_index);
            auto std_bit = getStd(dst_index);
            bitB ^= std_bit;
            bitA ^= std_bit;
            states.push_back(3 - bitA - bitB * 2);
        }
        iter += transaction_size;
    }
    assert(iter == N);
    assert(ack == N);
    return seqB;
}


int wordload() {
    SeqType expected;
    std::default_random_engine e(42);
    for (auto i: Range(0, N)) {
        bool bit = e() % 2;
        expected.push_back(bit);
    }
    auto[a, b] = generate_a(N, expected);
    auto vb = verify_b(N, a, expected);

    std::array<int, 4> statistics = {};
    if (b != vb) {
        throw std::runtime_error("fuck");
    }
    for (auto i: Range(N)) {
        auto a_yes = expected[i] == a[i];
        auto b_yes = expected[i] == b[i];
        auto state = a_yes + b_yes * 2;
        statistics[state]++;
    }

    for (auto i: Range(0, 4)) {
        double count = statistics[i];
        std::cout << "state " << i << "->" << count << "->" << count / N << std::endl;
    }
    return 0;
}


void test_conversion() {
    std::default_random_engine e(42);
    int N = precision / 2;
    for (int step = 0; step < 1000; ++step) {
        SeqType expected;
        for (auto i: Range(0, N)) {
            bool bit = e() % 2;
            expected.push_back(bit);
        }
        auto states = bits2states(expected);
        auto bits = states2bits(states, N);
        assert(bits == expected);
    }
};

int main() {
//    test_conversion();
    wordload();
    return 0;
}


