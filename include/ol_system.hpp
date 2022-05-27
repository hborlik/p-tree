/**
 * @file ol_system.hpp
 * @brief 
 * @date 2022-02-07
 * 
 */
#ifndef OL_SYSTEM_HPP
#define OL_SYSTEM_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <math.h>

/**
 * @brief Parametric symbol
 * 
 * @tparam T 
 */
template<typename T>
struct Symbol {
    uint32_t RepSym = 0;
    T value{};

    constexpr Symbol(uint32_t Sym) : RepSym{Sym} {} // note not explicit
    constexpr Symbol(uint32_t Sym, T val) : RepSym{Sym}, value{std::move(val)} {}

    bool operator==(const Symbol<T>& o) noexcept {
        return RepSym == o.RepSym &&
            value == o.value;
    }
};

template<typename T>
using SymbolString = std::vector<Symbol<T>>;

/**
 * @brief Contextual neighborhood and symbol at a position in a SymbolString.
 * 
 *  neighborhood[1] is the symbol at the center of this context, should never be nullptr.
 *      neighborhood[0] and neighborhood[2] nullptr at boundary symbols.
 * 
 * @tparam T 
 */
template<typename T>
struct SymbolN {
    const Symbol<T>* neighborhood[3];

    const Symbol<T>* left() const noexcept {
        return neighborhood[0];
    } 

    const Symbol<T>* center() const noexcept {
        return neighborhood[1];
    }

    const Symbol<T>* right() const noexcept {
        return neighborhood[2];
    }

    bool matches(uint32_t a_l, uint32_t a, uint32_t a_r) const {
        // center is required
        return  (!neighborhood[0] || neighborhood[0]->RepSym == a_l || a_l == 0) &&
                !neighborhood[1] || neighborhood[1]->RepSym == a &&
                (!neighborhood[2] || neighborhood[2]->RepSym == a_r || a_r == 0);
    }
};

/**
 * @brief production rule for Symbol<T> with relative probability p
 * 
 * @tparam T 
 */
template<typename T>
struct Production {
    float p = 1.0f;
    const uint32_t A;

    Production(float p, uint32_t rule) : p{p}, A{rule} {}

    virtual bool matches(const SymbolN<T>& sym) const = 0;
    virtual SymbolString<T> translate(const SymbolN<T>& sym) const = 0;
};

/**
 * @brief Default Production rule, optionally context sensitive. Base rule consumes matching symbol
 *  mapping a_l < a > a_r → null
 * 
 * @tparam T 
 */
template<typename T>
struct DefaultProduction : public Production<T> {

    const uint32_t A_l, A_r;
    const uint32_t N;
    const Symbol<T> R;

    DefaultProduction(float p, uint32_t a, uint32_t N, Symbol<T> R) : Production<T>{p, a}, A_l{0}, A_r{0}, N{N}, R{R} {}

    /**
     * @brief Construct a new Context Sensitive Production object
     * 
     * @param a_l   left symbol id, 0 for none
     * @param a     strict predecessor id
     * @param a_r   right symbol id, 0 for none
     */
    DefaultProduction(float p, uint32_t a_l, uint32_t a, uint32_t a_r, uint32_t N, Symbol<T> R) : Production<T>{p, a}, A_l{a_l}, A_r{a_r}, N{N}, R{R} {}

    bool matches(const SymbolN<T>& sym) const override {
        return sym.matches(A_l, this->A, A_r);
    }

    SymbolString<T> translate(const SymbolN<T>& sym) const override {
        SymbolString<T> ret{};
        if (matches(sym)) {
            ret = SymbolString<T>(N, R);
        }
        return ret;
    }

};

template<typename T>
class LSystemTr {

    std::vector<const Production<T>*> productions{};
    std::unordered_map<uint32_t, std::vector<const Production<T>*>> symbol_mapper{};

public:
    void add_rule(const Production<T>* rule) {
        productions.push_back(rule);
        auto it = symbol_mapper.find(rule->A);
        if (it == symbol_mapper.end()) {
            symbol_mapper.insert({rule->A, {}});
        }
        symbol_mapper.at(rule->A).push_back(rule);
    }

    /**
     * @brief Evaluate a single step using system grammar
     * 
     * @param axiom 
     * @return SymbolString<T> 
     */
    SymbolString<T> evaluate(const SymbolString<T>& axiom) {
        SymbolString<T> B;
        const SymbolString<T>& A = axiom;
        std::vector<const Production<T>*> valid_rules{};
        SymbolString<T> tr{};
        for (int i = 0; i < A.size(); ++i) {
            const SymbolN<T> n {
                {i == 0 ? nullptr : &(A[i - 1]),
                &(A[i]),
                i < A.size() - 1 ? &(A[i + 1]) : nullptr}
            };

            auto itr = symbol_mapper.find(n.neighborhood[1]->RepSym);
            float p_total = 0.0f;
            if (itr != symbol_mapper.end()) {
                auto rules_for_sym = (*itr).second;
                const Production<T>* selected_rule = rules_for_sym[0];
                // have a rule set for this symbol, random weighted pick based on valid rules
                if (rules_for_sym.size() > 1) {
                    float w = rand() / float(RAND_MAX); // [0.0, 1.0]
                    valid_rules.clear();
                    for (const auto rule : rules_for_sym) {
                        if (rule->matches(n)) {
                            valid_rules.push_back(rule);
                            p_total += rule->p;
                        }
                    }
                    w *= p_total; // [0.0, p_total]
                    for (const auto rule : valid_rules) {
                        if (w <= rule->p) { // inclusive
                            selected_rule = rule;
                            break;
                        }
                        w -= rule->p;
                    }
                }
                tr = selected_rule->translate(n);
                B.insert(B.end(), tr.begin(), tr.end());
            } else { // no matching productions, symbol does not change
                B.push_back(*n.neighborhood[1]);
            }
        }
        return B;
    }
};

template<typename T>
void PrintSymbolString(const char* dict[], const SymbolString<T>& str) {
    for (const auto& s : str) {
        std::ostringstream out;
        if (s.value != 0.0f) {
            out.precision(2);
            out << "(";
            out << s.value;
            out << ")";
        }

        std::cout << dict[s.RepSym] << out.str() << " ";
    }
    std::cout << std::endl;
}

#endif // OL_SYSTEM_HPP