#pragma once
#ifndef DSSC_LATEX_HXX
#define DSSC_LATEX_HXX

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>

namespace dssc {
namespace latex {

inline
std::string
color(
    size_t const j
) {
    switch(j) {
        case 0: return "red";
        case 1: return "green";
        case 2: return "blue";
        case 3: return "cyan";
        case 4: return "magenta";
        case 5: return "brown";
        default: return "gray";
    }
}

template<class PROBLEM, class PROBLEM_ADAPTOR>
void
printEdgeCostDistribution(
        std::ostream & out,
        PROBLEM const & problem,
        PROBLEM_ADAPTOR const & problemAdaptor
) {
    typedef typename PROBLEM::Rational R;

    // collect costs
    std::vector<R> costsTrueJoins;
    std::vector<R> costsTrueCuts;
    for(size_t j = 0; j < problemAdaptor.numberOfElements(); ++j)
        for(size_t k = 0; k < j; ++k) {
            auto const cost = problemAdaptor.costOfEdge({j, k});
            if (problem.labels()[j] == problem.labels()[k]) {
                costsTrueJoins.push_back(cost);
            } else {
                costsTrueCuts.push_back(cost);
            }
        }

    std::sort(costsTrueJoins.begin(), costsTrueJoins.end());
    std::sort(costsTrueCuts.begin(), costsTrueCuts.end());
    size_t const numberOfCoefficients = costsTrueJoins.size() + costsTrueCuts.size();

    // print
    out << R"(
    \documentclass[tikz]{standalone}
    \usepackage{pgfplots}
        \pgfplotsset{compat=newest}
    \begin{document}
    \begin{tikzpicture}
    \begin{axis}[xlabel={Cost}, ylabel={Fraction of coefficients}, legend pos={north west}]
        \addplot[green] table[x index=0, y index=1, col sep=comma] {
)";

    for(size_t j = 0; j < costsTrueJoins.size(); ++j) {
        out << costsTrueJoins[j]
            << ", " << static_cast<double>(j+1) / numberOfCoefficients << std::endl;
    }

    out << R"(
    };
    \addplot[red] table[x index=0, y index=1, col sep=comma] {
)";

    for(size_t j = 0; j < costsTrueCuts.size(); ++j) {
        out << costsTrueCuts[j]
            << ", " << static_cast<double>(j+1) / numberOfCoefficients << std::endl;
    }

    out << R"(
    };
        \legend {true join, true cut}
    \end{axis}
    \end{tikzpicture}
    \end{document}
)";
}

template<class PROBLEM, class PROBLEM_ADAPTOR>
void
printEdgeCostDistribution(
        std::string const & fileName,
        PROBLEM const & problem,
        PROBLEM_ADAPTOR const & problemAdaptor
) {
    std::ofstream fileStream;
    fileStream.open(fileName, std::ios::out);
    if(fileStream) {
        printEdgeCostDistribution(fileStream, problem, problemAdaptor);
    }
    fileStream.close();
}

template<class PROBLEM, class PROBLEM_ADAPTOR>
void
printCostDistribution(
    std::ostream & out,
    PROBLEM const & problem,
    PROBLEM_ADAPTOR const & problemAdaptor
) {
    typedef typename PROBLEM::Rational R;

    // collect costs
    std::vector<R> costsTrueJoins;
    std::vector<R> costsTrueCuts;
    for(size_t j = 0; j < problemAdaptor.numberOfElements(); ++j)
    for(size_t k = 0; k < j; ++k)
    for(size_t l = 0; l < k; ++l) {
        auto const cost = problemAdaptor.costOfTriple({j, k, l});
        if(problem.labels()[j] == problem.labels()[k]
        && problem.labels()[k] == problem.labels()[l]) {
            costsTrueJoins.push_back(cost);
        }
        else{
            costsTrueCuts.push_back(cost);
        }
    }
    std::sort(costsTrueJoins.begin(), costsTrueJoins.end());
    std::sort(costsTrueCuts.begin(), costsTrueCuts.end());
    size_t const numberOfCoefficients = costsTrueJoins.size() + costsTrueCuts.size();

    // print
    out << R"(
        \documentclass[tikz]{standalone}
        \usepackage{pgfplots}
            \pgfplotsset{compat=newest}
        \begin{document}
        \begin{tikzpicture}
        \begin{axis}[xlabel={Cost}, ylabel={Fraction of coefficients}, legend pos={north west}]
            \addplot[green] table[x index=0, y index=1, col sep=comma] {
    )";

    for(size_t j = 0; j < costsTrueJoins.size(); ++j) {
        out << costsTrueJoins[j]
            << ", " << static_cast<double>(j+1) / numberOfCoefficients << std::endl;
    }

    out << R"(
        };
        \addplot[red] table[x index=0, y index=1, col sep=comma] {
    )";

    for(size_t j = 0; j < costsTrueCuts.size(); ++j) {
        out << costsTrueCuts[j]
            << ", " << static_cast<double>(j+1) / numberOfCoefficients << std::endl;
    }

    out << R"(
        };
            \legend {true join, true cut}
        \end{axis}
        \end{tikzpicture}
        \end{document}
    )";
}

template<class PROBLEM, class PROBLEM_ADAPTOR>
void
printCostDistribution(
    std::string const & fileName,
    PROBLEM const & problem,
    PROBLEM_ADAPTOR const & problemAdaptor
) {
    std::ofstream fileStream;
    fileStream.open(fileName, std::ios::out);
    if(fileStream) {
        printCostDistribution(fileStream, problem, problemAdaptor);
    }
    fileStream.close();
}

}
}

#endif
