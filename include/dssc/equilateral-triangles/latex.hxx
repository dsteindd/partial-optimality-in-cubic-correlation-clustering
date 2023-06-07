#pragma once
#ifndef DSSC_EQUILATERAL_LATEX_HXX
#define DSSC_EQUILATERAL_LATEX_HXX

#include "dssc/math.hxx"
#include "dssc/latex.hxx"

#include "problem.hxx"

namespace dssc {
    namespace latex {

        template<class R>
        void
        printProblem(
                std::ostream &out,
                ProblemEquilateralTriangles<R> const &problem,
                bool const printWithTruth
        ) {
            typedef Point<R, 2> Point2;

            out << R"(
        \documentclass[tikz]{standalone}
        \usepackage{tikz}

        \begin{document}
        \begin{tikzpicture}[scale=2]
            % \draw[->] (0, 0) -- (1.5, 0) node[anchor=east] {$x$};
            % \draw[->] (0, 0) -- (0, 1.5) node[anchor=north] {$y$};
    )";

            if (printWithTruth) {
                size_t const &numTriangles = problem.triangles().size();
                size_t const &numPoints = problem.points().size();

                // print triangles
                for (size_t j = 0; j < numTriangles; ++j) {
                    // ground truth triangle
                    auto const &triangle = problem.triangles()[j];

                    Point2 point1 = triangle.points()[0];
                    Point2 point2 = triangle.points()[1];
                    Point2 point3 = triangle.points()[2];

                    // triangle ground truth
                    out << "\\draw[" << color(j) << "!40] ("
                        << point1[0] << ", " << point1[1]
                        << ") -- ("
                        << point2[0] << ", " << point2[1]
                        << ") -- ("
                        << point3[0] << ", " << point3[1]
                        << ") -- cycle;"
                        << std::endl;
                }
                for (size_t j = 0; j < numPoints; ++j) {
                    Point2 point1 = problem.points()[j];
                    size_t label = problem.labels()[j];


                    out << "      \\draw plot[mark=*, mark size=0.1ex, mark options={draw=" << color(label) << ", fill="
                        << color(label) << "}] coordinates {("
                        << point1[0] << ", " << point1[1]
                        << ")};" << std::endl;
                }


            } else {
                // print points
                for (size_t j = 0; j < problem.points().size(); ++j) {
                    auto const &point = problem.points()[j];
//                    auto const label = problem.labels()[j];
                    out
                            << "      \\draw plot[mark=*, mark size=0.1ex, mark options={draw=gray, fill=gray}] coordinates {("
                            << point[0] << ", " << point[1]
                            << ")};" << std::endl;
                }
            }

            out << R"(
        \end{tikzpicture}
        \end{document}
    )";
        }

        template<class R>
        inline
        void
        printProblem(
                std::string const &fileName,
                ProblemEquilateralTriangles<R> const &problem,
                bool const printWithTruth
        ) {
            std::ofstream fileStream;
            fileStream.open(fileName, std::ios::out);
            if (fileStream) {
                printProblem(fileStream, problem, printWithTruth);
            }
            fileStream.close();
        }

    }
}

#endif
