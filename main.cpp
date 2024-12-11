#include "TSP.hpp"

#include <iostream>


int main() {
    cost_matrix_t cm1 = {{INF, 10, 8,   19, 12},
                      {10, INF, 20,  6,  3},
                      {8,   20, INF, 4,  2},
                      {19,  6,  4, INF,  7},
                      {12,  3,  2,   7, INF}};

    /* Rozwiązania:
     * 32 : 2 3 4 1 0
     * 32 : 1 4 3 2 0
     */

    cost_matrix_t cm2 {
            {INF, 12,   3,  45,   6},
            {78, INF,  90,  21,   3},
            { 5,  56, INF,  23,  98},
            {12,   6,   8, INF,  34},
            { 3,  98,   3,   2, INF}
    };

    /* Rozwiązanie:
     * 30 : 4 3 2 0 1
    */

    cost_matrix_t cm3 {
            {INF,  3,  4,  2,  7},
            {3,  INF,  4,  6,  3},
            {4,  4,  INF,  5,  8},
            {2,  6,  5,  INF,  6},
            {7,  3,  8,  6,  INF},
    };

    /* Rozwiązania:
     * 19 : 4 3 0 2 1
     * 19 : 1 2 0 3 4
    */


//    std::cout << "petla do rozwiązywania" << std::endl;
//
//    tsp_solutions_t wynik = solve_tsp(cm1);
//
//    for(const auto& solution : wynik){
//        for (const auto& path : solution.path){
//            std::cout << " " << path  << ", ";
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << "koejny przykład" << std::endl;




    std::vector<std::vector<int>> mat1 {
            {INF, 12,   3,  45,   6},
            {78, INF,  90,  21,   3},
            { 5,  56, INF,  23,  98},
            {12,   6,   8, INF,  34},
            { 3,  98,   3,   2, INF}
    };
// rozw.: 2, 5, 4, 3, 1, 2

    std::vector<std::vector<int>> mat2 {
            {INF,  3,  4,  2,  7},
            {3,  INF,  4,  6,  3},
            {4,  4,  INF,  5,  8},
            {2,  6,  5,  INF,  6},
            {7,  3,  8,  6,  INF},
    };
// rozw.: 1, 3, 2, 5, 4, 1

    tsp_solutions_t wynik2 = solve_tsp(mat2);

    for(const auto& solution : wynik2){
        for (const auto& path : solution.path){
            std::cout << " " << path + 1  << ", ";
        }
        std::cout << "\n koszt: " << solution.lower_bound;
        std::cout << std::endl;
    }
    return EXIT_SUCCESS;
}
