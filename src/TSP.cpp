#include "TSP.hpp"
#include <algorithm>
#include <stack>
#include <optional>
#include <iostream>
#include <set>


std::ostream& operator<<(std::ostream& os, const CostMatrix& cm) { //obciązenie operatora aby łatwo wyświetlać macierz
    for (std::size_t r = 0; r < cm.size(); ++r) {
        for (std::size_t c = 0; c < cm.size(); ++c) {
            const auto& elem = cm[r][c];
            os << (is_inf(elem) ? "INF" : std::to_string(elem)) << " ";
        }
        os << "\n";
    }
    os << std::endl;

    return os;
}

/* PART 1 */

/**
 * Create path from unsorted path and last 2x2 cost matrix.
 * @return The vector of consecutive vertex.
 */
path_t StageState::get_path() {
    //pierwsza operacja na macierzy 2 na 2
    this->reduce_cost_matrix(); // redukowanie macierzy kosztów
    NewVertex new_vertex = this->choose_new_vertex(); // tworzenie nowej ścieżki
    this->update_cost_matrix(new_vertex.coordinates); //odnawianie macierzy na podstawi nowej ściezki
    this->append_to_path(new_vertex.coordinates); //dodanie współrzednych do unsorted path

    // druga operacja n amacierzy 1 na 1
    this->reduce_cost_matrix(); // redukcja macierzy 1 na 1
    NewVertex new_vertex2 = this->choose_new_vertex(); // towrzenie nowej ścieżki (ostatniej)
    this->update_cost_matrix(new_vertex2.coordinates); // aktualizacja macierzy kosztów juz powinna byc cała w infach
    this->append_to_path(new_vertex2.coordinates); // dodanie ostatnich kordynatów do unsorted path

    //sortowanie unsorted_paht_
    std::vector<std::pair<std::size_t, std::size_t>> path_of_pairs; // wektor par
    for (const auto& vertex : unsorted_path_){
        // dodajemy jeden aby operowac na indeksach mateamtycznych nie komputerowych
        path_of_pairs.emplace_back(vertex.row , vertex.col ); //emplace back towrzy element w miejscu wywołania, w tym wypadku std::pair
    }

    path_t path_result {path_of_pairs[0].first, path_of_pairs[0].second}; // pierwsze elementy trafiają do wektora wyniku

    path_of_pairs.erase(path_of_pairs.begin()); // elementy są juz w wektorze wynikowym wiec usuwamy go z wektora par
    while (path_result.size() != 5) {
        for (auto& pair: path_of_pairs) {
            // poczatek petli
            auto it1 = std::find(path_result.begin(), path_result.end(), pair.first); // tworzenie iteratora który odpowiada miejscu występowania lub nie elementu
            if (it1 == path_result.end()){ // nie znalazł 1.1
                auto it2 = std::find(path_result.begin(), path_result.end(), pair.second);
                if (it2 != path_result.end()){ // znalazł 2.2
                    path_result.insert(it2, pair.first); // wiec wrzucamy first przed znaleziony second
                    path_of_pairs.erase(std::remove(path_of_pairs.begin(), path_of_pairs.end(), pair), path_of_pairs.end());//usunieci nie potrzebnej juz pary
                }
            }else{ // znalazł 1
                auto it2 = std::find(path_result.begin(), path_result.end(), pair.second);
                if (it2 == path_result.end()) { // ale nie znalazł 2
                    path_result.insert(it1 + 1, pair.second); // wiec wstawiamy second za znaleziony first
                    path_of_pairs.erase(std::remove(path_of_pairs.begin(), path_of_pairs.end(), pair),
                                        path_of_pairs.end());//usunieci nie potrzebnej juz pary
                }
            }
        }
    }
    return path_result;
}

/**
 * Get minimum values from each row and returns them.
 * @return Vector of minimum values in row.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_rows() const {
    std::vector<cost_t> min_val(matrix_.size()); // wektor minimalnych wartości w danej kolumnie
    std::size_t counter = 0;
    for (const auto& row : matrix_) {  // pętla do iterowania po wierszach macierzy
        min_val[counter] = row[0]; // przyjmujemy ze na poczatku najmniejsza bedzie pierwsza liczba
        for (const auto& elem : row) {
            if((elem < min_val[counter]) && !(is_inf(elem))) {// jęsli element jest mniejszy od aktualnego na tej pozycji przepisujemy
                min_val[counter] = elem;
            }
        }
        if (is_inf(min_val[counter])) {
            min_val[counter] = 0;
        }
        counter++;
    }
    return min_val;
}

/**
 * Reduce rows so that in each row at least one zero value is present.
 * @return Sum of values reduced in rows.
 */
cost_t CostMatrix::reduce_rows() {
    std::vector<cost_t> min_vec_rows = get_min_values_in_rows();
    std::size_t counter = 0;
    cost_t sum = 0;
    for (auto& row : matrix_){
        for(auto& elem : row){
            if (!is_inf(elem)){
                elem -= min_vec_rows[counter];
                sum += min_vec_rows[counter];
            }
        }
        counter++;
    }
    return sum;
}

/**
 * Get minimum values from each column and returns them.
 * @return Vector of minimum values in columns.
 */
std::vector<cost_t> CostMatrix::get_min_values_in_cols() const {
    std::vector<cost_t> min_vec_cols(matrix_.size());
    std::size_t counter = 0;
    for(std::size_t cols = 0; cols < matrix_[0].size(); cols ++) {
        min_vec_cols[counter] = matrix_[0][cols];
        for(std::size_t rows = 0; rows < matrix_.size(); rows++) {
            if ((matrix_[rows][cols] < min_vec_cols[counter]) && !is_inf(matrix_[rows][cols])) {
                min_vec_cols[counter] = matrix_[rows][cols];
            }
        }
        if (is_inf(min_vec_cols[counter])) {
            min_vec_cols[counter] = 0;
        }
        counter ++;
    }
    return min_vec_cols;
}

/**
 * Reduces rows so that in each column at least one zero value is present.
 * @return Sum of values reduced in columns.
 */
cost_t CostMatrix::reduce_cols() {
    std::vector<cost_t> vec = get_min_values_in_cols();
    cost_t sum = 0;
    std::size_t fuse = 0;
    for (auto& row : matrix_){
        for(std::size_t c = 0; c < vec.size(); c++){
            if (!is_inf(row[c])) {
                row[c] -= vec[c];
                if (fuse == 0){
                    sum += vec[c];
                }
            }
        }
        fuse++;
    }
    return sum;
}

/**
 * Get the cost of not visiting the vertex_t (@see: get_new_vertex())
 * @param row
 * @param col
 * @return The sum of minimal values in row and col, excluding the intersection value.
 */
cost_t CostMatrix::get_vertex_cost(std::size_t row, std::size_t col) const {
    // Sprawdzenie rozmiaru macierzy (żeby upewnić się, że row i col są w granicach)
    if (row >= matrix_.size() || col >= matrix_[0].size()) {
        throw std::out_of_range("Indeks poza zakresem macierzy");
    }

    // Szukamy minimalnej wartości w kolumnie (pomijając punkt (row, col))
    cost_t acumin_col = std::numeric_limits<cost_t>::max();  // Inicjalizujemy do bardzo dużej liczby
    for (std::size_t r = 0; r < matrix_.size(); r++) {
        if (r != row) {  // Pomijamy wiersz `row`
            acumin_col = std::min(acumin_col, matrix_[r][col]);
        }
    }

    // Szukamy minimalnej wartości w wierszu (pomijając punkt (row, col))
    cost_t acumin_row = std::numeric_limits<cost_t>::max();  // Inicjalizujemy do bardzo dużej liczby
    for (std::size_t c = 0; c < matrix_[row].size(); c++) {
        if (c != col) {  // Pomijamy kolumnę `col`
            acumin_row = std::min(acumin_row, matrix_[row][c]);
        }
    }
    return acumin_col + acumin_row;
}

/* PART 2 */

/**
 * Choose next vertex to visit:
 * - Look for vertex_t (pair row and column) with value 0 in the current cost matrix.
 * - Get the vertex_t cost (calls get_vertex_cost()).
 * - Choose the vertex_t with maximum cost and returns it.
 * @param cm
 * @return The coordinates of the next vertex.
 */
NewVertex StageState::choose_new_vertex() {
    //szukanie wspólrzednych gdzie wartość = 0
    std::vector<path_t> pairs_of_0 ; //wektor wektorów z std::size_t
    for(std::size_t r = 0; r < matrix_.size(); r++){
        for(std::size_t c = 0; c < matrix_[0].size(); c++){
            if (matrix_[r][c] == 0){
                pairs_of_0.push_back({r, c});
            }
        }
    }
    // obliczanie kosztów par współrzednych
    path_t cost_of_pairs;
    cost_of_pairs.reserve(pairs_of_0.size()); // definicja z rezerwacją miejsca
    for (const auto& pair : pairs_of_0){
        cost_of_pairs.push_back(matrix_.get_vertex_cost(pair[0], pair[1]));//użycie funkcji get_vertex_cost która liczy koszt dla pary współrzednych
    }
    //szukanie największego kosztu i indexu które posłuży za identyfikacje dla których współrzednych ten koszt jest największy
    cost_t thelargest = cost_of_pairs[0];
    std::size_t index_of_large = 0; //definicja przed petlą aby później odnależć wspólrzędne
    for (std::size_t index = 0; index < cost_of_pairs.size(); index++){
        if (cost_of_pairs[index] > thelargest){
            thelargest = cost_of_pairs[index];
            index_of_large = index;
        }
    }
    //inicjalizacja obiketu struktury przechowującej współrzędne

    vertex_t vertex(pairs_of_0[index_of_large][0], pairs_of_0[index_of_large][1]);
    //dodanie vertexa do unsorted path

    //unsorted_path_.push_back(vertex);  nie trzeba bo w głównej pętli to juz jest

    //inicjalizacja pbiketu klasy przechowującej współrzędne i koszt
    NewVertex new_vertex(vertex, thelargest);
    return new_vertex;
}

//krok 3

/**
 * Update the cost matrix with the new vertex.
 * @param new_vertex
 */
void StageState::update_cost_matrix(vertex_t new_vertex) {
    //usuwanie wiersza
    //cost_matrix_t macierz = matrix_.get_matrix();
//    if (new_vertex.row < matrix_.size()){
//        matrix_.set_matrix().erase(matrix_.set_matrix().begin() + new_vertex.row);
//    }
//    //usuwanie kolumny
//    if (new_vertex.col < matrix_[0].size()) {
//        for (auto& row : matrix_.set_matrix()){
//            row.erase(row.begin() + new_vertex.col);
//        }
//    }
//    //blokowanie przejścia odwrotnego
    matrix_[new_vertex.col][new_vertex.row] = INF;

    if (new_vertex.row < matrix_.size() && new_vertex.col < matrix_[0].size()) {
        for (std::size_t c = 0; c < matrix_[0].size(); c++) {
            matrix_[new_vertex.row][c] = INF;
        }
        for (std::size_t r = 0; r < matrix_.size(); r++){
            matrix_[r][new_vertex.col] = INF;
        }
    }
    std::set<std::size_t> set_col {0,1,2,3,4};
    std::set<std::size_t> set_row {0,1,2,3,4};
    set_row.erase(new_vertex.row);
    set_col.erase(new_vertex.col);
    for (const auto& vertex : unsorted_path_){
        set_col.erase(vertex.col);
        set_row.erase(vertex.row);
    }

    for (std::size_t r = 0; r < matrix_.size(); r++){
        for (std::size_t c = 0; c < matrix_[0].size(); c++){
            if(!(is_inf(matrix_[r][c]))){
                if ((set_col.find(c) == set_col.end()) && (set_row.find(r) == set_row.end())){
                    matrix_[r][c] = INF;
                }
            }
        }
    }
}

/**
 * Reduce the cost matrix.
 * @return The sum of reduced values.
 */
cost_t StageState::reduce_cost_matrix() {
    cost_t sum = 0;
    for (const auto& elem : matrix_.get_min_values_in_cols()) {
        if (!( elem == 0 && is_inf(elem))){
            sum += matrix_.reduce_cols();
            break;
        }
    }
    for (const auto& elem : matrix_.get_min_values_in_rows()) {
        if (!( elem == 0 && is_inf(elem))){
            sum += matrix_.reduce_rows();
            break;
        }
    }
    return sum;
}

/**
 * Given the optimal path, return the optimal cost.
 * @param optimal_path
 * @param m
 * @return Cost of the path.
 */
cost_t get_optimal_cost(const path_t& optimal_path, const cost_matrix_t& m) {
    cost_t cost = 0;
    for (std::size_t idx = 1; idx < optimal_path.size(); ++idx) {
        cost += m[optimal_path[idx - 1]][optimal_path[idx]];
    }

    // Add the cost of returning from the last city to the initial one.
    cost += m[optimal_path[optimal_path.size() - 1]][optimal_path[0]];

    return cost;
}

/**
 * Create the right branch matrix with the chosen vertex forbidden and the new lower bound.
 * @param m
 * @param v
 * @param lb
 * @return New branch.
 */
StageState create_right_branch_matrix(cost_matrix_t m, vertex_t v, cost_t lb) {
    CostMatrix cm(m);
    cm[v.row][v.col] = INF;
    return StageState(cm, {}, lb);
}

/**
 * Retain only optimal ones (from all possible ones).
 * @param solutions
 * @return Vector of optimal solutions.
 */
tsp_solutions_t filter_solutions(tsp_solutions_t solutions) {
    cost_t optimal_cost = INF;
    for (const auto& s : solutions) {
        optimal_cost = (s.lower_bound < optimal_cost) ? s.lower_bound : optimal_cost;
    }

    tsp_solutions_t optimal_solutions;
    std::copy_if(solutions.begin(), solutions.end(),
                 std::back_inserter(optimal_solutions),
                 [&optimal_cost](const tsp_solution_t& s) { return s.lower_bound == optimal_cost; }
    );

    return optimal_solutions;
}

/**
 * Solve the TSP.
 * @param cm The cost matrix.
 * @return A list of optimal solutions.
 */
tsp_solutions_t solve_tsp(const cost_matrix_t& cm) {

    StageState left_branch(cm);

    // The branch & bound tree.
    std::stack<StageState> tree_lifo;

    // The number of levels determines the number of steps before obtaining
    // a 2x2 matrix.
    std::size_t n_levels = cm.size() - 2;

    tree_lifo.push(left_branch);   // Use the first cost matrix as the root.

    cost_t best_lb = INF;
    tsp_solutions_t solutions;

    while (!tree_lifo.empty()) { // pętla dziala dopóki stos jest pełny

        left_branch = tree_lifo.top(); //leftbranch nadpisywany górnym elmentem stosu
        tree_lifo.pop(); //usuwanie dodanego wyżej szczytu stosu do left branchu

        while (left_branch.get_level() != n_levels && left_branch.get_lower_bound() <= best_lb) {
            // Repeat until a 2x2 matrix is obtained or the lower bound is too high...

            if (left_branch.get_level() == 0) {
                left_branch.reset_lower_bound();
            }

            // 1. Reduce the matrix in rows and columns.
            cost_t new_cost = left_branch.reduce_cost_matrix();

            // 2. Update the lower bound and check the break condition.
            left_branch.update_lower_bound(new_cost);
            if (left_branch.get_lower_bound() > best_lb) {
                break;
            }

            // 3. Get new vertex and the cost of not choosing it.
            // @TODO (KROK 2)
            NewVertex new_vertex = left_branch.choose_new_vertex();

            // 4. @TODO Update the path - use append_to_path method.
            left_branch.append_to_path(new_vertex.coordinates);

            // 5. @TODO (KROK 3) Update the cost matrix of the left branch.
            left_branch.update_cost_matrix(new_vertex.coordinates);

            // 6. Update the right branch and push it to the LIFO.
            cost_t new_lower_bound = left_branch.get_lower_bound() + new_vertex.cost;
            tree_lifo.push(create_right_branch_matrix(cm, new_vertex.coordinates,
                                                      new_lower_bound));
        }

        if (left_branch.get_lower_bound() <= best_lb) {
            // If the new solution is at least as good as the previous one,
            // save its lower bound and its path.
            best_lb = left_branch.get_lower_bound();
            path_t new_path = left_branch.get_path();
            solutions.push_back({get_optimal_cost(new_path, cm), new_path});
        }
    }

    return filter_solutions(solutions); // Filter solutions to find only optimal ones.
}