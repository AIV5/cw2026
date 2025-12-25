#define GET(X, K) (((X) >> (K)) & 1u)
#define SET(X, K, V) (X) |= (V) << (K)
#define NEXT(K, X, B) (X ^= B[__builtin_ctz(++K)])

row_t bin_gauss(row_t * restrict augmented_matrix, int nrows,
        row_t * restrict nsolutions, row_t * restrict kernel_basis) {
    int column_order[COLS];
    int leading_columns = 0;
    for (int k = 0; k < COLS; ++k) {
        if (leading_columns < nrows && !GET(augmented_matrix[leading_columns], k)) {
            for (int i = leading_columns + 1; i < nrows; ++i) {
                if (GET(augmented_matrix[i], k)) {
                    augmented_matrix[leading_columns] ^= augmented_matrix[i];
                    break;
                }
            }
        }
        if (leading_columns == nrows || !GET(augmented_matrix[leading_columns], k)) {
            column_order[COLS-1-k + leading_columns] = k;
            continue;
        }
        column_order[leading_columns] = k;
        for (int i = leading_columns + 1; i < nrows; ++i) {
            if (GET(augmented_matrix[i], k)) {
                augmented_matrix[i] ^= augmented_matrix[leading_columns];
            }
        }
        ++leading_columns;
    }
    for (int i = leading_columns; i < nrows; ++i) {
        if (GET(augmented_matrix[i], COLS)) {
            *nsolutions = 0;
            return 0;
        }
    }
    int kernel_dimension = COLS - leading_columns;
    *nsolutions = 1u << kernel_dimension;
    for (int i = leading_columns - 1; i > 0; --i) {
        int k = column_order[i];
        for (int j = 0; j < i; ++j) {
            if (GET(augmented_matrix[j], k)) {
                augmented_matrix[j] ^= augmented_matrix[i];
            }
        }
    }
    row_t specific_solution = 0;
    for (int i = 0; i < leading_columns; ++i) {
        int k = column_order[i];
        SET(specific_solution, k, GET(augmented_matrix[i], COLS));
    }
    for (int i = 0; i < kernel_dimension; ++i) {
        int k = column_order[leading_columns + i];
        SET(kernel_basis[i], k, 1u);
        for (int j = 0; j < leading_columns; ++j) {
            SET(kernel_basis[i], column_order[j], GET(augmented_matrix[j], k));
        }
    }
    return specific_solution;
}

