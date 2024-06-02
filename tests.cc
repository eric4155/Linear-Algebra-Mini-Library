#ifdef CATCH_CONFIG_MAIN
#  undef CATCH_CONFIG_MAIN
#endif
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "catch.hpp"
#include "matrix.hpp"


TEST_CASE("12x2 matrix", "[12x12 inverse]") {
    std::vector<double> data_non_identity = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1,
    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3,
    5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4,
    6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5,
    7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6,
    8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7,
    9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8,
    10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
    };

    Matrix test(12, 12, data_non_identity);
    Matrix identity(12,12);
    identity.SetToIdentity();
    REQUIRE(test.Invert() * test == identity);
}

TEST_CASE("12x12 matrix det", "[12x2 det]") {
    std::vector<double> data_non_identity = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1,
    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3,
    5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4,
    6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5,
    7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6,
    8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7,
    9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8,
    10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
    };

    Matrix test(12, 12, data_non_identity);
    REQUIRE(test.Determinant() == 4829554409472);
}


TEST_CASE("Matrix Inversion: 3x3 Matrices", "[matrix-inversion-3x3]") {
    std::vector<double> data1 = {
        2, -1, 0,
        -1, 2, -1,
        0, -1, 2
    };
    Matrix test1(3, 3, data1);
    Matrix test1_inverse = test1.Invert();
    Matrix identity1(3, 3);
    identity1.SetToIdentity();
    REQUIRE(test1 * test1_inverse == identity1);

    std::vector<double> data2 = {
        4, -2, 1,
        -2, 4, -2,
        1, -2, 4
    };
    Matrix test2(3, 3, data2);
    Matrix test2_inverse = test2.Invert();
    Matrix identity2(3, 3);
    identity2.SetToIdentity();
    REQUIRE(test2 * test2_inverse == identity2);

    std::vector<double> data3 = {
        3, 0, 2,
        2, 0, -2,
        0, 1, 1
    };
    Matrix test3(3, 3, data3);
    Matrix test3_inverse = test3.Invert();
    Matrix identity3(3, 3);
    identity3.SetToIdentity();
    REQUIRE(test3 * test3_inverse == identity3);
}

TEST_CASE("Matrix Inversion: 4x4 Matrices", "[matrix-inversion-4x4]") {
    std::vector<double> data1 = {
        4, 1, 2, 3,
        1, 3, 0, 1,
        2, 0, 2, 0,
        3, 1, 0, 3
    };
    Matrix test1(4, 4, data1);
    Matrix test1_inverse = test1.Invert();
    Matrix identity1(4, 4);
    identity1.SetToIdentity();
    REQUIRE(test1 * test1_inverse == identity1);

    std::vector<double> data2 = {
        2, 0, 1, 1,
        0, 3, 1, 2,
        1, 1, 4, 1,
        1, 2, 1, 3
    };
    Matrix test2(4, 4, data2);
    Matrix test2_inverse = test2.Invert();
    Matrix identity2(4, 4);
    identity2.SetToIdentity();
    REQUIRE(test2 * test2_inverse == identity2);

    std::vector<double> data3 = {
        1, 2, 0, 1,
        2, 4, 1, 2,
        0, 1, 2, 1,
        1, 2, 1, 3
    };
    Matrix test3(4, 4, data3);
    Matrix test3_inverse = test3.Invert();
    Matrix identity3(4, 4);
    identity3.SetToIdentity();
    REQUIRE(test3 * test3_inverse == identity3);
}

TEST_CASE("Matrix Inversion: 5x5 Matrices", "[matrix-inversion-5x5]") {
    std::vector<double> data2 = {
        2, 1, 0, 0, 0,
        1, 2, 1, 0, 0,
        0, 1, 2, 1, 0,
        0, 0, 1, 2, 1,
        0, 0, 0, 1, 2
    };
    Matrix test2(5, 5, data2);
    Matrix test2_inverse = test2.Invert();
    Matrix identity2(5, 5);
    identity2.SetToIdentity();
    REQUIRE(test2 * test2_inverse == identity2);
}

TEST_CASE("Matrix Inversion: 6x6 Matrices", "[matrix-inversion-6x6]") {
    std::vector<double> data2 = {
        3, 0, 2, 0, 0, 1,
        0, 2, 0, 1, 0, 0,
        2, 0, 1, 0, 1, 0,
        0, 1, 0, 2, 0, 0,
        0, 0, 1, 0, 3, 2,
        1, 0, 0, 0, 2, 3
    };
    Matrix test2(6, 6, data2);
    Matrix test2_inverse = test2.Invert();
    Matrix identity2(6, 6);
    identity2.SetToIdentity();
    REQUIRE(test2 * test2_inverse == identity2);
}

TEST_CASE("Matrix Determinant: 3x3 Matrices", "[matrix-determinant-3x3]") {
    std::vector<double> data1 = {
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    };
    Matrix test1(3, 3, data1);
    REQUIRE(test1.Determinant() == 0);

    std::vector<double> data2 = {
        3, 0, 2,
        2, 0, -2,
        0, 1, 1
    };
    Matrix test2(3, 3, data2);
    REQUIRE(test2.Determinant() == 10);

    std::vector<double> data3 = {
        6, 1, 1,
        4, -2, 5,
        2, 8, 7
    };
    Matrix test3(3, 3, data3);
    REQUIRE(test3.Determinant() == -306);
}

TEST_CASE("Matrix Determinant: 4x4 Matrices", "[matrix-determinant-4x4]") {
    std::vector<double> data1 = {
        1, 0, 2, -1,
        3, 0, 0, 5,
        2, 1, 4, -3,
        1, 0, 5, 0
    };
    Matrix test1(4, 4, data1);
    REQUIRE(test1.Determinant() == 30);

    std::vector<double> data2 = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };
    Matrix test2(4, 4, data2);
    REQUIRE(test2.Determinant() == 0);

    std::vector<double> data3 = {
        2, 0, 1, 1,
        0, 3, 0, 1,
        1, 0, 4, 1,
        1, 1, 1, 5
    };
    Matrix test3(4, 4, data3);
    // REQUIRE(test3.Determinant() == static_cast<double>(86));
}

TEST_CASE("Matrix Determinant: 5x5 Matrices", "[matrix-determinant-5x5]") {
    std::vector<double> data1 = {
        1, 2, 3, 4, 5,
        5, 4, 3, 2, 1,
        1, 0, 0, 1, 0,
        0, 1, 0, 0, 1,
        1, 1, 1, 1, 1
    };
    Matrix test1(5, 5, data1);
    REQUIRE(test1.Determinant() == 0);

    std::vector<double> data2 = {
        2, 0, 0, 1, 2,
        0, 3, 1, 0, 2,
        0, 1, 4, 0, 1,
        1, 0, 0, 5, 1,
        2, 2, 1, 1, 3
    };
    Matrix test2(5, 5, data2);
    // REQUIRE(test2.Determinant() == static_cast<double>(-36));
}

TEST_CASE("Matrix Determinant: 6x6 Matrices", "[matrix-determinant-6x6]") {
    std::vector<double> data1 = {
        1, 0, 0, 0, 0, 1,
        0, 1, 0, 0, 1, 0,
        0, 0, 1, 1, 0, 0,
        0, 0, 1, 1, 0, 0,
        0, 1, 0, 0, 1, 0,
        1, 0, 0, 0, 0, 1
    };
    Matrix test1(6, 6, data1);
    REQUIRE(test1.Determinant() == 0);

    std::vector<double> data2 = {
        3, 0, 2, 0, 0, 1,
        0, 2, 0, 1, 0, 0,
        2, 0, 1, 0, 1, 0,
        0, 1, 0, 2, 0, 0,
        0, 0, 1, 0, 3, 2,
        1, 0, 0, 0, 2, 3
    };
    Matrix test2(6, 6, data2);
    // REQUIRE(test2.Determinant() == static_cast<double>(-72));
} 


TEST_CASE("Matrix-Vector Multiplication: General cases", "[matrix-vector-multiplication]") {
    SECTION("Case 1: Identity matrix") {
        std::vector<double> mat_data = {
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        };
        Matrix identity_matrix(3, 3, mat_data);
        std::vector<double> vec_data = {1, 2, 3};
        Vector vector(vec_data);
        Vector result = identity_matrix * vector;
        std::vector<double> expected_result = {1, 2, 3};
        REQUIRE(result.GetData() == expected_result);
    }

    SECTION("Case 2: Zero matrix") {
        std::vector<double> mat_data(9, 0);
        Matrix zero_matrix(3, 3, mat_data);
        std::vector<double> vec_data = {1, 2, 3};
        Vector vector(vec_data);
        Vector result = zero_matrix * vector;
        std::vector<double> expected_result(3, 0);
        REQUIRE(result.GetData() == expected_result);
    }

    SECTION("Case 3: Non-square matrix") {
        std::vector<double> mat_data = {
            1, 2,
            3, 4,
            5, 6
        };
        Matrix matrix(3, 2, mat_data);
        std::vector<double> vec_data = {1, 1};
        Vector vector(vec_data);
        Vector result = matrix * vector;
        std::vector<double> expected_result = {3, 7, 11};
        REQUIRE(result.GetData() == expected_result);
    }

    SECTION("Case 4: Different dimensions") {
        std::vector<double> mat_data = {
            2, 4, 1,
            3, 1, 0
        };
        Matrix matrix(2, 3, mat_data);
        std::vector<double> vec_data = {1, 2, 3};
        Vector vector(vec_data);
        Vector result = matrix * vector;
        std::vector<double> expected_result = {13, 5};
        REQUIRE(result.GetData() == expected_result);
    }

    SECTION("Case 5: Large numbers") {
        std::vector<double> mat_data = {
            1e9, 2e9,
            3e9, 4e9
        };
        Matrix matrix(2, 2, mat_data);
        std::vector<double> vec_data = {1, 1};
        Vector vector(vec_data);
        Vector result = matrix * vector;
        std::vector<double> expected_result = {3e9, 7e9};
        REQUIRE(result.GetData() == expected_result);
    }

    SECTION("Case 6: Negative numbers") {
        std::vector<double> mat_data = {
            -1, -2,
            -3, -4
        };
        Matrix matrix(2, 2, mat_data);
        std::vector<double> vec_data = {1, -1};
        Vector vector(vec_data);
        Vector result = matrix * vector;
        std::vector<double> expected_result = {1, 1};
        REQUIRE(result.GetData() == expected_result);
    }

    SECTION("Case 7: Single element") {
        std::vector<double> mat_data = {4};
        Matrix matrix(1, 1, mat_data);
        std::vector<double> vec_data = {2};
        Vector vector(vec_data);
        Vector result = matrix * vector;
        std::vector<double> expected_result = {8};
        REQUIRE(result.GetData() == expected_result);
    }

    SECTION("Case 8: Mixed values") {
        std::vector<double> mat_data = {
            0.5, 1.5,
            2.5, 3.5
        };
        Matrix matrix(2, 2, mat_data);
        std::vector<double> vec_data = {2, 3};
        Vector vector(vec_data);
        Vector result = matrix * vector;
        std::vector<double> expected_result = {5.5, 15.5};
        REQUIRE(result.GetData() == expected_result);
    }
}


TEST_CASE("Matrix Transpose: General cases", "[matrix-transpose]") {
    SECTION("Case 1: Identity matrix") {
        std::vector<double> mat_data = {
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        };
        Matrix identity_matrix(3, 3, mat_data);
        Matrix result = identity_matrix.Transpose();
        REQUIRE(result == identity_matrix);
    }

    SECTION("Case 2: Zero matrix") {
        std::vector<double> mat_data(9, 0);
        Matrix zero_matrix(3, 3, mat_data);
        Matrix result = zero_matrix.Transpose();
        REQUIRE(result == zero_matrix);
    }

    SECTION("Case 3: Non-square matrix (3x2)") {
        std::vector<double> mat_data = {
            1, 2,
            3, 4,
            5, 6
        };
        Matrix matrix(3, 2, mat_data);
        Matrix result = matrix.Transpose();
        std::vector<double> expected_data = {
            1, 3, 5,
            2, 4, 6
        };
        Matrix expected_matrix(2, 3, expected_data);
        REQUIRE(result == expected_matrix);
    }

    SECTION("Case 4: Non-square matrix (2x3)") {
        std::vector<double> mat_data = {
            1, 2, 3,
            4, 5, 6
        };
        Matrix matrix(2, 3, mat_data);
        Matrix result = matrix.Transpose();
        std::vector<double> expected_data = {
            1, 4,
            2, 5,
            3, 6
        };
        Matrix expected_matrix(3, 2, expected_data);
        REQUIRE(result == expected_matrix);
    }

    SECTION("Case 5: Single element matrix") {
        std::vector<double> mat_data = {4};
        Matrix matrix(1, 1, mat_data);
        Matrix result = matrix.Transpose();
        REQUIRE(result == matrix);
    }

    SECTION("Case 6: Negative numbers") {
        std::vector<double> mat_data = {
            -1, -2,
            -3, -4
        };
        Matrix matrix(2, 2, mat_data);
        Matrix result = matrix.Transpose();
        std::vector<double> expected_data = {
            -1, -3,
            -2, -4
        };
        Matrix expected_matrix(2, 2, expected_data);
        REQUIRE(result == expected_matrix);
    }

    SECTION("Case 7: Mixed values (2x2)") {
        std::vector<double> mat_data = {
            0.5, 1.5,
            2.5, 3.5
        };
        Matrix matrix(2, 2, mat_data);
        Matrix result = matrix.Transpose();
        std::vector<double> expected_data = {
            0.5, 2.5,
            1.5, 3.5
        };
        Matrix expected_matrix(2, 2, expected_data);
        REQUIRE(result == expected_matrix);
    }

    SECTION("Case 8: Large matrix") {
        std::vector<double> mat_data = {
            1, 2, 3, 4,
            5, 6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15, 16
        };
        Matrix matrix(4, 4, mat_data);
        Matrix result = matrix.Transpose();
        std::vector<double> expected_data = {
            1, 5, 9, 13,
            2, 6, 10, 14,
            3, 7, 11, 15,
            4, 8, 12, 16
        };
        Matrix expected_matrix(4, 4, expected_data);
        REQUIRE(result == expected_matrix);
    }
}


TEST_CASE("Matrix Solve: Large systems", "[matrix-solve]") {
    SECTION("Case 1: 4x4 system") {
        std::vector<double> mat_data = {
            1, 0, 2, -1,
            3, 0, 0, 5,
            2, 1, 4, -3,
            1, 0, 5, 0
        };
        Matrix matrix(4, 4, mat_data);
        std::vector<double> vec_data = {1, 2, 3, 4};
        Vector vector(vec_data);
        std::vector<double> expected_solution = {static_cast<double>(-1) / 6, 1.5, static_cast<double>(5) / 6, 0.5};
        Vector result = matrix.Solve(vector);
        for (size_t i = 0; i < expected_solution.size(); ++i) {
            REQUIRE(result.GetData()[i] == Approx(expected_solution[i]).epsilon(1e-5));
        }
    }
}

TEST_CASE("Matrix Solve: Unique solution cases", "[matrix-solve]") {
    SECTION("Case 1: Identity matrix") {
        std::vector<double> mat_data = {
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        };
        Matrix identity_matrix(3, 3, mat_data);
        std::vector<double> vec_data = {1, 2, 3};
        Vector vector(vec_data);
        Vector result = identity_matrix.Solve(vector);
        REQUIRE(result.GetData() == vec_data);
    }

    SECTION("Case 2: Simple 2x2 system") {
        std::vector<double> mat_data = {
            2, 1,
            1, 3
        };
        Matrix matrix(2, 2, mat_data);
        std::vector<double> vec_data = {3, 4};
        Vector vector(vec_data);
        std::vector<double> expected_solution = {1, 1};
        Vector result = matrix.Solve(vector);
        REQUIRE(result.GetData() == expected_solution);
    }

    SECTION("Case 3: 3x3 system with negative numbers") {
        std::vector<double> mat_data = {
            -2, 1, 2,
            1, -2, 1,
            2, 1, -2
        };
        Matrix matrix(3, 3, mat_data);
        std::vector<double> vec_data = {-1, -1, 1};
        Vector vector(vec_data);
        std::vector<double> expected_solution = {static_cast<double>(-1) / 4, 0.0 , static_cast<double>(-3) / 4};
        Vector result = matrix.Solve(vector);
        REQUIRE(result.GetData() == expected_solution);
    }
}