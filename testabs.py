import numc as nc

def test_add():
    mat1 = nc.Matrix([[1, -2, 3,], [4, -5, 6], [7, -8, 9]]);
    absolute = nc.Matrix(3, 3)
    sol = nc.Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    absolute = abs(mat1)
    for i in range(3):
        for j in range(3):
            assert sol[i][j] == absolute[i][j], "Should be equal"


if __name__ == "__main__":
    test_add()
    print("Everything passed")
