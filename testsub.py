import numc as nc

def test_add():
    mat1 = nc.Matrix(3, 3, 1)
    mat2 = nc.Matrix(3, 3, 2)
    sub = nc.Matrix(3, 3)
    sol = nc.Matrix(3, 3, 1)
    sub = mat2 - mat1
    for i in range(3):
        for j in range(3):
            assert sol[i][j] == sub[i][j], "Should be equal"


if __name__ == "__main__":
    test_add()
    print("Everything passed")
