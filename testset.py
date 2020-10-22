import numc as nc

def test_add():
    mat1 = nc.Matrix([[1, 1], [1, 1]]);
    mat1.set(0, 0, 3.1)
    assert mat1[0][0] == 3.1, "Should be equal"
    mat1.set(1,1, 3.455)
    assert mat1[1][1] == 3.455, "Should be equal"



if __name__ == "__main__":
    test_add()
    print("Everything passed")
