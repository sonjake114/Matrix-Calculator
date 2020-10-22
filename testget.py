import numc as nc

def test_add():
    mat1 = nc.Matrix([[1.2234, 1.6777], [133434, 1253535]]);
    val = mat1.get(0, 0)
    assert val == 1.2234, "Should be equal"
    val2 = mat1.get(1,1)
    assert val2 == 1253535, "Should be equal"



if __name__ == "__main__":
    test_add()
    print("Everything passed")
