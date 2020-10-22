import numc as nc

def test_add():
    mat1 = nc.Matrix([[1.1, 2.2 , 3.3  , 4.4 ], [5.5, 6.6 , 7.7 ,8.8], [ 9.9 , 10.10 ,11.11  ,12.12 ], [13.13  , 14.14  ,  15.15 , 16.16 ]]);
    absolute = nc.Matrix(4, 4)
    sol = nc.Matrix([[103.752, 112.486 , 123.893  , 135.3 ], [234.124, 257.862 , 287.837 ,317.812 ], [ 335.5646 , 372.0278 ,417.4901  ,462.9524  ], [454.3788 , 503.7274  , 565.3475 , 626.9676 ]]);
    absolute = mat1**2
    print(absolute)

def test_ezy():
    mat2 = nc.Matrix([[1, 2, 3,], [1, 2, 3], [1, 2, 3]]);
    solmat = nc.Matrix(3,3)
    solmat = mat2**2
    for i in range(3):
        for j in range(3):
            print(solmat[i][j])
        print("")

def test_ezy1():
    mat2 = nc.Matrix([[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3 ,4], [1, 2, 3, 4]]);
    solmat = nc.Matrix(4,4)
    solmat = mat2**2
    for i in range(4):
        for j in range(4):
            print(solmat[i][j])
        print("")

def test_ezy2():
    mat2 = nc.Matrix([[1, 2, 3, 4 ,5 ,6 ,7 ,8], [1, 2, 3, 4 ,5, 6, 7, 8], [1, 2, 3 ,4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3 ,4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8]]);
    solmat = nc.Matrix(8,8)
    solmat = mat2**2
    for i in range(8):
        for j in range(8):
            print(solmat[i][j])
        print("")



if __name__ == "__main__":
    test_ezy2()
    print("Everything passed")
