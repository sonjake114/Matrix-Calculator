# numc

Here's what I did in project 4:
I used basic row-major matrix. so the row*width+col of index would give (row, col) of a certain 2d matrix. For simple add, subtract, neg, abs functions I just integrated the data portion of matrix using the simple iteration and preformed operations each element one by one. For more complicated operation such as mul and pow. I used the vector dot product to multiply each row and col. Also, I used transpose of matrix to improve the speed because the cache blocking will improve the speed of algorithm. In terms of speed up, I majorly used unrolled loop. So I unrolled every major loop I can find using the method we were taught in lab. And I used OpenMP Parallel computation using 4 corse with hyperthreads.

I tried implementing the SIMD instructions using loop unrolling and just normal AIVX. However, I had extremely hard time debugging this because I was missing rows. Sometimes the values would simply just not computed. I tried my best but I think I ran out of time even after using all three slip days.

I was really suprised by the transposing of matrix can give to cache hitting ration. I remember in lab 9 that transposing using matrix block caheing can improve the performance by significant margin. However, I was really suprised by how many extra points I got just from implementing the cache block.  

Also, I was really suprised by how little efficient algorithms for matrix calculations are there. I tried my absolute best to look for different sources that mention about the way they improved their matrix calculation speed. But basically the most efficient way was cache blocking with different cache block sizes. Although I did come across this Strassen algorithm that claims it could do the matrix multiplication in best cast at O(n^2.374555). This is not necessarily significant improvement from the previous O(n^3), so I did not bother thinking about implementing it using that algorithm.

If I had time, I would love to build a website or some sort of more of a user friendly ways of presenting the matrix operations. I am really interested in app developments as well, such as unity and android studio, really gets my attention. I would love to make an app or website that allows user to input their on row and cols to make matrix and does various operations at it's own.

Also, I think there is a lot of cases where parallel programming would benefit a lot of perforamce. For example, the project that we did in CS61B using signpost is a great example. We could build the board very efficiently using the parallel optimization because everytime, board does the same thing of just getting n*n board. Also, the 3rd project that we did. tablut project using the knights and capture methods. I think it would be great in terms of performace to use parallelization in those projects. Thank you for the great project! It was nice to play with various performance enhancement but also very challenging.


-
