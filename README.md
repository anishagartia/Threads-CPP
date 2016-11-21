# Threads-CPP

In this project, we are going to compute a two–dimensional Discrete Fourier Transform of a given input image using threads. For this, we first compute the one–dimensional transforms on eacy row, followed by the one–dimensional transforms on each column. Following are the specifications of the implementations:

1. We will use the efficient Danielson–Lanczos approach for the one–dimensional transforms. This approach has a running time proportional to Nlog2N as opposed to the N2 running time of the non-threaded algorithm.
2. We are going to use 16 threads, plus the “main” thread, on a single computing platform to cooperate to perform the two–dimensional transform. All threads will have access to the same memory locations (ie. the original 2d image).
3. The main thread will create each of the 16 helper threads and pass the thread id as the argument to the thread function. Each thread will perform a 1D transform on a set of rows (identical to the way we did this with MPI). They should then barrier and then do the second 1D transform on the columns. There are a few ways to implement this barrier. here we use send and wait function of PThread.

## Steps of the Danielson–Lanczos algorithm: 
1. Reordering the initial h array to bit reversed order.
2. Precomputing the Weights for each stage of butterfly diagram.
3. Performing the transofrmation. 

We also implement the inverse fourier transform of the image.
