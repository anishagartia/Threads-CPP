// Threaded two-dimensional Discrete FFT transform
// YOUR NAME HERE
// ECE8893 Project 2


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <math.h>
#include <cmath>
#include <sstream>
#include <fstream>
#include <pthread.h>
#include <sys/time.h>
#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;

//GLOBAL VARIABLES visible to all threads
int N;
Complex *W;
Complex *h;
int ImageHeight;
int ImageWidth;
int nThreads = 16;

//pthread_t pt[16];

pthread_mutex_t startCountMutex;
pthread_mutex_t exitMutex;
pthread_mutex_t coutMutex;
pthread_cond_t exitCond;
int startCount;

// -----------------------------------------------------------
// Function: TransposeIm
// Parameters: 
// h - Reference to Input image matrix of type Complex
// w - width of input image referenced by H
// ht - height of input image referenced by H
// -----------------------------------------------------------
 
void TransposeIm(Complex* h, int w, int ht)
{
  for (int r = 0; r < ht; r++){
    for (int c = 0; c < r; c++){
      Complex temp = h[(r * w) + c];
      h[(r * w) + c]  = h[(c * w) + r];
      h[(c * w) + r]  = temp;
    }
  }
}

// -----------------------------------------------------------
// Function: ReverseBits - to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the number of points in the 1D transform.
// Parameters: 
// v - number whose bits is to be reversed
// returned value: number with reversed bits
// -----------------------------------------------------------
 
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// -----------------------------------------------------------
// Function: MyBarrierInit
// Initialises global mutexes used for the barrier
// -----------------------------------------------------------
 
// Call MyBarrier_Init once in main
void MyBarrier_Init()// you will likely need some parameters)
{ 
  // All mutex and condition variables must be "initialized"
  pthread_mutex_init(&exitMutex,0);
  pthread_mutex_init(&coutMutex,0);
  pthread_mutex_init(&startCountMutex,0);
  pthread_cond_init(&exitCond, 0);
}

// -----------------------------------------------------------
// Function: MyBarrier
// Implements barrier using pthread condition wait, which wats for the exit conditon, and exit mutex.
// -----------------------------------------------------------
 
// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier() // Again likely need parameters
{
  pthread_cond_wait(&exitCond, &exitMutex);
}
                    
// -----------------------------------------------------------
// Function: TransformID
// Performs 1D transform of the row of length N pointed by 
// complex pointer h
// -----------------------------------------------------------
 
void Transform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)

  for(int i=2; i<= N; i = i*2){
    //cout << "here 4_i_" << i << endl;
    for(int k=0; k<(i/2); k++){
      //cout << "here 4_K_" << k << endl;
      Complex W_this = W[k*(N/i)];
      //cout << "W " <<  k*(N/i) << endl;
      for(int j=0; j< (N/i); j++){
        Complex ans_p = h[(j*i)+k] + (W_this * h[(j*i)+k+(i/2)]);  
        //cout << "here 4_4_1"<< endl;
        Complex ans_n = h[(j*i)+k] - (W_this * h[(j*i)+k+(i/2)]);
        //cout << "here 4_4_2" << endl;
        h[(j*i)+k] = ans_p;
        //cout << "here 4_4_3"  << endl;
        h[(j*i)+k+(i/2)] = ans_n;
      }
    }
  }
}



// -----------------------------------------------------------
// Function: InverseTransformID
// Performs 1D Inversetransform of the row of length N pointed by 
// complex pointer h
// -----------------------------------------------------------
 
void InverseTransform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)

  for(int i=2; i<= N; i = i*2){
    //cout << "here 4_i_" << i << endl;
    for(int k=0; k<(i/2); k++){
      //cout << "here 4_K_" << k << endl;
      Complex W_this = W[k*(N/i)];
      //cout << "W " <<  k*(N/i) << endl;
      for(int j=0; j< (N/i); j++){
        Complex ans_p = h[(j*i)+k] + (W_this * h[(j*i)+k+(i/2)]);  
        //cout << "here 4_4_1"<< endl;
        Complex ans_n = h[(j*i)+k] - (W_this * h[(j*i)+k+(i/2)]);
        //cout << "here 4_4_2" << endl;
        h[(j*i)+k] = ans_p;
        //cout << "here 4_4_3"  << endl;
        h[(j*i)+k+(i/2)] = ans_n;
      }
    }
  }
  Complex byN((1.0/N),0);
  for (int d = 0; d<N; d++){
    h[d] = byN * h[d];
  }
}

// -----------------------------------------------------------
// Function: Transform2DTHread
// Starting point for thread
// Parameter: v is the thread number
// This function calls Transform1D to perform 1D transform of 
// ImageWidth/nThreads rows. 
// -----------------------------------------------------------
 
void* Transform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete

  unsigned long myID = (unsigned long)v;
  int rowsPerThread = ImageHeight/nThreads;
  int startingRow = myID * rowsPerThread;
  
  pthread_mutex_lock(&coutMutex);
  cout << "Executing thread id_"<< myID << endl; 
  pthread_mutex_unlock(&coutMutex);

  // Perform 1D of each row
  for (int r=0; r < rowsPerThread; r++){
    int thisRow = startingRow + r;

    pthread_mutex_lock(&coutMutex);
    //cout << "thread id_"<< myID << " thisRow_" << thisRow << endl; 
    pthread_mutex_unlock(&coutMutex);

    Transform1D( h+(thisRow*N), N );    
  }    
  
  // Send conditional signals
  pthread_mutex_lock(&startCountMutex);
  startCount--;
  if (startCount ==0){
    // Last to exit, notify main
    pthread_mutex_unlock(&startCountMutex);
    pthread_mutex_lock(&exitMutex);
    pthread_cond_signal(&exitCond);
    pthread_mutex_unlock(&exitMutex);
    pthread_mutex_lock(&coutMutex);
    cout << "Exiting all threads normally" << endl;
    pthread_mutex_unlock(&coutMutex);
  }
  else {
    pthread_mutex_unlock(&startCountMutex);
    pthread_mutex_lock(&coutMutex);
    cout << "Exiting thread id " << myID << " normally" << endl;
    pthread_mutex_unlock(&coutMutex);
  }
  return 0;
}



// -----------------------------------------------------------
// Function: InverseTransform2DTHread
// Starting point for thread
// Parameter: v is the thread number
// This function calls Transform1D to perform 1D transform of 
// ImageWidth/nThreads rows. 
// -----------------------------------------------------------
 
void* InverseTransform2DTHread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete

  unsigned long myID = (unsigned long)v;
  int rowsPerThread = ImageHeight/nThreads;
  int startingRow = myID * rowsPerThread;
  
  pthread_mutex_lock(&coutMutex);
  cout << "Executing thread id_"<< myID << endl; 
  pthread_mutex_unlock(&coutMutex);

  // Perform 1D of each row
  for (int r=0; r < rowsPerThread; r++){
    int thisRow = startingRow + r;

    pthread_mutex_lock(&coutMutex);
    //cout << "thread id_"<< myID << " thisRow_" << thisRow << endl; 
    pthread_mutex_unlock(&coutMutex);

    InverseTransform1D( h+(thisRow*N), N );    
  }    
  
  // Send conditional signals
  pthread_mutex_lock(&startCountMutex);
  startCount--;
  if (startCount ==0){
    // Last to exit, notify main
    pthread_mutex_unlock(&startCountMutex);
    pthread_mutex_lock(&exitMutex);
    pthread_cond_signal(&exitCond);
    pthread_mutex_unlock(&exitMutex);
    pthread_mutex_lock(&coutMutex);
    cout << "Exiting all threads normally" << endl;
    pthread_mutex_unlock(&coutMutex);
  }
  else {
    pthread_mutex_unlock(&startCountMutex);
    pthread_mutex_lock(&coutMutex);
    cout << "Exiting thread id " << myID << " normally" << endl;
    pthread_mutex_unlock(&coutMutex);
  }
  return 0;
}


// -----------------------------------------------------------
// Function: Transform2D
// - Reads image into array
// - Reorders the row eements to reversedBits order,
// - Creates the threads(to perform 1D trasnform of rows), 
// - calls barrier to wait till all threads have executed,
// - transposes the matrix,
// - creates threads again(to perform 1D transform of rows),
// - calls barrier to wait till all threads ahve executed,
// - transposes matrix,
// - saves it to file.
// -----------------------------------------------------------
 

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Create the global pointer to the image array data
  // Create 16 threads
  // Create ordered array - bit reversed in place
  // Wait for all threads complete
  // Write the transformed data
  int w = image.GetWidth();
  int ht = image.GetHeight();
  h = image.GetImageData();
  N = w;
  ImageWidth = w;
  ImageHeight = ht;
  W = new Complex[N/2];
  unsigned rev_ind[N];
  cout <<"w " <<  w << endl;
  cout << "ht " << ht << endl;  

  // Create array rev_ind with correct order of input indices
  for (int i =0; i < w; i++){
    rev_ind[i] = ReverseBits(i);
    //cout << rev_ind[i] << endl;
  }
  //cout << "here 1" << endl;

  // Run through each row and rearrnge elements to correct order (bit reversed)
  for (int r = 0; r < ht; r++){
    for (int c = 0; c < w; c++){
      if (rev_ind[c] <= (unsigned)c){
        continue;
      }
      else {
        Complex temp;
        temp = h[(r*w)+c];
        h[(r*w) + c] = h[(r*w)+rev_ind[c]];
        h[(r*w)+rev_ind[c]] = temp;
      }
    }
  }
  cout << "here 2" << endl;
  // Create array of weights for size w/2
  for (int n =0; n < (w/2); n++){
    W[n].real = cos(2*M_PI*n/w);
    W[n].imag = -sin(2*M_PI*n/w);
  }
  //cout << "here 3" << endl;

  // Create Threads (FIRST)
  startCount = nThreads;
  for (int i =0; i<nThreads; i++){
    pthread_t pt;

    //pthread_mutex_lock(&coutMutex);
    //cout << "here 3_in_"<< (void*)i << endl;
    //pthread_mutex_unlock(&coutMutex);

    pthread_create(&pt, NULL, Transform2DTHread, (void*)i);

    //pthread_mutex_lock(&coutMutex);
    //cout << "here 3_out_" << i << endl;
    //pthread_mutex_unlock(&coutMutex);
  }
  
  // Call Barrier
  MyBarrier();

  //cout << "here 4" << endl; 
  // Barrier, wait for all threads to complete 1D transform
  // Once all threads complete 1D Transform
  // Transpose h
  // Do rowise 1D FFT again
  // Wait for all threads to exit.
  // return to main/function


  cout << " Writing to MyAfter1D.txt " << endl;
  string out1d("MyAfter1D.txt");
  image.SaveImageData(out1d.c_str(), h, w, ht);

  // transpose image of 1D FFTs
  TransposeIm(h,w,ht); 

  MyBarrier_Init();

  // Run through entire row and rearrange them in bir reversed order
  for (int r = 0; r < ht; r++){
    for (int c = 0; c < w; c++){
      if (rev_ind[c] <= (unsigned)c){
        continue;
      }
      else {
        Complex temp;
        temp = h[(r*w)+c];
        h[(r*w) + c] = h[(r*w)+rev_ind[c]];
        h[(r*w)+rev_ind[c]] = temp;
      }
    }
  }

  // Create Threads (SECOND)
  //cout << "here 5" << endl; 
  startCount = nThreads;
  for (int i =0; i<nThreads; i++){
    pthread_t pt;

    //pthread_mutex_lock(&coutMutex);
    //cout << "here 5_in_"<< (void*)i << endl;
    //pthread_mutex_unlock(&coutMutex);

    pthread_create(&pt, NULL, Transform2DTHread, (void*)i);

    //pthread_mutex_lock(&coutMutex);
    //cout << "here 5_out_"<< i << endl;
    //pthread_mutex_unlock(&coutMutex);
  }

  // Call Barrier
  MyBarrier();

  // Transpose Matrix
  //cout << "here 6" << endl; 
  TransposeIm(h,w,ht); 
  
  //Write result to file
  //cout << "here 7" << endl; 
  string out2d("MyAfter2D.txt");
  cout << " Writing to " << out2d << endl;
  image.SaveImageData(out2d.c_str(), h, w, ht);
  //cout << "here 8" << endl; 

/* 
  // Uncomment for main to perfom all rows 1D FFT.
  // Perform 1D of each row
  for (int r=0; r < ht; r++){
    Transform1D( h+(r*w), w );    
    //cout << "here 4_row_" << r << endl;
  }    

  cout << "here 5" << endl;
  string out1d("myAfter1D.txt");
  image.SaveImageData(out1d.c_str(), h, w, ht);
  cout << "here 6" << endl;
*/
}




// -----------------------------------------------------------
// Function: InverseTransform2D
// - Reads image into array
// - Reorders the row eements to reversedBits order,
// - Calculates Weights
// - Creates the threads(to perform 1D trasnform of rows), 
// - calls barrier to wait till all threads have executed,
// - transposes the matrix,
// - creates threads again(to perform 1D transform of rows),
// - calls barrier to wait till all threads ahve executed,
// - transposes matrix,
// - saves it to file.
// -----------------------------------------------------------
 

void InverseTransform2D(const char* inputFN) 
{ // Do the 2D transform here.
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Create the global pointer to the image array data
  // Create 16 threads
  // Create ordered array - bit reversed in place
  // Wait for all threads complete
  // Write the transformed data
  int w = ImageWidth;;
  int ht = ImageHeight;
  W = new Complex[N/2];
  unsigned rev_ind[N];
  cout <<"w " <<  w << endl;
  cout << "ht " << ht << endl;  

  // Create array rev_ind with correct order of input indices
  for (int i =0; i < w; i++){
    rev_ind[i] = ReverseBits(i);
    //cout << rev_ind[i] << endl;
  }
  //cout << "here 1" << endl;

  // Run through each row and rearrnge elements to correct order (bit reversed)
  for (int r = 0; r < ht; r++){
    for (int c = 0; c < w; c++){
      if (rev_ind[c] <= (unsigned)c){
        continue;
      }
      else {
        Complex temp;
        temp = h[(r*w)+c];
        h[(r*w) + c] = h[(r*w)+rev_ind[c]];
        h[(r*w)+rev_ind[c]] = temp;
      }
    }
  }
  //cout << "here 2" << endl;
  // Create array of weights for size w/2
  for (int n =0; n < (w/2); n++){
    W[n].real = cos(2*M_PI*n/w);
    W[n].imag = sin(2*M_PI*n/w);
  }
  //cout << "here 3" << endl;

  // Create Threads (FIRST)
  startCount = nThreads;
  for (int i =0; i<nThreads; i++){
    pthread_t pt;

    //pthread_mutex_lock(&coutMutex);
    //cout << "here 3_in_"<< (void*)i << endl;
    //pthread_mutex_unlock(&coutMutex);

    pthread_create(&pt, NULL, InverseTransform2DTHread, (void*)i);

    //pthread_mutex_lock(&coutMutex);
    //cout << "here 3_out_" << i << endl;
    //pthread_mutex_unlock(&coutMutex);
  }
  
  // Call Barrier
  MyBarrier();

  //cout << "here 4" << endl; 
  // Barrier, wait for all threads to complete 1D transform
  // Once all threads complete 1D Transform
  // Transpose h
  // Do rowise 1D FFT again
  // Wait for all threads to exit.
  // return to main/function


  cout << " Writing to MyAfter1DInverse.txt " << endl;
  string out1d("MyAfter1DInverse.txt");
  image.SaveImageData(out1d.c_str(), h, w, ht);

  // transpose image of 1D FFTs
  TransposeIm(h,w,ht); 

  MyBarrier_Init();

  // Run through entire row and rearrange them in bir reversed order
  for (int r = 0; r < ht; r++){
    for (int c = 0; c < w; c++){
      if (rev_ind[c] <= (unsigned)c){
        continue;
      }
      else {
        Complex temp;
        temp = h[(r*w)+c];
        h[(r*w) + c] = h[(r*w)+rev_ind[c]];
        h[(r*w)+rev_ind[c]] = temp;
      }
    }
  }

  // Create Threads (SECOND)
  //cout << "here 5" << endl; 
  startCount = nThreads;
  for (int i =0; i<nThreads; i++){
    pthread_t pt;

    //pthread_mutex_lock(&coutMutex);
    //cout << "here 5_in_"<< (void*)i << endl;
    //pthread_mutex_unlock(&coutMutex);

    pthread_create(&pt, NULL, InverseTransform2DTHread, (void*)i);

    //pthread_mutex_lock(&coutMutex);
    //cout << "here 5_out_"<< i << endl;
    //pthread_mutex_unlock(&coutMutex);
  }

  // Call Barrier
  MyBarrier();

  // Transpose Matrix
  //cout << "here 6" << endl; 
  TransposeIm(h,w,ht); 
  
  //Write result to file
  //cout << "here 7" << endl; 
  string out2d("MyAfterInverse.txt");
  cout << " Writing to " << out2d << endl;
  image.SaveImageDataReal(out2d.c_str(), h, w, ht);
  //cout << "here 8" << endl; 
}



int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  
  MyBarrier_Init();
  
  // Main holds the exit mutex until waiting for exitCond condition
  pthread_mutex_lock(&exitMutex);
  Transform2D(fn.c_str()); // Perform the transform.
 
  //pthread_mutex_lock(&exitMutex);
  InverseTransform2D(fn.c_str());

}  
