//#include <mpi.h>
//#include <iostream>
//#include <math.h>
//using namespace std;
//#define n 4
//
//int main(int argc, char** argv) 
//{
//	int rank, size;
//	int matrix[n][n];
//	int* myElements;
//	int numOfMyElements;
//	int localMax = 0, globalMax = 0;
//	MPI_Datatype evenRows;
//	MPI_Datatype extendedEvenRows{};
//	
//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//	if (floor((double)n / (size - 1)) != ((double)n / (size - 1)))
//	{
//		cout << "Dimension isn't divisible by number of processes" << endl;
//		exit(-1);
//	}
//
//	if (rank == 0)
//	{
//		cout << "Input matrix elements" << endl;
//
//		for (int i = 0; i < n; i++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				cin >> matrix[i][j];
//			}
//		}
//		
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	int count = n / 2; // broj parnih vrsta
//	int blockLength = n / size;
//	int stride = 2 * n;
//
//	MPI_Type_vector(count, blockLength, stride, MPI_INT, &evenRows);
//	MPI_Type_commit(&evenRows);
//	MPI_Type_create_resized(evenRows, 0, (n / size) * sizeof(int), &extendedEvenRows);
//	MPI_Type_commit(&extendedEvenRows);
//
//	numOfMyElements = blockLength * count;
//	myElements = new int[numOfMyElements];
//
//	MPI_Scatter(&matrix[0][0], 1, extendedEvenRows, myElements, numOfMyElements, MPI_INT, 0, MPI_COMM_WORLD);
//
//	localMax = myElements[0];
//	for (int i = 1; i < numOfMyElements; i++)
//	{
//		if (myElements[i] > localMax)
//			localMax = myElements[i];
//	}
//
//	MPI_Reduce(&localMax, &globalMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
//
//	if (rank == 0)
//	{
//		cout << "Global maximum is : " << globalMax << endl;
//	}
//
//	delete[] myElements;
//	MPI_Finalize();
//}