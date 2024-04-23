#include <mpi.h>
#include <iostream>
using namespace std;
#define m 4
#define l 2
#define n 4

struct 
{
    double value;
    int rank;
}in, out;

//void z1(int argc, char** argv)
//{
//    int n;
//    double* arr;
//    int rank, size;
//    double min;
//
//    MPI_Status stat;
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//
//    if (rank == 0)
//    {
//        cout << "Input array length" << endl;
//        cin >> n;
//        cout << "Input array elements" << endl;
//        arr = new double[n];
//
//        for (int i = 0; i < n; i++)
//        {
//            cin >> arr[i];
//        }
//
//        int numOfElToSend = n / size;
//        for (int i = 0; i < size - 1; i++)
//        {
//            MPI_Send(&numOfElToSend, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
//            MPI_Send(arr + (i + 1) * numOfElToSend, numOfElToSend, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
//        }
//
//        double* minArr = new double[size - 1];
//        for (int i = 0; i < size - 1; i++)
//        {
//            MPI_Recv(minArr + i, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &stat);
//        }
//
//        min = arr[0];
//        for (int i = 0; i < numOfElToSend; i++)
//        {
//            if (arr[i] < min)
//                min = arr[i];
//        }
//
//        for (int i = 1; i < size - 1; i++)
//        {
//            if (minArr[i] < min)
//                min = minArr[i];
//        }
//
//        cout << "Minimum is " << min << endl;
//
//    }
//    else
//    {
//        MPI_Recv(&n, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
//        arr = new double[n];
//        MPI_Recv(arr, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
//
//        min = arr[0];
//        for (int i = 1; i < n; i++)
//        {
//            if (arr[i] < min)
//                min = arr[i];
//        }
//
//        MPI_Send(&min, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//    }
//
//
//    MPI_Finalize();
//}

//void z2(int argc, char** argv)
//{
//    int a[k][m];
//    int b[m][l];
//
//	int c[k][l];
//
//	int rank, size;
//	MPI_Status stat;
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);  
//
//    if (m % size != 0)
//    {
//		cout << "Error: m % size != 0" << endl;
//		return;
//	}
//
//    if (rank == 0)
//    {
//        for (int i = 0; i < k; i++)
//        {
//            for (int j = 0; j < m; j++)
//            {
//				a[i][j] = i * m + j;
//			}
//		}
//        for (int i = 0; i < m; i++)
//        {
//            for (int j = 0; j < l; j++)
//            {
//                b[i][j] = i * l + j;
//            }
//        }
//    }
//
//    MPI_Datatype column, resizedColumn;
//    MPI_Type_vector((m / size) * k, 1, size, MPI_INT, &column);
//    MPI_Type_create_resized(column, 0, sizeof(int), &resizedColumn);
//    MPI_Type_commit(&resizedColumn);
//
//    int* localA = new int[k * m / size];
//    MPI_Scatter(a, 1, resizedColumn, localA, (m / size) * k, MPI_INT, 0, MPI_COMM_WORLD);
//
//  //  for (int i = 0; i < k; i++)
//  //  {
//  //      for (int j = 0; j < m / p; j++)
//  //      {
//		//	cout << "Process " << rank << " " << localA[i][j] << endl;
//		//}
//  //  }
//
//    MPI_Datatype row, resizedRow;
//    MPI_Type_vector(m / size, l, l * size, MPI_INT, &row);
//    MPI_Type_create_resized(row, 0, l * sizeof(int), &resizedRow);
//    MPI_Type_commit(&resizedRow);
//
//    int* localB = new int[(m / size) * l];
//    MPI_Scatter(b, 1, resizedRow, localB, (m / size) * l, MPI_INT, 0, MPI_COMM_WORLD);
//
// //   for (int i = 0; i < m / p; i++)
// //   {
// //       for (int j = 0; j < l; j++)
// //       {
//	//		cout << "Process " << rank << " " << localB[i][j] << endl;
//	//	}
//	//}
//
//    int localC[k][l];
//    for (int i = 0; i < k; i++)
//    {
//        for (int j = 0; j < l; j++)
//        {
//            localC[i][j] = 0;
//            for (int z = 0; z < m / size; z++)
//            {
//				localC[i][j] += localA[i * m / size + z] * localB[z * l + j];
//			}   
//        }
//    }
//
//    //for (int i = 0; i < k; i++)
//    //{
//    //   for (int j = 0; j < l; j++)
//    //   {
//   	//	    cout << "Process " << rank << " " << localC[i][j] << endl;
//   	//   }
//    //}
//
//    int localMin = localA[0];
//    for (int i = 0; i < k; i++)
//    {
//        for (int j = 0; j < m / size; j++)
//        {
//			if (localA[i * m / size +j] < localMin)
//				localMin = localA[i * m / size + j];
//		}
//	}
//    in.rank = rank;
//    in.value = localMin;
//
//    MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
//
//    int globalMinRank;
//
//    if (rank == 0)
//    {
//		globalMinRank = out.rank;
//        cout << "Global min rank is " << globalMinRank << endl;
//	}
//    MPI_Bcast(&globalMinRank, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//    MPI_Reduce(localC, c, k * l, MPI_INT, MPI_SUM, globalMinRank, MPI_COMM_WORLD);
//
//    if (rank == globalMinRank)
//    {
//		cout << "Result matrix" << endl;
//        for (int i = 0; i < k; i++)
//        {
//            for (int j = 0; j < l; j++)
//            {
//                cout << c[i][j] << " "; 
//            }
//            cout << endl;
//        }
//    }
//
//    MPI_Finalize();
//}

//void z3(int argc, char** argv)
//{
//    int a[n][n];
//    int b[n];
//    int c[n];
//    MPI_Init(&argc, &argv);
//
//    int rank, size;
//    MPI_Status stat;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    if (sqrt(size) != floor(sqrt(size)))
//    {
//        cout << "Number of processes must be square" << endl;
//        exit(1);
//    }
//
//    if (rank == 0)
//    {
//        for (int i = 0; i < n; i++)
//        {
//			b[i] = i;
//            for (int j = 0; j < n; j++)
//            {
//				a[i][j] = i * n + j;
//			}
//		}   
//    }
//
//    int q = sqrt(size);
//    MPI_Datatype vector, resizedVector;
//    MPI_Type_vector(n / q, 1, q, MPI_INT, &vector);
//    MPI_Type_create_resized(vector, 0, sizeof(int), &resizedVector);
//    MPI_Type_commit(&resizedVector);
//
//    MPI_Comm rowComm, columnComm;
//    int rowRank, columnRank;
//    MPI_Comm_split(MPI_COMM_WORLD, rank / q, 0, &rowComm);
//    MPI_Comm_split(MPI_COMM_WORLD, rank % q, 0, &columnComm);
//
//    MPI_Comm_rank(rowComm, &rowRank);
//    MPI_Comm_rank(columnComm, &columnRank);
//
//    int* localB = new int[n / q];
//    if (columnRank == 0)
//    {
//        MPI_Scatter(b, 1, resizedVector, localB, n / q, MPI_INT, 0, rowComm);
//    }
//    MPI_Bcast(localB, n / q, MPI_INT, 0, columnComm);
//
//    //for (int i = 0; i < n / q; i++)
//    //{
//    //  	cout << "Process " << rank << " " << localB[i] << endl;
//    //}
//
//    MPI_Datatype row;
//    MPI_Type_vector(n / q * n / q, 1, q, MPI_INT, &row);
//    MPI_Type_commit(&row);
//
//    int* localA = new int[n / q * n / q];
//
//    if (rank == 0)
//    {
//        for (int i = 0; i < n / q; i++)
//            for (int j = 0; j < n; j += q)
//            {
//				localA[i * n / q + j / q] = a[i][j];
//			}
//        for (int i = 1; i < size; i++)
//        {
//            MPI_Send(&a[i / q * n / q][i % q], 1, row, i, 0, MPI_COMM_WORLD);
//        }
//    }
//    else
//    {
//        MPI_Recv(localA, n / q * n / q, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
//    }
//
//    int* localC = new int[n / q];
//    for (int i = 0; i < n / q; i++)
//    {
//		localC[i] = 0;
//		for (int j = 0; j < n / q; j++)
//		{ 
//            localC[i] += localA[i * n / q + j] * localB[j];
//        }
//    }
//
//    // sum by row
//    int* localCSummed = new int[n / q];
//    MPI_Reduce(localC, localCSummed, n / q, MPI_INT, MPI_SUM, 0, rowComm);
//
//    // gather by column
//    if (rowRank == 0)
//    {
//        MPI_Gather(localCSummed, n / q, MPI_INT, c, n / q, MPI_INT, 0, columnComm);
//
//        if (rank == 0)
//        {
//			cout << "result vector" << endl;
//            for (int i = 0; i < n; i++)
//            {
//				cout << c[i] << " ";
//			}
//			cout << endl;
//		}   
//    }
//
//    MPI_Finalize();
//}

void z4(int argc, char** argv)
{
    int rank, size, k = 2;
    int a[n][n], b[n][n], c[n][n];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    if (size != (n / k) * (n / k))
    {
        cout << "Error: size != (n / k) * (n / k) = " << (n / k) * (n / k) << endl;
        exit(1);
    }
    else if (k > n)
    {
        cout << "Error: k > n" << endl;
		exit(1);
    }

    if (rank == 0)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                a[i][j] = i * n + j;
                b[i][j] = i * n + j;
            }
        }
    }

    MPI_Datatype rowBlock, columnBlock;
    MPI_Type_contiguous(k * n, MPI_INT, &rowBlock);
    MPI_Type_commit(&rowBlock);
    MPI_Type_vector(n, k, n, MPI_INT, &columnBlock);
    MPI_Type_commit(&columnBlock);

    int* localA = new int[k * n];
    int* localB = new int[n * k];
    if (rank == 0)
    {
        // self inicialization
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < n; j++)
            {
                localA[i * n + j] = a[i][j];
            }
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
            {
                localB[i * k + j] = b[i][j];
            }
        }

        // distribution of data
        for (int i = 1; i < size; i++)
        {
            MPI_Send(&a[(i / (n / k)) * k][0], 1, rowBlock, i, 0, MPI_COMM_WORLD);
            MPI_Send(&b[0][(i % (n / k))* k], 1, columnBlock, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
		MPI_Recv(localA, k * n, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(localB, n * k, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

    int* localC = new int[k * k];
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            localC[i * k + j] = 0;

            for (int z = 0; z < n; z++)
                localC[i * k + j] += localA[i * n + z] * localB[z * k + j];
        }
    }

    MPI_Comm rowComm, columnComm;
    int rowRank, columnRank;
    MPI_Comm_split(MPI_COMM_WORLD, rank / (n / k), 0, &rowComm);
    MPI_Comm_split(MPI_COMM_WORLD, rank % (n / k), 0, &columnComm);
    MPI_Comm_rank(rowComm, &rowRank);
    MPI_Comm_rank(columnComm, &columnRank);

    int* localCGathered = new int[k * n];
    for (int i = 0; i < k; i++)
        MPI_Gather(&localC[i * k], k, MPI_INT, &localCGathered[i * n], k, MPI_INT, 0, rowComm);

    if (rowRank == 0)
        MPI_Gather(localCGathered, k * n, MPI_INT, c, k * n, MPI_INT, 0, columnComm);


    if (rank == 0)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
				cout << c[i][j] << " ";
			}
			cout << endl;
		}   
    }

    MPI_Finalize();
    delete[] localA;
    delete[] localB;
    delete[] localC;
}

int main(int argc, char** argv)
{
	z4(argc, argv);
	return 0;
}

