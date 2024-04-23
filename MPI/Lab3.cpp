//#include <mpi.h>
//#include <iostream>
//#define n 3
//#define m 4
//
//using namespace std;
//
///*1. Napisati MPI program koji kreira komunikator koga čine svi procesi sa identifikatorima deljivim sa 5. 
//Master proces (P0) svim procesima ove grupe šalje po jednu kolone matrice A. 
//Odštampati identifikatore procesa koji pripadaju novom komunikatoru, a čija suma elemenata primljene kolone matrice A nije manja od zadate vrednosti v. */
//
//void z1(int argc, char ** argv)
//{
//    int a[n][n] = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
//    int v = 10;
//    MPI_Init(&argc, &argv);
//
//    int rank = 0;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    int size = 0;
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    if (size < 10 || size >= 15)
//    {
//        cout << "Number of processes must be [10, 14] and yours was " << size << endl;
//        exit(1);
//    }
//
//    MPI_Group group, newGroup;
//    MPI_Comm_group(MPI_COMM_WORLD, &group);
//
//    int count = (size / 5) + 1;
//    int* members = new int[count];
//    for (int i = 0; i < size; i += 5)
//    {
//        members[i / 5] = i;
//    }
//
//    MPI_Group_incl(group, count, members, &newGroup);
//
//    MPI_Comm newComm;
//    MPI_Comm_create(MPI_COMM_WORLD, newGroup, &newComm);
//
//    int groupRank, groupSize = 0;
//    MPI_Group_rank(newGroup, &groupRank);
//    MPI_Group_size(newGroup, &groupSize);
//
//    MPI_Datatype column, columnResized;
//    MPI_Type_vector(n, 1, n, MPI_INT, &column);
//    MPI_Type_create_resized(column, 0, 1 * sizeof(int), &columnResized);
//    MPI_Type_commit(&columnResized);
//
//    if (groupRank >= 0)
//    {
//        int* vector = new int[n];
//        MPI_Scatter(a, 1, columnResized, vector, 3, MPI_INT, 0, newComm);
//
//        cout << "My group rank " << rank << " and size " << size << endl << "My new rank " << groupRank << " and size " << groupSize << endl;
//        cout << "My vector is: ";
//        int sum = 0;
//        for (int i = 0; i < n; i++)
//        {
//            cout << vector[i] << " ";
//            sum += vector[i];
//        }
//        cout << endl;
//        if (sum >= v)
//        {
//            cout << "My rank is " << rank << " and my group rank is " << groupRank << " and my sum is " << sum << endl;
//        }
//    }
//
//    MPI_Finalize();
//}
//
///*2. Napisati MPI program kojim se kreira dvodimenzionalna Cartesian struktura sa n vrsta i n kolona. 
//Podeliti procese u grupe koje odgovaraju gornjoj i donjoj trougaonoj matrici kreirane strukture. 
//Procese na dijagonali proizvoljno dodeliti jednoj od grupa. 
//U okviru svake grupe sumirati vrednosti identifikatora svih procesa koji pripadaju datoj grupi. 
//Master procesu komunikatora koji obuhvata sve procese dostaviti ove vrednosti i odštampati ih. 
//Ilustrovati raspored procesa i program testirati za različite dimenzije Cartesian strukture. */
//
//void z2(int argc, char** argv)
//{
//    MPI_Init(&argc, &argv);
//
//    int rank = 0;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    int size = 0;
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    if (size != n * n)
//    {
//        cout << "Number of processes must be " << n * n << " and yours was " << size << endl;
//        exit(1);
//    }
//
//    MPI_Comm cartesianComm;
//
//    int ndims = 2;
//    int dims[2] = { n, n };
//    int periods[2] = { false, true };
//    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, false, &cartesianComm);
//
//    int cartesianRank = 0;
//    MPI_Comm_rank(cartesianComm, &cartesianRank);
//
//    int coords[2] = { 0, 0 };
//    MPI_Cart_coords(cartesianComm, cartesianRank, ndims, coords);
//
//    MPI_Comm triangleComm;
//    int color = 0;
//    if (coords[0] > coords[1])
//    {
//		color = 0; // lower triangle
//	}
//	else if (coords[0] <= coords[1])
//	{
//        color = 1; // upper triangle and diagonal
//	}
//
//    MPI_Comm_split(cartesianComm, color, 0, &triangleComm);
//
//    int groupSum = 0;
//    MPI_Reduce(&cartesianRank, &groupSum, 1, MPI_INT, MPI_SUM, 0, triangleComm);
//    int groupRank = 0;
//    MPI_Comm_rank(triangleComm, &groupRank);
//
//    if (groupRank == 0)
//    {
//		cout << "My rank is " << rank << " and my group rank is " << groupRank << " and my group sum is " << groupSum << endl;
//	}
//
//    MPI_Finalize();
//}
//
///*
//Napisati MPI program kojim se kreira dvodimenzionalna Cartesian struktura sa n vrsta i m kolona. 
//U svakom procesu odštampati identifikatore i koordinate njegovog levog i desnog suseda na udaljenosti 3. 
//Ilustrovati raspored procesa i diskutovati dobijeno rešenje u zavisnosti od periodičnosti dimenzija.
//Program testirati za različite vrednosti n i m*/
//
//void z3(int argc, char** argv)
//{
//    MPI_Init(&argc, &argv);
//
//    int size = 0;
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    if (size != n * m)
//    {
//		cout << "Number of processes must be " << n * m << " and yours was " << size << endl;
//		exit(1);
//	}   
//
//    MPI_Comm cartesianComm;
//
//    int ndims = 2;
//    int dims[2] = { n, m };
//    int periods[2] = { false, true };
//    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, false, &cartesianComm);
//
//    int rank = 0;
//    MPI_Comm_rank(cartesianComm, &rank);
//
//    int coords[2] = { 0, 0 };
//    MPI_Cart_coords(cartesianComm, rank, ndims, coords);
// 
//    int rankDown  = 0, rankUp = 0;
//    MPI_Cart_shift(cartesianComm, 0, 3, &rankUp, &rankDown );
//    int rankLeft = 0, rankRight = 0;
//    MPI_Cart_shift(cartesianComm, 1, 3, &rankLeft, &rankRight);
//
//    if (rankRight >= 0)
//    {
//		int rightCoords[2] = { 0, 0 };
//		MPI_Cart_coords(cartesianComm, rankRight, 2, rightCoords);
//        cout << "My rank is " << rank << " and my right rank is " << rankRight << endl;
//        cout << "My coordinates are " << coords[0] << " " << coords[1] << " and my right coordinates are " << rightCoords[0] << " " << rightCoords[1] << endl;
//    }
//    if (rankLeft >= 0)
//    {
//        int leftCoords[2] = { 0, 0 };
//        MPI_Cart_coords(cartesianComm, rankLeft, 2, leftCoords);
//        cout << "My rank is " << rank << " and my left rank is " << rankLeft << endl;
//        cout << "My coordinates are " << coords[0] << " " << coords[1] << " and my left coordinates are " << leftCoords[0] << " " << leftCoords[1] << endl;
//    }
//
//    if (rankDown  >= 0)
//    {
//        int upperCoords[2] = { 0, 0 };
//        MPI_Cart_coords(cartesianComm, rankDown, 2, upperCoords);
//		cout << "My rank is " << rank << " and my upper rank is " << rankDown  << endl;
//        cout << "My coordinates are " << coords[0] << " " << coords[1] << " and my upper coordinates are " << upperCoords[0] << " " << upperCoords[1] << endl;
//	}
//    if (rankUp >= 0)
//    {
//		int lowerCoords[2] = { 0, 0 };
//		MPI_Cart_coords(cartesianComm, rankUp, 2, lowerCoords);
//		cout << "My rank is " << rank << " and my lower rank is " << rankUp << endl;
//		cout << "My coordinates are " << coords[0] << " " << coords[1] << " and my lower coordinates are " << lowerCoords[0] << " " << lowerCoords[1] << endl;  
//	}
//
//	MPI_Finalize();
//}
//
//void z4(int argc, char** argv)
//{
//    int rank, size;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    if (size != n * m)
//    {
//        cout << "Number of processes must be " << n * m << " and yours was " << size << endl;
//        exit(1);
//    }
//
//    int k = 0;
//    MPI_Comm cartesianComm;
//
//    int ndims = 2;
//    int dims[2] = { n, m };
//    int periods[2] = { false, true };
//    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, false, &cartesianComm);
//
//    int rankLeft, rankRight;
//    MPI_Cart_shift(cartesianComm, 1, 2, &rankLeft, &rankRight);
//
//    if (rankLeft >= 0)
//    {
//		cout << "My rank is " << rank << " and my left rank is " << rankLeft << endl;
//	}
//    if (rankRight >= 0)
//    {
//		cout << "My rank is " << rank << " and my right rank is " << rankRight << endl;
//	}
//
//    k += rank;
//    int newK = 0;
//    MPI_Status status;
//
//    MPI_Sendrecv(&k, 1, MPI_INT, rankRight, 0, &newK, 1, MPI_INT, rankLeft, 0, cartesianComm, &status);
//
//    if (rank % 2 == 0)
//    {
//        cout << "My rank is " << rank << " and my old k is " << k << " and my new k is " << newK << endl;
//    }
//
//    MPI_Finalize();
//}
//
//int main(int argc, char** argv)
//{
//    z4(argc, argv);
//}
//
