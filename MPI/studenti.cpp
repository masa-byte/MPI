//#include <mpi.h>
//#include <iostream>
//
//using namespace std;
//
//#define n 10
//
//typedef struct Student {
//	int index;
//	char name[n];
//	char surname[n];
//	float avgGrade;
//	int year;
//};
//
//
//void students(int argc, char** argv)
//{
//	Student s;
//	int rank, size;
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//	if (rank == 0)
//	{
//		cout << "Input student's index, name, surname, average grade and year of study" << endl;
//		cin >> s.index >> s.name >> s.surname >> s.avgGrade >> s.year;
//	}
//
//	MPI_Datatype studentType;
//	int count = 5;
//	MPI_Datatype types[5] = { MPI_INT, MPI_CHAR, MPI_CHAR, MPI_FLOAT, MPI_INT };
//	int arrOfBlockLengths[5] = { 1, n, n, 1, 1};
//	MPI_Aint offsets[5];
//	offsets[0] = offsetof(Student, index);
//	offsets[1] = offsetof(Student, name);
//	offsets[2] = offsetof(Student, surname);
//	offsets[3] = offsetof(Student, avgGrade);
//	offsets[4] = offsetof(Student, year);
//
//	MPI_Type_create_struct(count, arrOfBlockLengths, offsets, types, &studentType);
//	MPI_Type_commit(&studentType);
//
//	MPI_Bcast(&s, 1, studentType, 0, MPI_COMM_WORLD);
//
//	cout << "Process " << rank << " received student with index " << s.index << ", name " << s.name << ", surname " << s.surname << ", average grade " << s.avgGrade << " and year of study " << s.year << endl;
//
//	MPI_Finalize();
//}
//
//int main(int argc, char** argv)
//{
//	students(argc, argv);
//	return 0;
//}