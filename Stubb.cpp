/*************************************************************
* Author:		Justin.Urbany
* Filename:		Stub.cpp
* Date Created:	1/23/17
* Modifications:	1/24/17-added insertion and selection
1/25/17 fixed insertion
1/29/17 added heap and quick
1/30/17 added merge and shell
1/31/17 fixed shell
**************************************************************/
#define _CRTDBG_MAP_ALLOC
#include <iostream>
#include <algorithm>
#include<fstream>
#include <vector>
#include<time.h>
#include "Array.h"
#include <crtdbg.h>
using std::vector;
using std::fstream;
using std::ios;
using std::swap;
using std::cout;
using std::endl;
/*************************************************************
*
* Lab/Assignment: Lab 2- Sorts
*
* Overview:
*	This program is designed to test sort methods on c arrays,
*	vectors, and my array class with integers to sort the data
*	structures from smallest to greatest integer value
*
* Input:
*   The program doesn't take any input from the user takes
*	specified command line arguments that will be the given
*   value of integers that are wanted to sort
*
* Output:
*   In this program it will output the sort method what data
*	structure it is prefroming the sort on and how many inters
*	that data structure is holding
************************************************************/
template <typename T>
void BruteForceBubble(T & unsorted, int n);
template <typename T>
void FlaggedBubble(T & unsorted, int n);
template <typename T>
void Selection(T & unsorted, int n);
template <typename T>
void Insertion(T & unsorted, int n);
template<typename T>
void Shell(T & unsorted, int n);
template<typename T>
void MoveDown(T & unsorted, int first_node, int last_index);
template<typename T>
void Heap(T & unsorted, int n);
template<typename T>
void Merge(T & unsorted, int n);
template<typename T>
void MergeSort(T & unsorted, int * temp, int leftindex, int rightindex);
template<typename T>
void MergeSorts(T & unsorted, int * temp, int left, int right, int right_end);
template<typename T>
void Quick(T & unsorted, int n);
template<typename T>
void QuickSort(T & unsorted, int first, int last);
template<typename T>
void Randomize(T & sorted, int n);
template<typename T>
void DeepCopy(T & sorted, T & unsorted, int n);


int main(int argc, char *argv[])
{
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); //memory leak checker
	fstream sortstime("Stimes.txt", ios::out | ios::app); //opens file for output in append mode
	int n = atoi(argv[1]); //takes the 2nd command line arguement and converts it to an int
	int * unsortedarray = new int[n]; //make a pointer to ints and set the length to passed in parameter
	vector<int> unsorted_vector(n); //make a vector with passed in length
	Array<int> my_unsortedarray(n); //make an array with passed in length

	int * sortedarray = new int[n];
	vector<int> sorted_vector(n);
	Array<int> my_sortedarray(n);

	clock_t time = 0; //make a clock to hold time when the sort started
	double total_time = 0.00000; //will hold the total time it took to complete sort


	Randomize(unsortedarray, n); //randomize the values for the c style array
	Randomize(unsorted_vector, n); //randomize the values for the vector
	Randomize(my_unsortedarray, n); //randomize the values for my array class array

		//generic C array

		DeepCopy(sortedarray, unsortedarray, n);
		time=clock();
		BruteForceBubble(sortedarray, n);
		total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
		cout << "BruteForceBubble sort on C array for " << n << " elements took: " << total_time << endl;
		//sortstime << "BruteForceBubble sort on C array for " << n << " elements took: " << total_time << endl;

		DeepCopy(sortedarray, unsortedarray, n);
		time = clock();
		FlaggedBubble(sortedarray, n);
		total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
		cout << "FlaggedBubble sort on C array for "<<n<<" elements took: " << total_time << endl;
		//sortstime<< "FlaggedBubble sort on C array for " << n << " elements took: " << total_time << endl;

		DeepCopy(sortedarray, unsortedarray, n);
		time = clock();
		Selection(sortedarray, n);
		total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
		cout << "Selection sort on C array for " << n << " elements took: " << total_time << endl;
		//sortstime<< "Selection sort on C array for " << n << " elements took: " << total_time << endl;

		DeepCopy(sortedarray, unsortedarray, n);
		time = clock();
		Insertion(sortedarray, n);
		total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
		cout << "Insertion sort on C array for " << n << " elements took: " << total_time << endl;
		//sortstime<< "Insertion sort on C array for " << n << " elements took: " << total_time << endl;

		DeepCopy(sortedarray, unsortedarray, n);
		time = clock();
		Shell(sortedarray, n);
		total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
		cout << "Shell sort on C array for " << n << " elements took: " << total_time << endl;
		//sortstime << "Shell sort on C array for " << n << " elements took: " << total_time << endl;

		DeepCopy(sortedarray, unsortedarray, n);
		time = clock();
		Heap(sortedarray, n);
		total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
		cout << "Heap sort on C array for " << n << " elements took: " << total_time << endl;
		//sortstime << "Heap sort on C array for " << n << " elements took: " << total_time << endl;

	DeepCopy(sortedarray, unsortedarray, n);
	time = clock();
	Merge(sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Merge sort on C array for " << n << " elements took: " << total_time << endl;
	//sortstime << "Merge sort on C array for " << n << " elements took: " << total_time << endl;

	time = clock();
	DeepCopy(sortedarray, unsortedarray, n);
	Quick(sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Quick sort on C array for " << n << " elements took: " << total_time << endl;
	//sortstime << "Quick sort on C array for " << n << " elements took: " << total_time << endl;


	//Vectors

	DeepCopy(sorted_vector, unsorted_vector, n);
	time = clock();
	BruteForceBubble(sorted_vector, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "BruteForceBubble sort on vector for " << n << " elements took: " << total_time << endl;
	//sortstime << "BruteForceBubble sort on vector for " << n << " elements took: " << total_time << endl;

	DeepCopy(sorted_vector, unsorted_vector, n);
	time = clock();
	FlaggedBubble(sorted_vector, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "FlaggedBubble sort on vector for " << n << " elements took: " << total_time << endl;
	//sortstime << "FlaggedBubble sort on vector for " << n << " elements took: " << total_time << endl;

	DeepCopy(sorted_vector, unsorted_vector, n);
	time = clock();
	Selection(sorted_vector, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Selection sort on vector for " << n << " elements took: " << total_time << endl;
	//sortstime << "Selection sort on vector for " << n << " elements took: " << total_time << endl;

	DeepCopy(sorted_vector, unsorted_vector, n);
	time = clock();
	Insertion(sorted_vector, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Insertion sort on vector for " << n << " elements took: " << total_time << endl;
	//sortstime << "Insertion sort on vector for " << n << " elements took: " << total_time << endl;

	DeepCopy(sorted_vector, unsorted_vector, n);
	time = clock();
	Shell(sorted_vector, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Shell sort on vector for " << n << " elements took: " << total_time << endl;
	//sortstime << "Shell sort on vector for " << n << " elements took: " << total_time << endl;

	DeepCopy(sorted_vector, unsorted_vector, n);
	time = clock();
	Heap(sorted_vector, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Heap sort on vector for " << n << " elements took: " << total_time << endl;
	//sortstime << "Heap sort on vector for " << n << " elements took: " << total_time << endl;

	DeepCopy(sorted_vector, unsorted_vector, n);
	time = clock();
	//Merge(sorted_vector, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Merge sort on vector for " << n << " elements took: " << total_time << endl;
	//sortstime<< "Merge sort on vector for " << n << " elements took: " << total_time << endl;

	DeepCopy(sorted_vector, unsorted_vector, n);
	time = clock();
	Quick(sorted_vector, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Quick sort on vector for " << n << " elements took: " << total_time << endl;
	//sortstime<< "Quick sort on vector for " << n << " elements took: " << total_time << endl;

	//Myarray Class


	DeepCopy(my_sortedarray, my_unsortedarray, n);
	time = clock();
	BruteForceBubble(my_sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "BruteForceBubble sort for my Array class for " << n << " elements took: " << total_time << endl;
	//sortstime<< "BruteForceBubble sort for my Array class for " << n << " elements took: " << total_time << endl;

	DeepCopy(my_sortedarray, my_unsortedarray, n);
	time = clock();
	FlaggedBubble(my_sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "FlaggedBubble sort for my Array class for " << n << " elements took: " << total_time << endl;
	//sortstime << "FlaggedBubble sort for my Array class for " << n << " elements took: " << total_time << endl;

	DeepCopy(my_sortedarray, my_unsortedarray, n);
	time = clock();
	Selection(my_sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Selection sort for my Array class for " << n << " elements took: " << total_time << endl;
	//sortstime << "Selection sort for my Array class for " << n << " elements took: " << total_time << endl;

	DeepCopy(my_sortedarray, my_unsortedarray, n);
	time = clock();
	Insertion(my_sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Insertion sort for my Array class for " << n << " elements took: " << total_time << endl;
	//sortstime << "Insertion sort for my Array class for " << n << " elements took: " << total_time << endl;

	DeepCopy(my_sortedarray, my_unsortedarray, n);
	time = clock();
	Shell(my_sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Shell sort for my Array class for " << n << " elements took: " << total_time << endl;
	//sortstime << "Shell sort for my Array class for " << n << " elements took: " << total_time << endl;

	DeepCopy(my_sortedarray, my_unsortedarray, n);
	time = clock();
	Heap(my_sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Heap sort for my Array class for " << n << " elements took: " << total_time << endl;
	//sortstime << "Heap sort for my Array class for " << n << " elements took: " << total_time << endl;

	DeepCopy(my_sortedarray, my_unsortedarray, n);
	time = clock();
	Merge(my_sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Merge sort for my Array class for " << n << " elements took: " << total_time << endl;
	//sortstime << "Merge sort for my Array class for " << n << " elements took: " << total_time << endl;

	DeepCopy(my_sortedarray, my_unsortedarray, n);
	time = clock();
	Quick(my_sortedarray, n);
	total_time = ((clock() - time) / (double)CLOCKS_PER_SEC);
	cout << "Quick sort for my Array class for " << n << " elements took: " << total_time << endl;
	//sortstime << "Quick sort for my Array class for " << n << " elements took: " << total_time << endl;

	sortstime.close(); //close the files

	delete[] sortedarray; //delete the sorted array
	delete[] unsortedarray; //delete the unsorted array

	return 0;
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the integers in
*			least to greatest value by moving the greatest value
*			to the most significant position with every pass
*
* Precondition: Takes in a data structure and the length of the data
structure
* Postcondition:
*				The passed in data structure is sorted from descending
*				to ascending order
************************************************************************/
template <typename T>
void BruteForceBubble(T & unsorted, int n)
{
	for (int i = 0; i < n; ++i) //loop for every value in the array
	{
		for (int j = 0; j < (n - 1); ++j) //loop for every value in the array
		{
			if (unsorted[j] > unsorted[j + 1]) //if current greater then next
			{
				swap(unsorted[j], unsorted[j + 1]); //swap there position
			}
		}
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the integers in
*			least to greatest value by moving the greatest value
*			to the most significant position with every pass also checks
*			to see if the data structure is already in order if it is
*			quit
*
* Precondition: Takes in a data structure and the length of the data
structure
* Postcondition:
*				The passed in data structure is sorted from descending
*				to ascending order
************************************************************************/
template<typename T>
void FlaggedBubble(T & unsorted, int n)
{
	bool sorted = false; //create a check to see if sorted
	int pass = 0; //can decrement number of loops since greatest value will
				  //be in most significant position after each of the second loops
	for (int i = 0; i < n && sorted == false; ++i)
	{
		sorted = true;
		for (int j = 0; j < (n - pass - 1); ++j)
		{
			if (unsorted[j] > unsorted[j + 1])
			{
				sorted = false;
				swap(unsorted[j], unsorted[j + 1]);
			}
		}
		pass++; //increment the pass
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the data structure
*			from least to greatest value
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure sorted from least to most significant
*				value
************************************************************************/
template<typename T>
void Selection(T & unsorted, int n)
{
	int max = 0; //make variable max
	for (int i = 0; i < (n - 1); i++)
	{
		max = i;
		for (int j = (i + 1); j < n; ++j)
		{
			if (unsorted[j] < unsorted[max])
			{
				max = j;
			}
		}
		swap(unsorted[max], unsorted[i]);
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the data structure
*			from least to greatest value
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure sorted from least to most significant
*				value
************************************************************************/
template<typename T>
void Insertion(T & unsorted, int n)
{
	int temp;
	int j = 0;
	for (int i = 1; i < n; ++i)
	{
		temp = unsorted[i];
		for (j = i; j > 0 && temp < unsorted[j - 1]; --j)
		{
			unsorted[j] = unsorted[j - 1];
		}
		unsorted[j] = temp;
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the data structure
*			from least to greatest value
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure sorted from least to most significant
*				value
************************************************************************/
template<typename T>
void Shell(T & unsorted, int n)
{
	int i = 0;
	int * increments = new int[n];
	int h = 0;
	for (h = 1; h < n; ++i)
	{
		increments[i] = h;
		h = 3 * h + 1;
	}
	i;
	for (i; i >= 0; --i)
	{
		h = increments[i];
		for (int hCnt = h; hCnt < (2 * h); ++hCnt)
		{
			for (int j = hCnt; j < n; j++)
			{
				int temp = unsorted[j];
				int k = j;
				while (k - h >= 0 && temp < unsorted[k - h])
				{
					unsorted[k] = unsorted[k - h];
					k = k - h;
				}
				unsorted[k] = temp;
			}
		}
	}
	delete[]increments;
}
/**********************************************************************
* Purpose: The purpose of this function is to create a heap of data or give
*			data a form of organization to where it is somewhat sorted
* Precondition: Take in the data structure take in the position of the first
*				element that needs sorting and the last element that needs sorting
* Postcondition: Have the data structure sorted in chunks that are from least
*				to greatest value
************************************************************************/
template<typename T>
void MoveDown(T & unsorted, int first_node, int last_index)
{
	int largest = first_node * 2 + 1;
	while (largest <= last_index)
	{
		if (largest < last_index && unsorted[largest] < unsorted[largest + 1])
		{
			largest++;
		}
		if (unsorted[first_node] < unsorted[largest])
		{
			swap(unsorted[first_node], unsorted[largest]);
			first_node = largest;
			largest = first_node * 2 + 1;
		}
		else
		{
			largest = last_index + 1;
		}
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the data structure
*			from least to greatest value
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure sorted from least to most significant
*				value
************************************************************************/
template<typename T>
void Heap(T & unsorted, int n)
{
	for (int i = n / 2; i > 0; --i)
	{
		MoveDown(unsorted, i, n - 1);
	}
	for (int i = n - 1; i >= 0; --i)
	{
		swap(unsorted[0], unsorted[i]);
		MoveDown(unsorted, 0, i - 1);
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the data structure
*			from least to greatest value
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure sorted from least to most significant
*				value
************************************************************************/
template<typename T>
void Merge(T & unsorted, int n)
{
	int * temp = new int[n];
	MergeSort(unsorted, temp, 0, n - 1);

	delete[]temp;
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the data structure
*			from least to greatest value
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure sorted from least to most significant
*				value
************************************************************************/
template<typename T>
void MergeSort(T & unsorted, int * temp, int leftindex, int rightindex)
{
	int mid;
	if (leftindex < rightindex)
	{
		mid = (leftindex + rightindex) / 2;
		MergeSort(unsorted, temp, leftindex, mid);
		MergeSort(unsorted, temp, (mid + 1), rightindex);
		MergeSorts(unsorted, temp, leftindex, (mid + 1), rightindex);
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to sort the data structure
*			from least to greatest value
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure sorted from least to most significant
*				value
************************************************************************/
template<typename T>
void MergeSorts(T & unsorted, int * temp, int left, int right, int right_end)
{
	int left_end = right - 1;
	int temp_pos = left;
	int right_beg = 0;
	int num_elements = left;//right_end - left +1;

	while (temp_pos <= left_end && right <= right_end)
	{
		if (unsorted[left] <= unsorted[right])
		{
			temp[temp_pos] = unsorted[left];
			temp_pos++;
			left++;
		}
		else
		{
			temp[temp_pos] = unsorted[right];
			temp_pos++;
			right++;
		}
	}
	while (temp_pos <= left_end)
	{
		temp[temp_pos] = unsorted[left];
		left++;
		temp_pos++;
	}
	while (temp_pos <= right_end)
	{
		temp[temp_pos] = unsorted[right];
		temp_pos++;
		right++;
	}
	for (int i = left; i <= right_end; i++)
	{
		unsorted[i] = temp[i];
	}

}
/**********************************************************************
* Purpose: The purpose of this function is to sort the data structure
*			from least to greatest value
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure sorted from least to most significant
*				value
************************************************************************/
template<typename T>
void Quick(T & unsorted, int n)
{
	if (n >= 2)
	{
		int max = 0;
		for (int i = 1; i < n; ++i)
		{
			if (unsorted[max] < unsorted[i])
			{
				max = i;
			}
		}
		swap(unsorted[n - 1], unsorted[max]);
		QuickSort(unsorted, 0, n - 2);
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to sort a chunck of the data
*			structure from least to most siginifcant value
* Precondition: Take in the data structure the compare value and the position
*				of the last value that needs to be sorted
* Postcondition: Have the chunk of data sorted
************************************************************************/
template<typename T>
void QuickSort(T & unsorted, int first, int last)
{
	int small = first + 1;
	int large = last;
	int pivot = unsorted[first];
	while (small <= large)
	{
		while (unsorted[small] < pivot)
		{
			small = small + 1;
		}
		while (unsorted[large] > pivot)
		{
			large = large - 1;
		}
		if (small < large)
		{
			swap(unsorted[small], unsorted[large]);
			small++;
			large--;
		}
		else
		{
			small = small + 1;
		}
	}
	swap(unsorted[large], unsorted[first]);
	if (first < (large - 1))
	{
		QuickSort(unsorted, first, large - 1);

	}
	if (last > large + 1)
	{
		QuickSort(unsorted, large + 1, last);
	}
}
/**********************************************************************
* Purpose: The purpose of this function is to give a data structure random
*			integer values for the length of the data structure
* Precondition: Take in the data structure and length of the data structure
* Postcondition: Have the data structure set with random integer values
************************************************************************/
template<typename T>
void Randomize(T & sorted, int n)
{
	srand(time(NULL));
	for (int i = 0; i < n; ++i)
	{
		sorted[i] = rand();
	}
}
/**********************************************************************
* Purpose: The purpose of this function is transfer all the data contained
*			in a data structure to another data structure
* Precondition: Take in two data structures and their lengths
* Postcondition: Have all the values put into the data structure in the
*				same order
************************************************************************/
template<typename T>
void DeepCopy(T & sorted, T & unsorted, int n)
{
	for (int i = 0; i < n; i++)
	{
		sorted[i] = unsorted[i];
	}
}