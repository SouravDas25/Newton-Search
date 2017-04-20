#ifndef CMATRIX_H
#define CMATRIX_H

#include <stdio.h>

template <class T>
class CMatrix
{
public:
    CMatrix( int rows, int cols)
    {
        row = rows;
        col = cols;

        data = new T*[rows]; // replaced "int" for "T"

        for (int i = 0; i < row; i++) {
            data[i] = new T [cols]; // replaced "int" for "T"
        }

        for(int i = 0; i < row; i++) {
            for(int j = 0; j < cols; j++) {
                data[i][j] = (T) 0; // replaced "int" for "T"
            }
        }
    }

    void print();
    void setRow(int col,T * array)
	{
		for(int i = 0; i < this->row; i++)
		{
			data[i][col] = array[i];
		}
	}
    void setCol(int row,T*array)
	{
		for(int j = 0; j < this->col; j++)
		{
			data[row][j] = array[j];
		}
	}
    T& operator()(int row, int col);

    T **data;
    int row,col;
};

template <class T>
void CMatrix<T>::print ()
{
    int i,j;

    for (i=0;i < row;i++) // Here you used to have row hard coded to 4
    {
        for(j=0;j < col;j++) // Here you used to have col hard coded to 4
        {
            printf("%10g    ",(float) data[i][j]);
        }
        printf("\n");
    }
}

// Recently added
template<class T> T& CMatrix<T>::operator()(int row, int col)
{
    return data[row][col];
}

#endif // CMATRIX_H
