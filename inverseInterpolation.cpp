

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "Vectors\cppmat.h"

using namespace std;

int factorial(int n)
{
	if(n<=1)return 1;
	return n*factorial(n-1);
}

template <typename T>
CMatrix<T>& getDiffTable(T * fxs ,int n)
{
	CMatrix <T> * pm = new CMatrix<T>(n,n);
	CMatrix<T>& m = *pm;
	m.setRow(0,fxs);
	
	int row,col;
	for(col=1;col<n;col++)
	{
		for(row=0;row<n-col;row++)
		{
			m(row,col) = m(row+1,col-1) - m(row,col-1);
		}
	}
	return m;
}

template <typename T>
void getDeltaFXsFwd(T* arr,CMatrix<T>& m)
{
	for(int i=0; i < m.row ;i++)
	{
		arr[i] = m(0,i);
	}
}

template <typename T>
T ufactFwd(int n,T uk)
{
	int i;T fl = uk;
	for(i=0;i<n;i++)
	{
		fl *= (uk-(T)i);
	}
	return fl;
}

template <typename T>
T fnd_UthFwd(int k,T uk ,T dny0,T d0y0)
{
	return uk - ( ufactFwd(k,uk) / factorial(k+1) ) * ( dny0 / d0y0 ) ;
}

template <typename T>
T inv_interpolFwd(T y,T * dfxs,T* xs,int n)
{
	int i; T u1;
	if(n>2)u1 = (y-dfxs[0])/dfxs[1];
	else return xs[0];
	T uk = u1, 
	u_1 = u1, 
	h = xs[1] - xs[0];
	for(i=1;i<n;i++)
	{
		uk = fnd_UthFwd(i,uk,dfxs[i+1],dfxs[1]);
		if( uk == u_1 )
		{
			return xs[0] + uk*h;
		}
		u_1 = uk;
	}
	return xs[0] + uk*h;
}

template <typename T>
void getDeltaFXsBck(T* arr,CMatrix<T>& m)
{
	for(int i=0; i < m.row ;i++)
	{
		arr[i] = m(m.row-1-i,i);
	}
}

template <typename T>
T ufactBck(int n,T uk)
{
	int i;T fl = uk;
	for(i=0;i<n;i++)
	{
		fl *= (uk+(T)i);
	}
	return fl;
}

template <typename T>
T fnd_UthBck(int k,T uk ,T dny0,T d0y0)
{
	return uk - ( ufactBck(k,uk) / factorial(k+1) ) * ( dny0 / d0y0 ) ;
}

template <typename T>
T inv_interpolBck(T y,T * dfxs,T* xs,int n)
{
	int i; T u1;
	if(n>2)u1 = (y-dfxs[0])/dfxs[1];
	else return xs[n-1];
	T uk = u1, 
	u_1 = u1, 
	h = xs[1] - xs[0];
	for(i=1;i<n;i++)
	{
		uk = fnd_UthBck(i,uk,dfxs[i+1],dfxs[1]);
		if( uk == u_1 )
		{
			return xs[n-1] + uk*h;
		}
		u_1 = uk;
	}
	return xs[n-1] + uk*h;
}

template <typename T>
void getDeltaFXsMid(T* arr,CMatrix<T>& m)
{
	int mid = ((m.row+1)>>1)%m.row;
	for(int i=0; i < m.row ;i++)
	{
		arr[i] = m(mid,i);
		mid = ((mid+1)>>1)%m.row;//m.col
	}
}

template <typename T>
T ufactMid(int n,T uk)
{
	int i;T fl = uk;
	for(i=0;i<n;i++)
	{
		fl *= (uk-(T)i);
	}
	if(n>1)fl*=(uk+(T)1);
	return fl;
}

template <typename T>
T fnd_UthMid(int k,T uk ,T dny0,T d0y0)
{
	return uk - ( ufactMid(k,uk) / factorial(k+1) ) * ( dny0 / d0y0 ) ;
}

template <typename T>
T inv_interpolMid(T y,T * dfxs,T* xs,int n)
{
	int i;T u1;
	if(n>2) u1  = (y-dfxs[0])/dfxs[1];
	else return xs[0];
	T uk = u1, 
	u_1 = u1, 
	h = xs[1] - xs[0];
	for(i=1;i<n;i++)
	{
		uk = fnd_UthMid(i,uk,dfxs[i+1],dfxs[1]);
		if( uk == u_1 )
		{
			return xs[((n+1)>>1)] + uk*h;
		}
		u_1 = uk;
	}
	return xs[((n+1)>>1)] + uk*h;
}

int count=0;
int bs(int *a,int low,int high,int number);

int NewtonSearch(int *arr,int n, int value)
{
	int high = n-1;
	int low = 0;
	while(high-low > 10 ) 
	{
		count++;
		int qua = (high-low)/3;
		int mid1 = low+qua;
		int mid2 = low+qua*2;
		
		
		double s = (high-low)/10;
		int nx = (double)high/s;
		double dfxs[nx+1];
		memset(dfxs,0,sizeof(double)*(nx+1));
		double x[nx+1] ;
		double fx[nx+1] ;
		memset(x,0,sizeof(double)*(nx+1));
		memset(fx,0,sizeof(double)*(nx+1));
		double i;int cnt = 0;
		for(i=low;i<=high&&cnt<nx;i+=s)
		{
			x[cnt] = i;
			fx[cnt++] = (double)arr[(int)i];
		}
		x[cnt] = i;
		fx[cnt++] = (double)arr[(int)i];
		CMatrix <double> m = getDiffTable(fx,cnt);
		int index;
		if( value < arr[mid1] )
		{
			getDeltaFXsFwd(dfxs,m);
			index = (int)inv_interpolFwd((double)value,dfxs,x,cnt);
		}
		else if( value < arr[mid2] )
		{
			getDeltaFXsMid(dfxs,m);
			index = (int)inv_interpolMid((double)value,dfxs,x,cnt);
		}
		else
		{
			getDeltaFXsBck(dfxs,m);
			index = (int)inv_interpolBck((double)value,dfxs,x,cnt);
		}
		if(arr[index] == value) return index;
		if(arr[index] < value) low = index + 1;
		if(arr[index] > value) high = index -1;
	}
	return bs(arr,low,high,value);
}

void sort(int* a,int size)
{
    int i,p,k;
    for(i=0;i<=size;i++)
    {
        for(p=0;p<=size;p++)
        {
            if(a[i]<a[p])
            {
                k=a[i];
                a[i]=a[p];
                a[p]=k;
            }
        }
    }
}

int bs(int *a,int low,int high,int number)
{
    int mid;
    if(number>a[high]||number<a[low]) return count;
    if(a[low] == number) return count;
    if(a[high]== number) return count;
    while(low<=high)
    {
    	count++;
        mid = (high+low)>>1;
        if(a[mid] == number) return mid;
        else if(a[mid]> number) high = mid-1;
        else low = mid+1;
	}
    return -1;
}

int main()
{
	
	int n = 6;
	double fx[n] = {1,4,9,16,25,36}  , x[n] = {1,2,3,4,5,6}, dfxs[n];
	//double fx[n] = {1.1471202,1.1476828,1.1482466,1.1488115,1.1493776,1.1499448} , x[n] = {0.536,0.537,0.538,0.539,0.540,0.541},dfxs[n];
	int i;
	for(i=0;i<n;i++)
	{
		//printf("Enter the value of f(%d) : ",i);
		//scanf("%d",&arr[i]);
	}
	
	CMatrix < double> m = getDiffTable(fx,n);
	
	m.print();
	
	getDeltaFXsMid(dfxs,m);
	double v;
	printf("Enter the value to Find : ");scanf("%lf",&v);
	
	printf("fi =  %G ",inv_interpolMid(( double)v,dfxs,x,n));
	
	/*int n = 1000;
	int b ,i;
    int a[n];
    srand(time(NULL));
    for(i=0;i<n;i++)
    {
        a[i] = (rand()+rand())%(n*10);
        if(a[i]<0)a[i] = a[i]*-1;
    }
    sort(a,n);
    for(i=0;i<n;i++)
    {
        //printf(" %d:%d, ",i,a[i]);
    }
    for(i=0;i<20;i++)
    {
    	/*int v;
    	printf("\n\nEnter the value to Find : ");scanf("%d",&v);
    	printf("\ncount = %d,fi =  %d",count,NewtonSearch(a,n,v));*/
    	/*int pos = rand()%n;
    	b = a[pos];
    	count = 0;
        printf("\n\nNumber %d Found by bs:ns in %d ",b, count , bs(a,0,n,b));
        count = 0;
        printf(": %d iterations at index %d",count,NewtonSearch(a,n,b));
    }*/
	
}
