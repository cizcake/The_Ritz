// Метод Ритца.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "fstream"
#include "iostream" 
#include "math.h"
using namespace std;
const int n = 34;// число узлов на каждой стороне 
double h = 1. / (n - 1); // шаг метода
const int m = (n-2)*(n-2); // 
float A[m][m];/*матрица
			   m строк 
			   m столбцов */
float R[m];//правая часть
float x[m]; // вектор решений 
void gauss(float a[m][m], float y[m],float x[m]) // Метод Гаусса с выбором Главного элемента 
{ 
	
	float max;
	int k, index;
	const double eps = 0.0001;  // точность
	k = 0;
	while (k < m)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < m; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
		}
		for (int j = 0; j < m; j++)
		{
			float temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		float temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < m; i++)
		{
			float temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < m; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = m - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
}
void Check(float a[m][m], float x[m], float fv[m])// проверка решения 
{
	int i, j;
	float sum;
	for (i = 0; i<m; i++)
	{
		sum = 0.0;
		for (j = 0; j<m; j++) sum = sum + a[i][j] * x[j];
		if (fabs(fv[i] - sum) >= 1e-4) { printf("Error") ;   }
	}
}

int main()
{
	system("chcp 1251");
	system("cls");
	ofstream f("input.txt");
	ofstream g("right.txt");
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			  A[i][j] = 0;
		}
	}
		for(int i=1;i<m;i++)
			for (int j = 1; j < m; j++)
			{ /* Обработка Угловых Ячеек*/
				if (i == 0 && j == 0)
				{
					A[i][j] += 2 / 3.;
				}
				else if (i == n  && j == 0)
				{
					A[i][j] += 2 / 3.;
				}

				else if (i == n  && j == n)
				{
					A[i][j] += 2 / 3.;
				}

				else if (i == 0 && j == n)
				{
					A[i][j] += 2 / 3.;
				}
				/* Обработка Граничных Ячеек*/
				else if (i == 0 && j > 0 && j < n)
				{
					A[(i - 1)*(n - 2) + j - 1][(i - 1)*(n - 2) + j - 1] += 2 / 3.;
					A[(i - 2)*(n - 2) + j - 1][(i - 2)*(n - 2) + j - 1] += 2 / 3.;
					A[(i - 1)*(n - 2) + j - 1][(i - 2)*(n - 2) + j - 1] += -1 / 6.;
					A[(i - 2)*(n - 2) + j - 1][(i - 1)*(n - 2) + j - 1] += -1 / 6.;

				}
				else  if (j == 0 && i > 0 && i < n)
				{
					A[(i - 2)*(n - 2) + j][(i - 2)*(n - 2) + j] += 2 / 3.;
					A[(i - 2)*(n - 2) + j][(i - 2)*(n - 2) + j - 1] += -1 / 6.;
					A[(i - 2)*(n - 2) + j - 1][(i - 2)*(n - 2) + j] += -1 / 6.;
					A[(i - 2)*(n - 2) + j - 1][(i - 2)*(n - 2) + j - 1] += 2 / 3.;

				}


				else if (i == n && j > 0 && j<n)
				{
					A[(i - 1)*(n - 2) + j][(i - 1)*(n - 2) + j] += 2 / 3.;
					A[(i - 1)*(n - 2) + j][(i - 2)*(n - 2) + j] += -1 / 6.;
					A[(i - 2)*(n - 2) + j][(i - 1)*(n - 2) + j] += -1 / 6.;
					A[(i - 2)*(n - 2) + j][(i - 2)*(n - 2) + j] += 2 / 3.;

				}

				else if (j == n && i>0 && i < n)

				{
					A[(i - 1)*(n - 2) + j][(i - 1)*(n - 2) + j] += 2 / 3.;
					A[(i - 1)*(n - 2) + j][(i - 1)*(n - 2) + j - 1] += -1 / 6.;
					A[(i - 1)*(n - 2) + j - 1][(i - 1)*(n - 2) + j] += -1 / 6.;
					A[(i - 1)*(n - 2) + j - 1][(i - 1)*(n - 2) + j - 1] += 2 / 3.;
				}
				/*Обработка Внутренних Ячеек*/
				else if (i > 1 && i < n - 1 && j>1 && j < n - 1)
				{
					A[(i - 1)*(n - 2) + j][(i - 1)*(n - 2) + j] += 2/ 3.;
					A[(i - 1)*(n - 2) + j][(i - 2)*(n - 2) + j] += - 1 / 6.;
					A[(i - 1)*(n - 2) + j][(i - 2)*(n - 2) + j - 1] += - 1 / 3.;   // ФИ 1
					A[(i - 1)*(n - 2) + j][(i - 1)*(n - 2) + j - 1] += -1 / 6.;

					A[(i - 2)*(n - 2) + j][(i - 1)*(n - 2) + j] += - 1 / 6.;
					A[(i - 2)*(n - 2) + j][(i - 2)*(n - 2) + j] += 2 / 3.;
					A[(i - 2)*(n - 2) + j][(i - 2)*(n - 2) + j - 1] += - 1 / 6.;   // Фи 2
					A[(i - 2)*(n - 2) + j][(i - 1)*(n - 2) + j - 1] += - 1 / 3.;

					A[(i - 2)*(n - 2) + j - 1][(i - 1)*(n - 2) + j] += - 1 / 3.;
					A[(i - 2)*(n - 2) + j - 1][(i - 2)*(n - 2) + j] += - 1 / 6.;
					A[(i - 2)*(n - 2) + j - 1][(i - 2)*(n - 2) + j - 1] += 2 / 3.;   // Фи 3
					A[(i - 2)*(n - 2) + j - 1][(i - 1)*(n - 2) + j - 1] += - 1 / 6.;

					A[(i - 1)*(n - 2) + j - 1][(i - 1)*(n - 2) + j] += -1 / 6. ;
					A[(i - 1)*(n - 2) + j - 1][(i - 2)*(n - 2) + j] += -1 / 3.;
					A[(i - 1)*(n - 2) + j - 1][(i - 2)*(n - 2) + j - 1] += -1 / 6.;    // Фи 4
					A[(i - 1)*(n - 2) + j - 1][(i - 1)*(n - 2) + j - 1] += 2 / 3.;
				}
			};

		
		for (int i = 1; i<m; i++)
		{
			for (int j = 1; j<m; j++)
			{
				f << A[i][j] << " ";     // Запись Матрицы элементов в файл 
			}
			f << '\n';
		}
		
	
		for (int i = 1; i <= m; i++)
		{
			R[i] = 0; //Формирование столбца правой части 
		};

		for (int i = 2; i<m; i++)
			for (int j = 2; j<m; j++)
			{
				R[(i - 1)*(n - 2) + j] += (1 - 2 * h*(-2 + i + j) - 4 * h*h*(-1 + i - i*j + j)) / (4* h*h);
				R[(i - 2)*(n - 2) + j] += (-1 + 4 * h*h*(i - i*j) + 2 * h*(-1 + i + j)) / (4 * h*h);
				R[(i - 2)*(n - 2) + j - 1] += (1 + 4 * h*h*i*j - 2 * h*(i + j)) / (4 * h*h);
				R[(i - 1)*(n - 2) + j - 1] += (-1 - 4 * h*h*(i*j - j) + 2 * h*(-1 + i + j)) / (4 * h*h);

	} 
		for (int i = 0; i < m; i++)
		{
			 g << R[i] << "\n"; // Запись правой части системы в файл
		}
		gauss(A, R, x); // Использование метода 
		for (int i = 0; i < m; i++)
			cout << "x[" << i << "]=" << x[i] << endl;// вывод решения 
			
		Check(A, x, R); // проверка решения 
		g.close();
		f.close(); 
		getchar();
		getchar();
    return 0;
}

