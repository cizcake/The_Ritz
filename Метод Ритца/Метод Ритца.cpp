// ����� �����.cpp: ���������� ����� ����� ��� ����������� ����������.
//

#include "stdafx.h"
#include "fstream"
#include "iostream" 
#include "math.h"
using namespace std;
const int n = 34;// ����� ����� �� ������ ������� 
double h = 1. / (n - 1);
const int m = (n-2)*(n-2);
double A[m][m];//�������
double R[m];//������ �����
double x[m];
void gauss(double a[m][m], double y[m],double x[m])
{ 
	
	double max;
	int k, index;
	const double eps = 0.1;  // ��������
	k = 0;
	while (k < m)
	{
		// ����� ������ � ������������ a[i][k]
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
		// ������������ �����
		if (max < eps)
		{
			// ��� ��������� ������������ ���������
			cout << "������� �������� ���������� ��-�� �������� ������� ";
			cout << index << " ������� A" << endl;
		}
		for (int j = 0; j < m; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// ������������ ���������
		for (int i = k; i < m; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // ��� �������� ������������ ����������
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // ��������� �� �������� ���� �� ����
			for (int j = 0; j < m; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// �������� �����������
	for (k = m - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
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
			{ 
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
				//�������
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
				else if (i > 1 && i < n - 1 && j>1 && j < n - 1)
				{
					A[(i - 1)*(n - 2) + j][(i - 1)*(n - 2) + j] += 2/ 3.;
					A[(i - 1)*(n - 2) + j][(i - 2)*(n - 2) + j] += - 1 / 6.;
					A[(i - 1)*(n - 2) + j][(i - 2)*(n - 2) + j - 1] += - 1 / 3.;   // �� 1
					A[(i - 1)*(n - 2) + j][(i - 1)*(n - 2) + j - 1] += -1 / 6.;

					A[(i - 2)*(n - 2) + j][(i - 1)*(n - 2) + j] += - 1 / 6.;
					A[(i - 2)*(n - 2) + j][(i - 2)*(n - 2) + j] += 2 / 3.;
					A[(i - 2)*(n - 2) + j][(i - 2)*(n - 2) + j - 1] += - 1 / 6.;   // �� 2
					A[(i - 2)*(n - 2) + j][(i - 1)*(n - 2) + j - 1] += - 1 / 3.;

					A[(i - 2)*(n - 2) + j - 1][(i - 1)*(n - 2) + j] += - 1 / 3.;
					A[(i - 2)*(n - 2) + j - 1][(i - 2)*(n - 2) + j] += - 1 / 6.;
					A[(i - 2)*(n - 2) + j - 1][(i - 2)*(n - 2) + j - 1] += 2 / 3.;   // �� 3
					A[(i - 2)*(n - 2) + j - 1][(i - 1)*(n - 2) + j - 1] += - 1 / 6.;

					A[(i - 1)*(n - 2) + j - 1][(i - 1)*(n - 2) + j] += -1 / 6. ;
					A[(i - 1)*(n - 2) + j - 1][(i - 2)*(n - 2) + j] += -1 / 3.;
					A[(i - 1)*(n - 2) + j - 1][(i - 2)*(n - 2) + j - 1] += -1 / 6.;    // �� 4
					A[(i - 1)*(n - 2) + j - 1][(i - 1)*(n - 2) + j - 1] += 2 / 3.;
				}
			};


		for (int i = 1; i<m; i++)
		{
			for (int j = 1; j<m; j++)
			{
				f << A[i][j] << " ";
			}
			f << '\n';
		}
		//������ �����
	
		for (int i = 1; i <= m; i++)
		{
			R[i] = 0;
		};

		for (int i = 1; i<m; i++)
			for (int j = 1; j<m; j++)
			{
				R[(i - 1)*(n - 2) + j] += (1. - 2. * h*(-2. + i + j) - 4. * h*h*(-1 + i - i*j + j)) / (4. * h*h);
				R[(i - 2)*(n - 2) + j] += (-1 + 4 * h*h*(i - i*j) + 2 * h*(-1 + i + j)) / (4 * h*h);
				R[(i - 2)*(n - 2) + j - 1] += (1 + 4 * h*h*i*j - 2 * h*(i + j)) / (4 * h*h);
				R[(i - 1)*(n - 2) + j - 1] += (-1 - 4 * h*h*(i*j - j) + 2 * h*(-1 + i + j)) / (4 * h*h);

	}
		for (int i = 0; i < m; i++)
		{
			 g << R[i] << "\n"; // �� ����� ����� ������� ��������� ��� �������� �������
		}
		gauss(A, R,x);
		for (int i = 0; i < m; i++)
			cout << "x[" << i << "]=" << x[i] << endl;// ����� ������� 
			
		g.close();
		f.close(); 
		getchar();
		getchar();
    return 0;
}

