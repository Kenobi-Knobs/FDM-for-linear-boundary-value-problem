#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

//define param
double a = 0.4, b = 0.7; // [a, b]
double alfa1 = 1, alfa2 = 0, A = 2; // alfa1 * y(a) + alfa2 * y'(a) = A
double beta1 = 1, beta2 = 2, B = 0.7;// beta1 * y(b) + beta2 * y'(b) = B

//p(x)
double p(double x){
	return -3;
}
//q(x)
double q(double x){
	return 1/x;
}
//f(x)
double f(double x){
	return 1;
}

double ai(int n, double h, double pi, int i){
	if(i < n){
		return 1 - ((h*pi)/2);
	}else{
		return -beta2;
	}
}

double bi(int n, double h, double qi, int i){
	if(i < n && i > 0){
		return h*h * qi - 2;
	}else if(i == 0){
		return h * alfa1 - alfa2;
	}else if(i == n){
		return h * beta1 + beta2; 
	}
		
}

double ci(int n, double h, double pi, int i){
	if(i < n && i > 0){
		return 1 + ((h*pi)/2);
	}else if(i == 0){
		return alfa2;
	}
		
}

double di(int n, double h, double fi, int i){
	if(i < n && i > 0){
		return h*h * fi;
	}else if(i == 0){
		return h * A;
	}else if(i == n){
		return h * B; 
	}
}

 
void forward(int n, double **a, double *b)
{
    double v;
    for(int k = 0,i,j,im; k < n - 1; k++)
    {
        im = k;
        for(i = k + 1; i < n; i++)
        {
            if(fabs(a[im][k]) < fabs(a[i][k]))
            {
                im = i;
            }
        }
        if(im != k)
    	{
            for(j = 0; j < n; j++)
            {
                v = a[im][j];
                a[im][j] = a[k][j];
                a[k][j]  = v;
            }
            v     = b[im];
            b[im] = b[k];
            b[k]  = v;
        }
        for(i = k + 1; i < n; i++)
        {
            v = 1.0*a[i][k]/a[k][k];
            a[i][k] = 0;
            b[i]    = b[i] - v*b[k];
            if(v != 0)
            for(j = k + 1; j < n; j++)
            {
                a[i][j] = a[i][j] - v*a[k][j];
            }
        }
    }
}
 
void reverse(int n, double **a, double *b, double *x)
{
    double s = 0;
    x[n - 1] = 1.0*b[n - 1]/a[n - 1][n - 1];
    for(int i = n - 2, j; 0 <= i; i--)
    {
        s = 0;
        for(j = i + 1; j < n; j++)
        {
            s = s+a[i][j]*x[j];
        }
        x[i] = 1.0*(b[i] - s)/a[i][i];
    }
}

double *finite(int n){
	//define matrix
    double **matrA = new double*[n+1];
    for(int i=0; i<n+1; i++){
            matrA[i]=new double[n+1];
    }
	double *matrB = new double[n+1];	
    double *X = new double[n+1]; //vector X[n+1]
	//eval h
	double h = (b - a) / n;
	//cout << "h = " << h << endl;
	
	
	//eval pi, qi, fi
	double mp[n];
	double mq[n];
	double mf[n];
	
	double x = a;
	for(int i = 0; i < n; i++){
		mp[i] = p(x);
		mq[i] = q(x);
		mf[i] = f(x);
		x += h;
	}
	
	// zeroing
	for(int i = 0; i <= n; i++){
		for(int j = 0; j <= n; j++){
			matrA[i][j] = 0;
		}
	}
	
	// fill matrix free element
	for(int i = 0; i <= n; i++){
		matrB[i] = di(n, h, mf[i], i);
	}
	// fill 0 row
	matrA[0][0] = bi(n, h, mq[0], 0);
	matrA[0][1] = ci(n, h, mp[0], 0);
	
	// fill i row
	for(int i = 1; i <= n-1; i++){
		matrA[i][i-1] =	ai(n, h,mp[i],i);
		matrA[i][i] = bi(n, h, mq[i], i);
		matrA[i][i+1] =	ci(n, h, mp[i], i);
	}
	// fill n row
	matrA[n][n-1] = ai(n, h,mp[n],n);
	matrA[n][n] = bi(n, h, mq[n], n);
/*	
	cout << "System of linear equations: \n";
	
	for(int i = 0; i <= n; i++){
		for(int j = 0; j <= n; j++){
			cout << fixed << setprecision(3) <<  matrA[i][j] << "\t";
		}
		cout << " |= "<< matrB[i] << endl;
	}
*/	

    forward(n+1, matrA, matrB);
    reverse(n+1, matrA, matrB, X);

	return X;
}

int main(){
	int n;
	int flag = 0;
	double eps = 0.0001;
	int counter = 0;
	int mode;
	
	cout << "select mode: \n1)default\n2)accuracy\n";
	cout << "mode: ";
	cin >> mode;
	cout << "n = ";
	cin >> n;
	
	if(mode > 1){
		while(flag == 0){
			int altN = n * 2;
			double *X = new double[n+1];
			double *X1 = new double[altN+1];
	
			X = finite(n);
			X1 = finite(altN);
			
			for(int i = 0; i < n; i++){
				if(((fabs(X[i]-X1[i*2]))/15) > eps){
					counter++;	
				}
			}
			
			if(counter == 0){
				cout << "_______________________________" << endl;
				cout << "Given accuracy is reached on:" << endl;
				cout << "h = " << (b - a) / n << endl;
				cout << "n = " << n << endl;
				cout << "accuracy = " << eps << endl;
				printf("Results :\r\n");
				cout << "X" << '\t' << "Y" << '\t' << "|y0 - y1| / 15" << endl;
				cout << "_______________________________" << endl;
			    double x = a;
				for(int i = 0; i <= n; i++){
					cout << fixed << setprecision(5) << x << "\t|" << X[i] << "|\t" << ((fabs(X[i]-X1[i*2]))/15) << "|" << endl;
					x += (b - a) / n;
				}
				cout << "_______________________________" << endl;
				flag = 1;
			}
			else{
				n = altN;
				counter = 0;
				delete X;
				delete X1;
			}
		}
	}
	if(mode == 1){
		double *X = new double[n+1];
		X = finite(n);
		printf("Results :\r\n");
		double x = a;
		cout << "  x\t  y" << endl;
		for(int i = 0; i <= n; i++){
			cout << fixed << setprecision(5) << x << '\t' << X[i] << endl;
			x += (b - a) / n;
		}	
	}
	
	system("pause");
	return 0;
}


